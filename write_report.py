from __future__ import print_function
import argparse
import gzip
import itertools
import json
import numpy as np
import os
import zipfile
from collections import defaultdict
from scipy import stats

class ResultLoader(object):
  def __init__(self, tree_summary_fn, mutation_list_fn, mutation_assignment_fn, min_ssms):
    self._tree_summary_fn = tree_summary_fn
    self._mutation_list_fn = mutation_list_fn
    self._mutation_assignment_fn = mutation_assignment_fn
    self._min_ssms = min_ssms

    self.mutlist = None
    self.tree_summary = None
    self._subclone_idx_map = defaultdict(dict)

    self._load_tree_data()

  def _convert_keys_to_ints(self, dic):
    keys = dic.keys()
    for key in dic.keys():
      dic[int(key)] = dic[key]
      del dic[key]

  def _load_tree_data(self):
    with gzip.GzipFile(self._tree_summary_fn) as treesummf:
      self.tree_summary = json.load(treesummf)['trees']
    self._convert_keys_to_ints(self.tree_summary)
    for tree_idx, tree_features in self.tree_summary.items():
      self._convert_keys_to_ints(tree_features['populations'])
      self._convert_keys_to_ints(tree_features['structure'])

    with gzip.GzipFile(self._mutation_list_fn) as mutlistf:
      self.mutlist = json.load(mutlistf)
    self.num_ssms = len(self.mutlist['ssms'])

    for tree_idx, tree_features in self.tree_summary.items():
      small_nodes = self._remove_small_nodes(tree_idx, tree_features['populations'])
      self._remove_nodes_from_tree_structure(small_nodes, tree_features['structure'])
      self._renumber_nodes(tree_idx)

      # Though we don't use mutass now, this step is necessary to correct the
      # number of SSMs and CNVs in each node.
      mutass = self.load_mut_assignments(tree_idx)
      for sidx, subclone in tree_features['populations'].items():
        # Note that only mutass entries for subclones with (> 0 SSMs or > 0
        # CNVs) will exist. Thus, no mutass entry will exist for node 0, as it
        # never has SSMs or CNVs.
        if sidx == 0:
          continue
        for mut_type in ('ssms', 'cnvs'):
          subclone['num_%s' % mut_type] = len(mutass[sidx][mut_type])

  def _renumber_nodes(self, tree_idx):
    subclone_idxs = self.tree_summary[tree_idx]['populations'].keys()

    num_removed = 0
    for subclone_idx in range(1, max(subclone_idxs) + 1):
      if subclone_idx not in self.tree_summary[tree_idx]['populations']:
        num_removed += 1
        continue
      if num_removed > 0:
        self._subclone_idx_map[tree_idx][subclone_idx] = subclone_idx - num_removed

    # By proceeding in sorted order, we guarantee we're not overwriting a
    # single element twice, which would give the wrong values.
    for subclone_idx in sorted(subclone_idxs):
      if subclone_idx not in self._subclone_idx_map[tree_idx]:
        continue
      new_idx = self._subclone_idx_map[tree_idx][subclone_idx]

      self.tree_summary[tree_idx]['populations'][new_idx] = self.tree_summary[tree_idx]['populations'][subclone_idx]
      del self.tree_summary[tree_idx]['populations'][subclone_idx]

      if subclone_idx in self.tree_summary[tree_idx]['structure']:
        self.tree_summary[tree_idx]['structure'][new_idx] = self.tree_summary[tree_idx]['structure'][subclone_idx]

    # We must also renumber children in the structure -- just renumbering
    # parents isn't enough.
    for subclone_idx in sorted(subclone_idxs):
      # Node has no children.
      if subclone_idx not in self.tree_summary[tree_idx]['structure']:
        continue
      children = self.tree_summary[tree_idx]['structure'][subclone_idx]
      children = [
        self._subclone_idx_map[tree_idx][c]
        if c in self._subclone_idx_map[tree_idx]
        else c
        for c in children
      ]
      self.tree_summary[tree_idx]['structure'][subclone_idx] = children

  def _remove_small_nodes(self, tree_idx, populations):
    small_nodes = set()

    subclone_idxs = sorted(populations.keys())
    last_phi = None
    last_idx = None

    for subclone_idx in subclone_idxs:
      for p, children in self.tree_summary[tree_idx]['structure'].items():
        if subclone_idx in children:
          parent = p
          break
      subclone = populations[subclone_idx]
      # Ensure this node's phi is <= the phi of its preceding sibling node, if any exists.
      if subclone_idx > 0 and last_idx in self.tree_summary[tree_idx]['structure'][parent]:
        assert subclone['phi'] <= last_phi
      last_phi = subclone['phi']
      last_idx = subclone_idx

      if subclone_idx == 0 or subclone['num_ssms'] >= self._min_ssms:
        continue
      small_nodes.add(subclone_idx)

    for node_idx in small_nodes:
      del populations[node_idx]
      # Mark node as removed.
      self._subclone_idx_map[tree_idx][node_idx] = None

    return small_nodes

  def _remove_nodes_from_tree_structure(self, removed, tree_structure):
    removed = set(removed)

    for parent, children in tree_structure.items():
      to_remove = []
      to_add = []

      for child in children:
        if child not in removed:
          continue
        if child in tree_structure:
          grandchildren = tree_structure[child]
          to_add += grandchildren
        to_remove.append(child)

      tree_structure[parent] = [
        node for node in (tree_structure[parent] + to_add)
        if node not in removed
      ]
      tree_structure[parent].sort()

    for rem in removed:
      if rem in tree_structure:
        del tree_structure[rem]

  def _move_muts_to_best_node(self, muts, mutass, populations):
    for mut_type in ('ssms', 'cnvs'):
      for mut in muts[mut_type]:
        mut_id = mut['id']
        mut_stats = self.mutlist[mut_type][mut_id]
        ref_reads = np.mean(mut_stats['ref_reads'])
        total_reads = np.mean(mut_stats['total_reads'])
        implied_phi = 2 * (total_reads - ref_reads) / total_reads
        implied_phi = min(implied_phi, 1.0)

        lowest_phi_delta = 1
        best_node = None
        for pidx, pop in populations.items():
          phi_delta = abs(pop['phi'] - implied_phi)
          # Don't allow assignments to the non-cancerous root node.
          if phi_delta < lowest_phi_delta and pidx != 0:
            lowest_phi_delta = phi_delta
            best_node = pidx

        mutass[best_node][mut_type].append(mut)

  def _reassign_muts(self, tree_idx, mutass):
    deleted_muts = []

    for sidx in sorted(self._subclone_idx_map[tree_idx].keys()):
      new_idx = self._subclone_idx_map[tree_idx][sidx]
      if new_idx is None: # Node was removed.
        deleted_muts.append(mutass[sidx])
      else:
        mutass[new_idx] = mutass[sidx]
      del mutass[sidx]

    for dm in deleted_muts:
      self._move_muts_to_best_node(dm, mutass, self.tree_summary[tree_idx]['populations'])

  def _renumber_mut_nodes(self, tree_idx, mutass):
    for sidx in sorted(self._subclone_idx_map[tree_idx].keys()):
      new_idx = self._subclone_idx_map[tree_idx][sidx]
      mutass[new_idx] = mutass[sidx]
      del mutass[sidx]

  def _load_assignments(self, mutf, tree_idx):
    mutass = json.loads(mutf.read('%s.json' % tree_idx))
    mutass = mutass['mut_assignments']
    self._convert_keys_to_ints(mutass)
    self._reassign_muts(tree_idx, mutass)
    return mutass

  def load_mut_assignments(self, tree_idx):
    with zipfile.ZipFile(self._mutation_assignment_fn) as mutf:
      return self._load_assignments(mutf, tree_idx)

  def load_all_mut_assignments(self):
    with zipfile.ZipFile(self._mutation_assignment_fn) as mutf:
      for zinfo in mutf.infolist():
        tree_idx = int(zinfo.filename.split('.')[0])
        yield (tree_idx, self._load_assignments(mutf, tree_idx))

  def get_ssm_names(mutass):
    def _ssm_key(name):
      chrom, pos = name.split('_')
      chrom = chrom.lower()
      pos = int(pos)

      if chrom == 'x':
        chrom = 100
      elif chrom == 'y':
        chrom = 101
      else:
        chrom = int(chrom)
      return (chrom, pos)

    _, mut_assignments = next(load_mut_assignments(mutass))
    ssm_names = []
    for _, muts in mut_assignments.items():
      ssm_names += [m['name'] for m in muts['ssms']]
    ssm_names.sort(key = _ssm_key)

    idx_lookup = {name: i for (i, name) in enumerate(ssm_names)}
    return (ssm_names, idx_lookup)

class SubcloneStatsComputer(object):
  def __init__(self, tree_summary):
    self._tree_summary = tree_summary

  def calc(self):
    self._calc_global_stats()
    self._calc_pop_stats()

  def _find_clonal_node(self, pops):
    indices = [k for k in pops.keys() if k > 0]
    return min(indices)

  def _calc_global_stats(self):
    cancer_pop_counts = []
    cellularities = []

    for tree_idx, tree_features in self._tree_summary.items():
      pops = tree_features['populations']

      # Tree may not have any canceerous nodes left after removing nodes with <
      # 3 SSMs.  (Note that non-cancerous root node will always remain,
      # however.) In such cases, skip this tree.
      if len(pops) == 1:
        continue

      cancer_pop_counts.append(len(pops) - 1)
      clonal_idx = self._find_clonal_node(pops)
      cellularities.append(pops[clonal_idx]['phi'])

    self.cancer_pops = intmode(cancer_pop_counts)
    self.cellularity = np.mean(cellularities)

  def _calc_pop_stats(self):
    K = self.cancer_pops
    phis_sum, num_ssms_sum = np.zeros(K), np.zeros(K)
    trees_examined = 0

    for tree in self._tree_summary.values():
      pops = tree['populations']
      if len(pops) - 1 != K:
        continue
      vals = [(pops[pidx]['phi'], pops[pidx]['num_ssms']) for pidx in sorted(pops.keys()) if pidx != 0]
      phis, num_ssms = zip(*vals)
      phis_sum += phis
      num_ssms_sum += num_ssms
      trees_examined += 1

    self.phis = phis_sum / trees_examined
    self.num_ssms = num_ssms_sum / trees_examined
    self.num_ssms = self._round_to_int(self.num_ssms)

  def _round_to_int(self, arr, dtype=np.int64):
    '''Round each array element to the nearest integer, while distributing the
    rounding errors equitably amongst members.'''
    rounded = np.rint(arr)
    rerr = np.sum(rounded - arr)
    assert np.isclose(np.rint(rerr), rerr)

    # Correct rounding error by subtracting 1 from the largest node, then
    # moving on to the next and continuing until no error remains. Assuming
    # that elements in arr sum to an integer, the accumulated rounding error
    # will always be an integer, and will always be <= 0.5*len(arr).
    biggest_idxs = list(np.argsort(rounded))
    while not np.isclose(rerr, 0):
      biggest_idx = biggest_idxs.pop()
      if rerr > 0:
        rounded[biggest_idx] -= 1
      else:
        rounded[biggest_idx] += 1
      rerr = np.sum(rounded - arr)

    rounded = rounded.astype(dtype)
    assert np.isclose(np.sum(arr), np.sum(rounded))
    return rounded

class ClusterMembershipComputer(object):
  def __init__(self, loader, subclone_stats):
    self._loader = loader
    self._subclone_stats = subclone_stats

  def calc(self):
    assignments = {}

    # Binomial assignment won't work for areas of altered copy number, where
    # phi/2 is no longer the correct fraction of chromatids carrying the
    # mutation.
    for ssm_id, ssm in self._loader.mutlist['ssms'].items():
      d, a  = np.mean(ssm['total_reads']), np.mean(ssm['ref_reads'])
      ssm_idx = int(ssm_id[1:])

      best_pop = None
      best_prob = float('-inf')
      for cancer_pop_idx, phi in enumerate(self._subclone_stats.phis):
        prob = stats.binom.logpmf(d - a, d, phi / 2)
        if prob > best_prob:
          best_prob = prob
          best_pop = cancer_pop_idx
      assignments[ssm_idx] = best_pop

    idxs = sorted(assignments.keys())
    for ssm_idx in idxs:
      yield ('s%s' % ssm_idx, assignments[ssm_idx] + 1)

class CoassignmentComputer(object):
  def __init__(self, loader):
    self._loader = loader

  def compute_coassignments(self):
    num_ssms = self._loader.num_ssms
    coass = np.zeros((num_ssms, num_ssms))
    num_trees = 0

    for tree_idx, mut_assignments in self._loader.load_all_mut_assignments():
      num_trees += 1
      for subclone_idx, muts in mut_assignments.items():
        ssms = muts['ssms']
        for ssm1, ssm2 in itertools.combinations(ssms, 2):
          ssm1_idx = int(ssm1['id'][1:])
          ssm2_idx = int(ssm2['id'][1:])
          coass[ssm1_idx, ssm2_idx] += 1.0
          coass[ssm2_idx, ssm1_idx] += 1.0

    coass /= num_trees
    # Set diagonal to 1
    coass += np.diag(np.ones(num_ssms))
    return coass

class SsmRelationComputer(object):
  def __init__(self, loader):
    self._loader = loader

  class VertexRelation(object):
    ancestor_desc = 1
    desc_ancestor = 2
    sibling       = 3

  def _determine_vertex_relations(self, tree_structure):
    relations = defaultdict(lambda: self.VertexRelation.sibling)

    def _mark_desc(par, desc):
      for descendant in desc:
        relations[(par, descendant)] = self.VertexRelation.ancestor_desc
        relations[(descendant, par)] = self.VertexRelation.desc_ancestor
        if descendant in tree_structure:
          _mark_desc(par, tree_structure[descendant])

    for parent, children in tree_structure.items():
      _mark_desc(parent, children)

    return relations

  def compute_ancestor_desc(self):
    num_ssms = self._loader.num_ssms
    ancestor_desc = np.zeros((self._loader.num_ssms, num_ssms))
    num_trees = 0

    for tree_idx, mut_assignments in self._loader.load_all_mut_assignments():
      num_trees += 1
      vert_relations = self._determine_vertex_relations(self._loader.tree_summary[tree_idx]['structure'])

      ssm_assignment_map = {}
      for subclone_idx, muts in mut_assignments.items():
        for ssm in muts['ssms']:
          ssm_assignment_map[ssm['id']] = subclone_idx

      # Working with a list is faster than a dict, despite the former having O(n)
      # lookup and the latter O(1). The Python list is also faster than a NumPy
      # array. See
      # http://stackoverflow.com/questions/22239199/numpy-array-indexing-and-or-addition-seems-slow
      ssm_assignments = []
      last_idx = -1
      for ssm_idx in sorted([int(s[1:]) for s in ssm_assignment_map.keys()]):
        assert ssm_idx == last_idx + 1
        last_idx = ssm_idx
        ssm_assignments.append(ssm_assignment_map['s%s' % ssm_idx])
      assert len(ssm_assignments) == num_ssms

      ssm_combos = itertools.combinations(range(num_ssms), 2)
      for ssm1_idx, ssm2_idx in ssm_combos:
        ssm1_subclone = ssm_assignments[ssm1_idx]
        ssm2_subclone = ssm_assignments[ssm2_idx]

        if ssm1_subclone == ssm2_subclone:
          continue
        if vert_relations[(ssm1_subclone, ssm2_subclone)] == self.VertexRelation.ancestor_desc:
          ancestor_desc[ssm1_idx, ssm2_idx] += 1.0
        if vert_relations[(ssm1_subclone, ssm2_subclone)] == self.VertexRelation.desc_ancestor:
          ancestor_desc[ssm2_idx, ssm1_idx] += 1.0

    ancestor_desc /= num_trees
    return ancestor_desc

class NodeRelationComputer(object):
  def __init__(self, loader, num_cancer_pops):
    self._loader = loader
    self._num_cancer_pops = num_cancer_pops

  def compute_relations(self):
    adj_matrix = np.zeros((self._num_cancer_pops + 1, self._num_cancer_pops + 1))

    for tree_idx, tree_features in self._loader.tree_summary.items():
      structure = tree_features['structure']
      # Only examine populations with mode number of nodes.
      if len(tree_features['populations']) - 1 != self._num_cancer_pops:
        continue
      for parent, children in tree_features['structure'].items():
        for child in children:
          adj_matrix[parent, child] += 1.0

    most_common_parents = adj_matrix.argmax(axis=0)
    return most_common_parents

def intmode(iter):
  return int(stats.mode(iter)[0][0])

def main():
  parser = argparse.ArgumentParser(
    description='Write SMC-Het Challenge outputs',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--min-ssms', dest='min_ssms', type=int, default=3,
    help='Minimum number of SSMs to retain a subclone')
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  parser.add_argument('output_dir',
    help='Directory in which to save Challenge outputs')
  args = parser.parse_args()

  loader = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment, args.min_ssms)

  ssc = SubcloneStatsComputer(loader.tree_summary)
  ssc.calc()
  with open(os.path.join(args.output_dir, '1A.txt'), 'w') as outf:
    print(ssc.cellularity, file=outf)
  with open(os.path.join(args.output_dir, '1B.txt'), 'w') as outf:
    print(ssc.cancer_pops, file=outf)
  with open(os.path.join(args.output_dir, '1C.txt'), 'w') as outf:
    for cluster_num, (num_ssms, phi) in enumerate(zip(ssc.num_ssms, ssc.phis)):
      print(cluster_num + 1, num_ssms, phi, sep='\t', file=outf)

  cmc = ClusterMembershipComputer(loader, ssc)
  with open(os.path.join(args.output_dir, '2A.txt'), 'w') as outf:
    # These will be in sorted order, so nth row refers to SSM with ID "s{n - 1}".
    for ssm_id, cluster in cmc.calc():
      print(cluster, file=outf)

  coassc = CoassignmentComputer(loader)
  coass_matrix = coassc.compute_coassignments()
  with open(os.path.join(args.output_dir, '2B.txt.gz'), 'w') as outf:
    np.savetxt(outf, coass_matrix, newline='\n')

  nrc = NodeRelationComputer(loader, ssc.cancer_pops)
  parents = nrc.compute_relations()
  with open(os.path.join(args.output_dir, '3A.txt'), 'w') as outf:
    for child, parent in enumerate(parents):
      # Root node doesn't have a parent, so this value will be meaningless.
      if child == 0:
        continue
      print(child, parent, sep='\t', file=outf)

  ssmrc = SsmRelationComputer(loader)
  anc_desc = ssmrc.compute_ancestor_desc()
  with open(os.path.join(args.output_dir, '3B.txt.gz'), 'w') as outf:
    np.savetxt(outf, anc_desc, newline='\n')

if __name__ == '__main__':
  main()
