from __future__ import print_function
import argparse
import numpy as np
import os
from scipy import stats
from pwgsresults.result_loader import ResultLoader

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
      # clonal_idx should always be 1, given the renumbering I do to remove
      # nonexistent nodes.
      assert clonal_idx == 1
      cellularities.append(pops[clonal_idx]['cellular_prevalence'])

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
      vals = [(pops[pidx]['cellular_prevalence'], pops[pidx]['num_ssms']) for pidx in sorted(pops.keys()) if pidx != 0]
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
    assert np.isclose(np.rint(rerr), rerr), '%s %s' % (np.rint(rerr), rerr)

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

class SsmAssignmentComputer(object):
  def __init__(self, loader):
    self._loader = loader

  def compute_ssm_assignments(self):
    num_ssms = self._loader.num_ssms
    for tree_idx, mut_assignments in self._loader.load_all_mut_assignments():
      num_pops = len(mut_assignments)
      ssm_ass = np.zeros((num_ssms, num_pops))
      for subclone_idx, muts in mut_assignments.items():
        ssm_ids = [int(ssm['id'][1:]) for ssm in muts['ssms']]
        ssm_ass[ssm_ids, subclone_idx - 1] = 1.0
        yield (tree_idx, ssm_ass)

class CoassignmentComputer(object):
  def __init__(self, loader):
    self._loader = loader

  def compute_coassignments(self):
    num_ssms = self._loader.num_ssms
    coass = np.zeros((num_ssms, num_ssms))
    num_trees = 0
    ssm_ass = SsmAssignmentComputer(self._loader)

    for tree_idx, ssm_ass in ssm_ass.compute_ssm_assignments():
      num_trees += 1
      coass += np.dot(ssm_ass, ssm_ass.T)
    coass /= num_trees
    return coass

class SsmRelationComputer(object):
  def __init__(self, loader):
    self._loader = loader

  def _determine_node_ancestry(self, tree_structure, num_pops):
    node_ancestry = np.zeros((num_pops, num_pops))
    def _mark_desc(par, desc):
      for descendant in desc:
        node_ancestry[par - 1, descendant - 1] = 1
        if descendant in tree_structure:
          _mark_desc(par, tree_structure[descendant])
    for parent, children in tree_structure.items():
      _mark_desc(parent, children)
    return node_ancestry

  def compute_ancestor_desc(self):
    ssm_ass = SsmAssignmentComputer(self._loader)
    num_ssms = self._loader.num_ssms
    ancestor_desc = np.zeros((num_ssms, num_ssms))
    num_trees = 0

    for tree_idx, ssm_ass in ssm_ass.compute_ssm_assignments():
      num_trees += 1
      tree_summ = self._loader.tree_summary[tree_idx]
      structure = tree_summ['structure']
      num_pops = ssm_ass.shape[1]
      node_ancestry = self._determine_node_ancestry(structure, num_pops)

      ssm_ancestry = np.dot(ssm_ass, node_ancestry)
      # ADM: ancestor-descendant matrix
      tree_adm = np.dot(ssm_ancestry, ssm_ass.T)
      ancestor_desc += tree_adm

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
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  parser.add_argument('output_dir',
    help='Directory in which to save Challenge outputs')
  args = parser.parse_args()

  loader = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment)
  outputs_to_write = set(('1A', '1B', '1C', '2A', '2B', '3A', '3B'))

  # ssc is used for outputs 1A, 1B, 1C, 2A, and 3A, so always create it, since
  # it will most likely be used by something.
  ssc = SubcloneStatsComputer(loader.tree_summary)
  ssc.calc()

  if '1A' in outputs_to_write:
    with open(os.path.join(args.output_dir, '1A.txt'), 'w') as outf:
      print(ssc.cellularity, file=outf)
  if '1B' in outputs_to_write:
    with open(os.path.join(args.output_dir, '1B.txt'), 'w') as outf:
      print(ssc.cancer_pops, file=outf)
  if '1C' in outputs_to_write:
    with open(os.path.join(args.output_dir, '1C.txt'), 'w') as outf:
      for cluster_num, (num_ssms, phi) in enumerate(zip(ssc.num_ssms, ssc.phis)):
        print(cluster_num + 1, num_ssms, phi, sep='\t', file=outf)

  if '2A' in outputs_to_write:
    cmc = ClusterMembershipComputer(loader, ssc)
    with open(os.path.join(args.output_dir, '2A.txt'), 'w') as outf:
      # These will be in sorted order, so nth row refers to SSM with ID "s{n - 1}".
      for ssm_id, cluster in cmc.calc():
        print(cluster, file=outf)

  if '2B' in outputs_to_write:
    coassc = CoassignmentComputer(loader)
    coass_matrix = coassc.compute_coassignments()
    np.savetxt(os.path.join(args.output_dir, '2B.txt.gz'), coass_matrix)

  if '3A' in outputs_to_write:
    nrc = NodeRelationComputer(loader, ssc.cancer_pops)
    parents = nrc.compute_relations()
    with open(os.path.join(args.output_dir, '3A.txt'), 'w') as outf:
      for child, parent in enumerate(parents):
        # Root node doesn't have a parent, so this value will be meaningless.
        if child == 0:
          continue
        print(child, parent, sep='\t', file=outf)

  if '3B' in outputs_to_write:
    ssmrc = SsmRelationComputer(loader)
    anc_desc = ssmrc.compute_ancestor_desc()
    np.savetxt(os.path.join(args.output_dir, '3B.txt.gz'), anc_desc)

if __name__ == '__main__':
  main()
