#!/usr/bin/env python2
import argparse
import gzip
import json
import numpy as np
import os
import sys
import util2
import zipfile
from collections import defaultdict

def summarize_subclones(tree):
  subclones = {}
  structure = defaultdict(list)
  mut_assignments = defaultdict(lambda: {'cnvs': [], 'ssms': []})
  idx = [0]

  def _traverse_r(vertex, parent):
    mutations = vertex.get_data()
    # vertex.params represents phis (i.e., population freqs) associated with
    # each sample.
    phi = np.mean(vertex.params)
    current_idx = idx[0]

    num_ssms = 0
    num_cnvs = 0
    for mut in mutations:
      if mut.id.startswith('s'):
        mut_assignments[current_idx]['ssms'].append({'id': mut.id, 'name': mut.name})
        num_ssms += 1
      elif mut.id.startswith('c'):
        mut_assignments[current_idx]['cnvs'].append({'id': mut.id})
        num_cnvs += 1
      else:
        raise Exception('Unknown mutation ID type: %s' % mut.id)

    # Preorder traversal is consistent with printo_latex.py, meaning index
    # values should correspond to same vertices.
    subclones[current_idx] = {
      'phi': phi,
      'num_ssms': num_ssms,
      'num_cnvs': num_cnvs,
    }

    for child in vertex.children():
      idx[0] += 1
      structure[current_idx].append(idx[0])
      _traverse_r(child, current_idx)

  _traverse_r(tree.root['node'], None)
  return (subclones, mut_assignments, structure)

def parse_json(fin):
  with open(fin) as fh:
    return json.load(fh)

def write_json(dataset_name, tree_file, num_trees, summaries_output, mutation_assignment_output):
  summaries = {
    'dataset_name': dataset_name,
    'trees': {},
  }

  reader = util2.TreeReader(tree_file)
  with zipfile.ZipFile(mutation_assignment_output, 'w', compression=zipfile.ZIP_DEFLATED) as muts_file:
    for idx, llh, tree in reader.load_trees_and_metadata(num_trees = num_trees, remove_empty_vertices = True):
      subclones, muts, structure = summarize_subclones(tree)

      summaries['trees'][idx] = {
        'llh': llh,
        'structure': structure,
        'subclones': subclones,
      }

      muts['dataset_name'] = dataset_name
      muts_file.writestr('%s.json' % idx, json.dumps(muts))
  reader.close()

  with gzip.GzipFile(summaries_output, 'w') as summf:
    json.dump(summaries, summf)

def main():
  parser = argparse.ArgumentParser(description='Write JSON files describing trees')
  parser.add_argument('--num-trees', '-n', dest='num_trees', type=int,
    help='Only examine given number of trees')
  parser.add_argument('dataset_name',
    help='Name identifying dataset')
  parser.add_argument('tree_file',
    help='File containing sampled trees')
  parser.add_argument('tree_summary_output',
    help='Output file for JSON-formatted tree summaries')
  parser.add_argument('mutation_assignment_output',
    help='Output file for JSON-formatted list of SSMs and CNVs assigned to each subclone')
  args = parser.parse_args()

  write_json(args.dataset_name, args.tree_file, args.num_trees, args.tree_summary_output, args.mutation_assignment_output)

if __name__ == '__main__':
  main()
