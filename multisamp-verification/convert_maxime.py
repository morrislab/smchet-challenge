import argparse
import csv
import numpy as np
from collections import defaultdict
import clustermaker
import inputparser
import json

def _extract_nums(S, dtype):
  return np.array(S.split(','), dtype=dtype)

def load_ssms(ssmfn, max_ssms=None):
  variants = {}

  with open(ssmfn) as F:
    reader = csv.DictReader(F, delimiter='\t')
    for row in reader:
      if max_ssms is not None and len(variants) >= max_ssms:
        break
      variant = {
        'id': row['id'],
        'name': row['name'],
        'var_reads': _extract_nums(row['var_reads'], np.int),
        'total_reads': _extract_nums(row['total_reads'], np.int),
        'omega_v': _extract_nums(row['var_read_prob'], np.float),
        'cluster': int(row['cluster']),
      }

      assert np.all(variant['total_reads'] >= variant['var_reads'])
      if np.any(np.isnan(variant['omega_v'])) or np.any(np.isinf(variant['omega_v'])):
        continue
      assert np.all(0 <= variant['omega_v']) and np.all(variant['omega_v'] <= 1)
      variant['omega_v'] = np.maximum(variant['omega_v'], 1e-5)

      variant['id'] = 's%s' % len(variants)
      assert variant['id'] not in variants
      variants[variant['id']] = variant

  S = len(next(iter(variants.values()))['var_reads'])
  for vid, V in variants.items():
    for K in ('var_reads', 'total_reads', 'omega_v'):
      assert len(V[K]) == S, '%s for %s is of length %s, but %s expected' % (K, vid, len(V[K]), S)

  return variants

def extract_clusters(variants):
  clusters = defaultdict(list)
  for vid, V in variants.items():
    clusters[V['cluster']].append(vid)
  # Don't permit empty clusters.
  clusters = [clusters[C] for C in sorted(clusters.keys())]
  return clusters

def sort_clusters_by_vaf(clusters, variants):
  supervars = clustermaker.make_cluster_supervars(clusters, variants)
  supervars = [supervars['S%s' % idx] for idx in range(len(supervars))]
  sv_vaf = np.array([S['vaf'] for S in supervars])
  mean_vaf = np.mean(sv_vaf, axis=1)
  order = np.argsort(-mean_vaf)
  return [clusters[idx] for idx in order]

def make_sampnames(variants):
  S = len(next(iter(variants.values()))['var_reads'])
  return ['Sample %s' % (sidx + 1) for sidx in range(S)]

def write_params(clusters, garbage, sampnames, paramsfn):
  params = {
    'samples': sampnames,
    'clusters': clusters,
    'garbage': garbage,
  }
  with open(paramsfn, 'w') as outf:
    json.dump(params, outf)

def remove_small_clusters(clusters):
  garbage = []
  filtered = []

  num_vars = np.sum([len(C) for C in clusters])
  threshold = np.round(0.01 * num_vars)
  for C in clusters:
    if len(C) >= threshold:
      filtered.append(C)
    else:
      garbage += C

  return (filtered, garbage)

def main():
  parser = argparse.ArgumentParser(
    description='LOL HI THERE',
    formatter_class=argparse.ArgumentDefaultsHelpFormatter
  )
  parser.add_argument('--remove-small-clusters', action='store_true')
  parser.add_argument('inssm_fn')
  parser.add_argument('outssm_fn')
  parser.add_argument('outparams_fn')
  args = parser.parse_args()

  variants = load_ssms(args.inssm_fn)
  clusters = extract_clusters(variants)
  clusters = sort_clusters_by_vaf(clusters, variants)
  inputparser.write_ssms(variants, args.outssm_fn)

  if args.remove_small_clusters:
    clusters, garbage = remove_small_clusters(clusters)
  else:
    garbage = []
  sampnames = make_sampnames(variants)
  write_params(clusters, garbage, sampnames, args.outparams_fn)

if __name__ == '__main__':
  main()
