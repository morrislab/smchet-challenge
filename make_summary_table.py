#!/usr/bin/env python2
import argparse
import csv
import gzip
import numpy as np
import ujson as json
import sys
import zipfile
from collections import defaultdict
from scipy import stats

def intmode(iter):
  return int(stats.mode(iter)[0][0])

def filter_subclones(subclones, min_ssms):
  return {
    int(idx): subclone for (idx, subclone) in subclones.items()
    if int(idx) != 0 and subclone['num_ssms'] >= min_ssms
  }

def find_clonal_node(subclones):
  indices = [k for k in subclones.keys() if k > 0]
  return min(indices)

def calc_subclone_stats(tree_summary, min_ssms):
  cancer_pop_counts = []
  cellularities = []

  for tree_idx, tree in tree_summary.items():
    try:
      subclones = tree['subclones']
      subclones = filter_subclones(subclones, min_ssms)
      cancer_pop_counts.append(len(subclones))
      clonal_idx = find_clonal_node(subclones)
      cellularities.append(subclones[clonal_idx]['phi'])
    except:
      from IPython import embed
      embed()
  return (intmode(cancer_pop_counts) - 1, np.mean(cellularities))

def calc_frac_clonal(mut_assignments, tree_summary, min_ssms):
  clonal_counts = defaultdict(int)
  subclonal_counts = defaultdict(int)
  ssm_ids = []

  with zipfile.ZipFile(mut_assignments) as mutf:
    for zinfo in mutf.infolist():
      tree_idx = zinfo.filename.split('.')[0]
      subclones = filter_subclones(tree_summary['trees'][tree_idx]['subclones'], min_ssms)
      clonal_idx = find_clonal_node(subclones)

      tree_assignments_json = mutf.read(zinfo)
      tree_assignments = json.loads(tree_assignments_json)

      for subclone_idx, muts in tree_assignments.items():
        subclone_idx = int(subclone_idx)
        is_clonal = (subclone_idx == 1)

        for ssm in muts['ssms']:
          ssm_id = ssm['id']
          ssm_ids.append(ssm_id)
          if is_clonal:
            clonal_counts[ssm_id] += 1
          else:
            subclonal_counts[ssm_id] += 1

  num_clonal_ssms = 0
  num_subclonal_ssms = 0

  for ssm_id in ssm_ids:
    if clonal_counts[ssm_id] >= subclonal_counts[ssm_id]:
      num_clonal_ssms += 1
    else:
      num_subclonal_ssms += 1

  return float(num_clonal_ssms) / (num_clonal_ssms + num_subclonal_ssms)

def get_cancer_type(dataset_name, pancancer_tsv_details):
  with open(pancancer_tsv_details) as tsv_details:
    reader = csv.DictReader(tsv_details, delimiter='\t')
    for row in reader:
      if row['Tumour WGS aliquot ID(s)'] == dataset_name:
        cancer_type = row['Project code'].split('-')[0]
        return cancer_type
  raise Exception('Could not find cancer type for dataset %s' % dataset_name)

def count_ssms(tree_summary):
  total_ssms = sum([s['num_ssms'] for s in tree_summary['trees']['0']['subclones'].values()])
  return total_ssms

def main():
  parser = argparse.ArgumentParser(description='Hooray')
  parser.add_argument('--min-ssms', dest='min_ssms', type=int, default=3,
    help='Minimum number of SSMs to retain a subclone')
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  # For this file, see
  # https://wiki.oicr.on.ca/display/PANCANCER/Available+PCAWG+Data#AvailablePCAWGData-SantaCruzPilotDataset
  # -- specifically, http://pancancer.info/santa_cruz_pilot/santa_cruz_pilot.v2.2015_0504.tsv.
  parser.add_argument('pancancer_tsv_details',
    help='TSV file listing sample details')
  args = parser.parse_args()

  with gzip.GzipFile(args.tree_summary) as json_fd:
    tree_summary = json.load(json_fd)

  total_ssms = count_ssms(tree_summary)
  if total_ssms < 10:
    print >> sys.stderr, 'Only %s SSMs. Exiting.' % total_ssms
    return

  dataset_name = tree_summary['dataset_name']
  cancer_type = get_cancer_type(dataset_name, args.pancancer_tsv_details)
  num_subclones, cellularity = calc_subclone_stats(tree_summary['trees'], args.min_ssms)
  ploidy = 'NA'
  frac_clonal = calc_frac_clonal(args.mutation_assignment, tree_summary, args.min_ssms)
  branching = 'NA'

  print('\t'.join([str(v) for v in (cancer_type, dataset_name, num_subclones, cellularity, ploidy, frac_clonal, branching)]))

if __name__ == '__main__':
  main()
