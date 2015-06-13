#!/usr/bin/env python2
import argparse
import csv
import gzip
import json
import numpy as np
from collections import defaultdict
from scipy import stats

def intmode(iter):
  return int(stats.mode(iter)[0][0])

def filter_subclones(subclones):
  min_ssms = 0
  return {int(idx): subclone for (idx, subclone) in subclones.items() if int(idx) != 0 and subclone['num_ssms'] >= min_ssms}

def parse_json(json_fn):
  with gzip.GzipFile(json_fn) as json_fd:
    return json.load(json_fd)

def calc_subclone_stats(tree_summary):
  subclone_counts = []
  cellularities = []

  for tree in tree_summary['trees']:
    subclones = tree['subclones']
    subclones = filter_subclones(subclones)
    subclone_counts.append(len(subclones))
    cellularities.append(subclones[1]['phi'])

  return (intmode(subclone_counts), np.mean(cellularities))

def calc_frac_clonal(mut_assignments):
  clonal_counts = defaultdict(int)
  subclonal_counts = defaultdict(int)
  ssm_ids = []

  for tree in mut_assignments['trees']:
    for subclone_idx, features in tree.items():
      is_clonal = (subclone_idx == '1')
      for ssm in features['ssms']:
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

def main():
  parser = argparse.ArgumentParser(description='Hooray')
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

  tree_summary = parse_json(args.tree_summary)
  mut_assignments = parse_json(args.mutation_assignment)

  dataset_name = tree_summary['dataset_name']
  cancer_type = get_cancer_type(dataset_name, args.pancancer_tsv_details)
  num_subclones, cellularity = calc_subclone_stats(tree_summary)
  ploidy = 'NA'
  frac_clonal = calc_frac_clonal(mut_assignments)
  branching = 'NA'

  print('\t'.join([str(v) for v in (cancer_type, dataset_name, num_subclones, cellularity, ploidy, frac_clonal, branching)]))

if __name__ == '__main__':
  main()
