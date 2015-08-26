import argparse
import gzip
import json
import zipfile
from write_report import ResultLoader

def munge_tree_summaries(loader, munged_tree_summ_fn):
  with gzip.GzipFile(munged_tree_summ_fn, 'w') as outf:
    json.dump({
      'dataset_name': loader.dataset_name,
      'trees': loader.tree_summary
    }, outf)

def munge_mutass(loader, munged_mutass_fn):
  with zipfile.ZipFile(munged_mutass_fn, 'w', compression=zipfile.ZIP_DEFLATED) as muts_file:
    for idx, ass in loader.load_all_mut_assignments():
      ass = {
	'dataset_name': loader.dataset_name,
	'mut_assignments': ass,
      }
      muts_file.writestr('%s.json' % idx, json.dumps(ass))

def main():
  parser = argparse.ArgumentParser(
    description='Modify JSON tree summaries to exclude small nodes',
		formatter_class=argparse.ArgumentDefaultsHelpFormatter,
  )
  parser.add_argument('--min-ssms', dest='min_ssms', type=float, default=0.01,
    help='Minimum number of SSMs to retain a subclone')
  parser.add_argument('tree_summary',
    help='JSON-formatted tree summaries')
  parser.add_argument('mutation_list',
    help='JSON-formatted list of mutations')
  parser.add_argument('mutation_assignment',
    help='JSON-formatted list of SSMs and CNVs assigned to each subclone')
  parser.add_argument('munged_tree_summary',
    help='Munged tree summary output filename')
  parser.add_argument('munged_mutation_assignment',
    help='Munged mutation assignment output filename')
  args = parser.parse_args()

  loader = ResultLoader(args.tree_summary, args.mutation_list, args.mutation_assignment, args.min_ssms)
  munge_tree_summaries(loader, args.munged_tree_summary)
  munge_mutass(loader, args.munged_mutation_assignment)

main()
