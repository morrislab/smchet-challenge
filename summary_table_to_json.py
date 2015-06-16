import csv
import json
import sys

def main():
  table_fn = sys.argv[1]
  stats = {}
  with open(table_fn) as tablefd:
    reader = csv.DictReader(tablefd, delimiter='\t')
    for row in reader:
      stats[row['samplename']] = {
        'num_subclones': int(row['num_subclones']),
        'purity': float(row['purity']),
        'frac_clonal': float(row['frac_clonal']),
      }

  print(json.dumps(stats))

main()
