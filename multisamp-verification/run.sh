#!/bin/bash
#set -euo pipefail

BASEDIR=$(dirname "$(readlink -f "$0")")
INDIR=$SCRATCH/pairtree/inputs/smchet-multisamp
RESULTSDIR=$SCRATCH/pairtree/results/smchet-multisamp
PAIRTREEDIR=~/work/pairtree

PARALLEL=80
TREE_CHAINS=40
TREES_PER_CHAIN=1000
PHI_ITERATIONS=10000

function convert_inputs {
  module load r gnu-parallel
  mkdir -p $INDIR

  for maximefn in $INDIR/*.Rda; do
    sampid=$(basename $maximefn .Rda)
    echo "Rscript --vanilla $BASEDIR/convert_maxime.R" \
      "$maximefn" \
      "$INDIR/$(echo $sampid | tr . _).ssm"
  done | parallel -j$PARALLEL --halt 1
}

function run_pairtree {
  mkdir -p $RESULTSDIR
  #rm -f $RESULTSDIR/SJ*.{html,json,csv,stdout}

  for ssmfn in $INDIR/*.ssm; do
    runid=$(basename $ssmfn .ssm)
    jobname="smchet_pairtree_${runid}"
    jobfn=$(mktemp)
    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=23:59:00"
      echo "#SBATCH --job-name $jobname"
      echo "#SBATCH --output=$JOBDIR/slurm_${jobname}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      echo "cd $RESULTSDIR && " \
        "python3 $PAIRTREEDIR/basic.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $TREE_CHAINS" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "$INDIR/$runid.ssm" \
        "$RESULTSDIR/$runid.results.npz" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function main {
    #convert_inputs
    run_pairtree
}

main
