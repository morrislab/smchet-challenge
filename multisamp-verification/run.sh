#!/bin/bash
#set -euo pipefail

BASEDIR=$(dirname "$(readlink -f "$0")")
JOBDIR=$SCRATCH/jobs
INDIR=$SCRATCH/pairtree/inputs/smchet-multisamp.dpclust_allclust
RESULTSDIR=$SCRATCH/pairtree/results/smchet-multisamp.dpclust_allclust
PAIRTREEDIR=~/work/pairtree

PARALLEL=80
TREE_CHAINS=40
TREES_PER_CHAIN=10000
PHI_ITERATIONS=100000

function convert_inputs {
  module load r gnu-parallel
  mkdir -p $INDIR

  rm -f $INDIR/*.{ssm,params.json}
  for maximefn in $INDIR/*.Rda; do
    sampid=$(basename $maximefn .Rda)
    echo "Rscript --vanilla $BASEDIR/convert_maxime.R" \
      "$maximefn" \
      "$INDIR/$(echo $sampid | tr . _).ssm.tmp &&" \
      "PYTHONPATH=$PAIRTREEDIR:$PYTHONPATH python3 $BASEDIR/convert_maxime.py" \
      "--remove-small-clusters" \
      "$INDIR/$(echo $sampid | tr . _).{ssm.tmp,ssm,params.json} &&" \
      "rm $INDIR/$(echo $sampid | tr . _).ssm.tmp"
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
        "OMP_NUM_THREADS=1 OMP_THREAD_LIMIT=1 MKL_NUM_THREADS=1 python3 $PAIRTREEDIR/basic.py" \
        "--seed 1" \
        "--parallel $PARALLEL" \
        "--tree-chains $TREE_CHAINS" \
        "--trees-per-chain $TREES_PER_CHAIN" \
        "--phi-iterations $PHI_ITERATIONS" \
        "--params $INDIR/$runid.params.json" \
        "$INDIR/$runid.ssm" \
        "$RESULTSDIR/$runid.results.npz" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) > $jobfn
    sbatch $jobfn
    rm $jobfn
  done
}

function plot {
  module load gnu-parallel

  for resultsfn in $RESULTSDIR/*.results.npz; do
    runid=$(basename $resultsfn | cut -d. -f1)
    (
      echo "cd $RESULTSDIR && " \
        "python3 $PAIRTREEDIR/plotter.py" \
        "$runid" \
        "$INDIR/$runid.ssm" \
        "$INDIR/${runid}.params.json" \
        "$RESULTSDIR/$runid.results.npz" \
        "$RESULTSDIR/$runid.results.html" \
        ">>$runid.stdout" \
        "2>>$runid.stderr"
    )
  done | parallel -j$PARALLEL

  cd $RESULTSDIR
  (
    echo "<ul>"
    for foo in $(ls -v *.results.html); do
      echo "<li><a href=$foo>$(echo $foo | cut -d. -f1)</a> ($(ls -lh $foo | awk '{print $5}'))</li>"
    done
    echo "</ul>"
  ) > index.html
}

function main {
  #convert_inputs
  run_pairtree
  #plot
}

main
