#!/bin/bash
set -euo pipefail
BASEDIR=$SCRATCH/smchet-multi
INDIR=$BASEDIR/inputs
OUTBASE=$BASEDIR/results
JOBDIR=$SCRATCH/jobs
PARALLEL=40

function main {
  cd $BASEDIR
  for ssmfn in $BASEDIR/inputs/*.sampled.ssm; do
    runid=$(basename $ssmfn | cut -d. -f1)
    jobname="smchet_pwgs_${runid}"
    OUTDIR=$OUTBASE/$runid

    (
      echo "#!/bin/bash"
      echo "#SBATCH --nodes=1"
      echo "#SBATCH --ntasks=$PARALLEL"
      echo "#SBATCH --time=24:00:00"
      echo "#SBATCH --job-name $jobname"
      echo "#SBATCH --output=$JOBDIR/slurm_${jobname}_%j.txt"
      echo "#SBATCH --mail-type=NONE"

      # Must source ~/.bash_host_specific to get PATH set properly for
      # Miniconda.
      echo "source $HOME/.bash_host_specific && " \
        "module load gsl &&" \
        "mkdir -p $OUTDIR &&" \
        "cd $OUTDIR && " \
        "python2 $HOME/.apps/phylowgs/multievolve.py" \
        "--num-chains $PARALLEL" \
        "--ssms $INDIR/${runid}.sampled.ssm" \
        "--cnvs $INDIR/${runid}.sampled.cnv" \
        "--params $INDIR/${runid}.params.json" \
        ">$runid.stdout" \
        "2>$runid.stderr"
    ) > $JOBDIR/job.sh
    sbatch $JOBDIR/job.sh
    rm $JOBDIR/job.sh
  done
}

main
