#!/bin/bash
#
# Updated for slurm. Also merge as we go, so we can run smaller arrays of
# longer jobs
#
#SBATCH --partition=short
#SBATCH --time=3:00:00
#SBATCH --mem=50M

N=${1:-300}
B=${2:-50}
scratch=/scratch/rblythe3
temp=$scratch/cdistblock.${SLURM_JOB_ID}.XXXX

mkdir -p $scratch

prev=$(mktemp $temp)
this=$(mktemp $temp)
for b in $(seq $B); do
  echo Running sample $b of $B at $(date)
  ./crumble -L1000 -N$N -t0.001 -b0 -u1200000 -e50 -a10 -i0.2 clusters > $this
  next=$(mktemp $temp)
  ./merge $prev $this > $next
  rm $prev
  prev=$next
done

rm $this
mv $prev data/L1000-N$N-t0.001-coft-run$SLURM_ARRAY_TASK_ID.dat
echo Completed ${SLURM_JOB_ID}.${SLURM_ARRAY_TASK_ID} at $(date)
