#!/bin/bash
#BATCH -J PE_align
#SBATCH -N 1
#SBATCH --array=1-2

## --------------  change to 82 jobs
echo "starting task id $SLURM_ARRAY_TASK_ID"
./STAR_alignment_PE.sh $SLURM_ARRAY_TASK_ID
