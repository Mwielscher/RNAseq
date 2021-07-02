#!/bin/bash
#BATCH -J PE_align
#SBATCH -N 1
#SBATCH --array=1-4

##  82 jobs


echo "starting task id $SLURM_ARRAY_TASK_ID"

./STAR_alignment_SE.sh $SLURM_ARRAY_TASK_ID

