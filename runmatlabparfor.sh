#!/bin/bash

#SBATCH --partition=p_mlp195             # Partition (job queue)
#SBATCH --job-name=mparfor           # Assign an short name to your job
#SBATCH --nodes=1                    # Number of nodes you require
#SBATCH --ntasks=1                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=8            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=8000                   # Real memory (RAM) required (MB)
#SBATCH --time=00:05:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file

mkdir -p /scratch/$USER/$SLURM_JOB_ID

module purge
module load MATLAB/R2019a

matlab -nosplash -nodisplay < pmc.m

rm -rf /scratch/$USER/$SLURM_JOB_ID

sleep 6
sacct --format=MaxRSS,MaxDiskRead,MaxDiskWrite,Elapsed,Nodelist -j $SLURM_JOBID
sleep 2

