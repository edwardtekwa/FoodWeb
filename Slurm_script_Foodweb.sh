#!/bin/bash

#SBATCH --partition=p_mlp195             # Partition (job queue)
#SBATCH --job-name=foodweb         # Assign an short name to your job
#SBATCH --nodes=2                    # Number of nodes you require
#SBATCH --ntasks=2                   # Total # of tasks across all nodes
#SBATCH --cpus-per-task=30            # Cores per task (>1 if multithread tasks)
#SBATCH --mem=120000                   # Real memory (RAM) required (MB)
#SBATCH --time=336:00:00              # Total run time limit (HH:MM:SS)
#SBATCH --output=slurm.%N.%j.out     # STDOUT output file

mkdir -p /scratch/$USER/$SLURM_JOB_ID

module purge
module load MATLAB/R2019a

matlab -nosplash -nodisplay < Make_warming_MultPtsstats_Parallel0.m

rm -rf /scratch/$USER/$SLURM_JOB_ID

sleep 6
sacct --format=MaxRSS,MaxDiskRead,MaxDiskWrite,Elapsed,Nodelist -j $SLURM_JOBID
sleep 2