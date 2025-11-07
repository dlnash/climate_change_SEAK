#!/bin/bash
#SBATCH --account=cwp140
#SBATCH --partition=shared-128
#SBATCH -N 1 # Ensure that all cores are on one machine
#SBATCH --nodes=1 # specifies how many nodes to use (each node has 24 cores)
#SBATCH --mem=32GB
#SBATCH -t 5:00:00 # Runtime in HH:MM:SS
#SBATCH -o slurm_out/output_%a # File to which STDOUT will be written
#SBATCH -e slurm_out/output_err%a # File to which STDERR will be written
#SBATCH --job-name=compute_trends
#SBATCH --mail-type=ALL # Type of email notification- BEGIN,END,FAIL,ALL
#SBATCH --mail-user=dnash@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=1-15%15 # runs jobs x-y, 100 at a time

# Activate Anaconda work environment
source /home/dnash/miniconda3/envs/SEAK-impacts/etc/conda/activate.d/*
 
linevar=`sed $SLURM_ARRAY_TASK_ID'q;d' calls_snow.txt`
eval $linevar
