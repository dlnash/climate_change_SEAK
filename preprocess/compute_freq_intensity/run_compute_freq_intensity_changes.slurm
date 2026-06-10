#!/bin/bash
#SBATCH --account=cwp140
#SBATCH --partition=shared-128
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=compute_freq_intensity
#SBATCH --output=logs/out_%A_%a.txt
#SBATCH --error=logs/err_%A_%a.txt
#SBATCH --mail-user=dnash@ucsd.edu # Email to which notifications will be sent
#SBATCH --array=0-24%10

OFFSET=$1

/home/dnash/miniconda3/envs/SEAK-impacts/bin/python compute_frequency_intensity_changes.py \
    --task-id ${SLURM_ARRAY_TASK_ID} \
    --offset ${OFFSET} \

# python compute_frequency_intensity_changes.py --task-id 0 ## run this to test one job
# python compute_frequency_intensity_changes.py --task-id 0 --print-jobs ## run this to print list of jobs
sbatch run_compute_freq_intensity_changes.slurm 0 # to run all jobs
