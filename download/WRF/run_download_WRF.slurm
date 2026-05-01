#!/bin/bash
#SBATCH --account=cwp140
#SBATCH --partition=shared-128
#SBATCH --nodes=1
#SBATCH --mem=32GB
#SBATCH --time=04:00:00
#SBATCH --job-name=wrf_download
#SBATCH --output=logs/out_%A_%a.txt
#SBATCH --error=logs/err_%A_%a.txt
#SBATCH --array=0-999%100

OFFSET=$1

/home/dnash/miniconda3/envs/SEAK-impacts/bin/python download_wrf.py \
    --task-id ${SLURM_ARRAY_TASK_ID} \
    --offset ${OFFSET} \
    --nproc 4

# sbatch run_download_WRF.slurm 0
sbatch run_download_WRF.slurm 1000