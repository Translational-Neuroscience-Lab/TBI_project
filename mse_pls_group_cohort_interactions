#!/bin/bash
#SBATCH --account=def-briannek
#SBATCH --job-name=pls_group_session_mse
#SBATCH --output=logs/pls_within_session_mse_%j.out
#SBATCH --error=logs/pls_within_session_mse_%j.err
#SBATCH --time=0-05:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module load matlab/2023b.2
mkdir -p logs

matlab -nosplash -nodesktop -r "addpath('/home/vcarriqu/scratch/MSE'); base_path = '/home/vcarriqu/scratch/MSE/dataset'; mse_pls_group_cohort_interactions(base_path); exit;"
