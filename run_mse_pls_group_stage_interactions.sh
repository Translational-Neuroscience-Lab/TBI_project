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

matlab -nosplash -nodesktop -r "addpath('/home/eshams/scratch/test'); base_path = '/home/eshams/scratch/test/Dataset/NLF/NLF_Older_Cohort'; mse_pls_group_stage_interactions(base_path); exit;"
