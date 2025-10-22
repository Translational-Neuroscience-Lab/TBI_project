#!/bin/bash
#SBATCH --account=def-briannek
#SBATCH --job-name=mse_graph_per_mouse
#SBATCH --output=logs/mse_graph_per_mouse_%j.out
#SBATCH --error=logs/mse_graph_per_mouse_%j.err
#SBATCH --time=0-01:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module load matlab/2023b.2
mkdir -p logs

matlab -nosplash -nodesktop -r "addpath('/home/vcarriqu/scratch/MSE'); base_path = '/home/vcarriqu/scratch/MSE/dataset/MSE_NLF_12m_chronic'; conditions = {'TBI', 'Sham'}; plot_mse_comparison_average(base_path, conditions); exit;"
