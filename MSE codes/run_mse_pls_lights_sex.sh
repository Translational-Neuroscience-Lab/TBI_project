#!/bin/bash
#SBATCH --account=def-briannek
#SBATCH --job-name=plot_mse
#SBATCH --output=logs/mse_pls_%j.out
#SBATCH --error=logs/mse_pls_%j.err
#SBATCH --time=0-00:50:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G

module load matlab/2023b.2
mkdir -p logs


matlab -nosplash -nodesktop -r "addpath('/home/vcarriqu/scratch/MSE'); base_path_Adult = '/home/vcarriqu/scratch/MSE'; mse_pls_lights_sex(base_path_Adult); exit;"
