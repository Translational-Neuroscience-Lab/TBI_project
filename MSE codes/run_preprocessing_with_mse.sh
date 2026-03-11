#!/bin/bash
#SBATCH --account=def-briannek
#SBATCH --job-name=preprocess_mouse
#SBATCH --output=logs/preprocess_mouse_%A_%a.out
#SBATCH --error=logs/preprocess_mouse_%A_%a.err
#SBATCH --array=1-18 #change this based on the number of edf files that you have
#SBATCH --time=0-02:30:00
#SBATCH --cpus-per-task=2
#SBATCH --mem=16G

module load matlab/2023b.2

# Find all EDF files
mapfile -t EDF_FILES < <(find /home/vcarriqu/scratch/MSE/dataset/MSE_NLF_6m_subacute -type f -name '*_export.edf' | sort) #change according to your base directory

# Get the specific EDF file for this array task
EDF=${EDF_FILES[$SLURM_ARRAY_TASK_ID-1]}

# Corresponding TSV path
TSV="${EDF/_export.edf/_scores.tsv}"

# Output base = parent of EDF + results folder
COND_FOLDER=$(dirname "$EDF")
OUTPUT_BASE="${COND_FOLDER}/Results"


echo "Running preprocessing for:"
echo "EDF: $EDF"
echo "TSV: $TSV"
echo "Output Base: $OUTPUT_BASE"

# Run MATLAB preprocessing
matlab -nosplash -nodesktop -nodisplay -r "addpath('/home/vcarriqu/scratch/MSE'); preprocessing_MSE_daylight('$EDF', '$TSV', '$OUTPUT_BASE'); exit;"
# find /home/vcarriqu/scratch/MSE/dataset/MSE_NLF_12m_chronic -type f -name "*.edf" | wc -l
