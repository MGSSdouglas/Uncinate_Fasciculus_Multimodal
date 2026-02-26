#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --account=def-naguibm
#SBATCH --output=sctmarkers_seurat_%j.out
#SBATCH --error=sctmarkers_seurat_%j.err
#SBATCH --cpus-per-task=4          # Number of CPU cores you want to use
#SBATCH --mem=256G                 # Adjust memory allocation

module load StdEnv/2023
module load r/4.4.0


Rscript 12-findmarkers.R
