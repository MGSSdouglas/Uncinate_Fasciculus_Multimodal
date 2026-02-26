#!/bin/bash
#SBATCH --time=6:00:00
#SBATCH --account=def-naguibm
#SBATCH --output=sctransform_ol_seurat_%j.out
#SBATCH --error=sctransform_ol_seurat_%j.err
#SBATCH --cpus-per-task=4          # Number of CPU cores you want to use
#SBATCH --mem=256G                 # Adjust memory allocation

module load StdEnv/2023
module load r/4.4.0

# Control threading for memory efficiency
export MKL_THREADING_LAYER=GNU
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_DISPLAY_ENV=TRUE
export OMP_THREAD_LIMIT=1


Rscript 23-scTransform_OLlineage.R
