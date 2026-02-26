#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-naguibm
#SBATCH --output=harmony_seurat_%j.out
#SBATCH --error=harmony_seurat_%j.err
#SBATCH --cpus-per-task=64
#SBATCH --mem=256G  

module load StdEnv/2023
module load r/4.4.0

export MKL_THREADING_LAYER=GNU
export MKL_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export OMP_DISPLAY_ENV=TRUE
export OMP_THREAD_LIMIT=1
export BLIS_NUM_THREADS=1


Rscript 10-integrate_harmony.R
