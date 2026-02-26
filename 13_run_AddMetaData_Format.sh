#!/bin/bash
#SBATCH --time=12:00:00
#SBATCH --account=def-naguibm
#SBATCH --output=conversions_seurat_%j.out
#SBATCH --error=conversions_seurat_%j.err
#SBATCH --cpus-per-task=4          # Number of CPU cores you want to use
#SBATCH --mem=256G                 # Adjust memory allocation

module load StdEnv/2023
module load r/4.4.0

module load python
python -m venv ~/envs/anndata_env
source ~/envs/anndata_env/bin/activate
pip install --upgrade pip
pip install anndata scanpy


Rscript 13-AddMetaData_Format.R
