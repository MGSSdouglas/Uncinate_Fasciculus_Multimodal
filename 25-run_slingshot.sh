#!/bin/bash
#SBATCH --time=150:00:00
#SBATCH --account=def-naguibm
#SBATCH --output=v2_slingshot_%j.out
#SBATCH --error=v2_slingshot_%j.err
#SBATCH --cpus-per-task=12          # Number of CPU cores you want to use
#SBATCH --mem=512G

module load StdEnv/2023
module load r/4.4.0

Rscript 25-slingshot.R
