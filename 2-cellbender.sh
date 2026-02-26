#!/bin/bash
#SBATCH --account=def-naguibm            # Replace with your actual SLURM account
#SBATCH --job-name=lib2_cellbender
#SBATCH -o cellbender_lib2.o        # Standard output log
#SBATCH -e cellbender_lib2.e        # Standard error log
#SBATCH --time=3:00:00                   # Max runtime (3 hours)
#SBATCH --mem=186G                       # Memory allocation
#SBATCH -N 1                             # Number of nodes
#SBATCH --gpus-per-node=1               # Request 1 GPU
#SBATCH --cpus-per-task=6               # Number of CPU threads
 

#Module loads
module load cuda/12.2
module load apptainer/1.3.5
 
 
#Run Cellbender Docker
#Parameters -https://cellbender.readthedocs.io/en/latest/reference/index.html
 
apptainer run -C \
        -B /home/kelperl/scratch/cellbenderpipeline2:/mnt \
        -B /home/kelperl/scratch/cellbenderpipeline2/tmpcellbender:/tmp \
        --env TMPDIR=/tmp \
        --pwd /mnt/cellbender_outs \
        --nv /home/kelperl/scratch/cellbenderpipeline2/cellbender.sif \
        cellbender remove-background --cuda \
        --posterior-batch-size 32 \
        --fpr 0.000005 \
        --learning-rate 0.00005 \
        --checkpoint /mnt/cellbender_outs/ckpt.tar.gz \
        --input /mnt/raw_feature_bc_matrix_lib2.h5 \
        --output /mnt/cellbender_outs/lib2/lib2_cellbender.h5
