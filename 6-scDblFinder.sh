#!/bin/bash
#SBATCH --time=2:00:00
#SBATCH --account=def-naguibm
#SBATCH --job-name=scDblFinder_lib1
#SBATCH --mem=20G
#SBATCH --output=%x-%j.out
#SBATCH --cpus-per-task=1


COUNTS=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/converted_cellranger_format
SCDBLFINDER_OUTDIR=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/20250707_scDblFinder_lib1



#load singularity

module load StdEnv/2023
module load apptainer/1.3.5
module load hdf5/1.14.2
module load r/4.4.0

apptainer exec --bind /home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/ /project/rrg-gturecki/kelperl/Demuxafy.sif scDblFinder.R -o $SCDBLFINDER_OUTDIR -t $COUNTS
