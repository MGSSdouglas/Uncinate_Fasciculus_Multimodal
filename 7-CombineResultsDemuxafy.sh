#!/bin/bash
#SBATCH --time=1:00:00
#SBATCH --account=def-naguibm
#SBATCH --mem=100G
#SBATCH --output=%x-%j.out
#SBATCH --job-name=lib1_combine


module load StdEnv/2023
module load apptainer/1.3.5


OUTDIR=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/20250714_combinedresults_lib1
DEMUXLET_OUTDIR=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/20250707_demuxlet_lib1
SCDBLFINDER_OUTDIR=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/20250707_scDblFinder_lib1

apptainer exec --bind /home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/ /project/rrg-gturecki/kelperl/Demuxafy.sif Combine_Results.R \
  -o $OUTDIR/combined_results.tsv \
  --scDblFinder $SCDBLFINDER_OUTDIR \
  --demuxlet $DEMUXLET_OUTDIR \
  --method "AnyDoublet" ## there are other methods that can also be used, please see the help message above for the other options
