#!/bin/bash
#SBATCH --time=16:00:00
#SBATCH --account=def-naguibm
#SBATCH --output=%x-%j.out
#SBATCH --mem=512G
#SBATCH --cpus-per-task=8
#SBATCH --job-name=filtBAM_lib1

module load samtools/1.20
module load StdEnv/2023


#reformat barcodes
awk '{print "CB:Z:" $1}' /home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/lib1_cellbender_cell_barcodes.csv > filter.txt


#rename BAM file
export BAM_FILE='/home/kelperl/scratch/seqbatch1_lib1/20250227_cellranger_pool1/20250226_UF_seqbatch1_lib1/outs/possorted_genome_bam.bam'

# Save the header lines
samtools view -H $BAM_FILE > SAM_header

# Filter alignments using filter.txt. Use LC_ALL=C to set C locale instead of UTF-8
samtools view $BAM_FILE | LC_ALL=C grep -F -f filter.txt > filtered_SAM_body

# Combine header and body
cat SAM_header filtered_SAM_body > filtered.sam

# Convert filtered.sam to BAM format
samtools view -b filtered.sam > filtered.bam
