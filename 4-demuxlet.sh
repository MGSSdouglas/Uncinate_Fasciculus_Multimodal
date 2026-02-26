#!/bin/bash
#SBATCH --time=96:00:00
#SBATCH --account=def-naguibm
#SBATCH --cpus-per-task=8
#SBATCH --mem=256G
#SBATCH --output=%x-%j.out
#SBATCH --mail-user=kelly.perlman@mail.mcgill.ca
#SBATCH --mail-type=END,FAIL
#SBATCH --job-name=lib1_demuxlet


module load StdEnv/2023
module load apptainer/1.3.5


VCF=/home/kelperl/scratch/genotype/mgss_genotypes_round_1_2_3_22_4_5_6_autosomes_R2_0.1_hwe_0.000001_MAF_0.05_no_CNV_in_genes_only.recode.vcf
BARCODES=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/lib1_cellbender_cell_barcodes.csv
BAM=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/filtered.bam
DEMUXLET_OUTDIR=/home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/20250707_demuxlet_lib1
FIELD='GT' 
INDS=/home/kelperl/scratch/seqbatch1_lib1/UF_lib1_ids.txt #This will use all the individuals in your reference SNP genotype $VCF. If your $VCF only has the individuals multiplexed in your pool, then the $INDS file is not required.


#popscle pileup

apptainer exec --bind /home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/ --bind /home/kelperl/scratch/genotype/ /project/rrg-gturecki/kelperl/Demuxafy.sif popscle_pileup.py \
        --sam $BAM \
        --vcf $VCF \
        --group-list $BARCODES \
        --out $DEMUXLET_OUTDIR/pileup \
        --sm-list $INDS



#popscle demuxlet

apptainer exec --bind /home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/ --bind /home/kelperl/scratch/genotype/ /project/rrg-gturecki/kelperl/Demuxafy.sif popscle demuxlet \
        --plp $DEMUXLET_OUTDIR/pileup \
        --vcf $VCF \
        --field $FIELD \
        --group-list $BARCODES \
        --geno-error-coeff 1.0 \
        --geno-error-offset 0.05 \
        --out $DEMUXLET_OUTDIR/demuxlet \
        --sm-list $INDS




#demuxlet summary
apptainer exec --bind /home/kelperl/scratch/cellbenderpipeline2/cellbender_outs/lib1/ --bind /home/kelperl/scratch/genotype/ /project/rrg-gturecki/kelperl/Demuxafy.sif bash Demuxlet_summary.sh $DEMUXLET_OUTDIR/demuxlet.best > $DEMUXLET_OUTDIR/demuxlet_summary.tsv

