#!/bin/bash
#SBATCH --time=72:00:00
#SBATCH --account=rrg-gturecki
#SBATCH --job-name=UF_pool7A
#SBATCH --output=/home/kelperl/projects/rrg-gturecki/kelperl/20241202_seqbatch2_libs7a-12b/%j.out
#SBATCH --cpus-per-task=4
#SBATCH --mem=100G
#SBATCH --mail-user=kelly.perlman@mail.mcgill.ca
#SBATCH --mail-type=END,FAIL

/project/rrg-gturecki/Software_Installations_and_Databases/cellranger_installation_8.0.1/cellranger-8.0.1/cellranger count --id 20241205_UF_seqbatch2_lib7A \
         --transcriptome /lustre06/project/6019280/Software_Installations_and_Databases/cellranger_references/refdata-gex-GRCh38-2024-A \
         --create-bam true \
         --fastqs /home/kelperl/projects/rrg-gturecki/kelperl/20241202_seqbatch2_libs7a-12b/rawdata/UF_pool7A \
         --sample UF_pool7A_MEC2489A1-1 \
         --jobmode local \
         --maxjobs 100 \
         --jobinterval 1000 \
         --localcores=8 \
         --localmem=64
