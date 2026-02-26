# Uncinate fasciculus Single-nucleus RNA sequencing (snRNA-seq) bioinformatics pipeline

Associated with publication: A multimodal characterization of the human uncinate fasciculus

Kelly Perlman, Sarah-Barnett Burns, John Kim, Valérie Pineau Noël, Armand Collin, Justine Major, Malosree Maitra, Anjali Chawla, Murielle Mardenli, Sébastien Jerczynski, Maria Antonietta Davoli, Gabriella Frosi, Julien Cohen-Adad, Daniel Côté, Richard Bazinet, Gustavo Turecki, Corina Nagy, Naguib Mechawar

An overview of the workflow can be found below:

Note: steps 1-7 run for each library 

1 - cellranger

2 - cellbender

3 - filter cellranger bam file by cellbender barcodes

4 - demuxlet
https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/Demuxlet.html

5 - convert h5 to 10x readable format for scdbl finder

6 - scDblfinder
https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/scDblFinder.html

7 - combine demuxafy results
https://demultiplexing-doublet-detecting-docs.readthedocs.io/en/latest/CombineResults.html

8 - assign cells, make seurat objects and merge seurat  objects

9 - normalize with scTransform

10 - batch correct with harmony integration

11 - clustering (join layers, find neighbours, find clusters, run umap, markers, etc..)

12 - find markers sct wilcoxin

13 - add subject info/coviates. convert to other subject formats

(Do map my cells online: https://knowledge.brain-map.org/mapmycells/process/)
14 - meta neighbour analyses and scclusteval 

15 - cluster stability analysis adjusted rand index

16 - QC pseudobulk (PCA)

17 - variance Partitioning analysis

18a - limma-voom DEG group - individual

18b - limma-voom DEG group - broad

19a - downstream group DEG analysis  - disgenet, GSEA, etc..

19b - violin plotting top genes

20a - limma-voom DEG age - individual

20b - limma-voom DEG age- broad

21a - Age DEG downstream analysis

22a - Cell proportions -individual

22b - Cell type proportions - broad 

22c - Age animal-human overlap analysis (with data from Ximerakis et al., 2019, Nat Neuro)

23 - OL-lineage cells only scTransform

24 - OL-lineage cells only harmony

25 - slingshot pseudotime analysis

26 - WGCNA on OL cluster


Please note: LLM AI tools were used to assist in the writing and editing of a subset of scripts.

