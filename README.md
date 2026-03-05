# H1H2_neuron_RNAseq_Batch1_DA
Pipeline for differential RNA-seq analysis comparing H1H1 vs. H2H2 iPSC-derived neurons (Batch 1)

## Pipeline Stages:
1. Prepocessing of the gene expression count matrix
   1) gene filtering by CPM and %samples
   2) CQN normalization
   code: 3_GeneExpr_Filtering_CQN_normalization.R
2. Build linear model
   1) correlation of biological and technical covariates with top PCs of the normalized gene expression matrix
   2) all VIF < 2.5
   code: 4_PCA_Build_linear_model.R
3. Differential analysis
   1) DESeq2
   code: 5_DEseq2.R

## File locations
data_provided/ contains the raw gene count matrix, gene effective length, picard metrics, and metadata.
code/ contains all the scripts
results/ contains some of the results from running the scripts
