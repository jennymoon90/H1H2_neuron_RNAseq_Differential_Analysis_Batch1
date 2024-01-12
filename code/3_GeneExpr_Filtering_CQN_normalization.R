rm(list = ls())
library(stringr)
library(tidyverse)
library(limma); library(edgeR); library(cqn); library(biomaRt)
library(ggplot2); library(gridExtra)

setwd("./results/")

##### Pipeline #####

### (0) Load RawData
### (1) Gene Filtering
### (2) Normalization
### (3) Outlier Removal

##### (0) Load RawData #####

rm(list = ls())
lnames = load("../data_provided/3_0_RawData.rda") # rsem_gene, rsem_gene_effLen, datMeta, datSeq_numeric

##### (1) Gene Filtering #####

## Only keep transcripts on chr1-Y
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gene_gc_content")

ensembl_db = listDatasets(useEnsembl(biomart="ENSEMBL_MART_ENSEMBL")) # search for "hsapiens_gene_ensembl" got GRCh38
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

geneAnno1 <- getBM(attributes = getinfo, filters = c("ensembl_gene_id"), values = substr(rownames(rsem_gene),1,15), mart = mart)
geneAnno2 = geneAnno1[geneAnno1$chromosome_name %in% c(1:22, "X", "Y"),]
idx = which(substr(rownames(rsem_gene),1,15) %in% geneAnno2$ensembl_gene_id)
rsem_gene = rsem_gene[idx,]
rsem_gene_effLen = rsem_gene_effLen[idx]

save(rsem_gene, rsem_gene_effLen, datMeta, datSeq_numeric, geneAnno1, file = "3_01_RawData_chr1toY.rda")

### Decide on cpm threshold and x% of samples for keeping genes
# Google: As a general rule, a good threshold can be chosen for a CPM value that corresponds to a count of 10.
rm(list = ls())
lnames = load("3_01_RawData_chr1toY.rda") # contains geneAnno1 as well
#range(datSeq_numeric$READ_PAIRS_EXAMINED) # 7.2-11 M

cpm <- apply(rsem_gene,2, function(x) (x/sum(x))*1000000)

df_idx = as.data.frame(matrix(ncol = 3, nrow = 0))
for (cpm_th in seq(0, 1, 0.1)) {
  for (x in seq(0, 0.5, 0.05)) {
    keep = apply(cpm > cpm_th, 1, sum)
    idx = which(keep > x * dim(cpm)[2]) ## cpm > cpm_th in x% of samples
    #cpm = cpm[idx,]

    df_idx = rbind(df_idx, c(cpm_th, x, length(idx)))
  }
}
colnames(df_idx) = c("cpm_th", "x_pctSamples", "n_genes_kept")

pdf("3_01_NumGenesKept_inXpctSamples_withDifferentCpmThreshold.pdf")
df_idx$cpm_th = factor(df_idx$cpm_th, levels = seq(0, 1, 0.1))
df_idx %>%
  ggplot(aes(x = x_pctSamples, y = n_genes_kept, col = cpm_th)) +
  geom_point() +
  geom_line()
dev.off()
# Observation: 30% of samples seem to be the general kneeing point. Line cpm_th = 0.3 seems to be stable for n_genes_kept. Also 18733 genes with CPM > 0.3 in at least 30% samples, consistent with Google's suggestion that human have 47k genes, 40% of which are expressed in neurons (~18800 genes).

### Keep genes with cpm > 0.3 in 30% of samples
keep = apply(cpm > 0.3, 1, sum)
idx = which(keep > 0.30*dim(cpm)[2]) ## cpm > 0.3 in 30% of samples
cpm = cpm[idx,]
## 18733 genes genes remain - about 31% of genes pass filter

### Remove any remaining genes with an effective gene length <= 50 bp
rsem_gene_effLen = rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
stopifnot(rownames(cpm) == names(rsem_gene_effLen))

geneAnno <- geneAnno1[match(substr(rownames(cpm),1,15),geneAnno1[,1]),]
geneAnno <- geneAnno[!is.na(geneAnno[,1]),] # 18733 genes have ensembl_gene_id

# Make sure to include key antisense genes at the 17q21 locus
geneAnno$ensembl_gene_id[geneAnno$hgnc_symbol == "MAPT-AS1"] # ENSG00000264589
geneAnno$ensembl_gene_id[geneAnno$hgnc_symbol == "KANSL1-AS1"] # ENSG00000214401
idx1 = which(substr(names(rsem_gene_effLen), 1, 15) == "ENSG00000264589")
idx2 = which(substr(names(rsem_gene_effLen), 1, 15) == "ENSG00000214401")
rsem_gene_effLen[idx1] # 2809.77
rsem_gene_effLen[idx2] # 98.61
hist(rsem_gene_effLen, xlim = c(0,8e3), breaks = 4000) # a bar every 200bp with breaks = 1000

idx = which(rsem_gene_effLen <= 50)
short_remove = rownames(cpm)[idx]
idx_remove = match(short_remove,rownames(cpm)) # 54 genes
cpm = cpm[-idx_remove,]
rsem_gene_effLen = rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
rsem_gene = rsem_gene[match(rownames(cpm), rownames(rsem_gene)),]

dim(rsem_gene)
## 18679 genes, 14 samples

### Remove "_PAR_Y" genes
# As they will incur duplicated rownames error later
idx_remove = which(grepl("_PAR_Y",rownames(rsem_gene))) # 13 genes
rsem_gene = rsem_gene[-idx_remove,]
rsem_gene_effLen = rsem_gene_effLen[-idx_remove]

dim(rsem_gene)
## 18666 genes, 14 samples

geneAnno <- geneAnno[match(substr(rownames(rsem_gene),1,15),geneAnno[,1]),]
save(geneAnno, file = "3_01_BioMart_geneAnno.rda") # save for later use

rm(cpm,keep,short_remove,df_idx,idx,idx1,idx2,idx_remove,lnames,x)
save(list = ls(), file = "3_01_GeneFiltered_Data.rda")

##### (2) Normalization #####
rm(list = ls())
lnames = load("3_01_GeneFiltered_Data.rda")

### Get CQN GC Content and Read Length Adjustment
geneAnno$effective_length = rsem_gene_effLen[match(geneAnno[,1],substr(names(rsem_gene_effLen),1,15))]
rownames(geneAnno) = geneAnno[,1]

cqn.dat <- cqn(rsem_gene, lengths = as.numeric(geneAnno$effective_length), x = as.numeric(geneAnno$percentage_gene_gc_content), lengthMethod = c("smooth"), sqn = FALSE) # run cqn with no quantile normalization
save(cqn.dat, file = "3_02_CQN_results.rda")

### limma-trend normalization, using the cqn offset values
dge2 <- DGEList(counts = rsem_gene)
dge2 <- calcNormFactors(dge2)
logCPM_offset <- cpm(dge2, log = TRUE, offset = cqn.dat$offset, prior.count = 0.5) # prior count helps to dampen variance of low count genes

datExpr = logCPM_offset
rsem_gene_effLen = rsem_gene_effLen[match(rownames(datExpr), names(rsem_gene_effLen))]

save(datExpr, rsem_gene_effLen, datMeta, datSeq_numeric, file="3_02_NormalizedGeneFilteredRNAseq.rda")

##### (3) Outlier Removal #####
rm(list = ls())
load("3_02_NormalizedGeneFilteredRNAseq.rda")

table(datMeta$differentiationbatch) # b1-4: 7 4 2 1
table(datMeta$haplotype) # 6 H1H1, 8 H2H2
table(datMeta[,c("differentiationbatch", "haplotype")])
# batch H1H1 H2H2
# b1    4    3
# b2    1    3
# b3    0    2
# b4    1    0
# some batches contain very few samples and has bias in haplotype

### PCA-based Outlier Removal
norm <- t(scale(t(datExpr),scale=F))
PCA_res <- prcomp(norm,center=FALSE)
varexp <- (PCA_res$sdev)^2 / sum(PCA_res$sdev^2)
sum(varexp[1:8]) # First 8 PCs explain > 80% of total variance
nPC = 8
PCs = as.data.frame(PCA_res$rotation[,1:nPC])

stopifnot(rownames(PCs) == datMeta$cell_line)
PCs$haplotype = datMeta$haplotype
PCs$batch = datMeta$differentiationbatch
PCs$apoe = datMeta$apoe
PCs$sex = datMeta$sex
PCs$rin = datMeta$rin
PCs$read_depth = datSeq_numeric$TOTAL_READS/1e6

p1 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = haplotype)) +
  geom_point(size = 2)
# there is some haplotype effect

p2 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = sex)) +
  geom_point(size = 2)
# there is sex effect

p3 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = apoe)) +
  geom_point(size = 2)
# there seems to be a gradient of E2-E3-E4 from right bottom corner to top left corner.

p4 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = batch)) +
  geom_point(size = 2)
# all mixed

p5 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = rin)) +
  geom_point(size = 2) +
  scale_color_distiller(palette = "RdYlBu")
# not really any rin effect

p6 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = read_depth)) +
  geom_point(size = 2)
  scale_color_distiller(palette = "RdYlBu") +
  labs(col = "read_depth\n(million)")
  # not really any read_depth effect

pdf("3_03_PCA.pdf", height = 4, width = 9)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, byrow = T, nrow = 2, align = "hv")
dev.off()

## No outliers were observed.
## In linear model, should include haplotype, sex, APOE2 and APOE4 dosage + necessary seqCovs

topPC = PCs[,1:nPC]
colnames(topPC) <- paste("PC",c(1:nPC),"_",(signif(varexp[c(1:nPC)],2)*100),"%",sep="")
save(topPC, file = "3_03_topPCs.rda")




