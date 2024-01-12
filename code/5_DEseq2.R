# 4_2_DEseq2_using_established_lm.R

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)

setwd("./results/")

##### Pipeline #####

### (1) Load linear model and data
### (2) Run DESeq2
### (3) Histogram, MA plot, and Dispersion plot

##### (1) Load linear model and data #####
## a. linear model
lnames = load("4_04_final_lm_expr.rda") # expression_lm

## b. count data
lnames = load("3_01_GeneFiltered_Data.rda") # datMeta, datSeq_numeric, rsem_gene, rsem_gene_effLen, geneAnno, geneAnno1, cpm_th

## c. meta data
lnames = load("4_01_Data_with_Scaled_datMeta_datSeq.rda") # datMeta, datSeq_numeric, Covariates, datExpr, rsem_gene_effLen

stopifnot(rownames(datExpr) == rownames(rsem_gene))
#rsem_gene = rsem_gene[match(rownames(datExpr), rownames(rsem_gene)),]
stopifnot(rownames(datMeta) == colnames(rsem_gene))
stopifnot(rownames(datSeq_numeric) == colnames(rsem_gene))

## c. load CQN results
lnames = load("3_02_CQN_results.rda") # cqn.dat
str(cqn.dat$glm.offset)
stopifnot(rownames(cqn.dat$glm.offset) == rownames(rsem_gene))

##### (2) Run DESeq2 #####
library(DESeq2)

## Input data to DESeq2
expression_lm
# [1] "lm(logRPKM ~ haplotypeH2H2 + sexF + APOE2 + APOE4 + PCT.UTR.BASES + GC.NC.60.79 + MEDIAN.5PRIME.TO.3PRIME.BIAS + ESTIMATED.LIBRARY.SIZE + batchb2, data = cur_data)"
Covariates$batchb2 = ifelse(Covariates$batch == "b2", 1, 0)

# Trial 1 - singularity issue
dds <- DESeqDataSetFromMatrix(countData = round(rsem_gene, 0), colData = Covariates, design = ~ haplotype + sex + APOE2 + APOE4 + PCT.UTR.BASES + GC.NC.60.79 + MEDIAN.5PRIME.TO.3PRIME.BIAS + ESTIMATED.LIBRARY.SIZE + batchb2)
dds <- DESeq(dds)
# Error in solve.default(xtwx + ridge) : system is computationally singular: reciprocal condition number = 3.5312e-20

# Trial 2 - singularity issue
dds <- DESeqDataSetFromMatrix(countData = round(rsem_gene, 0), colData = Covariates, design = ~ haplotype + sex + APOE2 + APOE4 + PCT.UTR.BASES + GC.NC.60.79 + MEDIAN.5PRIME.TO.3PRIME.BIAS + ESTIMATED.LIBRARY.SIZE)
dds <- DESeq(dds)
# Error in solve.default(xtwx + ridge) : system is computationally singular: reciprocal condition number = 3.5217e-20

# Trial 3 - successful
dds <- DESeqDataSetFromMatrix(countData = round(rsem_gene, 0), colData = Covariates, design = ~ haplotype + sex + APOE2 + APOE4 + PCT.UTR.BASES + GC.NC.60.79 + MEDIAN.5PRIME.TO.3PRIME.BIAS)
dds <- DESeq(dds)
# estimating size factors
# estimating dispersions
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# 1363 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

## Integrate CQN offset
cqnNormFactors <- exp(cqn.dat$glm.offset)
normFactors <- cqnNormFactors / exp(rowMeans(log(cqnNormFactors)))
normalizationFactors(dds) <- normFactors

## Run DESeq2
dds <- DESeq(dds)
# estimating size factors
# estimating dispersions
# found already estimated dispersions, replacing these
# gene-wise dispersion estimates
# mean-dispersion relationship
# final dispersion estimates
# fitting model and testing
# 1397 rows did not converge in beta, labelled in mcols(object)$betaConv. Use larger maxit argument with nbinomWaldTest

resultsNames(dds) # lists the coefficients: Intercept, haplotype_H2H2_vs_H1H1, sex_F_vs_M, APOE2, APOE4, PCT.UTR.BASES, GC.NC.60.79, MEDIAN.5PRIME.TO.3PRIME.BIAS
res <- results(dds, name="haplotype_H2H2_vs_H1H1")
save(dds, res, file = "5_02_DESeq2_resultHaplotype.rda") # 4_2_02_DESeq2_resultHaplotype_lm_HaplotypeSexAPOE2APOE4PCT.UTR.BASESGC.NC.60.79MEDIAN.5PRIME.TO.3PRIME.BIAS.rda

## Any significant genes?
formatC(range(res$padj, na.rm = T), digits = 2) # new 4.1e-23 to 1; old 2.8e-16 to 1

summary(res)
# out of 18666 with nonzero total read count
# adjusted p-value < 0.1
# LFC > 0 (up)       : 183, 0.98%
# LFC < 0 (down)     : 159, 0.85%
# outliers [1]       : 0, 0%
# low counts [2]     : 2895, 16%
# (mean count < 20)

# Order by pval
resOrdered <- res[order(res$pvalue),]

##### 3) Histogram, MA plot, and Dispersion plot #####
pdf("5_03_Histogram_wCQNglmoffset.pdf") # 4_2_03_Histogram_wCQNglmoffset.pdf
hist(res$pvalue) # More significant than lm()
dev.off()
# very significant at p-val < 0.05

pdf("5_03_MA_and_Dispersion_plots_wCQNglmoffset.pdf")  # 4_2_03_MA_and_Dispersion_plots_wCQNglmoffset.pdf
plotMA(res, ylim=c(-2,2)) # looks good
plotDispEsts(dds)
dev.off()
# looks good

