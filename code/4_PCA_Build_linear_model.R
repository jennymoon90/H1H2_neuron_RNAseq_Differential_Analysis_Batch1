# 4_1_PCA_build_linear_model_APOE33AsBaselineE2andE4asDosage_includeSex_VIF2p5.R

### Get top PCs that explain total 80% variance in the CQN normalized peak count (logRPKM), identify biological and technical covariates that significantly correlate with these top PCs.

# We include read depth into covariates.
# We use apoe33 as baseline, and use E2 and E4 dosage as covariates.
# Use VIF cutoff of 2.5

##### Pipeline #####

### 1) Load data and format
### 2) Load top PCs of datExpr
### 3) Calculate correlation between covariates and topPCs of datExpr
### 4) Select covariates for linear model

####################

rm(list = ls())
options(stringsAsFactors = F)
library(tidyverse)
library(ggplot2)
library(corrplot)
library(olsrr)

setwd("./results/")

##### (1) Load data and format #####
## Load limma.trend normalized logRPKM with GC correction by CQN, and metadata and seq covariates

load("3_02_NormalizedGeneFilteredRNAseq.rda")

datMeta$haplotype = factor(datMeta$haplotype, levels = c("H1H1", "H2H2"))
datMeta$sex = factor(datMeta$sex, levels = c("M", "F"))
colnames(datMeta) = gsub("differentiation", "", colnames(datMeta))
datMeta$batch = factor(datMeta$batch) # Levels: b1 b2 b3 b4

datMeta = as.data.frame(datMeta)
rownames(datMeta) = datMeta$cell_line
datMeta = datMeta[,-c(1:3,8)] # use APOE33 as baseline, and E2 and E4 dosage as covariates

rownames(datSeq_numeric) = datSeq_numeric$cell_line
datSeq_numeric = datSeq_numeric[,-1]
colnames(datSeq_numeric) = gsub("_", ".", colnames(datSeq_numeric))
# Note that I'm using the raw picard metrics values rather than scaled ones.
# no need to scale the datSeq_numeric

stopifnot(rownames(datMeta) == colnames(datExpr))
stopifnot(rownames(datSeq_numeric) == colnames(datExpr))

Covariates = cbind(datMeta[,-c(3,7)], datSeq_numeric) # 14 samples with 49 variables,

save(datMeta, datSeq_numeric, Covariates, datExpr, rsem_gene_effLen, file = "4_01_Data_with_Scaled_datMeta_datSeq.rda") # no need of rsem_gene_effLen later, but save it anyway

##### 2) Load top PCs of datExpr #####
lnames = load("3_03_topPCs.rda") # topPC

##### 3) Calculate correlation between covariates and topPCs of datExpr #####

## Calculate correlation of covariates with the top PCs
mod_mat_expr = paste(c(colnames(Covariates)), collapse = " + ")
mod_mat_expr = paste0("~ ", mod_mat_expr)
mod_mat = model.matrix(eval(parse(text = mod_mat_expr)), data = Covariates)[,-1] # Remove intercept
save(mod_mat, file = "4_03_mod_mat.rda")
mod_mat_withPC = cbind(topPC, mod_mat)

Cor = cor(mod_mat_withPC)
Cor_spearman = cor(mod_mat_withPC, method = "spearman")
colnames(Cor)
ncol(topPC) # 8, correct.
idx_spearman = ncol(topPC) + c(1:2,4:8)
Cor[idx_spearman,] = Cor_spearman[idx_spearman,]; Cor[,idx_spearman] = Cor_spearman[,idx_spearman]

## Find out which Covariates significantly correlate with the topPCs (like what Yuyan described in Feng 2022 BioRxiv)

# Correlation significance p-values
Cor_sig = matrix(nrow = nrow(Cor), ncol = ncol(Cor))
rownames(Cor_sig) = colnames(Cor_sig) = colnames(Cor)
for (i in 1:ncol(mod_mat_withPC)) {
  for (j in 1:ncol(mod_mat_withPC)) {
    tmp = cor.test(mod_mat_withPC[,i], mod_mat_withPC[,j])
    Cor_sig[i,j] = tmp$p.value
    # if (tmp$p.value < 0.05) {print(paste(i, colnames(Cor)[i], j, colnames(Cor)[j]))}
  }
}
# View(Cor_sig)
# Observation: Haplotype highly correlates with PC7 and PC3 (cor_pval = 0.0025 and 0.14 respectively)

# Correlation significance FDR
Cor_fdr = matrix(p.adjust(Cor_sig, method = "fdr"), nrow = nrow(Cor_sig))
rownames(Cor_fdr) = colnames(Cor_fdr) = colnames(Cor_sig)
for (i in 1:ncol(mod_mat_withPC)) {
  Cor_fdr[i,i] = 1
}
# View(Cor_fdr)
# Observation: Haplotype significantly correlates with PC7 (cor_fdr = 0.02)

### Get candidate covariates by Cor_FDR or Cor(_coef)
# Option a) Criteria: FDR < 0.05
Cov_pool_idx = which(apply(Cor_fdr[1:ncol(topPC),], 2, function(x) any(x < 0.05)))
(Cov_pool = colnames(Cor_fdr)[Cov_pool_idx]) # haplotypeH2H2, APOE4, MEDIAN.5PRIME.BIAS, PCT.CODING.BASES, PCT.INTRONIC.BASES, PCT.MRNA.BASES, PCT.RIBOSOMAL.BASES, PCT.UTR.BASES, RIBOSOMAL.BASES
# Option b) Criteria: |cor| > 0.4 [choose this in the end]
Cov_pool_idx = which(apply(Cor[1:ncol(topPC), (ncol(topPC) + 1) : ncol(Cor)], 2, function(x) any(abs(x) > 0.4)))
nPC = 8
(Cov_pool = colnames(Cor)[Cov_pool_idx + nPC]) # 49 columns, which means only 2 covariates (from mod_mat) were not included.
# Best to include as many as possible, exclude the ones that will make VIF > 2.5. Choose b) |cor| > 0.4

### For the easiness of selecting covariates, plot corrplot by ranking covs that show max sum of weighted correlation with the topPCs.

# Get the variance explained by each PC.
colnames(Cor)[1:nPC] # see the variance explained by each PC
VarExp = sapply(colnames(Cor)[1:nPC], function(x) unlist(str_split(x, "_"))[2])
VarExp = gsub("%", "", VarExp)
VarExp = as.numeric(VarExp)/100

Cov_pool
VarExp_idx = 1:nPC
Cov_select_df = data.frame(Cov_pool = Cov_pool, sumCorTimesVarExp = sapply(Cov_pool, function(x) t(abs(Cor[VarExp_idx,which(colnames(Cor) == x)])) %*% VarExp[VarExp_idx]))
Cov_select_df = Cov_select_df %>%
  arrange(desc(sumCorTimesVarExp))
# Note that haplotypeH2H2 doesn't come the first, but that's normal and fine.

# Corrplot of ranked candiate covariates (Cov_pool) with top PCs
pdf(paste0("4_03_Corrplot_topPCofLimmaTrendCQNnormedlogRPKM_Covariates_CovPoolOnly.pdf"), height = 35, width = 45)
corrplot(Cor[c(colnames(Cor)[1:ncol(topPC)], Cov_select_df$Cov_pool), c(colnames(Cor)[1:ncol(topPC)], Cov_select_df$Cov_pool)],method="ellipse",tl.pos = "lt",tl.col = "black", tl.srt = 45, addCoef.col = "darkgrey", tl.cex = 2, cl.cex = 2, number.cex = 2)
dev.off()
# Observation: See large blocks of colinear covariates

##### 4) Select covariates for linear model #####
Cov_pool
# [1] "haplotypeH2H2"                  "sexF"
# [3] "rin"                            "batchb2"
# [5] "batchb3"                        "batchb4"
# [7] "APOE2"                          "APOE4"
# [9] "ALIGNED.READS"                  "AT.DROPOUT"
# [11] "CODING.BASES"                   "CORRECT.STRAND.READS"
# [13] "ESTIMATED.LIBRARY.SIZE"         "GC.DROPOUT"
# [15] "GC.NC.0.19"                     "GC.NC.20.39"
# [17] "GC.NC.40.59"                    "GC.NC.60.79"
# [19] "GC.NC.80.100"                   "INCORRECT.STRAND.READS"
# [21] "INTERGENIC.BASES"               "INTRONIC.BASES"
# [23] "MEDIAN.3PRIME.BIAS"             "MEDIAN.5PRIME.BIAS"
# [25] "MEDIAN.5PRIME.TO.3PRIME.BIAS"   "MEDIAN.CV.COVERAGE"
# [27] "NUM.R1.TRANSCRIPT.STRAND.READS" "NUM.R2.TRANSCRIPT.STRAND.READS"
# [29] "NUM.UNEXPLAINED.READS"          "PCT.CODING.BASES"
# [31] "PCT.CORRECT.STRAND.READS"       "PCT.INTERGENIC.BASES"
# [33] "PCT.INTRONIC.BASES"             "PCT.MRNA.BASES"
# [35] "PCT.RIBOSOMAL.BASES"            "PCT.USABLE.BASES"
# [37] "PCT.UTR.BASES"                  "PERCENT.DUPLICATION"
# [39] "PF.ALIGNED.BASES"               "PF.BASES"
# [41] "PF.READS"                       "READ.PAIRS.EXAMINED"
# [43] "READ.PAIR.DUPLICATES"           "READ.PAIR.OPTICAL.DUPLICATES"
# [45] "RIBOSOMAL.BASES"                "TOTAL.CLUSTERS"
# [47] "TOTAL.READS"                    "UNMAPPED.READS"
# [49] "UTR.BASES"

### Test whether the full model violates rule "all VIF < 2.5"

## Full model
colnames(mod_mat)
logRPKM = datExpr[1,]
cur_data = as.data.frame(cbind(logRPKM, mod_mat))
expression_lm = paste0("lm(logRPKM ~ ", paste(Cov_pool, collapse = " + "), ", data = cur_data)")
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction))
# All covariates get VIF = Inf, too many covariates - overfit

## Based on the PCA plot (3_03_PCA.plot), include Haplotype, sex, APOE2 and APOE4 first, and then go down the columns to ensure all VIF < 2.5

## step1. set up base model
expression_lm = "lm(logRPKM ~ haplotypeH2H2 + sexF + APOE2 + APOE4 + PCT.UTR.BASES, data = cur_data)"
fit_infunction <- eval(parse(text = expression_lm))
(vif_df_infunction = ols_vif_tol(fit_infunction)) # Google search "vif cutoff": The default VIF cutoff value is 5
# all VIF < 2, great!

## step2. iterate to add more covariates
tmp = unlist(strsplit(expression_lm, split = ",|~"))[2]
tmp2 = gsub(" ", "", tmp)
cur_covs = unlist(strsplit(tmp2, "\\+"))
all(cur_covs %in% Cov_select_df$Cov_pool) # T

for (i in 1:nrow(Cov_select_df)) {
  if (! Cov_select_df$Cov_pool[i] %in% cur_covs) {
    print(paste0("i=", i, ", Try add covariate ", Cov_select_df$Cov_pool[i]))
    exp1 = unlist(strsplit(expression_lm, ","))[1]
    exp2 = unlist(strsplit(expression_lm, ","))[2]
    expA = paste0(" + ", Cov_select_df$Cov_pool[i], ",")
    test_lm = paste0(exp1, expA, exp2)

    fit_infunction <- eval(parse(text = test_lm))
    vif_df_infunction = ols_vif_tol(fit_infunction)

    if (all(vif_df_infunction$VIF < 2.5)) {
      expression_lm = test_lm
      print("covariate added")
      print(expression_lm)
      i = i+1
    } else {
      print("covariate denied")
      i = i+1
    }
  }
}

## Final model:
expression_lm
# [1] "lm(logRPKM ~ haplotypeH2H2 + sexF + APOE2 + APOE4 + PCT.UTR.BASES + GC.NC.60.79 + MEDIAN.5PRIME.TO.3PRIME.BIAS + ESTIMATED.LIBRARY.SIZE + batchb2, data = cur_data)"

save(expression_lm, file = "4_04_final_lm_expr.rda")
