##### RNA-seq Counts Processing (Gene-Level)
##### Based on 31_Jill_ASD_pancortical_RNAseq_datasets/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD-master/code/01_RNAseqProcessing/01_02_A_CountsProcessing.R from May 2020, Jillian Haney

rm(list = ls())
library(stringr)
library(tidyverse)
library(limma); library(edgeR); library(cqn); library(biomaRt)
#library(WGCNA); library(devtools); library(paralell)
# BiocManager::install("variancePartition")
library(ggplot2); library(gridExtra); library(variancePartition)
library(doParallel); library(lmerTest)

## Install custom 'earth' R package to allow infinite genes as input
Jill_DIR = "~/Documents/Documents/Geschwind_lab/LAB/Database_download/31_Jill_ASD_pancortical_RNAseq_datasets/Broad-transcriptomic-dysregulation-across-the-cerebral-cortex-in-ASD-master/"
#devtools::install(paste(Jill_DIR,"code/01_RNAseqProcessing/earth-infGenes/",sep=""))
#library(earth)

setwd("~/Documents/Documents/Geschwind_lab/LAB/Projects/Project_H1H2/1_Neurons/RNA-seq/Batch1/RNA_analysis/results/3_CountProcessing/")

##### Pipeline #####

### (0) Prepare RawData with RSEM counts, gene effective length, biological and technical covariates
### (1) Gene Filtering
### (2) Normalization
### (3) Outlier Removal
### (4) Prepare Metadata for Linear Model Building

##### (0) Prepare RawData #####
# Load RSEM counts
lnames = load("../1_v2_rsem/rsem_results.rda") # df_counts, df_tpm, gene_length, gene_effective_length
# Load metadata
lnames = load("../2_Covariates/Covariates.rda") # Cov_Biol, Cov_Tech
all(colnames(df_counts) %in% Cov_Biol$cell_line) # T
df_counts = df_counts[,Cov_Biol$cell_line] # Checked that df_counts rather than df_tpm looks more like Jill's rsem_gene.

rsem_gene = df_counts # 60708 genes
datMeta = Cov_Biol
datSeq_numeric = Cov_Tech
rsem_gene_effLen = gene_effective_length

save(rsem_gene, rsem_gene_effLen, datMeta, datSeq_numeric, file = "3_0_RawData.rda")

rm(list = ls())
lnames = load("3_0_RawData.rda") # df_counts, gene_effective_length, Cov_Biol, Cov_Tech

##### (1) Gene Filtering #####

## Only keep transcripts on chr1-Y
getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype")

ensembl_db = listDatasets(useEnsembl(biomart="ENSEMBL_MART_ENSEMBL")) # search for "hsapiens_gene_ensembl" got GRCh38
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")  

geneAnno1 <- getBM(attributes = getinfo, filters = c("ensembl_gene_id"), values = substr(rownames(rsem_gene),1,15), mart = mart) 
geneAnno2 = geneAnno1[geneAnno1$chromosome_name %in% c(1:22, "X", "Y"),]
idx = which(substr(rownames(rsem_gene),1,15) %in% geneAnno2$ensembl_gene_id)
rsem_gene = rsem_gene[idx,]
rsem_gene_effLen = rsem_gene_effLen[idx]

save(rsem_gene, rsem_gene_effLen, datMeta, datSeq_numeric, file = "3_v2_01_RawData_chr1toY.rda")

### Decide on cpm threshold and x% of samples for keeping genes
# Google: As a general rule, a good threshold can be chosen for a CPM value that corresponds to a count of 10.
rm(list = ls())
lnames = load("3_v2_01_RawData_chr1toY.rda")
range(datSeq_numeric$READ_PAIRS_EXAMINED) # 7.2-11 M

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

pdf("3_v2_01_NumGenesKept_inXpctSamples_withDifferentCpmThreshold.pdf")
df_idx$cpm_th = factor(df_idx$cpm_th, levels = seq(0, 1, 0.1))
df_idx %>%
  ggplot(aes(x = x_pctSamples, y = n_genes_kept, col = cpm_th)) +
  geom_point() +
  geom_line()
dev.off()
# Observation: 30% of samples seem to be the general kneeing point.

# For making ppt, just draw cpm > 1 line
pdf("3_v2_01_NumGenesKept_inXpctSamples_withCpmGT1.pdf", width = 4, height = 3)
df_idx[df_idx$cpm_th == 1, ] %>%
  ggplot(aes(x = x_pctSamples, y = n_genes_kept)) +
  geom_point(col = "dodgerblue") +
  geom_line(col = "dodgerblue") +
  theme_bw() +
  geom_vline(xintercept = 0.30, linetype = "dashed", color = "red") +
  geom_hline(yintercept = df_idx$n_genes_kept[which(df_idx$cpm_th == 1 & df_idx$x_pctSamples == 0.30)], linetype = "dashed", color = "red") +
  #annotate("text", x = 0.45, y = 18800, label = "cpm > 1") +
  ylab("Number of genes\nwith cpm > 1") +
  xlab("Percent of samples")
dev.off()

### Keep genes with cpm > 1 (based on the read depth) in 30% of samples
keep = apply(cpm > 1, 1, sum) 
idx = which(keep > 0.30*dim(cpm)[2]) ## cpm > 1 in 15% of samples
cpm = cpm[idx,]
## 16504 genes genes remain - about 27% of genes pass filter

### Remove any remaining genes with an effective gene length <= 50 bp
rsem_gene_effLen = rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
stopifnot(rownames(cpm) == names(rsem_gene_effLen))

lnames = load("3_02_BioMart_geneAnno.rda") # geneAnno
geneAnno$ensembl_gene_id[geneAnno$hgnc_symbol == "MAPT-AS1"] # ENSG00000264589
geneAnno$ensembl_gene_id[geneAnno$hgnc_symbol == "KANSL1-AS1"] # ENSG00000214401
idx1 = which(substr(names(rsem_gene_effLen), 1, 15) == "ENSG00000264589")
idx2 = which(substr(names(rsem_gene_effLen), 1, 15) == "ENSG00000214401")
rsem_gene_effLen[idx1] # 2809.77
rsem_gene_effLen[idx2] # 98.61

hist(rsem_gene_effLen, xlim = c(0,8e3), breaks = 4000) # a bar every 200bp with breaks = 1000

idx = which(rsem_gene_effLen <= 50)
short_remove = rownames(cpm)[idx]
idx_remove = match(short_remove,rownames(cpm)) # 4 genes
cpm = cpm[-idx_remove,]
rsem_gene_effLen = rsem_gene_effLen[match(rownames(cpm), names(rsem_gene_effLen))]
rsem_gene = rsem_gene[match(rownames(cpm), rownames(rsem_gene)),]

rm(cpm,keep,short_remove,idx,idx_remove,lnames)

dim(rsem_gene)
## 16474 genes, 14 samples

### Remove "_PAR_Y" genes
# As they will incur duplicated rownames error later
idx_remove = which(grepl("_PAR_Y",rownames(rsem_gene))) # 11 genes
rsem_gene = rsem_gene[-idx_remove,]
rsem_gene_effLen = rsem_gene_effLen[-idx_remove]

dim(rsem_gene)
## 16463 genes, 14 samples

save(list = ls(), file = "3_v2_01_GeneFiltered_Data.rda")

##### (2) Normalization #####

counts = rsem_gene; rm(rsem_gene)

### Get CQN GC Content and Read Length Adjustment

getinfo <- c("ensembl_gene_id","hgnc_symbol","chromosome_name","start_position",
             "end_position","strand","band","gene_biotype","percentage_gene_gc_content")

ensembl_db = listDatasets(useEnsembl(biomart="ENSEMBL_MART_ENSEMBL")) # search for "hsapiens_gene_ensembl" got GRCh38
mart <- useMart(biomart="ENSEMBL_MART_ENSEMBL",
                dataset="hsapiens_gene_ensembl"#,
                #host="grch38.ensembl.org" # default www.ensembl.org - which is currently grch38
                )  

geneAnno1 <- getBM(attributes = getinfo, filters = c("ensembl_gene_id"), values = substr(rownames(counts),1,15), mart = mart) 

geneAnno <- geneAnno1[match(substr(rownames(counts),1,15),geneAnno1[,1]),]
counts_cqn <- counts[!is.na(geneAnno[,1]),]
geneAnno <- geneAnno[!is.na(geneAnno[,1]),] # 17283 genes have ensembl_gene_id
save(geneAnno, file = "3_v2_02_BioMart_geneAnno.rda")

geneAnno$effective_length = rsem_gene_effLen[match(geneAnno[,1],substr(names(rsem_gene_effLen),1,15))] 
rownames(geneAnno) = geneAnno[,1] # Initially encountered Error "duplicate 'row.names' are not allowed". Run the following lines and solved the issue by removing "_PAR_Y" genes.

cqn.dat <- cqn(counts_cqn, lengths = as.numeric(geneAnno$effective_length), x = as.numeric(geneAnno$percentage_gene_gc_content), lengthMethod = c("smooth"), sqn = FALSE) ## Run cqn with no quantile normalization
save(cqn.dat, file = "3_v2_02_CQN_results.rda")

### limma-trend normalization, using the cqn offset values
dge2 <- DGEList(counts = counts_cqn)
dge2 <- calcNormFactors(dge2)
logCPM_offset <- cpm(dge2, log = TRUE, offset = cqn.dat$offset, prior.count = 0.5) # prior count helps to dampen variance of low count genes

datExpr = logCPM_offset
rsem_gene_effLen = rsem_gene_effLen[match(rownames(datExpr), names(rsem_gene_effLen))]

save(datExpr, rsem_gene_effLen, datMeta, datSeq_numeric, file="3_v2_02_NormalizedGeneFilteredRNAseq.rda")

##### (3) Outlier Removal #####
rm(list = ls())
load("3_v2_02_NormalizedGeneFilteredRNAseq.rda")

table(datMeta$differentiationbatch) # b1-4: 7 4 2 1
table(datMeta$haplotype) # 6 H1H1, 8 H2H2
table(datMeta[,c("differentiationbatch", "haplotype")]) 
# batch H1H1 H2H2
# b1    4    3
# b2    1    3
# b3    0    2
# b4    1    0

### PCA-based Outlier Removal
# Jill calculated which samples have a z-score greater than 3 on the first 10 PCs by batch and brain-region.
# In my case, the batch-haplotype groups are at max 4 samples, too small. So I'll just plot PCA and color by batch or haplotype
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
PCs$read_depth = datSeq_numeric$Read_depth/1e6

p1 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = haplotype)) +
  geom_point(size = 2) #+
  #theme(legend.position = "top")
# looks fine

p2 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = sex)) +
  geom_point(size = 2) #+
  #theme(legend.position = "top")
# looks fine

p3 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = apoe)) +
  geom_point(size = 2) #+
  #theme(legend.position = "top")
# there seems to be a gradient of E2-E3-E4 from right bottom corner to top left corner.

p4 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = batch)) +
  geom_point(size = 2) #+
  #theme(legend.position = "top")
# there seems to be sex effect

p5 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = rin)) +
  geom_point(size = 2) +
  #scale_color_gradient2(midpoint = mean(PCs$rin), low = "blue", mid = "yellow", high = "red") #+
  scale_color_distiller(palette = "RdYlBu") #+
#theme(legend.position = "top")
# not really any rin effect

p6 = PCs %>%
  ggplot(aes(x = PC1, y = PC2, col = read_depth)) +
  geom_point(size = 2) +
  #scale_color_gradient2(midpoint = mean(PCs$rin), low = "blue", mid = "yellow", high = "red") #+
  scale_color_distiller(palette = "RdYlBu") +
  labs(col = "read_depth\n(million)")

pdf("3_v2_03_PCA.pdf", height = 4, width = 9)
cowplot::plot_grid(p1, p2, p3, p4, p5, p6, byrow = T, nrow = 2, align = "hv")
dev.off()

## No outliers were observed.
## In linear model, should include haplotype, sex and APOE2 and APOE4 dosage + necessary seqCovs

topPC = PCs[,1:nPC]
colnames(topPC) <- paste("PC",c(1:nPC),"_",(signif(varexp[c(1:nPC)],2)*100),"%",sep="")
save(topPC, file = "3_v2_03_topPCs.rda")

##### (4) Prepare Metadata for Linear Model Building #####
rm(list=ls())
load("3_v2_02_NormalizedGeneFilteredRNAseq.rda")

### Format meta data for input into EARTH/MARS
datMeta$haplotype = factor(datMeta$haplotype, levels = c("H1H1", "H2H2"))
datMeta$sex = factor(datMeta$sex, levels = c("M", "F"))
colnames(datMeta) = gsub("differentiation", "", colnames(datMeta))
datMeta$batch = factor(datMeta$batch) # Levels: b1 b2 b3 b4

datMeta = as.data.frame(datMeta)
rownames(datMeta) = datMeta$cell_line
datMeta_model = datMeta[,c(4,5,9,7,10:12)] # Missing Age information for a lot of samples

### Scale the continuous datMeta_model variables
datMeta_model$rin = scale(datMeta_model$rin)
datMeta_model$rin <- sapply(datMeta_model$rin, function(x) { attributes(x) <- NULL; x }) # remove the attributes of scale or any other columns

### Scale the datSeq_numeric variables
#idx_no_scale = grep("N",datSeq[1,]) # not sure what this means in Jill's script. I did not run here.
#datSeq_scaled = apply(datSeq_numeric[,-idx_no_scale],2,scale)
rownames(datSeq_numeric) = datSeq_numeric$cell_line
datSeq_numeric = datSeq_numeric[,-1]
datSeq = as.data.frame(apply(datSeq_numeric, 2, scale))
rownames(datSeq) = rownames(datSeq_numeric)
datSeq_model = datSeq

save(datMeta_model, datSeq_model, datMeta, datSeq_numeric, datExpr, rsem_gene_effLen, file = "3_v2_04_Data_with_Scaled_datMeta_datSeq.rda")

# --- stop here ----

### Filter datSeq_model so that I am only evaluating terms that are not collinear with other covariates.
## For each covariate identify other covariates that have an adjusted R2 > 0.95.
## First, identify seq covariates (continuous) that are confounds for any meta covariates. 
allmat_meta = list()

for(i in c(1:dim(datMeta_model)[2])){
  tmp1 = datMeta_model[,i]
  allmat_meta[[colnames(datMeta_model)[i]]]=NA
  for(j in c(1:dim(datSeq_model)[2])){
    tmp2 = datSeq_model[,j]
    if(is.factor(tmp1)==TRUE & is.factor(tmp2)==FALSE){
      mod=summary(lm(tmp2~tmp1))
      r2=mod$adj.r.squared
      if(r2 >= 0.95){
        allmat_meta[[colnames(datMeta_model)[i]]] = c(allmat_meta[[colnames(datMeta_model)[i]]],colnames(datSeq_model)[j])
      }
    }else{
      mod=summary(lm(tmp1~tmp2))
      r2=mod$adj.r.squared
      if(r2 >= 0.95){
        allmat_meta[[colnames(datMeta_model)[i]]] = c(allmat_meta[[colnames(datMeta_model)[i]]],colnames(datSeq_model)[j])
      }
    }
  }
  allmat_meta[[colnames(datMeta_model)[i]]] = allmat_meta[[colnames(datMeta_model)[i]]][-1]
}
## No seq covariates show R2 > 0.95 with meta covariates

### Now identify confounds among the remaining Seq covariates
allmat_seq = list()

for(i in c(1:dim(datSeq_model)[2])){
  tmp1 = datSeq_model[,i]
  allmat_seq[[colnames(datSeq_model)[i]]]=NA
  for(j in c(1:dim(datSeq_model)[2])){
    tmp2 = datSeq_model[,j]
    if(i==j){
      next
    }else{
      mod=summary(lm(tmp1~tmp2))
      r2=mod$adj.r.squared
      if(r2 >= 0.95){
        allmat_seq[[colnames(datSeq_model)[i]]] = c(allmat_seq[[colnames(datSeq_model)[i]]],colnames(datSeq_model)[j])
      }
    }
  }
  allmat_seq[[colnames(datSeq_model)[i]]] = allmat_seq[[colnames(datSeq_model)[i]]][-1]
}
# Warning mesages:
# In summary.lm(lm(tmp1 ~ tmp2)) :
#   essentially perfect fit: summary may be unreliable
## Did see colinearity among seq covariates 

## Iterate through each covariate, determine confounds.
## For confounds, determine which of the covariates has the highest overall association with the top 5 expr PCs (adj R2).
## Then keep this covariate and remove the remaining from the loop.

## Get top expression PCs.

norm <- t(scale(t(datExpr),scale=F))
PC <- prcomp(norm,center = FALSE)
varexp <- (PC$sdev)^2 / sum(PC$sdev^2)
sum(varexp[1:8])
nPC = 8
topPC <- PC$rotation[,1:nPC] ## these first 8 explain > 80% of the variance
colnames(topPC) <- paste("PC",c(1:nPC),"_",(signif(varexp[c(1:nPC)],2)*100),"%",sep="")

## Now begin filtering algorithm.

allmat_seq_org = allmat_seq
covs_iterate = covs_keep = names(allmat_seq_org)

for(cov in covs_iterate){
  print(cov)
  if((cov %in% covs_keep)==FALSE){
    next
  }
  if(length(allmat_seq[[cov]])==0){
    next
  }else{
    keep0 = which(allmat_seq[[cov]] %in% covs_keep)
    covs_test = c(cov,allmat_seq[[cov]][keep0])
    r2_mat = matrix(NA,nrow=nPC,ncol=length(covs_test))
    rownames(r2_mat)=colnames(topPC)[1:nPC]
    colnames(r2_mat)= covs_test
    for(i in c(1:nPC)){
      for(j in c(1:length(covs_test))){
        mod = summary(lm(topPC[,i] ~ datSeq_model[,which(colnames(datSeq_model)==covs_test[j])]))
        r2_mat[i,j]=mod$adj.r.squared
      }
    }
    r2_sum = apply(r2_mat,2,sum)
    max_val = max(r2_sum)
    max_cov = names(r2_sum)[which(r2_sum==max_val)]
    if(length(max_cov) > 1){
      keep = sample(max_cov,1)
    }else{
      keep = max_cov
    }
    allmat_seq[[keep]] = allmat_seq[[keep]][-c(1:length(allmat_seq[[keep]]))]
    remove = covs_test[-which(covs_test==keep)]
    covs_keep = covs_keep[-which(covs_keep %in% remove)]
    print(length(covs_keep))
  }
}

## Keeping 31 covariates (removing 12 covariates).
### Filter datSeq_model.

datSeq_model = datSeq_model[,which(colnames(datSeq_model) %in% covs_keep)]

save(datMeta_model, datSeq_model, datExpr, file="3_04_dat4MARS.rda")

# ------------------ Tried the later parts, but STOP here -----------------------
# ------------------ Use the next script to build linear model ------------------

##### (5) Run CV MARS #####
### Run EARTH/MARS.
### Use earth-inifite R package installation.
### It is recommended to run this section on a high-performance computing cluster.
# I failed to install package 'paralell' on hoffman2 R/4.0.2 or 4.2.2

lnames = load("3_04_dat4MARS.rda")

runMARS <- function(gene.expr, covars, n.predictors=NULL, 
                    n.replicates=1, n.cores=1, allow.interaction=FALSE, 
                    batch.interactions=FALSE,allow.nonlinear=FALSE,cr=10) {
  
  # uses the packages `earth` to determine an appropriate linear model for expression, given
  # technical covariates.
  # Inputs
  #  gene.expr     - The gene expression matrix (genes x samples)
  #  covars        - The sample covariate matrix (samples x covariates)
  #  n.predictors  - The (maximum) number of predictors to use in the forward earth  model
  #  no.cross      - vector of covariates to exclude from cross terms
  #
  # Returns:
  #  earth         - the fitted `earth` object
  #  formula       - formula giving the right-hand-size of the covariate-correction LM
  #  terms         - the terms included
  #  term.imp      - importance for each term: this is the 80th percentile of max(beta)
  #  model.matrix  - the model matrix resulting from model.matrix(formula, covars)
  
  gene.expr = t(gene.expr)
  
  binary.terms = apply(covars,2,function(x) { 
    return(length(levels(as.factor(x)))==2 && min(x)==0 && max(x)==1)})
  square.terms = covars[,!binary.terms]^2
  colnames(square.terms) = paste0( colnames(square.terms), "^2")
  covars = cbind(covars, square.terms)
  
  if(allow.interaction==TRUE) {
    degree=2;
  } else {
    degree=1;
  }
  
  allowed.fx <- function(degree, pred, parents, namesx, first) {
    if(batch.interactions==FALSE & degree > 1) {
      bad.idx <- grep("BATCH", toupper(namesx))
      
      return(!(pred %in% bad.idx || which(parents!=0) %in% bad.idx))
    } else {
      return(TRUE)
    }
    
  }
  
  if(allow.nonlinear==TRUE) {
    linear=FALSE;
  } else {
    linear=TRUE;
  }
  
  out = mclapply(1:n.replicates, function(r) {
    e = earth(x=covars, y=gene.expr, trace=2, degree=degree,linpreds=linear, allowed=allowed.fx,nfold=cr)
    med= apply(abs(e$coefficients),1,median)
    beta = data.frame(row.names=1:length(med), var=as.character(names(med)), abs.beta=med, replicate=r, 
                      modfit = paste0("RSS=",e$rss, " RSQ=",e$rsq, " GCV=",e$gcv, " GRSQ=",e$grsq))
    all = list("beta"=beta,"e"=e)
  }, mc.cores=n.cores)
  
  return(out)
}

datSeq=datSeq_model; rm(datSeq_model) 
datMeta=datMeta_model; rm(datMeta_model)
# Note that both are 'data.frame' by str()

form_datMeta=paste0(colnames(datMeta)[c(1:ncol(datMeta))],collapse=" + ")
form_datMeta=paste0("~ ",form_datMeta,collapse="")

form_datSeq=paste0(colnames(datSeq)[c(1:ncol(datSeq))],collapse=" + ")
form_datSeq=paste0("~ 0 + ",form_datSeq,collapse="")

Y = datExpr; rm(datExpr)

X.meta = data.frame(model.matrix(as.formula(form_datMeta),data=datMeta))
X.all = data.frame(model.matrix(as.formula(form_datSeq),data=datSeq))

Y = data.frame(Y[,match(rownames(X.meta), colnames(Y))]) # Though I know that everything is matched already

## add meta to the two seq model matrices

X.all=data.frame(X.meta,X.all)

#Remove covariates with 0 variance

to_keep = !(apply(X.all,2,sd)==0) # only intercept will be removed
X.all = X.all[,to_keep]

output="./"

if(TRUE | !file.exists(paste0(output,"/3_05_MARSOutput.RData")))  {
  e.all = runMARS(gene.expr = Y, covars = X.all, n.replicates = 1, allow.interaction = F,batch.interactions=F,allow.nonlinear=F)
  save(e.all, file=paste0(output,"/3_05_MARSOutput.RData"))
  
  beta_e_all=e.all[[1]]$beta
  coef_e_all=e.all[[1]]$e$coefficients
  
  save(beta_e_all,coef_e_all,file=paste0(output,"/3_05_MARSOutput_subset.RData"))
}
# 2023/02/03 6:08PM, done in 1min

### Now, compile cross-validated MARS output

cv_coef_list=list()
cv_r2_list=list()

for(i in c(1:10)){
  cv_coef_list[[i]]=e.all[[1]]$e$cv.list[[i]]$coefficients
  cv_r2_list[[i]]=e.all[[1]]$e$cv.list[[i]]$rsq
}

#tmp = e.all[[1]]$e$cv.rsq.tab
cv_rsq_all=e.all[[1]]$e$cv.rsq.tab[,"mean"] # Lots of -Inf, especially in fold 5-10

save(cv_coef_list,cv_r2_list,cv_rsq_all,file=paste0(output,"3_05_MARSOutput_CVCompiled.rda"))

# See if I can continue

##### (6) Build Linear Model #####

rm(list=ls())

### It is recommended to use the provided CV MARS data to be consistent with the publication.
lnames = load("3_03_topPCs.rda") # topPC
lnames = load("3_04_Data_with_Scaled_datMeta_datSeq.rda") # "datMeta_model" "datSeq_model" "datMeta" "datSeq_numeric" "datExpr" "rsem_gene_effLen"
lnames = load("3_05_MARSOutput_subset.RData") # beta_e_all, coef_e_all
lnames = load("3_05_MARSOutput_CVCompiled.rda") # cv_coef_list, cv_r2_list, cv_rsq_all

pdf(file="3_06_MARS_CV_plots.pdf", width=16, height=12)
for(i in c(1:10)){
  print(i)
  dat = rbind(cv_coef_list[[i]])
  med_dat=abs(apply(dat,1,median))
  plot_dat=data.frame("Coef"=rownames(dat),"Median_Beta"=med_dat)
  plot_dat=plot_dat[-grep("(Intercept)",plot_dat[,1]),] # there is only (Intercept), so now nothing left. - Stop here, and I shall write my own script to build the linear model. 
  
  unq_rsq=signif(cv_r2_list[[i]],3)
  cv_rsq=signif(cv_rsq_all[[i]],3)
  
  g1 <- ggplot(plot_dat, aes(x=reorder(Coef, Median_Beta), y=Median_Beta)) + 
    geom_bar(stat="identity")  +
    coord_flip() + 
    geom_label(aes(y=Inf,x=-Inf,hjust=1,vjust=-0.1,label=paste0("R^2=",unq_rsq,"; CV R^2=",cv_rsq)),
               size=8) +
    xlab("Absolute Median Beta") +
    ylab("Coefficient") + 
    ggtitle(paste("Median Beta Across All Genes: Iter ",i)) +
    theme(plot.title = element_text(hjust = 0.5,size=28),
          axis.text.x =element_text(size=22),
          axis.title=element_text(size=26))
  
  grid.arrange(g1)
}
dev.off()



