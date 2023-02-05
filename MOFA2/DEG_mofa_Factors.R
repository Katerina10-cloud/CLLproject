# Gene expression signatures for CLL-PD (F4)

#data("MOFA2_CLL", "gene", "rna")
load("~/CLLproject/MOFA2/MOFA2_CLL.RData")
save(dds, file = "dds_mofa.RData")
save(corRes.sig, file = "corRes.sigF6.RData")

library(MOFA2)
library(genefilter)
library(DESeq2)
library(tidyverse)


library(gridExtra)
library(ggrepel)
library(ggbeeswarm)
library(GEOquery)
library(pheatmap)
library(tidygraph)
library(igraph)
library(ggraph)
library(ComplexHeatmap)
library(grid)
library(cowplot)


#Process MOFA model
#get factor values
facTab <- get_factors(
  MOFAobject, 
  factors = "Factor4",
  as.data.frame = TRUE
) %>% as_tibble()

#get all factors
facTab.all <- get_factors(
  MOFAobject,
  factors = "all",
  as.data.frame = TRUE
) %>% as_tibble()

#get weights
weightTab <- get_weights(
  MOFAobject,
  factors = "all",
  as.data.frame = TRUE
) %>% as_tibble() %>%
  mutate(feature = ifelse(feature == "IGHV.status", "IGHV", feature))

detach("package:MOFA2", unload = TRUE) #detach MOFA because it's "predict" function marks the "predict" function in glmnet

## Pre-processing data
#Processing RNAseq data

dds<-estimateSizeFactors(rna)
dds$factor <- facTab[match(dds$PatID, facTab$sample),]$value
dds$IGHV <- gene[match(dds$PatID, rownames(gene)),]$IGHV
ddsSub <- dds[,!is.na(dds$factor)]

#vst
ddsSub.vst <- varianceStabilizingTransformation(ddsSub)

## Identify genes associated with CLL-PD
#Correlation test using DESeq2
factorMatrix <- facTab.all %>% spread(key = factor, value = value) %>%
  data.frame() %>% column_to_rownames("sample") 
factorMatrix <- factorMatrix[colnames(ddsSub),]

factorMatrix <- factorMatrix[,-1]
designMat <- model.matrix(~1 +.,factorMatrix)
deRes <- DESeq(ddsSub, full = designMat, betaPrior = FALSE)

corRes.rna <- results(deRes, name = "Factor6", tidy = TRUE) %>%
  mutate(symbol = rowData(ddsSub)[row,]$symbol) %>%
  dplyr::rename(logFC = log2FoldChange, t = stat, P.Value = pvalue, adj.P.Val = padj,id=row)

#Heatmap of significantly correlated genes (5% FDR)
corRes.sig <- filter(corRes.rna, adj.P.Val < 0.05)
exprMat <- assay(ddsSub.vst)
plotMat <- exprMat[corRes.sig$id,]
plotMat <- plotMat[order(corRes.sig$logFC, decreasing = TRUE),order(designMat[,"LF4"])]

colAnno <- data.frame(row.names = colnames(exprMat),
                      LF4 = designMat[,"LF4"])
colnames(colAnno) <- c("CLL-PD")
plotMat <- mscale(plotMat, censor = 4)

breaks <- seq(-4,4,length.out = 100)

pheatmap(plotMat, scale = "none", cluster_cols = FALSE, cluster_rows = FALSE,
         annotation_col = colAnno,
         color = colorRampPalette(c(colList[2],"white",colList[1]))(length(breaks)),
         breaks = breaks,
         show_rownames = FALSE, show_colnames = FALSE)

#How many genes show significant correlation?
nrow(corRes.sig)

#percentage
nrow(corRes.sig)/nrow(ddsSub)

#How many genes show up-regulation?
nrow(filter(corRes.sig, t>0))

#How many genes show down-regulation?
nrow(filter(corRes.sig, t<0))
