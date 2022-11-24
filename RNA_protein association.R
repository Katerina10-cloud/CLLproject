library(limma)
library(DESeq2)
library(cowplot)
library(proDA)
library(pheatmap)
library(SummarizedExperiment)
library(tidyverse)
library(dplyr)

#load datasets
load("CLL.methylation.RData")
load("CLLproteomics-main/data/ddsrna_enc.RData")
load("CLLproteomics-main/data/proteomic_explore_enc.RData")
source("CLLproteomics-main/code/utils.R")
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, dev = c("png","pdf"))

#RNA-protein associations
#Preprocess transcriptomic and proteomic data

#Preprocessing RNA sequencing data
dds <- estimateSizeFactors(dds)
sampleOverlap <- intersect(colnames(protCLL), colnames(dds))
geneOverlap <- intersect(rowData(protCLL)$ensembl_gene_id, rownames(dds))
ddsSub <- dds[geneOverlap, sampleOverlap]
protSub <- protCLL[match(geneOverlap, rowData(protCLL)$ensembl_gene_id), sampleOverlap]

#how many gene don't have RNA expression at all?
noExp <- rowSums(counts(ddsSub)) == 0

#remove those genes in both datasets
ddsSub <- ddsSub[!noExp,]
protSub <- protSub[!noExp,]

protSub <- protSub[!duplicated(rowData(protSub)$name)]

geneOverlap <- intersect(rowData(protSub)$ensembl_gene_id, rownames(ddsSub))

ddsSub.vst <- varianceStabilizingTransformation(ddsSub)

#Calculate correlations between protein abundance and RNA expression
rnaMat <- assay(ddsSub.vst)
proMat <- assays(protSub)[["count_combat"]]
rownames(proMat) <- rowData(protSub)$ensembl_gene_id

corTab <- lapply(geneOverlap, function(n) {
  rna <- rnaMat[n,]
  pro.raw <- proMat[n,]
  res.raw <- cor.test(rna, pro.raw, use = "pairwise.complete.obs")
  tibble(id = n,
         p = res.raw$p.value,
         coef = res.raw$estimate)
}) %>% bind_rows() %>%
  arrange(desc(coef)) %>% mutate(p.adj = p.adjust(p, method = "BH"),
                                 symbol = rowData(dds[id,])$symbol,
                                 chr = rowData(dds[id,])$chromosome)

plotCorScatter <- function(inputTab, x, y, x_lab = "X", y_lab = "Y", title = "",
                           col = NULL, shape = NULL, showR2 = TRUE, annoPos = "right", legendPos = "left",
                           dotCol = colList, dotShape = c(16,1,17,2), textCol="darkred") {
  
  #prepare table for plotting
  plotTab <- tibble(x = inputTab[[x]],y=inputTab[[y]])
  if (!is.null(col)) plotTab <- mutate(plotTab, status = inputTab[[col]])
  if (!is.null(shape)) plotTab <- mutate(plotTab, statusShape = inputTab[[shape]])
  
  plotTab <- filter(plotTab, !is.na(x), !is.na(y))
  
  #prepare annotation values
  corRes <- cor.test(plotTab$x, plotTab$y)
  pval <- formatNum(corRes$p.value, digits = 1, format = "e")
  Rval <- formatNum(corRes$estimate, digits = 2)
  R2val <- formatNum(corRes$estimate^2, digits = 2)
  Nval <- nrow(plotTab)
  annoP <- bquote(italic("P")~"="~.(pval))
  
  if (showR2) {
    annoCoef <-  bquote(R^2~"="~.(R2val))
  } else {
    annoCoef <- bquote(R~"="~.(Rval))
  }
  annoN <- bquote(N~"="~.(Nval))
  
  corPlot <- ggplot(plotTab, aes(x = x, y = y))
  
  if (!is.null(col) & is.null(shape)) {
    corPlot <- corPlot + geom_point(aes(fill = status), shape =21, size =3) +
      scale_fill_manual(values = dotCol)
  } else if (is.null(col) & !is.null(shape)) {
    corPlot <- corPlot + geom_point(aes(shape = statusShape), color = dotCol[1], size=3) +
      scale_shape_manual(values = dotShape)
  } else if (!is.null(col) & !is.null(shape)) {
    corPlot <- corPlot + geom_point(aes(shape = statusShape, color = status), size=3) +
      scale_shape_manual(values = dotShape) +
      scale_color_manual(values = dotCol)
  }
  else {
    corPlot <- corPlot + geom_point(fill = dotCol[1], shape =21, size=3)
  }
  
  corPlot <- corPlot +   geom_smooth(formula = y~x,method = "lm", se=FALSE, color = "grey50", linetype ="dashed" )
  
  if (annoPos == "right") {
    
    corPlot <- corPlot + annotate("text", x = max(plotTab$x), y = Inf, label = annoN,
                                  hjust=1, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$x), y = Inf, label = annoP,
               hjust=1, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$x), y = Inf, label = annoCoef,
               hjust=1, vjust =6, size = 5, parse = FALSE, col= textCol)
    
  } else if (annoPos== "left") {
    corPlot <- corPlot + annotate("text", x = min(plotTab$x), y = Inf, label = annoN,
                                  hjust=0, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$x), y = Inf, label = annoP,
               hjust=0, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$x), y = Inf, label = annoCoef,
               hjust=0, vjust =6, size = 5, parse = FALSE, col= textCol)
  }
  corPlot <- corPlot + ylab(y_lab) + xlab(x_lab) + ggtitle(title) +
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1)) +
    theme_full + theme(legend.position = legendPos,
                       plot.margin = margin(13,13,13,13))
  corPlot
}

#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}
geneList <- c("BIN1","JUND","SYK")
plotList <- lapply(geneList, function(n) {
  geneId <- rownames(dds)[match(n, rowData(dds)$symbol)]
  stopifnot(length(geneId) ==1)
  plotTab <- tibble(x=rnaMat[geneId,],y=proMat[geneId,], IGHV=protSub$IGHV.status)
  coef <- cor(plotTab$x, plotTab$y, use="pairwise.complete")
  annoPos <- ifelse (coef > 0, "left","right")
  plotCorScatter(plotTab, "x","y", showR2 = FALSE, annoPos = annoPos, x_lab = "RNA expression", shape = "IGHV",
                 y_lab ="Protein expression", title = n,dotCol = colList[4], textCol = colList[1], legendPos="none")
})
goodCorPlot <- cowplot::plot_grid(plotlist = plotList, ncol =3)
goodCorPlot
