#DESeq2 analysis

load("CLLdata/Kwok_NC_data/rc_Kwok.Rdata")
load("CLLdata/Kwok_NC_data/metadata_Kwok1.Rdata")

library(DESeq2)
library(dplyr)
library(tidyverse)
library(ggplot2)
library(ggrepel)
library(ggfortify)

#dataframe trasformations
rc_Kwok <- rc_Kwok %>% remove_rownames %>% column_to_rownames(var="Name")
rc_Kwok = rc_Kwok[,-1] #remove 2 column

#metadata_Kwok <- subset(metadata_Kwok, CLLtype %in% c('Regression','Progressive'))
#metadata_Kwok$CLLtype <- as.factor(metadata_Kwok$CLLtype)

### Check that sample names match in both files
all(colnames(rc_Kwok) %in% rownames(metadata_Kwok))
all(colnames(rc_Kwok) == rownames(metadata_Kwok))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = rc_Kwok, colData = metadata_Kwok, design = ~ CLLtype)

#Pre-filter the genes which have low counts
dds <- dds[rowSums(counts(dds)) >= 10,]


#Perform different expression analysis
dds <- DESeq2::DESeq(dds)

#Plotting dispersion estimates
DESeq2::plotDispEsts(dds, main="Dispersion Estimates")

#Using contrast
res1 <- results(dds, contrast=c("CLLtype","Indolent","Progressive"), pAdjustMethod = "BH", tidy = TRUE)
res2 <- results(dds, contrast=c("CLLtype","Regression","Indolent"), pAdjustMethod = "BH", tidy = TRUE)
res3 <- results(dds, contrast=c("CLLtype","Regression","Progressive"), pAdjustMethod = "BH", tidy = TRUE)


#Ordering results based on padj value
res1 <- res1 %>% remove_rownames %>% column_to_rownames(var="row")
res1 <- res1[order(res1$padj),]

res2 <- res2 %>% remove_rownames %>% column_to_rownames(var="row")
res2 <- res2[order(res2$padj),]

res3 <- res3 %>% remove_rownames %>% column_to_rownames(var="row")
res3 <- res3[order(res3$padj),]

sum(res1$padj < 0.01 & res1$log2FoldChange > 1, na.rm=TRUE)

#Getting DE genes

#up-regulated (Indolent vs Progressive)
sum(res1$padj < 0.01 & res1$log2FoldChange > 1, na.rm=TRUE)
#down-regulated (Indolent vs Progressive)
sum(res1$padj < 0.01 & res1$log2FoldChange < -1, na.rm=TRUE)

#up-regulated (Regression vs Indolent)
sum(res2$padj < 0.01 & res2$log2FoldChange > 1, na.rm=TRUE)
#down-regulated (Regression vs Indolent)
sum(res2$padj < 0.01 & res2$log2FoldChange < -1, na.rm=TRUE)

#up-regulated (Regression vs Progressive)
sum(res3$padj < 0.01 & res3$log2FoldChange > 1, na.rm=TRUE)
#down-regulated (Regression vs Progressive)
sum(res3$padj < 0.01 & res3$log2FoldChange < -1, na.rm=TRUE)

#Filtering DE genes (up- & down- regulated)

res1Sig_up <- subset(res1, res1$padj < 0.01 & res1$log2FoldChange > 1)
res1Sig_down <- subset(res1, res1$padj < 0.01 & res1$log2FoldChange < -1)

res2Sig_up <- subset(res2, res2$padj < 0.01 & res2$log2FoldChange > 1)
res2Sig_down <- subset(res2, res2$padj < 0.01 & res2$log2FoldChange < -1)

res3Sig_up <- subset(res3, res3$padj < 0.01 & res3$log2FoldChange > 1)
res3Sig_down <- subset(res3, res3$padj < 0.01 & res3$log2FoldChange < -1)

res3 <- subset(res3, res3$padj < 0.05)

#Datasets trasformation
# Create the new column
res1 <- res1 %>% mutate(sig=padj<0.01)
res2 <- res2 %>% mutate(sig=padj<0.01)
res3 <- res3 %>% mutate(sig=padj<0.01)

res1 = res1 %>% as.data.frame() %>% filter(padj < 0.01) %>% arrange(desc(log2FoldChange), desc(padj))
res2 = res2 %>% as.data.frame() %>% filter(padj < 0.01) %>% arrange(desc(log2FoldChange), desc(padj))
res3 = res3 %>% as.data.frame() %>% filter(padj < 0.01) %>% arrange(desc(log2FoldChange), desc(padj))

# distribution of adjusted p-values
hist(res1$padj, col="lightblue", main = "Adjusted p-value distribution")

#EXPLORATORY ANALYSIS AND VISUALIZATION
#Volcano plot
library(EnhancedVolcano)

#Indolent vs Progressive

EnhancedVolcano(toptable = res1,            
                x = "log2FoldChange",           
                y = "padj",                     
                lab = rownames(res1),
                pCutoff = 0.01,
                FCcutoff = 1,
                title = "Indolent vs Progressive"
)

#Regression vs Indolent
EnhancedVolcano(toptable = res2,            
                x = "log2FoldChange",           
                y = "padj",                     
                lab = rownames(res2),
                pCutoff = 0.01,
                FCcutoff = 1,
                title = "Regression vs Indolent"
)

#Regression vs Progressive
EnhancedVolcano(toptable = res3,            
                x = "log2FoldChange",           
                y = "padj",                     
                lab = rownames(res3),
                pCutoff = 0.01,
                FCcutoff = 1,
                title = "Regression vs Progressive"
)

#Volcano plot (standard)
#p = ggplot2::ggplot(res1, ggplot2::aes(log2FoldChange, -log10(pvalue))) +
  ggplot2::geom_point(ggplot2::aes(col = sig), size = 0.5) +
  ggplot2::scale_color_manual(values = c("red", "purple", "grey")) +
  labs(x="log2(fold change)",
       y="-log10 (adj.p-value)",
       title="Volcano Plot of DESeq2 analysis")  +
  theme_bw()+
  theme(plot.title = element_text(hjust = 0.5), 
        legend.position="right", 
        legend.title = element_blank()) 
  
#p + ggrepel::geom_text_repel(data=deseq2Res[1:10, ], ggplot2::aes(label=rownames(deseq2Res[1:10, ])))

# MA plot
deseq2Res$lfcSE <- NULL
deseq2Res$stat <- NULL
deseq2Res$pvalue <- NULL
deseq2Res$sig <- NULL
res2 <- as.data.frame(deseq2Res) %>% 
  mutate(padj = ifelse(padj <= 0.01, TRUE, FALSE))

DESeq2::plotMA(res2, main="MA Plot (regr_progr)", ylim=c(-5,5))
abline(h=c(-1,1),col = "blue")


#Generate normalized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

dds.vst<-DESeq2::varianceStabilizingTransformation(dds)
exp_Kwok <- SummarizedExperiment::assay(dds.vst)
exp_Kwok <- as.data.frame(exp_Kwok)

#Heatmap of the sample-to-sample distances

library("RColorBrewer")
library("pheatmap")
sampleDists <- dist(t(expMatrix_Kwok))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata_Kwok$CLLtype)
colnames(sampleDistMatrix) <- paste(metadata_Kwok$CLLtype)
colors <- colorRampPalette( rev(brewer.pal(9, "Blues")) )(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main = "Heatmap of sample distances")

#Selecting DE significat genes
DE_genes1 <- DE_genes
exp_Kwok1 <- exp_Kwok
DE_genes1 <- merge(exp_Kwok, DE_genes1, by = "row.names")
DE_genes1 <- DE_genes1 %>% remove_rownames %>% column_to_rownames(var="Row.names")
DE_genes1 <- merge(DE_genes1, metadata_Kwok, by = "row.names")
DE_genes1 <- as.data.frame(DE_genes1)

DE_genes <- merge(exp_Kwok, diff_genes, by = "row.names")
DE_genes = DE_genes %>% 
  as.data.frame() %>% 
  filter(padj < 0.01) %>% 
  arrange(desc(log2FoldChange), 
          desc(padj))

#Heatmap of the most significant genes

library("ComplexHeatmap")
library(RColorBrewer)
library(circlize)

#mat <- assay(exp_Kwok)[head(order(deseq2Res$padj),50), ]
#mat <- mat - rowMeans(mat)

metadata_Kwok <- as.matrix(metadata_Kwok)
metadata_Kwok <- SummarizedExperiment::assay(metadata_Kwok)
metadata_Kwok <- metadata_Kwok[,-3]

ann_col_info <- as.data.frame(metadata_Kwok)
anno_info_colors = list(
  CLLtype = c(Regression = "lightgrey", 
               Indolent = "black", Progressive = "red")
)

DE_genes1 <- as.matrix(DE_genes1)
pheatmap(DE_genes1, 
         cluster_rows = FALSE,                       
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         scale = "row",
         annotation_col = ann_col_info,
         annotation_colors = anno_info_colors,
         main = "Heatmap of different gene expression")

#PCA plot

#DE_genes1 <- DE_genes1%>% mutate_if(is.character, as.numeric)
DE_genes1 <- t(DE_genes1)
DE_genes1 <- as.data.frame(DE_genes1)
DE_genes1 <- merge(DE_genes1, metadata_Kwok, by = "row.names")
DE_genes1 <- DE_genes1 %>% remove_rownames %>% column_to_rownames(var="Row.names")

pca <- prcomp(DE_genes1[,c(1:836)])
summary(pca)

autoplot(pca, data = DE_genes1, colour = "CLLtype", main = "PCA of CLL subtypes")

