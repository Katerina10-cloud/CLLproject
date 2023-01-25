#Enrichment analysis using clusterPofiler

BiocManager::install("pathview")
BiocManager::install("GenomicState")
BiocManager::install("ReactomePA")
BiocManager::install("AnnotationDbi")

library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicState)
library(ReactomePA)


load("Kwok/res1.RData")
load("Kwok/res2.RData")
load("Kwok/res3.RData")


gene_list <- res1Sig_up
gene_list$Name <- row.names(gene_list)
gene_list <- gene_list %>% relocate(Name, .before = baseMean)
original_gene_list <- gene_list$log2FoldChange
names(original_gene_list) <- gene_list$Name
# sort the list in decreasing order (required for clusterProfiler)
original_gene_list = sort(original_gene_list, decreasing = TRUE)

#Prepare Input
#Sorting genes in decreasing order
res1Sig <- res1Sig[ order(res1Sig$log2FoldChange, decreasing = TRUE), ]

res2Sig <- res2Sig[ order(res2Sig$log2FoldChange, decreasing = TRUE), ]

res3Sig <- res3Sig[ order(res3Sig$log2FoldChange, decreasing = TRUE), ]

#res3Sig_down <- res3Sig_down[ order(res3Sig_down$log2FoldChange), ]


#Gene Annotations (Can be used to add e.g. ENTREZ ID, ENSEMBL ID, etc. to gene name)
res1$SYMBOL <- row.names(res1)
res2$SYMBOL <- row.names(res2)
res3$SYMBOL <- row.names(res3)

res3 <- res3 %>% relocate(SYMBOL, .before = baseMean)
row.names(res3) <- 1:nrow(res3)

res3 <- na.omit(res3)

columns(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
cols <- c("ENTREZID")

#deGenes <- res1$SYMBOL
#deGenes <- select(org.Hs.eg.db, keys=deGenes, columns=cols, keytype="SYMBOL")
#deGenes<- na.omit(deGenes)


geneList1 <- as.vector(res1$log2FoldChange)
names(geneList1) <- res1$SYMBOL
ids <- bitr(names(geneList3), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
res3 <- merge(res3, ids, by = "SYMBOL")
res3 <- res3 %>% relocate(ENTREZID, .before = baseMean)

save(res3, file = "Kwok/Enrichment/res3.RData")


geneList3 <- as.vector(res3$log2FoldChange)
names(geneList3) <- res3$ENTREZID
# sort the list in decreasing order (required for clusterProfiler)
geneList3 = sort(geneList3, decreasing = TRUE)

#GO enrichment analysis
res1_gseGO <- gseGO(geneList = geneList1, ont = "ALL", OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500)
res2_gseGO <- gseGO(geneList = geneList2, ont = "ALL", OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500)
res3_gseGO <- gseGO(geneList = geneList3, ont = "ALL", OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 500)
head(summary(res3_gseGO)[,-5])

save(res3_gseGO, file = "Kwok/Enrichment/res3_gseGO.RData")

