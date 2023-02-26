#Enrichment analysis using clusterPofiler

BiocManager::install("msigdbr")
BiocManager::install("qpcR")
BiocManager::install("reactome.db")
BiocManager::install("AnnotationDbi")
install("rowr")
library(tidyverse)
library(ggplot2)
library(clusterProfiler)
library(enrichplot)
library(pathview)
library(org.Hs.eg.db)
library(AnnotationDbi)
library(GenomicState)
library(ggupset)

#Kwok dataset

load("CLLdata/Kwok_data_DEG/Enrichment/Kwok/res1.RData")

#Preparing unput
gene_list <- res1Sig_up
gene_list$Name <- row.names(gene_list)
gene_list <- gene_list %>% relocate(Name, .before = baseMean)
original_gene_list <- gene_list$log2FoldChange
names(original_gene_list) <- gene_list$Name
# sort the list in decreasing order (required for clusterProfiler)
original_gene_list = sort(original_gene_list, decreasing = TRUE)

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
names(geneList1) <- res1$ENTREZID
ids <- bitr(names(geneList3), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
res3 <- merge(res3, ids, by = "SYMBOL")
res3 <- res3 %>% relocate(ENTREZID, .before = baseMean)

save(res3, file = "Kwok/Enrichment/res3.RData")

geneList3 <- as.vector(res3$log2FoldChange)
names(geneList3) <- res3$ENTREZID
# sort the list in decreasing order (required for clusterProfiler)
geneList2 = sort(geneList2, decreasing = TRUE)

#Merging DE gene lists
X1 <- names(geneList1)
X2 <- names(geneList2)
X3 <- names(geneList3)
library(qpcR)
gene_clusters <- qpcR:::cbind.na(X1, X2, X3)
gene_clusters <- as.data.frame(gene_clusters)
gene_clusters[is.na(gene_clusters)] <- ""

#create list of 3
gc <- list(X1=X1, X2=X2, X3=X3)
str(gc)


#GSEA GO enrichment analysis
res1_gseGO <- gseGO(geneList = geneList1, ont = "ALL", OrgDb = org.Hs.eg.db, 
                      keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)
res2_gseGO <- gseGO(geneList = geneList2, ont = "ALL", OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)
res3_gseGO <- gseGO(geneList = geneList3, ont = "ALL", OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)
head(summary(res3_gseGO)[,-10])

#Removing redundancy of enriched GO terms
res1_gseGO <- simplify(res1_gseGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
  measure = "Wang", semData = NULL)

res2_gseGO <- simplify(res2_gseGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

res3_gseGO <- simplify(res3_gseGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

#Visualization of Functional enrichment results
#Dotplot
dotplot(res1_gseGO, showCategory=30, split=".sign", label_format = 50, font.size = 8, title = "GSEA_GO Indolent vs Progressive CLL") + facet_grid(.~.sign)
dotplot(res2_gseGO, showCategory=25, split=".sign", label_format = 60, font.size = 8, title = "GSEA_GO Regressive vs Indolent CLL") + facet_grid(.~.sign)
dotplot(res3_gseGO, showCategory=25, split=".sign", label_format = 70, font.size = 8, title = "GSEA_GO Regressive vs Progressive CLL") + facet_grid(.~.sign)


#Encrichment Map
res1_gseGO <- pairwise_termsim(res1_gseGO)
emapplot(res1_gseGO, lshowCategory = 15, ayout="kk")+ ggtitle("Enrichment map for GSEA_GO Indolent vs Progressive CLL")

res2_gseGO <- pairwise_termsim(res2_gseGO)
emapplot(res2_gseGO, showCategory = 15, cex_category=1.5,layout="kk")+ ggtitle("Enrichment map for GSEA_GO Regression vs Indolent CLL")

res3_gseGO <- pairwise_termsim(res3_gseGO)
emapplot(res3_gseGO, showCategory = 15, cex_category=1.5, layout="kk")+ ggtitle("Enrichment map for GSEA_GO Regression vs Progressive CLL")

# convert gene ID to Symbol
enrNetwork <- setReadable(res2_gseGO, 'org.Hs.eg.db', 'ENTREZID')
#Gene concept network
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList2, showCategory = 3)+ ggtitle("Cnetplot for GSEA_GO Regression vs Indolent CLL")
p1

#p2 <- cnetplot(enrNetwork, foldChange=geneList1, showCategory = 2, circular = TRUE, colorEdge = TRUE) 

#Heatmap-like functional classification
p1 <- heatplot(res1_gseGO, foldChange=geneList1, showCategory=5)

#Ridgeplot
p1 <- ridgeplot(res2_gseGO, showCategory=15) + labs(x = "enrichment distribution")+ ggtitle("Density plot for GSEA_GO Regression vs Indolent CLL")
p1

#Tree plot
res1 <- pairwise_termsim(res1_gseKEGG)
p1 <- treeplot(res1)+ ggtitle(" Heatmap plot of enriched terms gsea_KEGG Indolent vs Progressive CLL")
p1

res2 <- pairwise_termsim(res2_gseKEGG)
p1 <- treeplot(res2)+ ggtitle(" Heatmap plot of enriched terms gsea_KEGG Regressive vs Indolent CLL")
p1

res3 <- pairwise_termsim(res3_gseKEGG)
p1 <- treeplot(res3)+ ggtitle(" Heatmap plot of enriched terms gsea_KEGG Regressive vs Progressive CLL")
p1

#KEGG Gene Set Enrichment Analysis

# Create a vector of the gene universe
kegg_geneList3 <- res3$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_geneList3) <- res3$ENTREZID
# sort the list in decreasing order (required for clusterProfiler)
kegg_geneList3 = sort(kegg_geneList3, decreasing = TRUE)

#Create gseKEGG object
kegg_organism = "hsa"
res1_gseKEGG <- gseKEGG(geneList = kegg_geneList1, organism = kegg_organism,
            minGSSize = 3, maxGSSize    = 350, pvalueCutoff = 0.05, pAdjustMethod = "none",
               keyType = "ncbi-geneid")
res2_gseKEGG <- gseKEGG(geneList = kegg_geneList2, organism = kegg_organism,
                        minGSSize = 3, maxGSSize = 350, pvalueCutoff = 0.05, pAdjustMethod = "none",
                        keyType = "ncbi-geneid")
res3_gseKEGG <- gseKEGG(geneList = kegg_geneList3, organism = kegg_organism,
                        minGSSize = 3, maxGSSize    = 350, pvalueCutoff = 0.05, pAdjustMethod = "none",
                        keyType = "ncbi-geneid")

save(res3_gseKEGG, file = "Kwok/Enrichment/res3_gseKEGG.RData")

#Visualization of KEGG enrichment analysis

dotplot(res1_gseKEGG, showCategory = 25, label_format = 65, font.size = 10, title = "KEGG enriched pathways Indolent vs Progressive CLL" , split=".sign") + facet_grid(.~.sign)
dotplot(res2_gseKEGG, showCategory = 30, label_format = 65, font.size = 8, title = "KEGG enriched pathways Regressive vs Indolent CLL" , split=".sign") + facet_grid(.~.sign)
dotplot(res3_gseKEGG, showCategory = 30, label_format = 65, font.size = 8, title = "KEGG enriched pathways Regressive vs Progressive CLL" , split=".sign") + facet_grid(.~.sign)

#Over-representation analysis

#Preparing unput
geneList1 <- as.vector(res1$log2FoldChange)
names(geneList1) <- res1$ENTREZID
geneList1 = sort(geneList1, decreasing = TRUE)
gene1 <- names(geneList1)

sig_genes = subset(res1, padj < 0.05)
genes1 <- sig_genes$log2FoldChange
names(genes1) <- sig_genes$ENTREZID
# filter on min log2fold change (log2FoldChange > 1)
genes1 <- names(genes1)[abs(genes1)  > 1]


#Create the object GO
res1_oraGO <- enrichGO(gene = gene1, ont = "ALL", OrgDb = org.Hs.eg.db, 
                    keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 3, maxGSSize = 350)
res2_oraGO <- enrichGO(gene = gene2, ont = "ALL", OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)
res3_oraGO <- enrichGO(gene = gene3, ont = "ALL", OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)

#Removing redundancy of enriched GO terms
res1_oraGO <- simplify(res1_oraGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

res2_oraGO <- simplify(res2_oraGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

res3_oraGO <- simplify(res3_oraGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

head(summary(res1_oraGO))

save(res3_oraGO, file = "Kwok/Enrichment/res3_oraGO.RData")

#Visualization of Functional enrichment results
#Dotplot
dotplot(res1_oraGO, showCategory=30, label_format = 60, font.size = 12, title = "ORA_GO Indolent vs Progressive CLL") 
dotplot(res2_oraGO, showCategory=30, label_format = 60, font.size = 12, title = "ORA_GO Regressive vs Indolent CLL") 
dotplot(res3_oraGO, showCategory=30, label_format = 60, font.size = 8, title = "ORA_GO Regressive vs Progressive CLL") 

# convert gene ID to Symbol
enrNetwork1 <- setReadable(res1_oraGO, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork1, categorySize="pvalue", foldChange=geneList1, cex_label_gene = 0.6, cex_label_category = 1, color_category='firebrick')+ ggtitle("ORA_GO Indolent vs Progressive")
p1

#Upsetplot
upsetplot(res1_oraGO)+ ggtitle("Upsetplot for ORA_GO Indolent vs Progressive CLL")
upsetplot(res2_oraGO)+ ggtitle("Upsetplot for ORA_GO Regressive vs Indolent CLL")
upsetplot(res3_oraGO)+ ggtitle("Upsetplot for ORA_GO Regressive vs Progressive CLL")

#KEGG ORA analysis

#Preparing input
# Create a vector of the gene unuiverse
kegg_gene_list3 <- res3$log2FoldChange
# Name vector with ENTREZ ids
names(kegg_gene_list3) <- res3$ENTREZID
# sort the list in decreasing order (required for clusterProfiler)
kegg_gene_list3 = sort(kegg_gene_list3, decreasing = TRUE)
# Exctract significant results from df2
kegg_sig_genes3 = subset(res3, padj < 0.05)
# From significant results, we want to filter on log2fold change
kegg_genes3 <- kegg_sig_genes3$log2FoldChange
# Name the vector with the CONVERTED ID!
names(kegg_genes3) <- kegg_sig_genes3$ENTREZID
# filter on log2fold change (PARAMETER)
kegg_genes3 <- names(kegg_genes3)[abs(kegg_genes3) > 1]

kegg_gene3 <- names(kegg_gene_list3)

#Create enrichKEGG object
kegg_organism = "hsa"
res1_oraKEGG <- enrichKEGG(gene = kegg_gene1, organism = kegg_organism,
                        minGSSize = 3, maxGSSize    = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
res2_oraKEGG <- enrichKEGG(gene = kegg_gene2, organism = kegg_organism,
                           minGSSize = 3, maxGSSize    = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
res3_oraKEGG <- enrichKEGG(gene = kegg_gene3, organism = kegg_organism,
                           minGSSize = 3, maxGSSize    = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

save(res3_oraKEGG, file = "res3_oraKEGG.RData")

head(summary(res1_oraKEGG))

#Dotplot
dotplot(res1_oraKEGG, showCategory=20, font.size = 10, title = "KEGG ORA Indolent vs Progressive CLL")
dotplot(res2_oraKEGG, showCategory=30, font.size = 10, title = "KEGG ORA Regressive vs Indolent CLL")
dotplot(res3_oraKEGG, showCategory=30, font.size = 12, title = "KEGG ORA Regressive vs Progressive CLL")

#MSigDb gene set enrichment analysis
library(msigdbr)
library(fgsea)

msigdbr_species()

m_hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
head(m_hallmark)

msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

#MSigDb over-presentaton analysis
#Preparing input
geneList3 <- as.vector(res3$log2FoldChange)
names(geneList3) <- res3$ENTREZID
geneList3 = sort(geneList3, decreasing = TRUE)
gene3<- names(geneList3)

msig_H1 <- enricher(gene1, TERM2GENE=msig_H)
head(msig_H1)

msig_H2 <- enricher(gene2, TERM2GENE=msig_H)
head(msig_H2)

msig_H3 <- enricher(gene3, TERM2GENE=msig_H)
head(msig_H3)

#Visualization
# convert gene ID to Symbol
enrNetwork1 <- setReadable(msig_H1, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork1, categorySize="pvalue", foldChange=geneList1, cex_label_gene = 0.6, cex_label_category = 0.7, color_category='firebrick')+ ggtitle("Over-expression analysis Hallmark Indolent vs Progressive")
p1

enrNetwork3 <- setReadable(msig_H3, 'org.Hs.eg.db', 'ENTREZID')
p2 <- cnetplot(enrNetwork3, categorySize="pvalue", foldChange=geneList3, cex_label_gene = 0.6, cex_label_category = 0.7, color_category='firebrick') + ggtitle("Over-expression analysis Hallmark Regressive vs Progressive")
p2

#MSigDb gene set enrichment analysis
msig_H1_gsea <- GSEA(geneList1, TERM2GENE = msig_H)
head(msig_H1_gsea)

msig_H2_gsea <- GSEA(geneList2, TERM2GENE = msig_H)
head(msig_H2_gsea)

msig_H3_gsea <- GSEA(geneList3, TERM2GENE = msig_H)
head(msig_H3_gsea)

#Visualization
# convert gene ID to Symbol
enrNetwork1 <- setReadable(msig_H1_gsea, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork1, categorySize="pvalue", foldChange=geneList1, cex_label_gene = 0.7, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("GSEA Hallmark Indolent vs Progressive")
p1

enrNetwork2 <- setReadable(msig_H2_gsea, 'org.Hs.eg.db', 'ENTREZID')
p2 <- cnetplot(enrNetwork2, categorySize="pvalue", foldChange=geneList2, showCategory=3, cex_label_gene = 0.6, cex_label_category = 0.7, color_category='firebrick')+ ggtitle("GSEA Hallmark Regressive vs Indolent")
p2

dotplot(msig_H2_gsea, title = "GSEA Hallmark Regressive vs Indolent CLL")


##Reactome enrichment analysis##

library(ReactomePA)

#Reactome pathway over-representation analysis
#Preparing unput
geneList3 <- as.vector(res3$log2FoldChange)
names(geneList3) <- res3$ENTREZID
geneList3 = sort(geneList3, decreasing = TRUE)
gene3 <- names(geneList3)
gene3 <- names(geneList3)[abs(geneList3) < 1]

res3_ReactomeORA <- enrichPathway(gene=gene3,pvalueCutoff=0.05, readable=T)
head(summary(res3_ReactomeORA))

#Visualization
enrNetwork1 <- setReadable(res2_ReactomeORA, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork1, categorySize="pvalue", foldChange=geneList2, showCategory=10, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("Reactome ORA theme network Regressive vs Indolent")
p1

dotplot(res3_ReactomeORA, showCategory=20, font.size = 10, label_format = 50, title = "Reactome ORA Regressive vs Progressive CLL")

##Reactome pathway GSEA
#Preparing input
geneList3 <- as.vector(res3$log2FoldChange)
names(geneList3) <- res3$ENTREZID
geneList3= sort(geneList3, decreasing = TRUE)

res3_ReactomeGSEA <- gsePathway(geneList=geneList3, pvalueCutoff=0.05, minGSSize = 10, maxGSSize = 350)
head(summary(res3_ReactomeGSEA))

#Visualization
dotplot(res3_ReactomeGSEA, showCategory=20, font.size = 10, label_format = 50, title = "Reactome GSEA Regressive vs Progressive CLL")

enrNetwork1 <- setReadable(res2_ReactomeGSEA, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork1, categorySize="pvalue", foldChange=geneList2, showCategory=8, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("Reactome GSEA theme network Regressive vs Indolent")
p1


#####Enrichment analysis Nature cancer dataset MOFA2#####

load("CLLproject/MOFA2/corRes.sigF4.RData")

corRes.sig <- corRes.sig %>% relocate(symbol, .before = id)
corRes.sig  <- corRes.sig[,-2]
colnames(corRes.sig)[1] <- "SYMBOL"

geneList <- as.vector(corRes.sig$logFC)
names(geneList) <- corRes.sig$SYMBOL
ids <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
corRes <- merge(corRes.sig, ids, by = "SYMBOL")
corRes <- corRes %>% relocate(ENTREZID, .before = baseMean)

save(corRes, file = "corRes.sigF1.RData")

#GSEA GO

geneList <- as.vector(corRes$logFC)
names(geneList) <- corRes$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

F1_gseGO <- gseGO(geneList = geneList, ont = "ALL", OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
                    keyType = "ENTREZID", minGSSize = 10, maxGSSize = 350)
head(F1_gseGO)

dotplot(F1_gseGO, showCategory=30, split=".sign", label_format = 30, font.size = 14, title = "GSEA GO F1") + facet_grid(.~.sign)    

enrNetwork <- setReadable(F4_gseGO, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=8, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick', vertex.label.font=6)+ ggtitle("GSEA GO theme network Factor4 MOFA2")
p1

#KEGG GSEA

kegg_geneList <- corRes$logFC
names(kegg_geneList) <- corRes$ENTREZID
kegg_geneList = sort(kegg_geneList, decreasing = TRUE)

kegg_organism = "hsa"
F1_gseKEGG <- gseKEGG(geneList = kegg_geneList, organism = kegg_organism,
                        minGSSize = 10, maxGSSize    = 350, pvalueCutoff = 0.05, pAdjustMethod = "none",
                        keyType = "ncbi-geneid")
head(F1_gseKEGG)

dotplot(F1_gseKEGG, showCategory=25, split=".sign", label_format = 60, font.size = 8, title = "GSEA KEGG F1") + facet_grid(.~.sign)     

#ORA GO
geneList <- as.vector(corRes$logFC)
names(geneList) <- corRes$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
gene <- names(geneList)
gene <- names(geneList)[abs(geneList) > 1]

F1_oraGO <- enrichGO(gene = gene, ont = "ALL", OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)
head(F1_oraGO)

#Removing redundancy of enriched GO terms
F1_oraGO <- simplify(F1_oraGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

dotplot(F1_oraGO, showCategory=25, label_format = 40, font.size = 14, title = "ORA GO F1 (logFC > 1)")


#KEGG ORA
kegg_gene <- names(kegg_geneList)
kegg_gene <- names(geneList)[abs(geneList) < 1]
F1_oraKEGG <- enrichKEGG(gene = kegg_gene, organism = kegg_organism,
                           minGSSize = 10, maxGSSize    = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
head(F1_oraKEGG)

dotplot(F1_oraKEGG, showCategory=25, label_format = 60, font.size = 10, title = "KEGG ORA F1")


#MSigDb ORA Hallmark

m_hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
gene <- names(geneList)

msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

msigH_ora_F1 <- enricher(gene, TERM2GENE=msig_H)
head(msigH_ora_F1)

dotplot(msigH_F4, showCategory=10, label_format = 30, font.size = 10, title = "Hallmarks ORA Factor4")

enrNetwork <- setReadable(msigH_F4, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=4, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick', node_label="category")+ ggtitle("Hallmarks theme network F4")
p1

#MSigDb GSEA
msigH_gsea_F1 <- GSEA(geneList, TERM2GENE = msig_H)
head(msigH_gsea_F1)

dotplot(msigH_gsea_F1, showCategory=10, label_format = 40, font.size = 12, title = "Hallmarks GSEA F1")

enrNetwork <- setReadable(msigH_gsea_F1, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=3, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick', node_label="category")+ ggtitle("Hallmarks GSEA theme network F1")
p1

#Reactome pathway ORA

gene <- names(geneList)[abs(geneList) < 1]
gene <- names(geneList)

F1_ReactomeORA <- enrichPathway(gene=gene, pvalueCutoff=0.05, readable=T)
head(F1_ReactomeORA)

enrNetwork <- setReadable(F4_ReactomeORA, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=3, cex_label_gene = 0.7, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("Reactome ORA F4 (logFC > 1)")
p1

dotplot(F1_ReactomeORA, showCategory=25, label_format = 60, font.size = 8, title = "Reactome ORA F1")

#Reactome pathway GSEA
F1_ReactomeGSEA <- gsePathway(geneList=geneList, pvalueCutoff=0.05, minGSSize = 3, maxGSSize = 350)
head(F1_ReactomeGSEA)

enrNetwork <- setReadable(F4_ReactomeGSEA, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=3, cex_label_gene = 0.7, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("Reactome GSEA F4")
p1

dotplot(F1_ReactomeGSEA, showCategory=25, label_format = 60, font.size = 8, title = "Reactome GSEA F1")


#Biological theme comparison Kwok DE genes groups
#create list of 3
X1 <- names(geneList1)
X2 <- names(geneList2)
X3 <- names(geneList3)
gc <- list(X1=X1, X2=X2, X3=X3)
str(gc)

#Over-expression gene clusters comparison
#enrichGO
gc_oraGO <- compareCluster(geneCluster = gc, fun = enrichGO, OrgDb = org.Hs.eg.db, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
gc_oraGO <- setReadable(gc_oraGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(gc_oraGO) 
dotplot(gc_oraGO, showCategory = NULL, label_format = 60, font.size = 12, title = "Biological theme comparison ORA GO")
cnetplot(gc_oraGO, showCategory = NULL, categorySize="pvalue", cex_label_category = 1.2, node_label="category", )+ ggtitle("Gene-concept network ORA GO")

#enrichKEGG
kegg_organism = "hsa"
gc_oraKEGG <- compareCluster(geneCluster = gc, fun = enrichKEGG, organism = kegg_organism,  minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
dotplot(gc_oraKEGG, showCategory = NULL, label_format = 60, font.size = 10, title = "Biological theme comparison ORA KEGG")

#enrichPathway
gc_oraReactomePA <- compareCluster(geneCluster = gc, fun = enrichPathway, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
dotplot(gc_oraReactomePA, showCategory = 25, label_format = 80, font.size = 8, title = "Biological theme comparison ORA ReactomePA")
cnetplot(gc_oraReactomePA, showCategory = 20, categorySize="pvalue", cex_label_category = 1.2, node_label="category", )+ ggtitle("Gene-concept network ORA ReactomePA")

#GSEA gene clusters comparison
#create list of 3
X1 <- geneList1
X2 <- geneList2
X3 <- geneList3
gc_gsea <- list(X1=X1, X2=X2, X3=X3)
str(gc_gsea)

#gseGO
gc_gseaGO <- compareCluster(geneCluster = gc_gsea, fun = gseGO, OrgDb = org.Hs.eg.db, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
gc_gseaGO <- simplify(gc_gseaGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                      measure = "Wang", semData = NULL)

#gseReactomePA
gc_gseReactomePA <- compareCluster(geneCluster = gc_gsea, fun = gsePathway, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
dotplot(gc_gseReactomePA, showCategory = 30, label_format = 80, font.size = 8, title = "Biological theme comparison GSEA ReactomePA")

#gseKEGG
gc_gseKEGG <- compareCluster(geneCluster = gc_gsea, fun = gseKEGG, organism = kegg_organism,  minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
dotplot(gc_gseKEGG, showCategory = NULL, label_format = 70, font.size = 10, title = "Biological theme comparison GSEA KEGG")


##Biological theme comparison Kwok DEG groups & Nature cancer DEG groups associated with F1, F4, F6
#Preparing input
geneList_F1 <- as.vector(corRes$logFC)
names(geneList_F1) <- corRes$ENTREZID
geneList_F1 = sort(geneList_F1, decreasing = TRUE)
F1 <- names(geneList_F1)

#create list of 6
X1 <- names(geneList1)
X2 <- names(geneList2)
X3 <- names(geneList3)
gc <- list(X1=X1, X2=X2, X3=X3, F1=F1, F4=F4, F6=F6)
str(gc)

#Over-expression gene clusters comparison
#enrichGO
gc_oraGO <- compareCluster(geneCluster = gc, fun = enrichGO, OrgDb = org.Hs.eg.db, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
gc_oraGO <- setReadable(gc_oraGO, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
head(gc_oraGO) 
gc_oraGO <- simplify(gc_oraGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)
dotplot(gc_oraGO, showCategory = NULL, label_format = 80, font.size = 8, title = "Biological theme comparison ORA GO")

#enrichKEGG
kegg_organism = "hsa"
gc_oraKEGG <- compareCluster(geneCluster = gc, fun = enrichKEGG, organism = kegg_organism,  minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")

dotplot(gc_oraKEGG, showCategory = 20, label_format = 80, font.size = 7, title = "Biological theme comparison ORA KEGG")

#enrichPathway
gc_oraReactomePA <- compareCluster(geneCluster = gc, fun = enrichPathway, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
dotplot(gc_oraReactomePA, showCategory = 20, label_format = 85, font.size = 7, title = "Biological theme comparison ORA ReactomePA")
cnetplot(gc_oraReactomePA, showCategory = 20, categorySize="pvalue", cex_label_category = 1.2, node_label="category", )+ ggtitle("Gene-concept network ORA ReactomePA")

#GSEA gene clusters comparison
#create list of 6
X1 <- geneList1
X2 <- geneList2
X3 <- geneList3
F1 <- geneList_F1
F4 <- geneList_F4
F6 <- geneList_F6
gc_gsea <- list(X2=X2, X3=X3, F1=F1)
str(gc_gsea)

#gseGO
gc_gseaGO <- compareCluster(geneCluster = gc_gsea, fun = gseGO, OrgDb = org.Hs.eg.db, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
gc_gseaGO <- simplify(gc_gseaGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                      measure = "Wang", semData = NULL)
dotplot(gc_gseaGO, showCategory = 20, label_format = 85, font.size = 7, title = "Biological theme comparison GSEA GO")

#gseReactomePA
gc_gseReactomePA <- compareCluster(geneCluster = gc_gsea, fun = gsePathway, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
dotplot(gc_gseReactomePA, showCategory = 30, label_format = 85, font.size = 7, title = "Biological theme comparison GSEA ReactomePA")
head(gc_gseReactomePA)

regr2 <- data.frame(Entrez=names(geneList2), FC=geneList2)
regr2 <- regr2[abs(regr2$FC) > 1,]
regr2$othergroup[regr2$FC > 1] <- "upregulated"
regr2$othergroup[regr2$FC < -1] <- "downregulated"
x2_gseReactomePA <- compareCluster(Entrez~othergroup, data=regr2, fun="gseKEGG", minGSSize = 10, maxGSSize = 350)

#gseKEGG
kegg_organism = "hsa"
gc_gseKEGG <- compareCluster(geneCluster = gc_gsea, fun = gseKEGG, organism = kegg_organism,  minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
dotplot(gc_gseKEGG, showCategory = NULL, label_format = 80, font.size = 10, title = "Biological theme comparison GSEA KEGG")

###Functional profile comparison###
#Filtering DE genes (up- & down- regulated)
res2Sig_up <- subset(res2, res2$padj < 0.05 & res2$log2FoldChange > 1)
res2Sig_down <- subset(res2, res2$padj < 0.05 & res2$log2FoldChange < 0)

res3Sig_up <- subset(res3, res3$padj < 0.05 & res3$log2FoldChange > 1)
res3Sig_down <- subset(res3, res3$padj < 0.05 & res3$log2FoldChange < 0)

corResF1_up <- subset(corRes, corRes$padj < 0.05 & corRes$logFC > 0)
corResF1_down <- subset(corRes, corRes$padj < 0.05 & corRes$logFC < 0)

geneList3_up<- as.vector(res3Sig_up$log2FoldChange)
names(geneList3_up) <- res3Sig_up$ENTREZID
geneList3_up <- sort(geneList3_up, decreasing = TRUE)

geneList3_down<- as.vector(res3Sig_down$log2FoldChange)
names(geneList3_down) <- res3Sig_down$ENTREZID
geneList3_down <- sort(geneList3_down, decreasing = TRUE)

geneListF1_up <- as.vector(corResF1_up$logFC)
names(geneListF1_up) <- corResF1_up$ENTREZID
geneListF1_up <- sort(geneListF1_up, decreasing = TRUE)

geneListF1_down <- as.vector(corResF1_down$logFC)
names(geneListF1_down) <- corResF1_down$ENTREZID
geneListF1_down <- sort(geneListF1_down, decreasing = TRUE)

#create list of 6
X2_up <- names(geneList2_up)
X2_down <- names(geneList2_down)
X3_up <- names(geneList3_up)
X3_down <- names(geneList3_down)
F1_up <- names(geneListF1_up)
F1_down <- names(geneListF1_down)

gc <- list(X2_up=X2_up, X2_down=X2_down, X3_up=X3_up, X3_down=X3_down, F1_up=F1_up, F1_down=F1_down)

#ReactomePA ORA
gc_oraReactomePA <- compareCluster(geneCluster = gc, fun = enrichPathway, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
dotplot(gc_oraReactomePA, showCategory = 20, label_format = 90, font.size = 7, title = "Biological theme comparison ORA ReactomePA")

#ReactomePA GSEA
X2_up <- geneList2_up
X2_down <- geneList2_down
X3_up <- geneList3_up
X3_down <- geneList3_down
F1_up <- geneListF1_up
F1_down <- geneListF1_down

gc_gsea <- list(X2_up=X2_up, X2_down=X2_down, X3_up=X3_up, X3_down=X3_down, F1_up=F1_up, F1_down=F1_down)
gc_gseReactomePA <- compareCluster(geneCluster = gc_gsea, fun = gsePathway, minGSSize = 10, maxGSSize = 350, pvalueCutoff = 0.05)
dotplot(gc_gseReactomePA, showCategory = 30, label_format = 90, font.size = 7, title = "Biological theme comparison GSEA ReactomePA")


dotplot(gseReactomePA_NC_Kwok, showCategory = 30, label_format = 80, font.size = 7, title = "Comparison GSE ReactomePA (NC & Kwok)")
dotplot(gc_oraReactomePA, showCategory = 30, label_format = 90, font.size = 7, title = "Comparison ORA ReactomePA (NC & Kwok)")
dotplot(gc_oraGO, showCategory = 30, label_format = 90, font.size = 8, title = "Comparison ORA GO (NC & Kwok)")

gc_oraGO <- dropGO(gc_oraGO, level = NULL, term = NULL)
gc_oraGO <- simplify(gc_oraGO,cutoff = 0.7,by = "p.adjust",select_fun = min, measure = "Wang",semData = NULL)

