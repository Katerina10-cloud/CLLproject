#Enrichment analysis using clusterPofiler

BiocManager::install("msigdbr")
BiocManager::install("GenomicState")
BiocManager::install("reactome.db")
BiocManager::install("AnnotationDbi")

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

geneList2 <- as.vector(res1$log2FoldChange)
names(geneList1) <- res1$ENTREZID
ids <- bitr(names(geneList3), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
res3 <- merge(res3, ids, by = "SYMBOL")
res3 <- res3 %>% relocate(ENTREZID, .before = baseMean)

save(res3, file = "Kwok/Enrichment/res3.RData")

geneList3 <- as.vector(res3$log2FoldChange)
names(geneList3) <- res3$ENTREZID
# sort the list in decreasing order (required for clusterProfiler)
geneList3 = sort(geneList3, decreasing = TRUE)

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

#Enrichment analysis Nature cancer dataset MOFA2

load("CLLproject/MOFA2/corRes.sigF4.RData")

corRes.sig <- corRes.sig %>% relocate(symbol, .before = id)
corRes.sig  <- corRes.sig[,-2]
colnames(corRes.sig)[1] <- "SYMBOL"

geneList <- as.vector(corRes.sig$logFC)
names(geneList) <- corRes.sig$SYMBOL
ids <- bitr(names(geneList), fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
corRes <- merge(corRes.sig, ids, by = "SYMBOL")
corRes <- corRes %>% relocate(ENTREZID, .before = baseMean)

save(corRes, file = "corRes.sigF4.RData")

#GSEA GO

geneList <- as.vector(corRes$logFC)
names(geneList) <- corRes$ENTREZID
geneList = sort(geneList, decreasing = TRUE)

F4_gseGO <- gseGO(geneList = geneList, ont = "ALL", OrgDb = org.Hs.eg.db, pvalueCutoff = 0.05,
                    keyType = "ENTREZID", minGSSize = 3, maxGSSize = 500)
head(summary(F4_gseGO)[,-10])

enrNetwork <- setReadable(F4_gseGO, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=8, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick', vertex.label.font=6)+ ggtitle("GSEA GO theme network Factor4 MOFA2")
p1

#KEGG GSEA

kegg_geneList <- corRes$logFC
names(kegg_geneList) <- corRes$ENTREZID
kegg_geneList = sort(kegg_geneList, decreasing = TRUE)

kegg_organism = "hsa"
F4_gseKEGG <- gseKEGG(geneList = kegg_geneList, organism = kegg_organism,
                        minGSSize = 10, maxGSSize    = 350, pvalueCutoff = 0.05, pAdjustMethod = "none",
                        keyType = "ncbi-geneid")
head(F4_gseKEGG)

dotplot(F4_gseKEGG, showCategory=30, split=".sign", label_format = 50, font.size = 10, title = "GSEA KEGG Factor4 MOFA2") + facet_grid(.~.sign)     

#ORA GO
geneList <- as.vector(corRes$logFC)
names(geneList) <- corRes$ENTREZID
geneList = sort(geneList, decreasing = TRUE)
gene <- names(geneList)
gene <- names(geneList)[abs(geneList) > 1]

F4_oraGO <- enrichGO(gene = gene, ont = "ALL", OrgDb = org.Hs.eg.db, 
                       keyType = "ENTREZID", pvalueCutoff = 0.05, minGSSize = 10, maxGSSize = 350)
head(F4_oraGO)

#Removing redundancy of enriched GO terms
F4_oraGO <- simplify(F4_oraGO, cutoff = 0.7, by = "p.adjust", select_fun = min,
                       measure = "Wang", semData = NULL)

dotplot(F4_oraGO, showCategory=30, label_format = 50, font.size = 12, title = "ORA GO up-regulated Factor4 MOFA2")

#KEGG ORA
kegg_gene <- names(kegg_geneList)
kegg_gene <- names(geneList)[abs(geneList) < 1]
F4_oraKEGG <- enrichKEGG(gene = kegg_gene, organism = kegg_organism,
                           minGSSize = 10, maxGSSize    = 350, pvalueCutoff = 0.05, keyType = "ncbi-geneid")
head(F4_oraKEGG)

dotplot(F4_oraKEGG, showCategory=30, label_format = 50, font.size = 10, title = "KEGG ORA Factor4 MOFA2")


#MSigDb ORA Hallmark

m_hallmark <- msigdbr(species = "Homo sapiens", category = "H") 
gene <- names(geneList)

msig_H <- msigdbr(species = "Homo sapiens", category = "H") %>% 
  dplyr::select(gs_name, entrez_gene)

msigH_ora_F4 <- enricher(gene, TERM2GENE=msig_H)
head(msigH_ora_F4)

dotplot(msigH_F4, showCategory=10, label_format = 30, font.size = 10, title = "Hallmarks ORA Factor4")

enrNetwork <- setReadable(msigH_F4, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=4, cex_label_gene = 0.6, cex_label_category = 0.8, color_category='firebrick', node_label="category")+ ggtitle("Hallmarks theme network F4")
p1

#MSigDb GSEA
msigH_gsea_F4 <- GSEA(geneList, TERM2GENE = msig_H)
head(msigH_gsea_F4)

dotplot(msigH_gsea_F4, showCategory=10, label_format = 40, font.size = 10, title = "Hallmarks GSEA Factor4")

#Reactome pathway ORA

gene <- names(geneList)[abs(geneList) > 1]
gene <- names(geneList)

F4_ReactomeORA <- enrichPathway(gene=gene, pvalueCutoff=0.05, readable=T)
head(summary(F4_ReactomeORA))

enrNetwork <- setReadable(F4_ReactomeORA, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=3, cex_label_gene = 0.7, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("Reactome ORA F4 (logFC > 1)")
p1

dotplot(F4_ReactomeORA, showCategory=20, label_format = 50, font.size = 10, title = "Reactome ORA Factor4")

#Reactome pathway GSEA
F4_ReactomeGSEA <- gsePathway(geneList=geneList, pvalueCutoff=0.05, minGSSize = 10, maxGSSize = 350)
head(summary(F4_ReactomeGSEA))
enrNetwork <- setReadable(F4_ReactomeGSEA, 'org.Hs.eg.db', 'ENTREZID')
p1 <- cnetplot(enrNetwork, categorySize="pvalue", foldChange=geneList, showCategory=3, cex_label_gene = 0.7, cex_label_category = 0.8, color_category='firebrick')+ ggtitle("Reactome GSEA F4")
p1
