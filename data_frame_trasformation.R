library(dplyr)
library(tidyverse)
library(magrittr)
library(SummarizedExperiment)
library(devtools)
library(usethis)
library("rtracklayer")
load("CLLdata/Kwok_data/rc_rna_Lu.RData")
load("CLLdata/Kwok_data/rc_Kwok.Rdata")
load("CLLdata/Kwok_data/exp_Kwok.Rdata")
load("CLLproject/Survival_analysis/rna2.Rdata")
load("CLLdata/geneMut_LRP8_Mean.Rdata")
load("CLLdata/Kwok_data/rc_Kwok3.Rdata")

nc_Kwok = (read.csv("CLLdata/Kwok_data/Normalized_counts.csv", header=FALSE, sep=",")[-1,]) #skip first row
rc_Kwok <- rc_Kwok[,-2] #delete column
colnames(metadata_Kwok)[2] <- "CLLtype" #change column name
save(rc_Kwok, file = "rc_Kwok.Rdata")

metadata_Kwok <- read_excel("CLLdata/Kwok_data/metadata_Kwok.xlsx")

#Export to Exel
install.packages("openxlsx", dependencies=TRUE)
library("openxlsx")
write.xlsx(GO_list, file = "GO_Kwok_down.xlsx", colnames = TRUE)


metadata_Kwok <- metadata_Kwok[,-7 : -8] 
save(metadata_Kwok, file = "metadata_Kwok.Rdata")

exprMat1 <- exprMat1[,-212] #remove column

#moving the last column to the start
exprMat2 <- exprMat1 %>%
  select(hgnc_symbol, everything()) 

gene$ZMYM3 = as.factor(gene$ZMYM3) #convert numeric variable to factor (categorical)

reactomeGS3 <- genes.table_SYMBOLS

reactomeGS3 %>%
  relocate(hgnc_symbol, .before = ensembl_gene_id)

reactomeGS3<-reactomeGS3[,c(1305, 1:1304)]

exprMat3 <- exprMat2
exprMat3$hgnc_symbol[is.na(exprMat3$hgnc_symbol)] <- exprMat3$ensembl_gene_id[is.na(exprMat3$hgnc_symbol)]
                                                
exprMat2 <- exprMat2[,-2]

#replace NAs
library(dplyr)
df <- df %>% replace(is.na(.), 0)
#In a Specific Column, Replace NA Values
df %>% mutate(position = ifelse(is.na(position), 0, position))
#Replace NA Values in One Column by Another Column
exprMat3 <- exprMat2
exprMat3$hgnc_symbol[is.na(exprMat3$hgnc_symbol)] <- exprMat3$ensembl_gene_id[is.na(exprMat3$hgnc_symbol)]
#replace blancks by NA
exprMat3$hgnc_symbol[exprMat3$hgnc_symbol== ""] <- NA  
exprMat3$hgnc_symbol[is.na(exprMat3$hgnc_symbol)] <- exprMat3$ensembl_gene_id[is.na(exprMat3$hgnc_symbol)]

#remove blanks
rc_Kwok_3 <- rc_Kwok_3[!(rc_Kwok_3$hgnc_symbol == ""), ]

#rename colname
library(plyr)
rna1 <- rename(rna1,c('hgnc_symbol'='name')) 

geneTab <- filter(expMat1, name == "LPL") #filter one row

colnames(riskTab) <- riskTab[1, ] #convert first row to header

exprMat4 <- exprMat3
exprMat4 <- exprMat4[,-2] #remove 2 column

# reassigning row names
exprMat4 <- exprMat4 %>% remove_rownames %>% column_to_rownames(var="hgnc_symbol")

save(exprMat4, file = "exprMat4.Rdata")

#reactome dataframe trasformation
reactomeGS3 <- genes.table_SYMBOLS
reactomeGS3 %>%
 relocate(hgnc_symbol, .before = ensembl_gene_id) #move column

reactomeGS3<-reactomeGS3[,c(1305, 1:1304)]

reactomeGS3$hgnc_symbol[reactomeGS3$hgnc_symbol== ""] <- NA  
reactomeGS3$hgnc_symbol[is.na(reactomeGS3$hgnc_symbol)] <- reactomeGS3$ensembl_gene_id[is.na(reactomeGS3$hgnc_symbol)]
reactomeGS3 <- reactomeGS3[,-2]
reactomeGS3 <- reactomeGS3 %>% remove_rownames %>% column_to_rownames(var="hgnc_symbol")
reactomeGS3 <- t(reactomeGS3)

save(reactomeGS3, file = "reactomeGS3.Rdata")


#writing csv files
write.csv(gene, file="geneCSV.csv")
write.csv(exprMat, file="rnaCSV.csv")
write.csv(methData, file="methCSV.csv")

rna1 <- assay(rna)
rna1 <- as.data.frame(rna1)
write.csv(rna1, file="rna1CSV.csv")
save(rna1, file = "raw_counts_RNA.Rdata")

#survTab <- survival %>% 
select(patientID, OS, died, TTT, treatedAfter) %>%
  filter(patientID %in% colnames(expMat))
##Prepare genomic annotations
IGHVstatus = gene
IGHVstatus <- subset(IGHVstatus, select = c(IGHV))
IGHVstatus<- replace(IGHVstatus, is.na(IGHVstatus), 0)
IGHVstatus$patientID <- row.names(IGHVstatus)
IGHVstatus <- IGHVstatus %>%
  select(patientID, everything())
#table of known risks
riskTab1 <- select(survTab1, patientID) %>% 
  left_join(IGHVstatus[c("patientID", "IGHV")],  by = c(patientID = "patientID")) 
riskTab1$IGHV <- as.factor(riskTab1$IGHV)

#Extracting mutation's subgroups
survplotdata$Name <- row.names(survplotdata)
row.names(survplotdata) <- 1:nrow(survplotdata)
survplotdata <- survplotdata %>%
  relocate(Name, .before = TTT)
surv_data <- survplotdata
survplotdata <- subset(survplotdata, LRP8 %in% c('Mid'))
survplotdata <- survplotdata[,-2:-3]
geneMut_low <- merge(survplotdata, gene, by = "row.names")
geneMut_low <- geneMut_low[,-1]
geneMut_LRP8_Mean<- geneMut_LRP8_Mean %>% remove_rownames %>% column_to_rownames(var="Row.names")
geneMut_high <- gene_mut
geneMut_mid <- geneMut_mid %>% replace(is.na(.), 0)
geneMut_midMean <- sapply(geneMut_mid, mean) 
geneMut_midMean <- as.data.frame(geneMut_midMean)
geneMut_LRP8_Mean <- merge(geneMut_highMean, geneMut_lowMean, by = "row.names")
geneMut_LRP8_Mean <- merge(geneMut_LRP8_Mean, geneMut_midMean, by = "row.names")
save(geneMut_LRP8_Mean, file = "geneMut_LRP8_Mean.Rdata")
plotData <-geneMut_LRP8_Mean[,c('KRAS', 'NOTCH1', 'BRAF')]
geneMut_LRP8_Mean <- t(geneMut_LRP8_Mean)
plotData <- as.matrix(plotData)
rownames(plotData) <- c("high","low","mid")
barplot(prop.table(plotData) * 100, main = "Frequency of gene mutations (%)",
        col = rainbow(3), legend.text = rownames(plotData), args.legend = list(x = "topright", inset = c(0.80, 0)))

#Merge multiple dataframes
library(tidyverse)
#put all data frames into list
DEG_list <- list(F4_up, F6_dowm)

#merge all data frames in list
DEG_list %>% reduce(full_join, by='SYMBOL')
DEG_merged <- Reduce(function(x, y) merge(x, y, all=TRUE), DEG_list) 

#Remove dublicates by single column
DEG <- DEG_merged[!duplicated(DEG_merged$SYMBOL), ]
DEG <- DEG %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

DEG_NC <- merge(DEG, exprMat_NC, by = "row.names")
DEG_NC <- DEG_NC[,-2:-9] #delete column
DEG_NC <- DEG_NC %>% remove_rownames %>% column_to_rownames(var="Row.names")

DEG_NC <- t(scale(t(DEG_NC)))
save(DEG_NC, file = "DEG_NC.RData")

F4_DEG <- corRes
F6_DEG <- corRes

F4_up <- subset(F4_DEG, F4_DEG$logFC > 0)
F6_dowm <- subset(F6_DEG, F6_DEG$logFC < 0)

rc_NC <- rc_NC[,-2] #delete column
rc_NC <- rc_NC[!duplicated(rc_NC$GeneID), ]
rc_NC <- na.omit(rc_NC)
library(DESeq2)
rc_NC <- rc_NC %>% remove_rownames %>% column_to_rownames(var="GeneID")
rc_NC <- as.matrix(rc_NC)
rc_NC <- DESeq2::varianceStabilizingTransformation(rc_NC)
exprMat_NC <- rc_NC
save(exprMat_NC, file = "exprMat_NC.Rdata")
DEG_F6F4 <- merge(exprMat_NC, DEG_NC, by = "row.names")
DEG_F6F4 <- DEG_F6F4[,-212:-219]
DEG_F6F4 <- DEG_F6F4 %>% remove_rownames %>% column_to_rownames(var="Row.names")
DEG_F6F4 <- t(scale(t(DEG_F6F4)))
              
#Sorting a dataframe
library(tidyverse)
metadata_Kwok <- metadata_Kwok %>% arrange(CLLtype)
save(metadata_Kwok, file = "metadata_Kwok1.Rdata")

# reassigning row names
res3 <- res3 %>% remove_rownames %>% column_to_rownames(var="SYMBOL")
res3 <- res3[!duplicated(res3$SYMBOL), ] #remove duplicates

#Volcano plot
library(EnhancedVolcano)
library(ggrepel)
library(ggplot2)
library(ggfortify)

#Indolent vs Progressive

EnhancedVolcano(toptable = res3,            
                x = "log2FoldChange",           
                y = "padj",                     
                lab = rownames(res3),
                pCutoff = 0.01,
                FCcutoff = 1,
                title = "Regressive vs Progressive"
)

expMatrix_Kwok <- expMatrix_Kwok%>%
  relocate(SRR9140544, .before = SRR9140505)

#Heatmap of sample-to-sample distances
sampleDists <- dist(t(expMatrix_Kwok))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(metadata_Kwok$CLLtype)
colnames(sampleDistMatrix) <- paste(metadata_Kwok$CLLtype)

library("RColorBrewer")
library("pheatmap")
colors <- colorRampPalette( rev(brewer.pal(9, "PuRd")) )(255)
color = colorRampPalette(c(colList[2],"white",colList[1]))(100)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         cluster_cols = TRUE,
         cluster_rows = TRUE,
         show_rownames = TRUE, show_colnames = TRUE,
         main = "Heatmap of sample correlation")

#Hierarchical clustering

library("ComplexHeatmap")
library("RColorBrewer")
library("circlize")
library("grid")
library("tidyverse")


levels(riskTab1$IGHV)
levels(riskTab1$IGHV) <- c("U", "M")
                                

DEG_Kwok_3groups <- DEG_Kwok_3groups[,-3575:-3576]
DEG_Kwok_3groups <- as.matrix(DEG_Kwok_3groups)

gene_Kwok <- as.matrix(gene_Kwok)
ann_col_info <- as.data.frame(metadata_Kwok)
anno_info_colors = list(
  CLLtype = c(Regression = "blue", Indolent = "red", Progressive = "gray")
)

gene <- merge(gene, riskTab1, by = "patientID")
riskTab1 <- riskTab1 %>% remove_rownames %>% column_to_rownames(var="patientID")
gene <- gene[,-3] #delete column
gene <- as.data.frame(gene)
colnames(gene)[1] <- "ITGA4"

gene$patientID <- row.names(gene)
riskTab1$patientID <- row.names(riskTab1)
gene <- t(gene)

pheatmap(gene_Kwok, 
         cluster_rows = FALSE,                       
         cluster_cols = TRUE, 
         show_rownames = FALSE, 
         show_colnames = FALSE,
         scale = "row",
         color = colorRampPalette(c("navy", "white", "firebrick3"))(50),
         annotation_col = ann_col_info,
         annotation_colors = anno_info_colors,
         main = "Hierarchical clustering of Blood samples (ITGA4)")

#PCA
library(ggrepel)
library(ggfortify)
library(tidyverse)

DEG_Kwok_3groups <- t(DEG_Kwok_3groups)
DEG_Kwok_3groups <- as.data.frame(DEG_Kwok_3groups)
DEG_Kwok_3groups <- merge(DEG_Kwok_3groups, metadata_Kwok, by = "row.names")
DEG_Kwok_3groups <- DEG_Kwok_3groups %>% remove_rownames %>% column_to_rownames(var="Row.names")

singn_genes <- t(singn_genes)
singn_genes <- as.data.frame(singn_genes)
singn_genes <- merge(singn_genes, metadata_Kwok, by = "row.names")

gene_IGHV <- merge(gene, riskTab1, by = "row.names")
gene_IGHV <- gene_IGHV %>% remove_rownames %>% column_to_rownames(var="Row.names")
gene <- t(gene)
gene <- as.data.frame(gene)

pca <- prcomp(gene_IGHV[,c(1)], center = TRUE)

summary(pca)
autoplot(pca, data = gene, frame = TRUE, frame.type = 'norm', main = "PCA of ITGA4 subtypes")

#Filtering genes
sign_genes <- subset(all_genes_results, all_genes_results$padj < 0.05)
regressive <- regressive %>% remove_rownames %>% column_to_rownames(var="Row.names")
regressive <- merge(regressive, expMatrix_Kwok, by = "row.names")
regressive <- regressive[,-1:-7] #delete column

save(regressive, file = "regressive.Rdata")

res2_down <- subset(res2, res2$log2FoldChange < 0)

#Filtering GO results
#extract a dataframe with results from object of type Large gseaResult
gc_Regr_down_gseGO <- Regr_down_gseGO@result

#subset columns
gc_Regr_down_gseGO <- gc_Regr_down_gseGO[c("ID", "core_enrichment")]

#format core_enrichment column
gc_Regr_down_gseGO$core_enrichment <- gsub("/", ",", gc_Regr_down_gseGO$core_enrichment)
gc_Regr_down_gseGO$core_enrichment <- character(gc_Regr_down_gseGO$core_enrichment)

write.xlsx(gc_Ind_down_gseGO, file = "gc_Ind_down_gseGO.xlsx", colnames = TRUE)

GO_list <- list(GO_1, GO_2, GO_3, GO_4, GO_5, GO_6)

#merge all data frames in list
GO_list <- Reduce(function(x, y) merge(x, y, all=TRUE), GO_list) 

GO_list <- GO_list %>% remove_rownames %>% column_to_rownames(var="ID")

GO_list <- t(GO_list)
GO_list <- as.data.frame(GO_list)

#Dataframe rows as a list of vectors
GO_list <- as.list(GO_list)

GO.select <- GO_list
save(GO.select, file = "GO.select_down_Kwok.Rdata")

GO_6 <- enrichment[[6]]
GO_6 <- GO_6[7,]

save(GO_list, file = "GOenr_F4F6")
