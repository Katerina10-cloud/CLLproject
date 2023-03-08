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
colnames(corRes)[8] <- "padj" #change column name
save(rc_Kwok, file = "rc_Kwok.Rdata")

metadata_Kwok <- read_excel("CLLdata/Kwok_data/metadata_Kwok.xlsx")

#Export to Exel
install.packages("openxlsx", dependencies=TRUE)
library("openxlsx")
write.xlsx(res3, file = "res3.sig_Regr_vs_Progr.xlsx", colnames = TRUE)


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
DEG_list <- list(res1, res2, res3)

#merge all data frames in list
DEG_list %>% reduce(full_join, by='SYMBOL')
DEG_Kwok_merged <- Reduce(function(x, y) merge(x, y, all=TRUE), DEG_list) 

#Remove dublicates by single column
DEG_Kwok <- DEG_Kwok_merged[!duplicated(DEG_Kwok_merged$SYMBOL), ]
DEG_Kwok1 <- DEG_Kwok %>% remove_rownames %>% column_to_rownames(var="SYMBOL")

DEG_3groups_Kwok <- merge(DEG_Kwok1, expMatrix_Kwok, by = "row.names")
DEG_3groups_Kwok <- DEG_3groups_Kwok[,-2:-8] #delete column
DEG_3groups_Kwok <- DEG_3groups_Kwok %>% remove_rownames %>% column_to_rownames(var="Row.names")

DEG_Kwok_3groups <- t(scale(t(DEG_3groups_Kwok)))
save(DEG_Kwok_3groups, file = "DEG_Kwok_3groups.RData")


