load("CLLdata/Kwok_data/rc_Kwok.Rdata")
load("CLLdata/Kwok_data/metadata_Kwok.Rdata")

library(DESeq2)
library(dplyr)
library(tidyverse)

metadata_Kwok <- metadata_Kwok %>% remove_rownames %>% column_to_rownames(var="Run")
metadata_Kwok = metadata_Kwok[,-1]
colnames(metadata_Kwok)[2] <- "CLLtype"
rc_Kwok <- rc_Kwok %>% remove_rownames %>% column_to_rownames(var="Name")
rc_Kwok = rc_Kwok[,-2] #remove 2 column

### Check that sample names match in both files
all(colnames(rc_Kwok) %in% rownames(metadata_Kwok))
all(colnames(rc_Kwok) == rownames(metadata_Kwok))

## Create DESeq2Dataset object
dds <- DESeqDataSetFromMatrix(countData = rc_Kwok, colData = metadata_Kwok, design = ~ CLLtype)

#Generate normilized counts
dds <- estimateSizeFactors(dds)
sizeFactors(dds)
normalized_counts <- counts(dds, normalized=TRUE)

dds.vst<-DESeq2::varianceStabilizingTransformation(dds)
exp_Kwok <- assay(dds.vst)

save(exp_Kwok, file = "exp_Kwok.Rdata")

     