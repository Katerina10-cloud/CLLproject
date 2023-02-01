load("CLLdata/mofaCLL-main/data/reactomeGS.Rdata")
load("CLLdata/Kwok_data/rc_Kwok3_D.Rdata")
load("CLLdata/Kwok_data/rc_gc_Kwok3.Rdata")
load("RC_merge.Rdata")

library(biomaRt)
library(org.Hs.eg.db)

rc_gc_Kwok3 <- merge(rc_gc_Kwok3, SRR9140515_Counts.Rmatrix, by.x = "Geneid")
rc <- rc_gc_Kwok3
rc <- rc %>% remove_rownames %>% column_to_rownames(var="Geneid")

rc <- rc %>%
  relocate(ensembl_gene_id, .before = SRR9140505)
colnames(rc)[1] <- "ensembl_gene_id"

#Generate data frame
exprMat1 <- exprMat
exprMat1<- as.data.frame(exprMat1)

#exprMat1.df <- mapIds(org.Hs.eg.db, keys = rownames(exprMat1.df), keytype="ENSEMBL", column = "SYMBOL")
#exprMat1.df$Symbol <- mapIds(org.Hs.eg.db, keys = rownames(exprMat1.df), keytype="ENSEMBL", column = "SYMBOL")

rc$ensembl_gene_id <- row.names(rc)

#Extract Ensembl ID
ensembl_IDs <- row.names(rc)
ensembl_IDs <- as.data.frame(ensembl_IDs)


ensembl_IDs <- row.names(reactomeGS1)
ensembl_IDs <- as.data.frame(ensembl_IDs)

genes.table = NULL

ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#Generating my file
genes.table <- getBM(mart=ensembl, attributes= c("ensembl_gene_id","hgnc_symbol"), filters= "ensembl_gene_id", values=ensembl_IDs)

length(unique(genes.table$ensembl_gene_id))

genes.table_SYMBOLS <- merge(x = rc_Kwok3, y = genes.table, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T, all.y = T )

save(genes.table_SYMBOLS, file = "exprMat1_rna.Rdata")

##extracting NA SYMBOLS
NA_SYMBOLS <- genes.table_SYMBOLS %>% filter(hgnc_symbol == "")
NA_SYMBOLS1 <- NA_SYMBOLS %>% filter(gene_biotype == "protein_coding")

genes.table_name_symbols <- getBM(filters = "ensembl_gene_id", attributes= c("ensembl_gene_id", "external_gene_name", "hgnc_symbol", "gene_biotype"), values=ensembl_IDs, mart= ensembl)

genes.table_name_symbols <- merge(x = rc_Kwok3, y = genes.table_name_symbols, by.x = "ensembl_gene_id", by.y = "ensembl_gene_id", all.x = T, all.y = T )

#removing dublicates
rawCountFinal <- rawCounts_Kwok %>% group_by(external_gene_name) %>% summarise_if(is.numeric, sum)
save(rawCounts_Kwok, file = "rawCounts_Kwok.Rdata")

