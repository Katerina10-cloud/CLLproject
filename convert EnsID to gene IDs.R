library(biomaRt)
#Generate data frame
exprMat1 <- exprMat

exprMat1 <- data.frame(exprMat1)

exprMat1$ensembl_gene_id_version <- row.names(exprMat1)

#Extract Ensembl ID
ensembl_IDs <- row.names(exprMat1)
view (ensembl_IDs)
genes.table = NULL
ensembl <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

#Generating my file
genes.table <- getBM(filters= "ensembl_gene_id_version",
                     attributes= c("ensembl_gene_id_version",
                                   "hgnc_symbol", "gene_biotype"), values= ensembl_IDs, mart= ensembl)
affyids <- c("202763_at","209310_s_at","207500_at")
getBM(attributes = c('affy_hg_u133_plus_2', 'entrezgene_id'),
      filters = 'affy_hg_u133_plus_2',
      values = affyids, 
      mart = ensembl)
#controllo quali valori si perdono con genemaRT
ensLOST <- ensembl_IDs[!ensembl_IDs %in% genes.table$ensembl_gene_id_version]
View(ensLOST) #44 valori persi tutti quelli con PAR.Y

length(unique(genes.table$ensembl_gene_id_version))

