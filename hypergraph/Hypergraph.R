###Hypergraph Dynamic Correlation###

BiocManager::install("coin")
BiocManager::install("dynamicTreeCut")
BiocManager::install("org.Sc.sgd.db")

library(tidyverse)
library(igraph)
library(coin)
library(survival)
library(ks)
library(dynamicTreeCut)
library(GOstats)
library(org.Hs.eg.db)
library(clusterProfiler)

#Preliminary process
load('spellman_73_filled.bin')
load('Hypergraph/GO_select_yeast.bin')
dyn.load("hypergraph/csupp.dll")

array <- DEG_2down_Kwok

array <- array[sample(nrow(array), 1000),]

source("main.r") #load the main code file
source("utils.r")
source("visualize.r")


###Unsupervised hypergraph construction###
hp_us <- hypergraph_unsup(array, fdr=0.2, save_glist=T)

save(DEG_F6F4, file = "DEG_F6F4.Rdata")

#Visualize the full hypergraph and top-connected sub-hypergraph
full_hypergraph_us <- plot_entrie_unsup(hp_us$net, hp_us$label, folds=2)
plot_top_unsup(full_hypergraph_us$g, full_hypergraph_us$elist, hp_us$label, full_hypergraph_us$folds)

#Extract gene list
glist_us <- hp_us[["glist"]]
ids <- bitr(glist_us, fromType = "SYMBOL", toType = "ENTREZID", OrgDb=org.Hs.eg.db)
entrezID <- as.vector(ids$ENTREZID)

hp_us <- append(hp_us,list(entrezID=entrezID))

#GO enrichment test
library(org.Hs.eg.db)

GOen<-function(entrezID,label){
  library(GOstats)
  all.entrez=entrezID
  res=list()
  k=ncol(label)-1
  for (i in 1:k){
    sel.entrez=entrezID[which(label[,i]==1)]
    params <- new("GOHyperGParams", geneIds=sel.entrez
                  ,universeGeneIds=all.entrez, ontology="BP"
                  ,pvalueCutoff=0.01,conditional=F 
                  ,testDirection="over"
                  #,annotation="org.Sc.sgd.db") # yeast
                  ,annotation="org.Hs.eg.db") # human/skin
    over.pres<-hyperGTest(params)
    sum = summary(over.pres)
    id = sum$Size %in% 5:350
    res[[i]] <- sum[id,]
    cat(i,"out of",k,'done.\n')
  }
  return(res)
}

enrichment <- GOen(hp_us$entrezID, hp_us$label)


save(enrichment, file = "hypergraph/GOenr.RData")




#Supervised hypergraph construction
hp_sp <- hypergraph_sup(array, GO.select, fdr=0.2, save_triplets=T, save_glist=T)

#Visualization
full_hypergraph_sp <- plot_entrie_sup(hp_sp$net, hp_sp$label, GO.select, folds=10)
plot_top_sup(full_hypergraph_sp$g, full_hypergraph_sp$elist, hp_sp$label, GO.select, full_hypergraph_sp$folds)
plot_one_sup(full_hypergraph_sp$g, full_hypergraph_sp$elist, hp_sp$label, GO.select, full_hypergraph_sp$folds, 
             "GO:0006090", max_num_display = 16)

GO <- names(GO.select)
module_names <- Term(GO)
which(GO == "GO:0006090")
which(GO == "GO:0045937")
which(GO == "GO:2000241")
hyperedge <- c(232, 141, 1)

plot_gene_level(hyperedge, GO, hp_sp$net, hp_sp$label, hp_sp$triplets, hp_sp$glist, folds=10)
