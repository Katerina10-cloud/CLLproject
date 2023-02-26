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

#Preliminary process
load('spellman_73_filled.bin')
load('Hypergraph/GO_select_yeast.bin')
dyn.load("hypergraph/csupp.dll")

array <- DEG_Kwok_3groups

array <- array[sample(nrow(array), 1000),]

source("main.r") #load the main code file
source("utils.r")
source("visualize.r")


###Unsupervised hypergraph construction###
hp_us <- hypergraph_unsup(array, fdr=0.2, save_glist=T)


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
    id = sum$Size %in% 5:500
    res[[i]] <- sum[id,]
    cat(i,"out of",k,'done.\n')
  }
  return(res)
}

enrichment <- GOen(hp_us$entrezID, hp_us$label)

save(enrichment, file = "hypergraph/GOenr.RData")




#Supervised hypergraph construction
