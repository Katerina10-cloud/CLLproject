#Load packages and datasets
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, message = FALSE, dev = c("png","pdf"))
library(limma)
library(pheatmap)
library(jyluMisc)
library(survival)
library(survminer)
library(maxstat)
library(igraph)
library(tidygraph)
library(ggraph)
library(glmnet)
library(SummarizedExperiment)
library(cowplot)
library(tidyverse)

library(DESeq2)
library(proDA)
library(readr)
library(purrr)
library(stringr)
library(ggplot2)

load("CLLproteomics-main/data/patMeta_enc.RData")
load("CLLproteomics-main/data/ddsrna_enc.RData")
load("CLLproteomics-main/data/proteomic_explore_enc.RData")
load("CLLproteomics-main/output/deResList.RData") #precalculated differential expression
load("CLLproteomics-main/data/survival_enc.RData")
load("CLLproteomics-main/data/screenData_enc.RData")

source("C:/Users/rifug/Documents/CLLproteomics-main/code/utils.R")

viabMat <- screenData %>% filter(!lowQuality, ! Drug %in% c("DMSO","PBS"), patientID %in% colnames(protCLL)) %>%
  group_by(patientID, Drug) %>% summarise(viab = mean(normVal.adj.cor_auc)) %>%
  spread(key = patientID, value = viab) %>%
  data.frame(stringsAsFactors = FALSE) %>% column_to_rownames("Drug") %>%
  as.matrix()

proMat <- assays(protCLL)[["count"]]
proMat <- proMat[,colnames(viabMat)]

#Remove proteins without much variance (to lower multi-testing burden)

#sds <- genefilter::rowSds(proMat,na.rm=TRUE)
#proMat <- proMat[sds > genefilter::shorth(sds),]

ncol(proMat)

#Association test using univariate test
resTab.auc <- lapply(rownames(viabMat),function(drugName) {
  viab <- viabMat[drugName, ]
  batch <- protCLL[,colnames(viabMat)]$batch
  designMat <- model.matrix(~1+viab+batch)
  fit <- lmFit(proMat,  designMat)
  fit2 <- eBayes(fit)
  corRes <- topTable(fit2, number ="all", adjust.method = "BH", coef = "viab") %>% rownames_to_column("id") %>%
    mutate(symbol = rowData(protCLL[id,])$hgnc_symbol, Drug = drugName) %>%
    mutate(adj.P.Val = p.adjust(P.Value, method = "BH"))
}) %>% bind_rows() %>% arrange(P.Value)

#Bar plot to show the number of significant associations (5% FDR)
#Select significant associations (10% FDR)
resTab.sig <- filter(resTab.auc, adj.P.Val <= 0.05) %>%   
  select(Drug, symbol, id,logFC, P.Value, adj.P.Val)
plotTab <- resTab.sig %>%
  group_by(Drug) %>%
  summarise(n = length(id)) %>% ungroup()
ordTab <- group_by(plotTab, Drug) %>% summarise(total = sum(n)) %>%
  arrange(desc(total))
plotTab <- mutate(plotTab, Drug = factor(Drug, levels = ordTab$Drug)) %>%
  filter(n>0)

colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
theme_half <- theme_bw() + theme(axis.text = element_text(size=14),
                                 axis.title = element_text(size=16),
                                 axis.line =element_line(size=0.8),
                                 panel.border = element_blank(),
                                 axis.ticks = element_line(size=1.5),
                                 plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())

drugBar <- ggplot(plotTab, aes(x=Drug, y = n)) + geom_bar(stat="identity",fill=colList[4]) + 
  geom_text(aes(label = paste0(n)),vjust=-1,col="black", size=6) +
  ylim(0,500)+ #annotate("text", label = "Number of associations (10% FDR)", x=Inf, y=Inf,hjust=1, vjust=1, size=6)+
  theme_half + theme(axis.text.x = element_text(angle = 90, hjust=1, vjust=0.5)) +
  ylab("Number of associations (5% FDR)") + xlab("")

#Table of significant associations (5% FDR)
resTab.sig %>% mutate_if(is.numeric, formatC, digits=2, format= "e") %>%
  DT::datatable()

#Correlation plot of selected protein-drug pairs
#Stratified by IGHV and trisomy12
proMat.combat <- assays(protCLL)[["count_combat"]]
proMat.combat <- proMat.combat[,colnames(viabMat)]

pairList <- list(c("Cobimetinib","STAT2"),c("Trametinib", "PTPN11"), c("Ibrutinib","LYN"),c("Ibrutinib","ANXA2"))
plotList <- lapply(pairList, function(pair) {
  textCol <- "darkred"
  drugName <- pair[1]
  proteinName <- pair[2]
  id <- rownames(protCLL)[match(proteinName, rowData(protCLL)$hgnc_symbol)]
  plotTab <- tibble(patID = colnames(viabMat), 
                    viab = viabMat[drugName,],
                    expr = proMat.combat[id,]) %>%
    mutate(IGHV = protCLL[,patID]$IGHV.status,
           trisomy12 = protCLL[,patID]$trisomy12) %>%
    mutate(trisomy12 = ifelse(trisomy12 ==1,"yes","no")) %>%
    filter(!is.na(viab),!is.na(expr))
  
  pval <- formatNum(filter(resTab.sig, Drug == drugName, symbol == proteinName)$P.Value, digits = 1, format = "e")
  Rval <- sprintf("%1.2f",cor(plotTab$viab, plotTab$expr))
  Nval <- nrow(plotTab)
  annoP <- bquote(italic("P")~"="~.(pval))
  annoN <- bquote(N~"="~.(Nval))
  annoCoef <- bquote(R~"="~.(Rval))
  
  corPlot <- ggplot(plotTab, aes(x = viab, y = expr)) + 
    geom_point(aes(col = trisomy12, shape = IGHV), size=5) +
    scale_shape_manual(values = c(M = 19, U = 1)) + 
    scale_color_manual(values = c(yes = colList[2], no = colList[3])) +
    geom_smooth(formula = y~x,method = "lm", se=FALSE, color = "grey50", linetype ="dashed" ) +
    ggtitle(sprintf("%s ~ %s", drugName, proteinName)) +
    ylab("Protein expression") + xlab("Viability after treatment") +
    theme_full +
    theme(legend.position = "bottom") 
  
  if (Rval < 0) annoPos <- "right" else annoPos <- "left"
  
  if (annoPos == "right") {
    
    corPlot <- corPlot + annotate("text", x = max(plotTab$viab), y = Inf, label = annoN,
                                  hjust=1, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$viab), y = Inf, label = annoP,
               hjust=1, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = max(plotTab$viab), y = Inf, label = annoCoef,
               hjust=1, vjust =6, size = 5, parse = FALSE, col= textCol)
    
  } else if (annoPos== "left") {
    corPlot <- corPlot + annotate("text", x = min(plotTab$viab), y = Inf, label = annoN,
                                  hjust=0, vjust =2, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$viab), y = Inf, label = annoP,
               hjust=0, vjust =4, size = 5, parse = FALSE, col= textCol) +
      annotate("text", x = min(plotTab$viab), y = Inf, label = annoCoef,
               hjust=0, vjust =6, size = 5, parse = FALSE, col= textCol)
  }
  
  corPlot <- corPlot + ylab("Protein expression") + xlab("Viability after treatment") + 
    scale_y_continuous(labels = scales::number_format(accuracy = 0.1)) +
    scale_x_continuous(labels = scales::number_format(accuracy = 0.1))
  
  corPlot
})
drugCor <- cowplot::plot_grid(plotlist = plotList, ncol =2)
drugCor

#Association test with blocking for IGHV and trisomy12
testList <- filter(resTab.auc, adj.P.Val <= 0.05)
resTab.auc.block <- lapply(seq(nrow(testList)),function(i) {
  pair <- testList[i,]
  expr <- proMat[pair$id,]
  viab <- viabMat[pair$Drug, ]
  ighv <- protCLL[,colnames(viabMat)]$IGHV.status
  tri12 <- protCLL[,colnames(viabMat)]$trisomy12
  batch <- protCLL[,colnames(viabMat)]$batch
  res <- anova(lm(viab~ighv+tri12+batch+expr))
  data.frame(id = pair$id, P.Value = res["expr",]$`Pr(>F)`, symbol = pair$symbol,
             Drug = pair$Drug,
             P.Value.IGHV = res["ighv",]$`Pr(>F)`,P.Value.trisomy12 = res["tri12",]$`Pr(>F)`,
             P.Value.noBlock = pair$P.Value,
             stringsAsFactors = FALSE)
  
}) %>% bind_rows() %>% mutate(adj.P.Val = p.adjust(P.Value, method = "BH")) %>% arrange(P.Value)

#Prepare genomic annotations
geneMat <-  patMeta[match(colnames(proMat), patMeta$Patient.ID),] %>%
  select(Patient.ID, IGHV.status, del11q:U1) %>%
  mutate_if(is.factor, as.character) %>% mutate(IGHV.status = ifelse(IGHV.status == "M", 1,0)) %>%
  mutate_at(vars(-Patient.ID), as.numeric) %>% #assign a few unknown mutated cases to wildtype 
  mutate_all(replace_na,0) %>%
  data.frame() %>% column_to_rownames("Patient.ID")

#Protein markers for clincial outcomes
#Uni-variate model to identify proteins associated with outcomes
protCLL.sub <- protCLL[!rowData(protCLL)$chromosome_name %in% c("X","Y"),]
protMat <- assays(protCLL.sub)[["count_combat"]]


survTab <- survT %>% 
  select(patID, OS, died, TTT, treatedAfter, age,sex) %>%
  dplyr::rename(patientID = patID) %>%
  filter(patientID %in% colnames(protMat))

#function for cox regression
com <- function(response, time, endpoint, scale =FALSE) {
  
  if (scale) {
    #calculate z-score
    response <- (response - mean(response, na.rm = TRUE))/sd(response, na.rm=TRUE)
  }
  surv <- coxph(Surv(time, endpoint) ~ response)
  
  
  tibble(p = summary(surv)[[7]][,5],
         HR = summary(surv)[[7]][,2],
         lower = summary(surv)[[8]][,3],
         higher = summary(surv)[[8]][,4])
}
uniRes.ttt <- lapply(rownames(protMat), function(n) {
  testTab <- mutate(survTab, expr = protMat[n, patientID])
  com(testTab$expr, testTab$TTT, testTab$treatedAfter, TRUE) %>%
    mutate(id = n)
}) %>% bind_rows() %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p) %>% mutate(name = rowData(protCLL[id,])$hgnc_symbol) %>%
  mutate(outcome = "TTT")

uniRes.os <- lapply(rownames(protMat), function(n) {
  testTab <- mutate(survTab, expr = protMat[n, patientID])
  com(testTab$expr, testTab$OS, testTab$died, TRUE) %>%
    mutate(id = n)
}) %>% bind_rows() %>% mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p) %>% mutate(name = rowData(protCLL[id,])$hgnc_symbol) %>%
  mutate(outcome = "OS")

uniRes <- bind_rows(uniRes.ttt, uniRes.os) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

#table of known risks
riskTab <- select(survTab, patientID, age, sex) %>%
  left_join(patMeta[,c("Patient.ID","IGHV.status","TP53","trisomy12","del17p")], by = c(patientID = "Patient.ID")) %>% 
  mutate(TP53 = as.numeric(as.character(TP53)),
         del17p = as.numeric(as.character(del17p))) %>%
  mutate(`TP53.del17p` = as.numeric(TP53 | del17p),
         IGHV = factor(ifelse(IGHV.status %in% "U",1,0))) %>%
  select(-TP53, -del17p,-IGHV.status) %>%
  mutate(age = age/10) 

#Function to run multivariate Cox model on test table
runCox <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- coxph(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#Multi-variate test
cTab.ttt <- lapply(filter(uniRes, outcome == "TTT")$id, function(n) {
  risk0 <- riskTab
  expr <- protMat[n,]
  expr <- (expr - mean(expr,na.rm=TRUE))/sd(expr,na.rm = TRUE)
  risk1 <- riskTab %>% mutate(protExpr = expr[patientID])
  res0 <- summary(runCox(survTab, risk0, "TTT","treatedAfter"))
  fullModel <- runCox(survTab, risk1, "TTT","treatedAfter")
  res1 <- summary(fullModel)
  tibble(id = n, c0 = res0$concordance[1], c1 = res1$concordance[1],
         se0 = res0$concordance[2],se1 = res1$concordance[2],
         ci0 = se0*1.96, ci1 = se1*1.96,
         p = res1$coefficients["protExpr",5],
         fullModel = list(fullModel))
}) %>% bind_rows() %>% mutate(diffC = c1-c0) %>%
  arrange(desc(diffC)) %>%
  mutate(name=rowData(protCLL[id,])$hgnc_symbol,
         outcome = "TTT")

cTab.os <- lapply(filter(uniRes, outcome == "OS")$id, function(n) {
  risk0 <- riskTab
  expr <- protMat[n,]
  expr <- (expr - mean(expr,na.rm=TRUE))/sd(expr,na.rm = TRUE)
  risk1 <- riskTab %>% mutate(protExpr = expr[patientID])
  res0 <- summary(runCox(survTab, risk0, "OS","died"))
  fullModel <- runCox(survTab, risk1, "OS","died")
  res1 <- summary(fullModel)
  
  tibble(id = n, c0 = res0$concordance[1], c1 = res1$concordance[1],
         se0 = res0$concordance[2],se1 = res1$concordance[2],
         ci0 = se0*1.96, ci1 = se1*1.96,
         p = res1$coefficients["protExpr",5],
         fullModel = list(fullModel))
}) %>% bind_rows() %>% mutate(diffC = c1-c0) %>%
  arrange(desc(diffC)) %>%
  mutate(name=rowData(protCLL[id,])$hgnc_symbol,
         outcome = "OS")

cTab <- bind_rows(cTab.ttt, cTab.os) %>%
  mutate(p.adj = p.adjust(p, method = "BH")) %>%
  arrange(p)

#Forest plot of several markers as examples 
plotHazard <- function(survRes, protName, title = "", xLim = c(0.2,6)) {
  sumTab <- summary(survRes)$coefficients
  confTab <- summary(survRes)$conf.int
  #correct feature name
  nameOri <- rownames(sumTab)
  nameMod <- substr(nameOri, 1, nchar(nameOri) -1)
  plotTab <- tibble(feature = rownames(sumTab),
                    nameMod = substr(nameOri, 1, nchar(nameOri) -1),
                    HR = sumTab[,2],
                    p = sumTab[,5],
                    Upper = confTab[,4],
                    Lower = confTab[,3]) %>%
    mutate(feature = ifelse(nameMod %in% names(survRes$xlevels), nameMod, feature)) %>%
    mutate(feature = str_replace(feature, "[.]","/")) %>%
    mutate(feature = str_replace(feature, "[_]","-")) %>%
    mutate(feature = str_replace(feature, "IGHV","IGHV-U")) %>%
    mutate(candidate = ifelse(feature == "protExpr", "yes","no")) %>%
    mutate(feature = ifelse(feature == "protExpr", protName, feature)) %>%
    #arrange(desc(abs(p))) %>% 
    mutate(feature = factor(feature, levels = feature)) #%>%
  #mutate(type = ifelse(HR >1 ,"up","down")) %>%
  # mutate(Upper = ifelse(Upper > 10, 10, Upper))
  
  p <- ggplot(plotTab, aes(x=feature, y = HR, color = candidate)) +
    geom_hline(yintercept = 1, linetype = "dotted") +
    geom_point(position = position_dodge(width=0.8), size=3) +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, size=1) +
    geom_text(position = position_nudge(x = 0.3),
              aes(y = HR, label =  sprintf("italic(P)~'='~'%s'",
                                           formatNum(p, digits = 1))),
              color = "black", size =5, parse = TRUE) +
    scale_color_manual(values = c(yes = "darkred", no = "black")) +
    ggtitle(title) + scale_y_log10(limits = xLim) +
    ylab("Hazard ratio") +
    coord_flip() +
    theme_full +
    theme(legend.position = "none", axis.title.y = element_blank())
  return(p)
}
#Function to format floats
formatNum <- function(i, limit = 0.01, digits =1, format="e") {
  r <- sapply(i, function(n) {
    if (n < limit) {
      formatC(n, digits = digits, format = format)
    } else {
      format(n, digits = digits)
    }
  })
  return(r)
}

#MTSS1
protName <- "MTSS1"
outcomeName <- "TTT"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.prmt5 <- plotHazard(survRes, protName, outcomeName)
hr.prmt5

#MTSS1
protName <- "MTSS1"
outcomeName <- "OS"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.mtss1 <- plotHazard(survRes, protName, outcomeName)
hr.mtss1

protName <- "BIN1"
outcomeName <- "TTT"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.bin1 <- plotHazard(survRes, protName, outcomeName)
hr.bin1

#BIN1
protName <- "BIN1"
outcomeName <- "OS"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.bin1 <- plotHazard(survRes, protName, outcomeName)
hr.bin1

#MYO1E
protName <- "MYO1E"
outcomeName <- "TTT"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.myo1E <- plotHazard(survRes, protName, outcomeName)
hr.myo1E

#FOXO1
protName <- "FOXO1"
outcomeName <- "OS"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.foxo1 <- plotHazard(survRes, protName, outcomeName)
hr.foxo1

#FOXO1
protName <- "FOXO1"
outcomeName <- "TTT"
survRes <- filter(cTab, outcome == outcomeName , name == protName)$fullModel[[1]]
hr.foxo1 <- plotHazard(survRes, protName, outcomeName)
hr.foxo1


