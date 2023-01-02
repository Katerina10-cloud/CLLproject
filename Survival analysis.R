#Survival analysis

load("CLLproject/Survival_analysis/survTab.RData")
load("CLLproject/Survival_analysis/riskTab.RData")

library(data.table)
library(survival)
library(survminer)
library(tidyverse)
library(maxstat)
library(gridExtra)

#dataframe trasformations
rna2 <- rna2[c("LPL", "ZAP70", "CD38", "LRP1", "LRP5", "LRP6", "LRP8", "LDLR", "VLDLR", "LILRB4", "AEBP1", "LAMA5", "COL9A3", "MTSS1", "BIN1", "MYO1E", "COTL1", "LSP1", "PIGR", "IL2RA", "CXCR5", "CD44", "IFNGR1", "CD6", "ITGB1", "CD84", "SELL", "CEMIP2", "RHOB", "TGFBR2", "ALOX5", "RGCC", "LILRB4", "JUND", "KLF2", "ZNF331", "ARID5B", "EGR1", "CCDC88A", "NR4A3", "ID3", "SREBF2", "ZBTB10", "FOXN3", "ZNF395", "NFATC1", "ZEB2", "ZBTB18", "KLF9", "TGIF1", "NCOA1", "BCL6", "ZNF83", "HHEX", "FOSL2", "ZNF529", "SPI1", "ZNF532", "ZNF91", "ZNF831", "PBX3", "ZNF821", "ZNF432", "FOXO1", "ZNF264", "ZNF350", "RXRA", "MAFF", "HIVEP3", "MAFG"),]
rna2$name <- row.names(rna2)#convert rownames to column
rna2 <- rna2 %>%
  select(name, everything()) # moving the last column to the start
#change row names to a list of integers
row.names(rna2) <- 1:nrow(rna2)#change row names to a list of integers
colnames(rna2)[1]  <- "patientID"
riskTab2 <- riskTab2[,-3]
riskTab <- merge(rna2, riskTab1, by.x = "patientID") #merging 2 datasets
riskTab <- riskTab %>%
  select(IGHV, everything()) # moving the last column to the start
riskTab <- riskTab %>% relocate(patientID, .before = IGHV)

survTab <- survival %>% 
  filter(patientID %in% colnames(rna2))
survTab<- replace(survTab, is.na(survTab), 0)

#Function to format floats
#' @export
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
#function for italicize gene names
#' @export
italicizeGene <- function(featureNames) {
  geneNameList <- c("SF3B1","NOTCH1","TP53")
  specialCase <- "TP53/del17p"
  formatedList <- sapply(featureNames, function(n) {
    if (n %in% geneNameList) {
      sprintf("italic('%s')",n)
    } else if (n == specialCase) {
      sprintf("italic('TP53')~'/del17p'")
    } else if (n == "CLLPD") {
      "CLL-PD"
    }
    else {
      sprintf("'%s'",n)
    }
  })
  return(formatedList)
}
#function for plot hazard ratio
#' @export
plotHazard <- function(survRes, title = "") {
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
    arrange(desc(abs(p))) %>%
    mutate(feature = italicizeGene(feature)) %>%
    mutate(feature = factor(feature, levels = feature)) %>%
    mutate(type = ifelse(HR >1 ,"up","down")) %>%
    mutate(Upper = ifelse(Upper > 10, 10, Upper))
  
  ggplot(plotTab, aes(x=feature, y = HR, color = type)) +
    geom_hline(yintercept = 1, linetype = "dotted", color = "grey50") +
    geom_point(position = position_dodge(width=0.8), size=3, color = "black") +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.3, size=1,color = "grey20") +
    geom_text(position = position_nudge(x = 0.3),
              aes(y = HR, label =  sprintf("italic(P)~'='~'%s'",
                                           formatNum(p, digits = 1))),
              color = "black", size =4.5, parse = TRUE) +
    expand_limits(y=c(-0.5,0))+
    #scale_color_manual(values = c(up = colList[1], down = colList[2])) +
    ggtitle(title) + scale_y_log10() +
    scale_x_discrete(labels = parse(text = levels(plotTab$feature))) +
    ylab("Hazard ratio") +
    coord_flip() +
    theme_full +
    theme(legend.position = "none", axis.title.y = element_blank())
}
#Function to run multivariate Cox model on test table
#' @export
runCox <- function(survTab, riskTab, time, endpoint) {
  survTab <- select(survTab, patientID, !!time, !!endpoint) %>%
    dplyr::rename(time = !!time, endpoint = !!endpoint) %>%
    dplyr::filter(!is.na(time), !is.na(endpoint))
  testTab <- right_join(survTab, riskTab, by = "patientID") %>%
    select(-patientID)
  surv1 <- coxph(
    Surv(time, endpoint) ~
      .,
    data = testTab)
  return(surv1)
}

#List of genes("LPL", "ZAP70", "CD38", "LRP1", "LRP5", "LRP6", "LRP8", "LDLR", "VLDLR", "LILRB4", 
"AEBP1", "LAMA5", "COL9A3", "MTSS1", "BIN1", "MYO1E", "COTL1", "LSP1", "PIGR", "IL2RA", "CXCR5", 
"CD44", "IFNGR1", "CD6", "ITGB1", "CD84", "SELL", "CEMIP2", "RHOB", "TGFBR2", "ALOX5", "RGCC", "LILRB4", 
"JUND", "KLF2", "ZNF331", "ARID5B", "EGR1", "CCDC88A", "NR4A3", "ID3", "SREBF2", "ZBTB10", "FOXN3", "ZNF395", 
"NFATC1", "ZEB2", "ZBTB18", "KLF9", "TGIF1", "NCOA1", "BCL6", "ZNF83", "HHEX", "FOSL2", "ZNF529", "SPI1", 
"ZNF532", "ZNF91", "ZNF831", "PBX3", "ZNF821", "ZNF432", "FOXO1", "ZNF264", "ZNF350", "RXRA", "MAFF", "HIVEP3", "MAFG"),]

riskTab3 <- subset(riskTab, select = c("patientID", "IGHV", "LPL", "AEBP1", "LAMA5"))
rna2$name <- row.names(rna2)#convert rownames to column))

#Multivariate cox regression for gene of interest
resTTT <- runCox(survTab, riskTab3, "TTT", "treatedAfter")

haTTT <- plotHazard(resTTT, title = "TTT") +
  scale_y_log10(limits=c(0.2,5))
haTTT

resOS <- runCox(survTab, riskTab3, "OS", "died")
#summary(surv1)
haOS <- plotHazard(resOS,"OS") + scale_y_log10(limits=c(0.05,5))
haOS


#Kaplan-Meiler plots
#Survival analysis with gene expression
#Kaplan-Meiler plots

BiocManager::install("RegParallel")
library(Biobase)
library(dplyr)
library(survival)
library(survminer)
library(RegParallel)
library(ggplot2)

load("CLLproject/Survival_analysis/res.RData")
load("CLLproject/Survival_analysis/coxdata.RData")

riskTab2 <- riskTab
# reassigning row names
riskTab2 <- riskTab2 %>% remove_rownames %>% column_to_rownames(var="patientID")
survTab <- survTab %>% remove_rownames %>% column_to_rownames(var="patientID")

# transform the expression data to Z scores
riskTab2 <- scale(riskTab2)

#riskTab1 <- t(scale(t(riskTab1)))

coxdata <- merge(riskTab2, survTab,
                 by = 'row.names', all = TRUE)
# reassigning row names
coxdata <- coxdata %>% remove_rownames %>% column_to_rownames(var="Row.names")

#fit the Cox model independently for each gene
res.ttt <- RegParallel(
  data = coxdata,
  formula = 'Surv(TTT, treatedAfter) ~ [*]',
  FUN = function(formula, data)
    coxph(formula = formula,
          data = data,
          ties = 'breslow',
          singular.ok = TRUE),
  FUNtype = 'coxph',
  variables = colnames(coxdata),
  blocksize = 76,
  conflevel = 95
)
res.ttt <- res.ttt[!is.na(res.ttt$P),]
res.ttt

# extract RFS and probe data for downstream analysis
survplotdata <- coxdata[,c('TTT', 'treatedAfter', 'LPL', 'ZAP70')]
colnames(survplotdata) <- c('TTT', 'treatedAfter', 'LPL', 'ZAP70')

# set Z-scale cut-offs for high and low expression
highExpr <- 1.0
lowExpr <- -1.0
survplotdata$LPL <- ifelse(survplotdata$LPL >= highExpr, 'High',
                           ifelse(survplotdata$LPL <= lowExpr, 'Low', 'Mid'))
survplotdata$ZAP70 <- ifelse(survplotdata$ZAP70 >= highExpr, 'High',
                             ifelse(survplotdata$ZAP70 <= lowExpr, 'Low', 'Mid'))
survplotdata$RXRA <- ifelse(survplotdata$RXRA >= highExpr, 'High',
                            ifelse(survplotdata$RXRA <= lowExpr, 'Low', 'Mid'))
survplotdata$MAFF <- ifelse(survplotdata$MAFF >= highExpr, 'High',
                            ifelse(survplotdata$MAFF <= lowExpr, 'Low', 'Mid'))
survplotdata$HIVEP3 <- ifelse(survplotdata$HIVEP3 >= highExpr, 'High',
                              ifelse(survplotdata$HIVEP3 <= lowExpr, 'Low', 'Mid'))
survplotdata$MAFG <- ifelse(survplotdata$MAFG >= highExpr, 'High',
                            ifelse(survplotdata$MAFG <= lowExpr, 'Low', 'Mid'))

# relevel the factors to have mid as the ref level
survplotdata$LPL <- factor(survplotdata$LPL,
                           levels = c('Mid', 'Low', 'High'))
survplotdata$ZAP70 <- factor(survplotdata$ZAP70,
                             levels = c('Mid', 'Low', 'High'))
survplotdata$RXRA <- factor(survplotdata$RXRA,
                            levels = c('Mid', 'Low', 'High'))
survplotdata$MAFF <- factor(survplotdata$MAFF,
                            levels = c('Mid', 'Low', 'High'))
survplotdata$HIVEP3 <- factor(survplotdata$HIVEP3,
                              levels = c('Mid', 'Low', 'High'))
survplotdata$MAFG <- factor(survplotdata$MAFG,
                            levels = c('Mid', 'Low', 'High'))

p <- ggsurvplot(survfit(Surv(TTT, treatedAfter) ~ ZAP70,
                        data = survplotdata),
                data = survplotdata,
                risk.table = TRUE,
                title = "ZAP70 vs TTT",
                pval = TRUE,
                ggtheme = theme_bw(), 
                tables.theme = theme_cleantable(),
                risk.table.y.text.col = TRUE,
                risk.table.y.text = FALSE)
p


