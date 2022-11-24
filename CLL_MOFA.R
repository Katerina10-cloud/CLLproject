knitr::opts_chunk$set(warning = FALSE, message = FALSE)

BiocManager::install("data.frame")
BiocManager::install("org.Hs.eg.db")
BiocManager::install("MOFA2")
BiocManager::install("biomaRt")

#Part1_Prepare MOFA model

library(MOFA2)
library(reticulate)
library(DESeq2)
library(sva)
library(MultiAssayExperiment)
library(tidyverse)
library(BiocParallel)
options(stringsAsFactors=FALSE)

#set the global ggplot theme
theme_set(theme_bw() + theme(axis.text = element_text(size=12), 
                             axis.title = element_text(size=14),
                             plot.title = element_text(size = 15, hjust =0.5, face="bold")))

#Prepare each single omic dataset
#Load datasets
load("CLLproject/mofaCLL-main/data/gene.Rdata")
load("CLLproject/mofaCLL-main/data/rna.Rdata")
load("CLLproject/mofaCLL-main/data/meth.Rdata")

#RNA sequencing
rna.vst<-DESeq2::varianceStabilizingTransformation(rna)
exprMat <- assay(rna.vst)
nTop = 5000
sds <- genefilter::rowSds(exprMat)
exprMat <- exprMat[order(sds, decreasing = T)[1:nTop],]

#DNA methylation array
methData <- assays(meth)[["beta"]]
nTop = 5000
sds <- genefilter::rowSds(methData)
methData <- methData[order(sds, decreasing = T)[1:nTop],]

#Genetics
gene <- t(gene)

#Assemble multiAssayExperiment object
#Create object
mofaData <- list(mRNA = exprMat,
                 Mutations = gene, 
                 Methylation = methData)

# Create MultiAssayExperiment object 
mofaData <- MultiAssayExperiment::MultiAssayExperiment(
  experiments = mofaData
)

#Only keep samples that have at least three assays
useSamples <- MultiAssayExperiment::sampleMap(mofaData) %>%
  as_tibble() %>% group_by(primary) %>% summarise(n= length(assay)) %>%
  filter(n >= 3) %>% pull(primary)
mofaData <- mofaData[,useSamples]

#Dimensions for each dataset
experiments(mofaData)

#How many samples have the complete datasets
table(table(sampleMap(mofaData)$primary))

#Build MOFA object
#Build MOFA object from multiAssayExperiment object
MOFAobject <- create_mofa(mofaData)
MOFAobject

plot_data_overview(MOFAobject)

#Setup MOFA training parameters

#Define data options
data_opts <- get_default_data_options(MOFAobject)
data_opts


#Define model options
model_opts <- get_default_model_options(MOFAobject)
model_opts$num_factors <- 30

model_opts

model_opts$likelihoods[["Mutations"]] <- "bernoulli"
model_opts#$likelihoods


#Define training options
train_opts <- get_default_training_options(MOFAobject)
train_opts$convergence_mode <- "slow"
train_opts$seed <- 42
train_opts$maxiter <- 5000
train_opts$drop_factor_threshold <- 0.02
train_opts


#Run MOFA model
MOFAobject <- prepare_mofa(
  object = MOFAobject,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

#reticulate::use_python
#use_condaenv("r-reticulate", required = TRUE)

MOFAobject <- run_mofa(MOFAobject, outfile="/Users/rifug/Documents/CLLproject/CLL-main/data/MOFA2_CLL.Rdata")
save(MOFAobject, file = "MOFA2_CLL.RData")
saveRDS(MOFAobject, "MOFA2_CLL.hdf5")


#Part 2
#Preliminary analysis of the results
#### Variance explained by MOFA for each omic data


load("CLLproject/mofaCLL-main/data/gene.Rdata")
load("CLLproject/mofaCLL-main/data/demographic.Rdata")
load("CLLproject/mofaCLL-main/data/rna.Rdata")
load("CLLproject/mofaCLL-main/data/survival.Rdata")
load("CLLproject/mofaCLL-main/data/doublingTime.Rdata")
load("CLLproject/mofaCLL-main/data/reactomeGS.Rdata")
load("CLLproject/mofaCLL-main/data/protein.Rdata")
load("CLLproject/mofaCLL-main/data/MOFA2_CLL.Rdata")

library(mltools)
library(data.table)
library(survival)
library(survminer)
library(maxstat)
library(gridExtra)
library(car)
library(cowplot)
library(egg)

plot_factor_cor(MOFAobject)

r2 <-calculate_variance_explained(MOFAobject)
r2$r2_total
r2$r2_per_factor

varExpPlot <- plot_variance_explained(MOFAobject, censor = 0.15 )
plot_variance_explained(MOFAobject, plot_total = T)[[2]]
varExpPlot
plot_variance_explained(MOFAobject, max_r2=15)

plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "Factor1"
)
plot_weights(MOFAobject,
             view = "Mutations",
             factor = 1,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
plot_factor(MOFAobject, 
            factors = 1, 
            color_by = "IGHV",
            add_violin = TRUE,
            dodge = TRUE
)

plot_weights(MOFAobject,
             view = "mRNA",
             factor = 1,
             nfeatures = 10,     
             scale = T           
)
plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "IGHV"
) + labs(y="RNA expression")

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "negative",
                  color_by = "IGHV"
) + labs(y="RNA expression")

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
#Methylation profile
plot_top_weights(MOFAobject,
                 view = "Methylation",
                 factor = 1,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
plot_weights(MOFAobject,
             view = "Methylation",
             factor = 1,
             nfeatures = 10,     
             scale = T    
)
plot_data_heatmap(MOFAobject, 
                  view = "Methylation",
                  factor = 1,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)


#Factor2 profile
plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Factor2"
)
plot_weights(MOFAobject,
             view = "Mutations",
             factor = 2,
             nfeatures = 10,     # Top number of features to highlight
             scale = T           # Scale weights from -1 to 1
)
plot_top_weights(MOFAobject,
                 view = "Mutations",
                 factor = 2,
                 nfeatures = 10,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "trisomy12",
            add_violin = TRUE,
            dodge = TRUE
)

plot_weights(MOFAobject,
             view = "mRNA",
             factor = 2,
             nfeatures = 10,     
             scale = T           
)
plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "positive",
                  color_by = "trisomy12"
) + labs(y="RNA expression")

plot_data_scatter(MOFAobject, 
                  view = "mRNA",
                  factor = 1,  
                  features = 4,
                  sign = "negative",
                  color_by = "trisomy12"
) + labs(y="RNA expression")

plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  factor = 2,  
                  features = 25,
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = FALSE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)
#Factor heatmap
library(pheatmap)
#gene annotation
facMat <- t(get_factors(MOFAobject)[[1]])

seleGenes <- c("IGHV", "trisomy12")

colAnno <- tibble(Name = colnames(gene),
                  IGHV = gene["IGHV",],
                  trisomy12 = gene["trisomy12",]) %>%
  column_to_rownames("Name") %>% data.frame()
pheatmap(facMat, clustering_method = "complete", annotation_col = colAnno)


#Meaning of the first 2 factors
allWeights <- get_weights(MOFAobject,
                         views = "all",
                         factors = "all",
                         as.data.frame = TRUE) %>% as_tibble() %>%
  mutate(feature = ifelse(feature == "IGHV.status","IGHV",feature),
         factor = gsub("LF","F",factor))

allFactors <- get_factors(
  MOFAobject, 
  factors = "all",
  as.data.frame = TRUE
) %>% as_tibble() %>%
  mutate(factor = gsub("LF","F", factor))

patAnno <- gene[c("IGHV", "trisomy12")] %>% data.frame() %>%
  rownames_to_column("sample") %>%
  dplyr::filter(!is.na(IGHV),!is.na(trisomy12)) %>%
  mutate(IGHV = ifelse(IGHV==1, "M","U"),
         trisomy12 = ifelse(trisomy12 ==1, "yes","no"))

allFactors <- left_join(allFactors, 
                        patAnno,
                        by = "sample")

#Associations between F1 and epigenetic subtypes
plotTab <- filter(allFactors, factor %in% c("Factor1", "Factor2")) %>%
  mutate(Epitype = methCluster[sample]) %>%
  spread(key = factor ,value = value)

#set the global ggplot theme
colList <- c("#BC3C29FF","#0072B5FF","#E18727FF","#20854EFF","#7876B1FF","#6F99ADFF","#FFDC91FF","#EE4C97FF")
theme_full <- theme_bw() + theme(axis.text = element_text(size=14),
                                 axis.title = element_text(size=16),
                                 axis.line = element_blank(),
                                 panel.border = element_rect(size=1.5),
                                 axis.ticks = element_line(size=1.5),
                                 plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
violinLF1_Meth <- ggplot(filter(plotTab, !is.na(Epitype)), aes(x=Epitype, y=Factor1, fill = Epitype)) +
  geom_violin() + geom_point() + scale_fill_manual(values = colList[4:length(colList)]) +
  theme_full + theme(legend.position = "none", axis.title.y = element_text(vjust=-2)) +
  xlab("") + ggtitle("Epigenetic subtypes")
violinLF1_Meth

#Plot the separation of the samples by the first two factors
plotTab <- filter(allFactors, factor %in% c("Factor1","Factor2")) %>%
  spread(key = factor, value = value) %>% mutate(trisomy12 = factor(trisomy12)) %>%
  filter(!is.na(IGHV), !is.na(trisomy12))

pcaLF1 <- ggplot(plotTab, aes(x=Factor1, y=Factor2, color = trisomy12, 
                              shape = IGHV, label = sample)) + 
  geom_point(size=3) +
  scale_shape_manual(values = c(M = 16, U =1)) +
  scale_color_manual(values = c(no = colList[1], yes = colList[2])) +
  theme_full
pcaLF1

#Forest plot showing values of variance explained for each view
R2list <- calculate_variance_explained(MOFAobject)
plotTab <- R2list$r2_per_factor%>%
  as_tibble(rownames = "factor") %>%
  tidyr::pivot_longer(-factor, names_to = "view", values_to = "r2") %>%
  mutate(factor = str_replace(factor,"LF","F"))

#clusters <-cluster_samples(MOFAobject, k=4, factors = 5)
#clusters

#plot_weights_heatmap(MOFAobject)
#plot_weights_scatter(MOFAobject, factors = 1:2)
#run_tsne(MOFAobject, factors = "all", groups = "all")

#Combination of F1 and F2
p <- plot_factors(MOFAobject, 
                  factors = c(1,2), 
                  color_by = "IGHV",
                  shape_by = "trisomy12",
                  dot_size = 2.5,
                  show_missing = T
)

p <- p + 
  geom_hline(yintercept=-1, linetype="dashed") +
  geom_vline(xintercept=(-0.5), linetype="dashed")

print(p)

#F3 profile
plot_factor(MOFAobject, 
            factors = 2, 
            color_by = "Factor3"
)
plot_data_heatmap(MOFAobject, 
                  view = "mRNA",
                  features=15,
                  factor = 3,  
                  denoise = TRUE,
                  cluster_rows = FALSE, cluster_cols = TRUE,
                  show_rownames = TRUE, show_colnames = FALSE,
                  scale = "row"
)

#Assocations with RNAseq batch
allFactors <- get_factors(
  MOFAobject, 
  factors = "all",
  as.data.frame = TRUE
) %>% as_tibble() %>%
  mutate(factor = gsub("LF","F", factor))

batchTab <- filter(allFactors , factor == "Factor3") %>%
  mutate(batch = rna[,match(sample,colnames(rna))]$batch) %>% 
  filter(!is.na(batch)) %>%
  mutate(batch = paste0("batch", batch+1))

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
pval <- car::Anova(lm(value ~ factor(batch), batchTab))$`Pr(>F)`[1]
pval <- formatNum(pval, digits = 2)
pAnno <- bquote(italic("P")~"="~.(pval))
colListNew <- colList[-4]
pL3 <- ggplot(batchTab, aes(x=batch, y = value, col = batch)) +
  geom_boxplot() + 
  ggbeeswarm::geom_beeswarm() + 
  scale_color_manual(values = colListNew) +
  annotate("text", x=Inf, y=Inf, label=pAnno, hjust=1.5, vjust=1.5)+
  theme_full +
  theme(legend.position = "none") +
  ylab("F3 value") + xlab("") + ggtitle("F3 ~ RNAseq batch")
pL3

#Factor 5 profile
plot_top_weights(MOFAobject,
                 view = "mRNA",
                 factor = 5,
                 nfeatures = 5,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
#Factor 6 profile
plot_top_weights(MOFAobject,
                 view = "mRNA",
                 factor = 6,
                 nfeatures = 5,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)
##Factor 7 profile
plot_top_weights(MOFAobject,
                 view = "mRNA",
                 factor = 7,
                 nfeatures = 5,     # Top number of features to highlight
                 scale = T           # Scale weights from -1 to 1
)


#Associations of latent factors to clinical behaviors
#Univariate Cox regression
testTab <- left_join(allFactors, survival, by = c(sample = "patientID"))

#function for cox regression
#' @export
comSurv <- function(response, time, endpoint, scale =FALSE) {
  
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
#for OS
resOS <- filter(testTab, !is.na(OS)) %>%
  group_by(factor) %>%
  do(comSurv(.$value, .$OS, .$died, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
  mutate(Endpoint = "OS")

#for TTT
resTTT <- filter(testTab, !is.na(TTT)) %>%
  group_by(factor) %>%
  do(comSurv(.$value, .$TTT, .$treatedAfter, TRUE)) %>% ungroup() %>%
  arrange(p) %>% mutate(p.adj = p.adjust(p, method = "bonferroni")) %>%
  mutate(Endpoint = "TTT")
resOS
resTTT

#Plot p values and hazard ratios
plotTab <- bind_rows(resOS, resTTT) %>%
  filter(factor %in% c("Factor1","Factor2","Factor4"))

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

haPlot <- ggplot(plotTab, aes(x=factor, y = HR, col = Endpoint, dodge = Endpoint)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width =0.8), 
                aes(ymin = lower, ymax = higher), width = 0.3, size=1) + 
  geom_text(position = position_dodge2(width = 0.8),
            aes(x=as.numeric(as.factor(factor))+0.15,
                label = sprintf("italic(P)~'='~'%s'",
                                formatNum(p))),
            color = "black",size =5, parse = TRUE) +
  xlab("Factor") + ylab("Hazard ratio") +
  scale_y_log10(limits = c(0.5,4)) +
  coord_flip() + theme_full + theme(legend.title = element_blank(),
                                    legend.position = c(0.2,0.1),
                                    legend.background = element_blank(),
                                    legend.key.size = unit(0.5,"cm"),
                                    legend.key.width = unit(0.6,"cm"),
                                    legend.text = element_text(size=rel(1.2))) +
  scale_color_manual(values = c(OS = colList[3], TTT = colList[5])) 

haPlot

#Kaplan-Meiler plots
# Function for Kaplan-Meier plot
#' @export
km <- function(response, time, endpoint, titlePlot = "KM plot", pval = NULL,
               stat = "median", maxTime =NULL, showP = TRUE, showTable = FALSE,
               ylab = "Fraction", xlab = "Time (years)",
               table_ratio = c(0.7,0.3), yLabelAdjust = 0) {
  #function for km plot
  survS <- tibble(time = time,
                  endpoint = endpoint)
  
  if (!is.null(maxTime))
    survS <- mutate(survS, endpoint = ifelse(time > maxTime, FALSE, endpoint),
                    time = ifelse(time > maxTime, maxTime, time))
  
  if (stat == "maxstat") {
    ms <- maxstat.test(Surv(time, endpoint)  ~ response,
                       data = survS,
                       smethod = "LogRank",
                       minprop = 0.2,
                       maxprop = 0.8,
                       alpha = NULL)
    
    survS$group <- factor(ifelse(response >= ms$estimate, "high", "low"))
    p <- comSurv(survS$group, survS$time, survS$endpoint)$p
    
  } else if (stat == "median") {
    med <- median(response, na.rm = TRUE)
    survS$group <- factor(ifelse(response >= med, "high", "low"))
    p <- comSurv(survS$group, survS$time, survS$endpoint)$p
    
  } else if (stat == "binary") {
    survS$group <- factor(response)
    if (nlevels(survS$group) > 2) {
      sdf <- survdiff(Surv(survS$time,survS$endpoint) ~ survS$group)
      p <- 1 - pchisq(sdf$chisq, length(sdf$n) - 1)
    } else {
      p <- comSurv(survS$group, survS$time, survS$endpoint)$p
    }
  }
  
  if (is.null(pval)) {
    if(p< 1e-16) {
      pAnno <- bquote(italic("P")~"< 1e-16")
    } else {
      pval <- formatNum(p, digits = 1)
      pAnno <- bquote(italic("P")~"="~.(pval))
    }
    
  } else {
    pval <- formatNum(pval, digits = 1)
    pAnno <- bquote(italic("P")~"="~.(pval))
  }
  
  if (!showP) pAnno <- ""
  
  colorPal <- colList[1:length(unique(survS$group))]
  p <- ggsurvplot(survfit(Surv(time, endpoint) ~ group, data = survS),
                  data = survS, pval = FALSE,  conf.int = FALSE, palette = colorPal,
                  legend = ifelse(showTable, "none","top"),
                  ylab = "Fraction", xlab = "Time (years)", title = titlePlot,
                  pval.coord = c(0,0.1), risk.table = showTable, legend.labs = sort(unique(survS$group)),
                  ggtheme = theme_half + theme(plot.title = element_text(hjust =0.5),
                                               panel.border = element_blank(),
                                               axis.title.y = element_text(vjust =yLabelAdjust)))
  if (!showTable) {
    p <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size =5)
    return(p)
  } else {
    #construct a gtable
    pp <- p$plot + annotate("text",label=pAnno, x = 0.1, y=0.1, hjust =0, size=5)
    pt <- p$table + ylab("") + xlab("") + theme(plot.title = element_text(hjust=0, size =10))
    p <- plot_grid(pp,pt, rel_heights = table_ratio, nrow =2, align = "v")
    return(p)
  }
}

#KM plot for overall survival (OS)
facList <- sort(filter(resOS, p.adj <=0.02)$factor)
osList <- lapply(facList, function(x) {
  eachTab <- filter(testTab, factor == x) %>%
    select(value, OS, died) %>% filter(!is.na(OS))
  pval <- filter(resOS, factor == x)$p
  km(eachTab$value, eachTab$OS, eachTab$died, sprintf("%s VS Overall survival time", x),
     stat = "maxstat", pval = pval, showTable = TRUE)
})

grid.arrange(grobs = osList, ncol = 2)


#KM plot for time to treatment (TTT)
facList <- sort(filter(resTTT, p.adj <=0.01)$factor)
tttList <- lapply(facList, function(x) {
  eachTab <- filter(testTab, factor == x) %>%
    select(value, TTT, treatedAfter) %>% filter(!is.na(TTT))
  pval <- filter(resTTT, factor == x)$p
  km(eachTab$value, eachTab$TTT, eachTab$treatedAfter, sprintf("%s VS Time to treatment", x), stat = "maxstat",
     maxTime = 7, pval = pval, showTable = TRUE)
})

grid.arrange(grobs = tttList, ncol = 2)

#KM plot for subgroup defined by IGHV status and median latent factor values
groupTab <- filter(testTab, factor == "Factor4", !is.na(value), !is.na(IGHV)) %>%
  mutate(subgroup = ifelse(value > median(value), paste0(IGHV,"_highFactor4"), paste0(IGHV,"_lowFactor4"))) 

plotList <- list()

# TTT
theme_half <- theme_bw() + theme(axis.text = element_text(size=14),
                                 axis.title = element_text(size=16),
                                 axis.line =element_line(size=0.8),
                                 panel.border = element_blank(),
                                 axis.ticks = element_line(size=1.5),
                                 plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
plotList[["TTT"]] <- km(groupTab$subgroup, groupTab$TTT, groupTab$treatedAfter, "Time to treatment", stat = "binary", maxTime = 7, showP = TRUE, showTable = TRUE, yLabelAdjust = -10)


# OS
plotList[["OS"]] <- km(groupTab$subgroup, groupTab$OS, groupTab$died, "Overall survival", stat = "binary", maxTime = 7, showP = TRUE, showTable = TRUE, yLabelAdjust = -10)

grid.arrange(grobs = plotList, ncol = 2)

#KM plot for subgroup defined by median latent factor values of F1 and F4
groupTab <- filter(testTab, factor %in% c("Factor1","Factor4"), !is.na(value)) %>%
  spread(key = factor, value = value) %>%
  mutate(LF1group = ifelse(Factor1 > median(Factor1), "highFactor1", "lowFactor1"),
         LF4group = ifelse(Factor4 > median(Factor4), "highFactor4","lowFactor4")) %>%
  mutate(subgroup = paste0(LF1group, "_",LF4group)) 

plotList1 <- list()
# TTT
plotList1[["TTT"]] <- km(groupTab$subgroup, groupTab$TTT, groupTab$treatedAfter, "Time to treatment", stat = "binary", maxTime = 7, showP = TRUE, showTable = TRUE, yLabelAdjust = -15)

# OS
plotList1[["OS"]] <- km(groupTab$subgroup, groupTab$OS, groupTab$died, "Overall survival", stat = "binary", maxTime = 7, showP = TRUE, showTable = TRUE, yLabelAdjust = -15)

grid.arrange(grobs = plotList1, ncol = 2)

#Multi-variate Cox regression for Factor 4 (F4)
#Prepare data for multi-variate Cox regression
survTab <- survival
facTab <- filter(allFactors, factor == "Factor4")
riskTab <- gene[,c("IGHV","TP53","NOTCH1","del17p","SF3B1")] %>% 
  data.frame() %>% rownames_to_column("patientID") %>%
  mutate(`TP53.del17p` = as.numeric(TP53 | del17p)) %>%
  select(-TP53, -del17p) %>%
  mutate_if(is.numeric, as.factor) %>%
  left_join(select(demographic, patientID, age, sex), by = "patientID") %>%
  mutate(age = age/10) %>%
  mutate(F4 = facTab[match(patientID, facTab$sample),]$value,
         IGHV = factor(IGHV, levels = c(0,1))) %>%
  dplyr::rename(IGHV_mutated = IGHV, Sex_male = sex, Age = age) %>%
  filter(!is.na(Factor4))

#Correlation between F4 and Lymphocyte doubling time
#Univariate test
#Pearson’s correlation test
LDT <- doublingTime %>% mutate(Factor4 = facTab[match(patID, facTab$sample),]$value) %>%
  filter(!is.na(Factor4)) %>%
  mutate(IGHV = as.factor(facTab[match(patID, facTab$sample),]$IGHV))
corRes <- cor.test(log10(LDT$doubling.time), LDT$Factor4)
#Scatter plot of correlations
pval <- formatNum(corRes$p.value, digits = 1, format = "e")
annoN <- sprintf("n = %s", nrow(LDT))
annoP <- bquote(italic("P")~"="~.(pval))
annoCoef <- sprintf("coefficient = %1.2f",corRes$estimate)

corPlot <- ggplot(LDT, aes(x = Factor4, y = doubling.time/30)) + 
  geom_point(fill =colList[5], shape =21, size=3) + 
  geom_smooth(method = "lm", se=FALSE, color = "grey50", linetype ="dashed" ) +
  annotate("text", x = max(LDT$Factor4), y = Inf, label = annoN,
           hjust=1, vjust =1.5, size = 5, parse = FALSE, col= colList[1]) +
  annotate("text", x = max(LDT$Factor4), y = Inf, label = annoP,
           hjust=1, vjust =3.5, size = 5, parse = FALSE, col= colList[1]) +
  annotate("text", x = max(LDT$Factor4), y = Inf, label = annoCoef,
           hjust=1, vjust =5.5, size = 5, parse = FALSE, col= colList[1]) +
  ylab(bquote("doubling time (months)")) + ggtitle("Lymphocyte doubling time") +
  scale_y_log10() +
  theme_full

corPlot

#Consider IGHV status
#IGHV as a covariate
#ANOVA test
LDT <- doublingTime %>% mutate(Factor4 = facTab[match(patID, facTab$sample),]$value,
                               IGHV = as.factor(facTab[match(patID, facTab$sample),]$IGHV)) %>%
  filter(!is.na(Factor4),!is.na(IGHV))

corRes <- car::Anova(lm(log10(doubling.time) ~ Factor4 + IGHV, data = LDT))
corRes

#Only M-CLL
LDT.M <- filter(LDT, IGHV == "M") 
corRes.M <- cor.test(log2(LDT.M$doubling.time), LDT.M$Factor4)
corRes.M

#Only U-CLL
LDT.U <- filter(LDT, IGHV == "U") 
corRes.U <- cor.test(log2(LDT.U$doubling.time), LDT.U$Factor4)
corRes.U

#Scatter plot of correlations, stratified by IGHV
annoM <- sprintf("'M-CLL: n = %s, coefficient = %1.2f,'~italic(P)~'= %s'",nrow(LDT.M), corRes.M$estimate, formatNum(corRes.M$p.value, digits = 1, format = "e"))
annoU <- sprintf("'U-CLL: n = %s, coefficient = %1.2f,'~italic(P)~'= %s'", nrow(LDT.U), corRes.U$estimate,formatNum(corRes.U$p.value, digits = 1, format = "e"))

corPlot.IGHV <- ggplot(LDT, aes(x = Factor4, y = doubling.time/30, fill = IGHV, col = IGHV)) + 
  geom_point(shape=21, size=3, col = "black") + 
  geom_smooth(method = "lm", se=FALSE, linetype ="dashed" ) + 
  annotate("text", x = Inf, y = Inf, label = annoM,
           hjust=1.05, vjust =1.2, size = 4, parse = TRUE, color = colList[1]) +
  annotate("text", x = Inf, y = Inf, label = annoU,
           hjust=1.05, vjust =2.5, size = 4, parse = TRUE, color = colList[2]) +
  ylab("doubling time (months)") + ggtitle("Lymphocyte doubling time") +
  scale_fill_manual(values = c(M = colList[1],U=colList[2])) +
  scale_color_manual(values = c(M = colList[1],U=colList[2])) +
  scale_y_log10() +
  theme_full + theme(legend.position = "none")

corPlot.IGHV

#Correlations in untreated patients only
LDT.untreat <- doublingTime %>% mutate(Factor4 = facTab[match(patID, facTab$sample),]$value, 
                                       pretreat = demographic[match(patID, demographic$patientID),]$pretreat) %>%
  filter(!is.na(Factor4),!is.na(pretreat)) %>% filter(pretreat == 0)

corRes <- cor.test(LDT.untreat$doubling.time, LDT.untreat$Factor4)

annoText <- sprintf("'coefficient = %1.2f, '~italic(P)~'= %s'",corRes$estimate,formatNum(corRes$p.value, digits = 1, format = "e"))
corPlot.untreat <- ggplot(LDT.untreat, aes(x = Factor4, y = doubling.time/30)) + 
  geom_point(fill =colList[5], shape=21, size=3) + 
  geom_smooth(method = "lm", se=FALSE, color = "grey50", linetype ="dashed" ) + 
  annotate("text", x = 2.2, y = Inf, label = annoText,
           hjust=1, vjust =2, size = 4, parse = TRUE, col= colList[1]) +
  ylab("doubling time (months)") + ggtitle(sprintf("Lymphocyte doubling time\n(untreated patients only, n=%s)",nrow(LDT.untreat))) +
  scale_y_log10() +
  theme_full

corPlot.untreat

#Variance explalined for lymphocyte doubling time
onlyIGHV <- summary(lm(log2(doubling.time) ~ IGHV, data = LDT))
onlyF4 <- summary(lm(log2(doubling.time) ~ Factor4, data = LDT))
combined <- summary(lm(log2(doubling.time) ~ Factor4 + IGHV, data = LDT))
plotTab <- tibble(model = c("IGHV only","Factor4 only", "IGHV + Factor4"),
                  R2 = c(onlyIGHV$adj.r.squared, onlyF4$r.squared, combined$adj.r.squared)) %>%
  mutate(model = factor(model, levels = model))
explainedDT <- ggplot(plotTab, aes(x=model, y = R2)) + geom_bar(stat = "identity", aes(fill = model), width = 0.8) +
  coord_flip(expand = FALSE, xlim = c(0.5,3.5)) + theme_half + scale_fill_manual(values = colList) +
  geom_text(aes(x = model, y =0.01, label = model), hjust =0, fontface = "bold", size =5) +
  theme(axis.ticks.y = element_blank(), axis.text.y = element_blank(), legend.position = "none",
        axis.title.y = element_text( size =13)) +
  xlab("Predictors") + ylab("Variance explained") 
explainedDT

#F1 enrichment analysis
#GSEA on positive weights
res.positive <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               factor = 1,
                               sign = "positive"
)

plot_enrichment(res.positive,factor = 1,max.pathways = 10)

enL1 <- plot_enrichment(MOFAobject, res.positive, factor = 1, max.pathways = 10) +
  ylab(bquote("-log"[10]*"(adjusted "*italic("P")~"value)")) +
  ggtitle("Pathways enriched for Factor4") + theme_half

plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 1, 
  max.pathways = 5
)

##GSEA on negative weights
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA", factor = 1,
                               sign = "negative"
)
plot_enrichment(res.negative, factor = 1, max.pathways = 10)

plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 1, 
  max.pathways = 5
)

names(res.positive)

#F6 enrichment analysis
#GSEA on positive weights
res.all <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA", factor = 6,
                               sign = "all"
)
plot_enrichment(res.all, factor = 6, max.pathways = 10, alpha = 0.1)

enL1 <- plot_enrichment(MOFAobject, res.positive, factor = 5, max.pathways = 10) +
  ylab(bquote("-log"[10]*"(adjusted "*italic("P")~"value)")) +
  ggtitle("Pathways enriched for Factor4") + theme_half

plot_enrichment_detailed(
  enrichment.results = res.positive,
  factor = 6, 
  max.pathways = 5
)

##GSEA on negative weights
res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA", factor = 5,
                               sign = "negative"
)
plot_enrichment(res.negative, factor = 5, max.pathways = 10)

plot_enrichment_detailed(
  enrichment.results = res.negative,
  factor = 1, 
  max.pathways = 5
)

res.negative <- run_enrichment(MOFAobject, 
                               feature.sets = reactomeGS, 
                               view = "mRNA",
                               sign = "negative"
)
plot_enrichment_heatmap(res.negative)

# Factor 5 Assocations with CD4, CD8 expression
rna <- estimateSizeFactors(rna)
exprTab <- counts(rna[rowData(rna)$symbol %in% c("CD4","CD8A"),],normalized = TRUE) %>%
  t() %>% as_tibble(rownames = "sample") %>% 
  pivot_longer(-sample, names_to="id", values_to = "count") %>%
  mutate(symbol = rowData(rna)[id,]$symbol)


facTab <- filter(allFactors , factor == "Factor5")
plotTab <- left_join(exprTab, facTab, by = "sample")

pAnno <- bquote(italic("P")~"<"~.(10e-13))

theme_half <- theme_bw() + theme(axis.text = element_text(size=14),
                                 axis.title = element_text(size=16),
                                 axis.line =element_line(size=0.8),
                                 panel.border = element_blank(),
                                 axis.ticks = element_line(size=1.5),
                                 plot.title = element_text(size = 16, hjust =0.5, face="bold"),
                                 panel.grid.major = element_blank(),
                                 panel.grid.minor = element_blank())
pL5 <- ggplot(plotTab, aes(y=log2(count), x= value, col = symbol)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  theme_half +
  scale_color_manual(values = colList) +
  annotate("text", x=-Inf, y=Inf, label=pAnno, hjust=-0.5, vjust=5, size=5)+
  xlab("F5 value") + ylab(bquote("log"[2]*"(RNAseq count)")) +
  ggtitle("T cell marker expressions ~ F5") +
  theme(legend.position = c(0.8,0.15),
        legend.text = element_text(size=15),
        legend.title = element_text(size=15)) +
  ylim(0,13)
plot(pL5)

#Factor 6
exprTab <- counts(rna[rowData(rna)$symbol %in% c("SOD1","GPX4"),],normalized = TRUE) %>%
  t() %>% as_tibble(rownames = "sample") %>% 
  pivot_longer(-sample, names_to="id", values_to = "count") %>%
  mutate(symbol = rowData(rna)[id,]$symbol)


facTab <- filter(allFactors , factor == "F6")
plotTab <- left_join(exprTab, facTab, by = "sample")

p1 <- cor.test(~ value + log2(count), filter(plotTab, symbol == "SOD1"))
p2 <- cor.test(~ value + log2(count), filter(plotTab, symbol == "GPX4"))

pAnno <- bquote(italic("P")~"<"~.(10e-13))
corL6 <- ggplot(plotTab, aes(y=log2(count), x= value, col = symbol)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  theme_half +
  scale_color_manual(values = colList) +
  annotate("text", x=-Inf, y=Inf, label=pAnno, hjust=-0.5, vjust=5, size=5)+
  xlab("F6 value") + ylab(bquote("log"[2]*"(RNAseq count)")) +
  ggtitle("SOD1 and GPX4 ~ F6") +
  theme(legend.position = c(0.8,0.15),
        legend.text = element_text(size=15),
        legend.title  = element_text(size=15))

#source table
outTab <- tibble(pathway = rownames(fsea.results$pval),
                 PValue = fsea.results$pval[,1],
                 adj.PValue = fsea.results$pval.adj[,1]) %>%
  arrange(PValue)
exprTab <- counts(rna[rowData(rna)$symbol %in% c("SOD1","GPX4"),],normalized = TRUE) %>%
  t() %>% as_tibble(rownames = "sample") %>% 
  pivot_longer(-sample, names_to="id", values_to = "count") %>%
  mutate(symbol = rowData(rna)[id,]$symbol)


facTab <- filter(allFactors , factor == "F6")
plotTab <- left_join(exprTab, facTab, by = "sample")

p1 <- cor.test(~ value + log2(count), filter(plotTab, symbol == "SOD1"))
p2 <- cor.test(~ value + log2(count), filter(plotTab, symbol == "GPX4"))

pAnno <- bquote(italic("P")~"<"~.(10e-13))
corL6 <- ggplot(plotTab, aes(y=log2(count), x= value, col = symbol)) +
  geom_point() + geom_smooth(method = "lm", se = FALSE) +
  theme_half +
  scale_color_manual(values = colList) +
  annotate("text", x=-Inf, y=Inf, label=pAnno, hjust=-0.5, vjust=5, size=5)+
  xlab("F6 value") + ylab(bquote("log"[2]*"(RNAseq count)")) +
  ggtitle("SOD1 and GPX4 ~ F6") +
  theme(legend.position = c(0.8,0.15),
        legend.text = element_text(size=15),
        legend.title  = element_text(size=15))

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

#Factor 7
fsea.results <- MOFA2::run_enrichment(MOFAobject,view = "mRNA", factor = 7 , alpha = 0.2, feature.sets = reactomeGS)
enL7 <- MOFA2::plot_enrichment(MOFAobject, fsea.results, factor = 7, max.pathways = 5) +
  ylab(bquote("-log"[10]*"(adjusted "*italic("P")~"value)")) +
  ggtitle("Pathways enriched for F7") + theme_half
