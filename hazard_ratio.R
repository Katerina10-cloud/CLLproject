load("CLLdata/gene_untr.RData")
load("CLLproject/Survival_analysis/rna2.RData")
load("CLLproject/Survival_analysis/uniRes.os.RData")

# Select Rows by list of names
expMat <- expMat[c("LPL", "LRP1", "LRP5", "LRP6", "LRP8", "LDLR", "VLDLR", "LILRB4", "AEBP1", "LAMA5", "COL9A3", "MTSS1", "BIN1", "MYO1E", "COTL1", "LSP1", "PIGR", "IL2RA", "CXCR5", "CD44", "IFNGR1", "CD6", "ITGB1", "CD84", "SELL", "CEMIP2", "RHOB", "TGFBR2", "ALOX5", "RGCC", "LILRB4", "JUND", "KLF2", "ZNF331", "ARID5B", "EGR1", "CCDC88A", "NR4A3", "ID3", "SREBF2", "ZBTB10", "FOXN3", "ZNF395", "NFATC1", "ZEB2", "ZBTB18", "KLF9", "TGIF1", "NCOA1", "BCL6", "ZNF83", "HHEX", "FOSL2", "ZNF529", "SPI1", "ZNF532", "ZNF91", "ZNF831", "PBX3", "ZNF821", "ZNF432", "FOXO1", "ZNF264", "ZNF350", "RXRA", "MAFF", "HIVEP3", "MAFG"),]

survTab_untr <- survival %>% 
  filter(patientID %in% colnames(riskTab2_untr))

riskTab <- t(riskTab)
riskTab <- as.data.frame(riskTab)
colnames(riskTab) <- riskTab[1, ] #convert first row to header
riskTab <- riskTab[-1,]

rna2 <- as.matrix(rna2)
expMat <- rna2
expMat <- as.matrix(expMat)

riskTab <- riskTab %>% mutate_if(is.character, as.numeric)
riskTab <- as.matrix(riskTab)

#Univariate cox regression
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

#Univariate cox regression test
uniRes.ttt <- lapply(rownames(expMat), function(n) {
  testTab <- mutate(survTab, expr = expMat[n, patientID])
  com(testTab$expr, testTab$TTT, testTab$treatedAfter, TRUE) 
}) %>% bind_rows() %>% mutate(p.adj = p.adjust(p, method = "BH")) %>% arrange(p) %>% mutate(name = rownames(expMat)) %>% mutate(outcome = "TTT")

uniRes.os <- lapply(rownames(expMat), function(n) {
  testTab <- mutate(survTab, expr = expMat[n, patientID])
  com(testTab$expr, testTab$OS, testTab$died, TRUE) 
}) %>% bind_rows() %>% mutate(p.adj = p.adjust(p, method = "BH")) %>% arrange(p) %>% mutate(name = rownames(expMat)) %>% mutate(outcome = "OS")

uniRes <- bind_rows(uniRes.ttt, uniRes.os) %>%
  mutate(p.adj = p.adjust(p, method = "BH"))

save(uniRes, file="uniRes_untr.Rdata")

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
#Forest plot of selected gene markers
plotTab <- uniRes %>%
  filter(name %in% c("FOXO1", "ZNF264", "ZNF350", "RXRA", "MAFF", "HIVEP3", "MAFG"))

haPlot <- ggplot(plotTab, aes(x=name, y = HR, col = outcome, dodge = outcome)) + 
  geom_hline(yintercept = 1, linetype = "dotted") +
  geom_point(position = position_dodge(width=0.8)) +
  geom_errorbar(position = position_dodge(width =0.8), 
                aes(ymin = lower, ymax = higher), width = 0.3, size=1) + 
  geom_text(position = position_dodge2(width = 0.8),
            aes(x=as.numeric(as.factor(name))+0.15,
                label = sprintf("italic(P)~'='~'%s'",
                                formatNum(p))),
            color = "black",size =3, parse = TRUE) +
  xlab("name") + ylab("HR") +
  scale_y_log10(limits = c(0.5,4)) +
  coord_flip() + theme_full + theme(legend.title = element_blank(),
                                    legend.position = c(0.9,0.6),
                                    legend.background = element_blank(),
                                    legend.key.size = unit(0.5,"cm"),
                                    legend.key.width = unit(0.6,"cm"),
                                    legend.text = element_text(size=rel(1))) +
  scale_color_manual(values = c(OS = colList[3], TTT = colList[5], TFT = colList[4] )) +
  ggtitle("HR for gene of interest")
haPlot

#p <- ggplot(plotTab, aes(x=name, y=HR, ymin=lower, ymax=higher)) + 
  geom_linerange(size=8, colour="#a6d8f0") +
  geom_hline(aes(x=0, yintercept=1), lty=1) +
  geom_point(size=3, shape=21, fill="#008fd5", colour = "white", stroke = 1) +
  scale_y_continuous(limits = c(0.5, 2)) +
  coord_flip() +
  ggtitle("HR for gene of interest TTT") +
  theme_minimal()
#p


