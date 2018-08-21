# Big script of preliminary analysis

rm(list = ls())
setwd('/Users/ereznik/Documents/randomprojects/vascularmimicry/analysis/')
library(ggplot2)
library(pls)
library(reshape2)
library(caret)

load('../data/imported_VascularMet.RData')

# Generate log versions
logsn = log2(snmet)
logcl = log2(clmet)

#######
# PCA
#######
pc_clmet = prcomp(t(clmet),center = TRUE,scale = TRUE)
pc_clmet_res = as.data.frame( pc_clmet$rotation )
pc_clmet_res$CellLine = cld[rownames(pc_clmet_res),'Cell.line']
pc_clmet_res$Cultivation = cld[rownames(pc_clmet_res),'Cultivation']

ggplot(pc_clmet_res,aes(PC1,PC2,color = CellLine,shape = Cultivation)) + geom_point(size = 4,alpha = 0.4) + 
  facet_wrap(~Cultivation)

pc_snmet = prcomp(t(snmet),center = TRUE,scale = TRUE)
pc_snmet_res = as.data.frame( pc_snmet$rotation )
pc_snmet_res$CellLine = snd[rownames(pc_snmet_res),'Cell.line']
pc_snmet_res$Cultivation = snd[rownames(pc_snmet_res),'cultivation']

ggplot(pc_snmet_res,aes(PC1,PC2,color = CellLine,shape = Cultivation)) + geom_point(size = 4,alpha = 0.4) + 
  facet_wrap(~Cultivation) + theme_minimal()

#####
# Fit a linear model with PLS
#####

# Start with the cell lines
conversions = c('K8484' = 0, 'TB32048' = 1, '2' = 0, '3D (matrigel)' = 1)
clresponse = cld[,c('Cell.line','Cultivation')]
clresponse$Cell.line = conversions[clresponse$Cell.line]
clresponse$Cultivation = conversions[clresponse$Cultivation]
clresponse = as.matrix(clresponse)

#clplsr = plsda(x = clresponse, y = as.matrix(clmet))
#clpls =  plsr(as.matrix(clmet) ~ clresponse,scale = TRUE)

#####
# Fit a linear model with lm
#####
rescl = data.frame()
for (metname in colnames(logcl)){
  # Fit the model
  templm = lm(logcl[,metname] ~ cld$Cell.line + cld$Cultivation)
  rescl[metname,'TB32Effect'] = summary(templm)$coefficients['cld$Cell.lineTB32048','Estimate']
  rescl[metname,'TB32P'] = summary(templm)$coefficients['cld$Cell.lineTB32048','Pr(>|t|)']
  
  rescl[metname,'MatrigelEffect'] = summary(templm)$coefficients['cld$Cultivation3D (matrigel)','Estimate']
  rescl[metname,'MatrigelP'] = summary(templm)$coefficients['cld$Cultivation3D (matrigel)','Pr(>|t|)']
}

rescl$TB32Q = p.adjust(rescl$TB32P,method = 'BH')
rescl$MatrigelQ = p.adjust(rescl$MatrigelP,method = 'BH')

rescl$TB32Sig = 0
rescl[which(rescl$TB32Q < 0.1),'TB32Sig'] = 1

rescl$MatrigelSig = 0
rescl[which(rescl$MatrigelQ < 0.1),'MatrigelSig'] = 1

# Make plots
sigcl = which(rescl$TB32Sig == 1 | rescl$MatrigelSig == 1)
for (ii in sigcl){
  metname = rownames(rescl)[ii]
  pdata = data.frame( logcl[,metname], cld[rownames(logcl),] )
  colnames(pdata) = c('Metabolite',colnames(cld))
  ggplot(pdata,aes(Cell.line,Metabolite,color = Cultivation,fill = Cultivation)) + geom_boxplot(alpha = 0.1) + 
    theme_minimal(base_size = 16) + geom_point(position = position_jitterdodge(jitter.width = .1)) + 
    ggtitle(metname) + 
    ggsave(paste0('../results/lm/',metname,'.pdf'),height =7,width = 7)
}
