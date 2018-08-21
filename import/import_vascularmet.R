# Script to import Miller lab metabolomics

rm(list = ls())
setwd('/Users/ereznik/Documents/randomprojects/vascularmimicry/import/')
library(XLConnect)
source('sourcefiles.R')

sn = readWorksheetFromFile(snfile,sheet = 'Imputed Data')
cl = readWorksheetFromFile(cellfile,sheet = 'Imputed Data')

rownames(sn) = sn$Sample.Identification
rownames(cl) = cl$Sample.Identification

snmet = sn[,7:dim(sn)[2]]
clmet = cl[,8:dim(cl)[2]]

snd = sn[,1:6]
cld = cl[,1:7]

save(snmet,clmet,snd,cld,file = '../data/imported_VascularMet.RData')
