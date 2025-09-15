setwd("D:/GitHub/MR/")
library(openxlsx)
library(dplyr)
library(tidyr)
library(readr)
Origin_query <- read.table("TYK2_myquery.txt",header = T)
TYK2_myquery_stroke <- read.csv("TYK2_SNP_asIV.csv") 
TYK2_myquery_CAD <- read.csv("TYK2_SNP_asIV_CAD.csv")
# 
TYK2_myquery_stroke2 <- Origin_query[c(Origin_query$SNP %in% TYK2_myquery_stroke$SNP),]
TYK2_myquery_CAD2 <- Origin_query[c(Origin_query$SNP %in% TYK2_myquery_CAD$SNP),]
write.table(TYK2_myquery_stroke2,"TYK2_Stroke_myquery.txt",row.names = F,quote = F)
write.table(TYK2_myquery_CAD2,"TYK2_CAD_myquery.txt",row.names = F,quote = F)
# SNP Chr BP effect_allele other_allele Freq Probe Probe_Chr Probe_bp Gene Orientation b SE p
