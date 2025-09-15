setwd("D:/GitHub/MR/")
library(dplyr)
library(tidyr)
library(TwoSampleMR)
library(data.table)
# IV
JAK1_stroke <- read.csv("JAK1_SNP_asIV.csv")
JAK2_stroke <- read.csv("JAK2_SNP_asIV.csv")
JAK3_stroke <- read.csv("JAK3_SNP_asIV.csv")
TYK2_stroke <- read.csv("TYK2_SNP_asIV.csv")
JAK1_CAD <- read.csv("JAK1_SNP_asIV_CAD.csv")
JAK2_CAD <- read.csv("JAK2_SNP_asIV_CAD.csv")
JAK3_CAD <- read.csv("JAK3_SNP_asIV_CAD.csv")
TYK2_CAD <- read.csv("TYK2_SNP_asIV_CAD.csv")
# outcome 
stroke <- data.table::fread("D:/Ubuntu/rootfs/home/songqy9/downloaddir/smr-1.3.1-linux-x86_64/data/stroke_mygwas.ma")
CAD <- data.table::fread("D:/Ubuntu/rootfs/home/songqy9/downloaddir/smr-1.3.1-linux-x86_64/data/CAD_mygwas.ma")

# extract SNP from outcome---------local---------- 
# outcome 
stroke <- data.table::fread("D:/Ubuntu/rootfs/home/songqy9/downloaddir/smr-1.3.1-linux-x86_64/data/stroke_mygwas.ma")
CAD <- data.table::fread("D:/Ubuntu/rootfs/home/songqy9/downloaddir/smr-1.3.1-linux-x86_64/data/CAD_mygwas.ma")
# stroke
JAK1_stroke_outcome <- stroke[(stroke$SNP %in% JAK1_stroke$SNP),]
JAK2_stroke_outcome <- stroke[(stroke$SNP %in% JAK2_stroke$SNP),]
JAK3_stroke_outcome <- stroke[(stroke$SNP %in% JAK3_stroke$SNP),]
TYK2_stroke_outcome <- stroke[(stroke$SNP %in% TYK2_stroke$SNP),]
write.csv(JAK1_stroke_outcome,"JAK1_stroke_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK2_stroke_outcome,"JAK2_stroke_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK3_stroke_outcome,"JAK3_stroke_IV_outcome.csv",row.names = F,quote = F)
write.csv(TYK2_stroke_outcome,"TYK2_stroke_IV_outcome.csv",row.names = F,quote = F)
# CAD
JAK1_CAD_outcome <- CAD[(CAD$SNP %in% JAK1_CAD$SNP),]
JAK2_CAD_outcome <- CAD[(CAD$SNP %in% JAK2_CAD$SNP),]
JAK3_CAD_outcome <- CAD[(CAD$SNP %in% JAK3_CAD$SNP),]
TYK2_CAD_outcome <- CAD[(CAD$SNP %in% TYK2_CAD$SNP),]
write.csv(JAK1_CAD_outcome,"JAK1_CAD_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK2_CAD_outcome,"JAK2_CAD_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK3_CAD_outcome,"JAK3_CAD_IV_outcome.csv",row.names = F,quote = F)
write.csv(TYK2_CAD_outcome,"TYK2_CAD_IV_outcome.csv",row.names = F,quote = F)
# extract SNP from outcome---------online--------
# stroke
JAK1_stroke_outcome <- extract_outcome_data(snps = JAK1_stroke$SNP, outcomes = 'ukb-b-8714')
JAK2_stroke_outcome <- extract_outcome_data(snps = JAK2_stroke$SNP, outcomes = 'ukb-b-8714')
JAK3_stroke_outcome <- extract_outcome_data(snps = JAK3_stroke$SNP, outcomes = 'ukb-b-8714')
write.csv(JAK1_stroke_outcome,"JAK1_stroke_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK2_stroke_outcome,"JAK2_stroke_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK3_stroke_outcome,"JAK3_stroke_IV_outcome.csv",row.names = F,quote = F)
# CAD
JAK1_CAD_outcome <- extract_outcome_data(snps = JAK1_CAD$SNP, outcomes = 'ukb-b-15748')
JAK2_CAD_outcome <- extract_outcome_data(snps = JAK2_CAD$SNP, outcomes = 'ukb-b-15748')
JAK3_CAD_outcome <- extract_outcome_data(snps = JAK3_CAD$SNP, outcomes = 'ukb-b-15748')
write.csv(JAK1_CAD_outcome,"JAK1_CAD_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK2_CAD_outcome,"JAK2_CAD_IV_outcome.csv",row.names = F,quote = F)
write.csv(JAK3_CAD_outcome,"JAK3_CAD_IV_outcome.csv",row.names = F,quote = F)
