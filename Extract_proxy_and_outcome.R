setwd("D:/GitHub/MR/")
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
# proxy SNP
JAK2 <- data.table::fread("JAK2_myquery.txt") %>% .[,c("SNP","Chr","BP","Freq","b","SE","p")]
JAK2_2 <- JAK2[,c("SNP","Chr","BP")]
names(JAK2)<-c("proxy_snp.outcome","proxy_Chr","proxy_BP","proxy_Freq","proxy_b","proxy_SE","proxy_p")
names(JAK2_2)<-c("SNP","Chr_exposure","BP_exposure")
TYK2 <- data.table::fread("TYK2_myquery.txt") %>% .[,c("SNP","Chr","BP","Freq","b","SE","p")]
TYK2_2 <- TYK2[,c("SNP","Chr","BP")]
names(TYK2)<-c("proxy_snp.outcome","proxy_Chr","proxy_BP","proxy_Freq","proxy_b","proxy_SE","proxy_p")
names(TYK2_2)<-c("SNP","Chr_exposure","BP_exposure")
# JAK1_stroke <- read.csv("AK1_SNP_asIV.csv")
JAK2_stroke <- read.csv("JAK2_SNP_asIV.csv") %>% .[which(.$`proxy.outcome`==TRUE),] %>% dplyr::mutate(trait = rep("stroke",nrow(.)))
TYK2_stroke <- read.csv("TYK2_SNP_asIV.csv") %>% .[which(.$`proxy.outcome`==TRUE),] %>% dplyr::mutate(trait = rep("stroke",nrow(.)))
# JAK1_CAD <- read.csv("JAK1_SNP_asIV_CAD.csv")
JAK2_CAD <- read.csv("JAK2_SNP_asIV_CAD.csv") %>% .[which(.$`proxy.outcome`==TRUE),] %>% dplyr::mutate(trait = rep("CAD",nrow(.)))
TYK2_CAD <- read.csv("TYK2_SNP_asIV_CAD.csv") %>% .[which(.$`proxy.outcome`==TRUE),] %>% dplyr::mutate(trait = rep("CAD",nrow(.)))
# 
proxy_JAK2 <- rbind(JAK2_stroke,JAK2_CAD) %>% 
  .[,c("trait","SNP","proxy_snp.outcome","chr","pos","effect_allele.exposure","other_allele.exposure","proxy_a1.outcome","proxy_a2.outcome","eaf.outcome",
       "beta.exposure","beta.outcome","se.exposure","se.outcome","pval.exposure","pval.outcome")] %>% 
  merge(.,JAK2,by = "proxy_snp.outcome") %>% 
  merge(.,JAK2_2,by = "SNP")
proxy_TYK2 <- rbind(TYK2_stroke,TYK2_CAD) %>% 
  .[,c("trait","SNP","proxy_snp.outcome","chr","pos","effect_allele.exposure","other_allele.exposure","proxy_a1.outcome","proxy_a2.outcome","eaf.outcome",
       "beta.exposure","beta.outcome","se.exposure","se.outcome","pval.exposure","pval.outcome")] %>% 
  merge(.,TYK2,by = "proxy_snp.outcome") %>% 
  merge(.,TYK2_2,by = "SNP")
proxy <- rbind(proxy_JAK2,proxy_TYK2)
proxy_save <- proxy[,c("trait","SNP","Chr_exposure","BP_exposure","effect_allele.exposure","other_allele.exposure","beta.exposure","se.exposure","pval.exposure",
                       "proxy_snp.outcome","chr","pos","proxy_a1.outcome","proxy_a2.outcome","beta.outcome","se.outcome","pval.outcome")] %>% 
  dplyr::arrange(.,desc(trait))
names(proxy_save)<-c("Trait","SNP","Chr_exposure","BP_exposure","effect_allele","other_allele","beta.exposure","se.exposure","pval.exposure",
                     "proxy_snp.outcome","proxy_chr","proxy_pos","proxy_a1.outcome","proxy_a2.outcome","beta.outcome","se.outcome","pval.outcome")
write.csv(proxy_save,"tidy_reasult/proxy.csv",row.names = F,quote = F)

