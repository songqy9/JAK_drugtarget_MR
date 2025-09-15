setwd("D:/GitHub/MR/药靶")
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
library(VariantAnnotation)
library(gwasvcf)
library(gwasglue)
#beta = log(OR)

stroke_exposure2 <- VariantAnnotation::readVcf("ukb-b-8714.vcf.gz") %>% gwasvcf_to_TwoSampleMR(.)
stroke <- dplyr::select(stroke_exposure2,SNP,effect_allele.exposure,other_allele.exposure,eaf.exposure,beta.exposure,se.exposure,pval.exposure,ncase.exposure) %>% 
  dplyr::rename(A1=effect_allele.exposure,A2=other_allele.exposure,freq=eaf.exposure,b=beta.exposure,se=se.exposure,p=pval.exposure,n=ncase.exposure)
write.table(stroke,"stroke_mygwas.ma",row.names = F,quote = F)# 获取smr所需的ma文件