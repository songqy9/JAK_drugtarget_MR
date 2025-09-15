setwd("D:/GitHub/MR/")
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
result <- read.csv("result_CADukb-b-15748.csv",header = T)
result_smr_jak1 <- read.delim("JAK1_CAD_Msmr.msmr")
result_smr_jak2 <- read.delim("JAK2_CAD_Msmr.msmr")
result_smr_jak3 <- read.delim("JAK3_CAD_Msmr.msmr")
result_smr_TYK2 <- read.delim("TYK2_Msmr.msmr")
result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3,result_smr_TYK2)
mr_presso <- read.csv("result_presso_CADukb-b-15748.csv",header = T) %>% dplyr::select(2)
# IVW
IVW_sensitivity <- dplyr::select(result,exposure,method,pleiotraphy_intercept,pleiotraphy_p,pleiotraphy_se,heterogenicity_Q,heterogenicity_p)
exposure_method <- IVW_sensitivity[c(2,8,14,20),c(1,2)]
sensit <- IVW_sensitivity[c(1,7,13,19),3:5]
hetero <- IVW_sensitivity[c(2,8,14,20),7]
sensitivity_table <- cbind(exposure_method,sensit,hetero,mr_presso) %>% 
  dplyr::mutate(pleiotraphy_p=round(pleiotraphy_p,digits = 3),hetero=round(hetero,digits = 3)) %>% 
  dplyr::rename(`Exposure traits`=exposure,
                `MR method`=method,
                `MR-egger intercept`=pleiotraphy_intercept,
                `MR-egger intercept standard error`=pleiotraphy_se,
                `MR-egger intercept p Value` = pleiotraphy_p,
                `Cochran Q test p Value` = hetero,
                `p value for MR-PRESSO Global test`=MR.PRESSO.results.Global.Test.Pvalue)
write.csv(sensitivity_table,"tidy_reasult/CAD_IVW_sensitivitytable.csv",row.names = F,quote = F)
# smr 手动改
write.csv(result_smr,"tidy_reasult/CAD_smr_sensitivitytable.csv",row.names = F,quote = F)

