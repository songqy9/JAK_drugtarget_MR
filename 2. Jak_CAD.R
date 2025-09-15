setwd("D:/GitHub/MR/")
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
result_smr_jak1 <- read.delim("JAK1_CAD_Msmr.msmr")
result_smr_jak2 <- read.delim("JAK2_CAD_Msmr.msmr")
result_smr_jak3 <- read.delim("JAK3_CAD_Msmr.msmr")
result_smr_TYK2 <- read.delim("TYK2_CAD_Msmr.msmr")
result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3,result_smr_TYK2)
write.csv(result_smr,"tidy_reasult/CAD_smr_sensitivitytable.csv",row.names = F,quote = F)
# MR_exposure
exposure_read <- function(name){
  JAK = read_exposure_data(filename = paste0(name,"_5e-0x_freq_Fstat.csv"),sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",eaf_col = "eaf",samplesize_col = "samplesize",clump = FALSE)
  JAK = clump_data(JAK,clump_kb = 100,
                   clump_r2 = 0.3,
                   pop = "EUR")
  return(JAK)
}
# JAK1_exposure <- exposure_read("JAK1")
# JAK2_exposure <- exposure_read("JAK2")
# JAK3_exposure <- exposure_read("JAK3")
# TYK2_exposure <- exposure_read("TYK2")
# saveRDS(JAK1_exposure,"JAK1_clump.rds")
# saveRDS(JAK2_exposure,"JAK2_clump.rds")
# saveRDS(JAK3_exposure,"JAK3_clump.rds")
# saveRDS(TYK2_exposure,"TYK2_clump.rds")
JAK1_exposure<-readRDS("JAK1_clump.rds")
JAK2_exposure<-readRDS("JAK2_clump.rds")
JAK3_exposure<-readRDS("JAK3_clump.rds")
TYK2_exposure<-readRDS("TYK2_clump.rds")

# MR_harmonise_dat and get IV
harmoize_dat <- function(JAK,outcome_dataset,name){
  outcome = extract_outcome_data(snps = JAK$SNP,
                                 outcomes = outcome_dataset)
  dat <- harmonise_data(JAK,outcome)
  dat2 <- dplyr::filter(dat,mr_keep=TRUE)
  write.csv(dat2,paste0(name,"_SNP_asIV_CAD.csv"),row.names = F,quote = F)
  return(dat)
}
# MR result
mr_result <- function(JAK_dat,JAK_smr,name){
  JAK=mr(JAK_dat, method_list = c("mr_egger_regression","mr_ivw","mr_weighted_median","mr_weighted_mode"))
  JAK2 = mr.raps::mr.raps(b_exp = JAK_dat$beta.exposure,
                          b_out = JAK_dat$beta.outcome,
                          se_exp = JAK_dat$se.exposure,
                          se_out = JAK_dat$se.outcome,
                          loss.function = "tukey",
                          over.dispersion = FALSE)#如果egger有多效性 就要设置TRUE
  raps_row = data.frame(id.exposure = JAK$id.exposure[1], #加上raps的结果
                        id.outcome = JAK$id.outcome[1],
                        outcome = JAK$outcome[1],
                        exposure = JAK$exposure[1],
                        method = "Robust adjusted profile score",
                        nsnp = JAK$nsnp[1],
                        b = JAK2[["beta.hat"]],
                        se = JAK2[["beta.se"]],
                        pval = JAK2[["beta.p.value"]])
  smr_row = data.frame(id.exposure = JAK$id.exposure[1],# 加上smr的结果
                       id.outcome = JAK$id.outcome[1],
                       outcome = JAK$outcome[1],
                       exposure = JAK$exposure[1],
                       method = "Summary-data-based Mendelian Randomization",
                       nsnp = JAK$nsnp[1],
                       b = JAK_smr$b_SMR,
                       se = JAK_smr$se_SMR,
                       pval = JAK_smr$p_SMR)
  JAK = rbind(JAK,raps_row,smr_row)
  JAK=generate_odds_ratios(JAK)
  JAK_plei = mr_pleiotropy_test(JAK_dat)
  JAK_hetero = mr_heterogeneity(JAK_dat)
  # 整理mr的result
  out = data.frame(exposure = rep(name,6),
                   outcome = JAK$outcome,
                   method = JAK$method,
                   nsnp = JAK$nsnp,
                   beta = JAK$b,
                   se = JAK$se,
                   p = JAK$pval,
                   OR = JAK$or,
                   OR_lci95 = JAK$or_lci95,
                   OR_uci95 = JAK$or_uci95,
                   beta_lci = JAK$lo_ci,
                   beta_uci = JAK$up_ci,
                   pleiotraphy_intercept = c(JAK_plei$egger_intercept,rep(NA,5)),
                   pleiotraphy_se = c(JAK_plei$se,rep(NA,5)),
                   pleiotraphy_p = c(JAK_plei$pval,rep(NA,5)),
                   heterogenicity_Q = c(JAK_hetero$Q,rep(NA,4)),
                   heterogenicity_p = c(JAK_hetero$Q_pval,rep(NA,4)))
  return(out)
}
out_dat <- data.frame(exposure = NULL,
                      outcome = NULL,
                      method = NULL,
                      nsnp = NULL,
                      beta = NULL,
                      se = NULL,
                      p = NULL,
                      OR = NULL,
                      OR_lci95 = NULL,
                      OR_uci95 = NULL,
                      beta_lci = NULL,
                      beta_uci = NULL,
                      pleiotraphy_intercept = NULL,
                      pleiotraphy_se = NULL,
                      pleiotraphy_p = NULL,
                      heterogenicity_Q = NULL,
                      heterogenicity_p = NULL)
JAK_result <- function(JAK1_exposure,JAK2_exposure,JAK3_exposure,TYK2_exposure,result_smr_jak1,result_smr_jak2,result_smr_jak3,result_smr_TYK2,outcome_dataset,out_dat){
  JAK1_dat = harmoize_dat(JAK1_exposure,outcome_dataset,"JAK1")
  JAK2_dat = harmoize_dat(JAK2_exposure,outcome_dataset,"JAK2")
  JAK3_dat = harmoize_dat(JAK3_exposure,outcome_dataset,"JAK3")
  TYK2_dat = harmoize_dat(TYK2_exposure,outcome_dataset,"TYK2")
  JAK1_mr = mr_result(JAK1_dat,result_smr_jak1,"JAK1") 
  JAK2_mr = mr_result(JAK2_dat,result_smr_jak2,"JAK2")
  JAK3_mr = mr_result(JAK3_dat,result_smr_jak3,"JAK3")
  TYK2_mr = mr_result(TYK2_dat,result_smr_TYK2,"TYK2")
  out_dat = rbind(JAK1_mr,JAK2_mr,JAK3_mr,TYK2_mr)
  write.csv(out_dat,paste0("result_CAD",outcome_dataset,".csv"),row.names = F,quote = F)
  return(out_dat)
}
# mr presso
presso <- function(JAK1_exposure,JAK2_exposure,JAK3_exposure,TYK2_exposure,outcome_dataset,out_dat){
  JAK1_dat = harmoize_dat(JAK1_exposure,"ukb-b-15748","JAK1")
  JAK2_dat = harmoize_dat(JAK2_exposure,"ukb-b-15748","JAK2")
  JAK3_dat = harmoize_dat(JAK3_exposure,"ukb-b-15748","JAK3")
  TYK2_dat = harmoize_dat(TYK2_exposure,"ukb-b-15748","TYK2")
  res_presso_JAK1 = run_mr_presso(JAK1_dat)
  res_presso_JAK2 = run_mr_presso(JAK2_dat)
  res_presso_JAK3 = run_mr_presso(JAK3_dat)
  res_presso_TYK2 = run_mr_presso(TYK2_dat)
  out = data.frame(exposure = c("JAK1","JAK2","JAK3","TYK2"),
                  `MR-PRESSO results Global Test Pvalue`=c(res_presso_JAK1[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                                           res_presso_JAK2[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                                           res_presso_JAK3[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]],
                                                           res_presso_TYK2[[1]][["MR-PRESSO results"]][["Global Test"]][["Pvalue"]]))
  write.csv(out,paste0("result_presso_CAD",outcome_dataset,".csv"),row.names = F,quote = F)
  return(out)
}
# Scatter plot add smr and raps
harmoize_dat2 <- function(JAK,outcome_dataset){
  outcome = extract_outcome_data(snps = JAK$SNP,
                                 outcomes = outcome_dataset)
  dat <- harmonise_data(JAK,outcome)
  # dat2 <- dplyr::filter(dat,mr_keep=TRUE)
  # write.csv(dat2,paste0(name,"_SNP_asIV.csv"),row.names = F,quote = F)
  return(dat)
}
my_scatterplot_add_smr_raps <- function(JAK_dat,JAK_smr,name){
  JAK_dat = dplyr::mutate(JAK_dat,outcome=rep("CAD",nrow(JAK_dat)))
  JAK_dat = dplyr::mutate(JAK_dat,exposure=rep(name,nrow(JAK_dat)))
  JAK=mr(JAK_dat, method_list = c("mr_egger_regression","mr_ivw","mr_weighted_median","mr_weighted_mode"))
  JAK2 = mr.raps::mr.raps(b_exp = JAK_dat$beta.exposure,
                          b_out = JAK_dat$beta.outcome,
                          se_exp = JAK_dat$se.exposure,
                          se_out = JAK_dat$se.outcome,
                          loss.function = "tukey",
                          over.dispersion = FALSE)#如果egger有多效性 就要设置TRUE
  raps_row = data.frame(id.exposure = JAK$id.exposure[1],
                        id.outcome = JAK$id.outcome[1],
                        outcome = JAK$outcome[1],
                        exposure = JAK$exposure[1],
                        method = "Robust adjusted profile score",
                        nsnp = JAK$nsnp[1],
                        b = JAK2[["beta.hat"]],
                        se = JAK2[["beta.se"]],
                        pval = JAK2[["beta.p.value"]])
  smr_row = data.frame(id.exposure = JAK$id.exposure[1],
                       id.outcome = JAK$id.outcome[1],
                       outcome = JAK$outcome[1],
                       exposure = JAK$exposure[1],
                       method = "Summary-data-based Mendelian Randomization",
                       nsnp = JAK$nsnp[1],
                       b = JAK_smr$b_SMR,
                       se = JAK_smr$se_SMR,
                       pval = JAK_smr$p_SMR)
  JAK = rbind(JAK,raps_row,smr_row)
  
  #scatterplot
  return(mr_scatter_plot(JAK,JAK_dat))
}


#coronary artery disease
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"ebi-a-GCST005195",out_dat) #有多效性 

result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"ebi-a-GCST003116",out_dat) # 不相关 
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"ukb-d-I9_CHD",out_dat) #不相关
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"finn-b-I9_CHD",out_dat) #不相关
# coronary athrosclerosis
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"ukb-d-I9_CORATHER",out_dat) #不相关
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"finn-b-I9_CORATHER",out_dat) #临近相关 有多效性

# CAD other
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,TYK2_exposure,result_smr_jak1,result_smr_jak2,result_smr_jak3,result_smr_TYK2,outcome_dataset = "ukb-b-15748",out_dat) # 无多效性*
result <- presso(JAK1_exposure,JAK2_exposure,JAK3_exposure,TYK2_exposure,"ukb-b-15748",out_dat)
#CAD plot
JAK1_dat <- harmoize_dat2(JAK1_exposure,"ukb-b-15748")
JAK2_dat <- harmoize_dat2(JAK2_exposure,"ukb-b-15748")
JAK3_dat <- harmoize_dat2(JAK3_exposure,"ukb-b-15748")
TYK2_dat <- harmoize_dat2(TYK2_exposure,"ukb-b-15748")
# 1
JAK1_plot<-my_scatterplot_add_smr_raps(JAK1_dat,result_smr_jak1,"JAK1")
tiff(filename = paste0("tidy_reasult/JAK1_CAD_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300 )
par(cex = 2);
par(mar = c(4,4,4,4))
JAK1_plot
dev.off()
# 2
JAK2_plot<-my_scatterplot_add_smr_raps(JAK2_dat,result_smr_jak2,"JAK2")
tiff(filename = paste0("tidy_reasult/JAK2_CAD_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
JAK2_plot
dev.off()
# 3
JAK3_plot<-my_scatterplot_add_smr_raps(JAK3_dat,result_smr_jak3,"JAK3")
tiff(filename = paste0("tidy_reasult/JAK3_cad_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
JAK3_plot
dev.off()
# 4
TYK2_plot <- my_scatterplot_add_smr_raps(TYK2_dat,result_smr_TYK2,"TYK2")
tiff(filename = paste0("tidy_reasult/TYK2_cad_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
TYK2_plot
dev.off()

result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"finn-b-I9_ANGIO",out_dat) #相关 有多效
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"ukb-b-7869",out_dat) # 相关 无多效

# stroke
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,result_smr_jak1,result_smr_jak2,result_smr_jak3,"ukb-b-8714",out_dat) #无水平多效性 ###



result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"ebi-a-GCST005838",out_dat) #不相关
result <- JAK_result(JAK1_exposure,JAK2_exposure,JAK3_exposure,"finn-b-C_STROKE",out_dat) #不相关

# finn-b-I9_REVACS_EXNONE emergencycoronaryrevascularization 正相关
JAK2_MR <- MR_dat(JAK2_exposure)
JAK2_MR$
# 'ukb-b-8714' stroke 负相关
# ukb-d-I9_STR_EXH stroke 负相关
# ebi-a-GCST005838 stroke 
# finn-b-I9_STR_SAH_EXNONE includingSAH 
# ukb-d-I50 heartfailure 临近_负相关
# ebi-a-GCST011364 myocardialinfarction 负相关
# ebi-a-GCST005195 coronaryheartdisease 负相关
# ebi-a-GCST003116 coronaryartherydisease 
# ukb-b-7869 coronaryangioplasty
# ukb-b-11064 coronarybypass 临近_负相关
# ukb-b-15748 coronaryangioplasty+bypass 负相关
# ukb-d-I9_CORATHER coronaryathrosclerosis 不相关
# ukb-d-I9_CHD majorcoronaryevent
# ebi-a-GCST005194 Coronaryarterydisease 负相关
# finn-b-I9_ATHSCLE athrosclerosis 正相关
# finn-b-I9_ANGIO Coronaryangiopasty 正相关
# finn-b-I9_CORATHER coronaryathrosclerosis 临近_正相关
# finn-b-I9_CHD majorcoronaryevent
# ukb-b-13391 transient attack 负相关
# finn-b-I9_PVT PVT 临近_负相关

# positive control
# finn-b-RHEUMA_NOS 正相关
# ukb-b-11874 负相关

JAK3_MR <- MR_dat(JAK3_exposure)
JAK3_MR
# ukb-b-453 myocardialinfarction_anterior 正相关
# ebi-a-GCST005195 coronaryheartdisease 临近_正相关
# ukb-d-I9_CORATHER coronaryathroselerosis 临近_正相关




























