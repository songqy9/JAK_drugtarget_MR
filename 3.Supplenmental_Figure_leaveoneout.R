setwd("D:/GitHub/MR/")
library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
library(ggplot2)
result_smr_jak1 <- read.delim("JAK1_Msmr.msmr")
result_smr_jak2 <- read.delim("JAK2_Msmr.msmr")
result_smr_jak3 <- read.delim("JAK3_Msmr.msmr")
result_smr_TYK2 <- read.delim("TYK2_Stroke_Msmr.msmr")
result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3,result_smr_TYK2)
# MR_exposure
exposure_read <- function(name){
  JAK = read_exposure_data(filename = paste0(name,"_5e-0x_freq_Fstat.csv"),sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",eaf_col = "eaf",samplesize_col = "samplesize",clump = FALSE)
  JAK = clump_data(JAK,clump_kb = 100,
                   clump_r2 = 0.3,
                   # clump_r2 = 0.5,
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
  # write.csv(dat2,paste0(name,"_SNP_asIV.csv"),row.names = F,quote = F)# 获取最终的工具变量IV
  return(dat)
}
# MR result
# leave_one_out
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK1)~"on stroke"))
                                            })
  res
}
# a leave group out analysis of 30% most influence SNP
# subset_snp <- function(dat,partion,filename){
#   # dat1 = mr_leaveoneout(dat = dat)
#   # dat1 = dat1[,c("SNP","b")]
#   dat1 = dplyr::mutate(dat,b = beta.outcome/beta.exposure)
#   dat1 = dat1[,c("SNP","b")]
#   dat2 = dplyr::arrange(dat1,b)
#   subset_number = nrow(dat2)-round(nrow(dat2)*partion)
#   dat_tmp = tail(dat2,subset_number)
#   SNP_tmp = dat_tmp$SNP
#   if("All" %in% SNP_tmp){
#     subset_number = subset_number + 1
#   }else{
#     subset_number = subset_number
#   }
#   dat3 = tail(dat2,subset_number)
#   dat3 = dat3[which(!(dat3$SNP == "All")),]
#   dat_har_sub = dat[which(dat$SNP %in% dat3$SNP),]
#   write.csv(dat_har_sub,paste0("30%subset_leaveoneout_",filename,".csv"),row.names = F,quote = F)
#   leaveoneout_subset = mr_leaveoneout(dat_har_sub)
#   leaveoneout_subset = leaveoneout_subset[which(leaveoneout_subset$SNP == "All"),]
#   leaveoneout_subset = dplyr::mutate(leaveoneout_subset,SNP=dplyr::if_else(SNP == "All","30% top SNP","mistake"))
#   return(leaveoneout_subset)
# }
# a <- subset_snp(JAK2_har_CAD,0.3,"JAK2_CAD")

leaveoneout_plotanddata <- function(dat,exposure_name,outcome_name){
  # dat = harmoize_dat(dat)
  subset_filename = paste0(exposure_name,"_",outcome_name)
  # subset_leaveoneout_result = subset_snp(dat,0.3,subset_filename)
  leave_one_out = mr_leaveoneout(dat = dat)
  # leave_one_out = rbind(leave_one_out,subset_leaveoneout_result)
  leave_one_out = dplyr::mutate(leave_one_out,exposure = rep(exposure_name,nrow(leave_one_out)))
  leave_one_out = dplyr::mutate(leave_one_out,outcome = rep(outcome_name,nrow(leave_one_out)))
  # leave_one_out_writeout1 = leave_one_out[which(!(leave_one_out$SNP %in% c("All","30% top SNP"))),] 
  # leave_one_out_writeout1 = dplyr::arrange(leave_one_out_writeout1,desc(b))
  # leave_one_out_writeout2 = leave_one_out[which(leave_one_out$SNP == "All"),] 
  # leave_one_out_writeout3 = leave_one_out[which(leave_one_out$SNP == "30% top SNP"),] 
  # leave_one_out_writeout = rbind(leave_one_out_writeout1,leave_one_out_writeout3,leave_one_out_writeout2)
  # write.csv(leave_one_out_writeout,paste0(exposure_name,"_",outcome_name,"leaveoneout.csv"),row.names = F,quote = F)
  res = leaveonout_plot_correct(leaveoneout_results = leave_one_out)  
  return(res)
}

JAK1_har <- harmoize_dat(JAK1_exposure,outcome_dataset = "ukb-b-8714",name = "JAK1")
JAK2_har <- harmoize_dat(JAK2_exposure,outcome_dataset = "ukb-b-8714",name = "JAK2")
JAK3_har <- harmoize_dat(JAK3_exposure,outcome_dataset = "ukb-b-8714",name = "JAK3")
TYK2_har <- harmoize_dat(TYK2_exposure,outcome_dataset = "ukb-b-8714",name = "TYK2")
JAK1_har_CAD <- harmoize_dat(JAK1_exposure,outcome_dataset = "ukb-b-15748",name = "JAK1")
JAK2_har_CAD <- harmoize_dat(JAK2_exposure,outcome_dataset = "ukb-b-15748",name = "JAK2")
JAK3_har_CAD <- harmoize_dat(JAK3_exposure,outcome_dataset = "ukb-b-15748",name = "JAK3")
TYK2_har_CAD <- harmoize_dat(TYK2_exposure,outcome_dataset = "ukb-b-15748",name = "TYK2")
mr_leaveoneout_plot()
# plot
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                           axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK1)~"on stroke"))
                                            })
  res
}
JAK1_leaveoneout<-leaveoneout_plotanddata(JAK1_har,exposure_name = "JAK1",outcome_name = "stroke")
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK2)~"on stroke"))
                                            })
  res
}
JAK2_leaveoneout<-leaveoneout_plotanddata(JAK2_har,exposure_name = "JAK2",outcome_name = "stroke")
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK3)~"on stroke"))
                                            })
  res
}
JAK3_leaveoneout<-leaveoneout_plotanddata(JAK3_har,exposure_name = "JAK3",outcome_name = "stroke")
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(TYK2)~"on stroke"))
                                            })
  res
}
TYK2_leaveoneout<-leaveoneout_plotanddata(TYK2_har,exposure_name = "TYK2",outcome_name = "stroke")############
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK1)~"on CAD"))
                                            })
  res
}
JAK1_leaveoneout_CAD<-leaveoneout_plotanddata(JAK1_har_CAD,exposure_name = "JAK1",outcome_name = "CAD")
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK2)~"on CAD"))
                                            })
  res
}
JAK2_leaveoneout_CAD<-leaveoneout_plotanddata(JAK2_har_CAD,exposure_name = "JAK2",outcome_name = "CAD")
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(JAK3)~"on CAD"))
                                            })
  res
}
JAK3_leaveoneout_CAD<-leaveoneout_plotanddata(JAK3_har_CAD,exposure_name = "JAK3",outcome_name = "CAD")
leaveonout_plot_correct <- function (leaveoneout_results) 
{
  requireNamespace("ggplot2", quietly = TRUE)
  requireNamespace("plyr", quietly = TRUE)
  res <- plyr::dlply(leaveoneout_results, c("id.exposure", 
                                            "id.outcome"), function(d) {
                                              d <- plyr::mutate(d)
                                              if (sum(!grepl("All", d$SNP)) < 3) {
                                                return(blank_plot("Insufficient number of SNPs"))
                                              }
                                              d$up <- d$b + 1.96 * d$se
                                              d$lo <- d$b - 1.96 * d$se
                                              d$tot <- 1
                                              # d$tot[d$SNP != "All"] <- 0.01
                                              # d$SNP <- as.character(d$SNP)
                                              # nom <- d$SNP[which(!(d$SNP %in% c("All","30% top SNP")))]
                                              # nom <- nom[order(d$b)]
                                              # d <- rbind(d, d[nrow(d), ])
                                              # d$SNP[nrow(d) - 1] <- ""
                                              # d$b[nrow(d) - 1] <- NA
                                              # d$up[nrow(d) - 1] <- NA
                                              # d$lo[nrow(d) - 1] <- NA
                                              # d$SNP <- ordered(d$SNP, levels = c("All","","30% top SNP", nom))
                                              d$tot[d$SNP != "All"] <- 0.01
                                              d$SNP <- as.character(d$SNP)
                                              nom <- d$SNP[d$SNP != "All"]
                                              nom <- nom[order(d$b)]
                                              d <- rbind(d, d[nrow(d), ])
                                              d$SNP[nrow(d) - 1] <- ""
                                              d$b[nrow(d) - 1] <- NA
                                              d$up[nrow(d) - 1] <- NA
                                              d$lo[nrow(d) - 1] <- NA
                                              d$SNP <- ordered(d$SNP, levels = c("All", "", nom))
                                              ggplot2::ggplot(d, ggplot2::aes(y = SNP, x = b)) + ggplot2::geom_vline(xintercept = 0, 
                                                                                                                     linetype = "dotted") + ggplot2::geom_errorbarh(ggplot2::aes(xmin = lo, 
                                                                                                                                                                                 xmax = up, size = as.factor(tot), colour = as.factor(tot)), 
                                                                                                                                                                    height = 0) + ggplot2::geom_point(ggplot2::aes(colour = as.factor(tot))) + 
                                                ggplot2::geom_hline(ggplot2::aes(yintercept = which(levels(SNP) %in% 
                                                                                                      "")), colour = "grey") + ggplot2::scale_colour_manual(values = c("#4169E1", 
                                                                                                                                                                       "red")) + ggplot2::scale_size_manual(values = c(0.4, 
                                                                                                                                                                                                                       1)) + ggplot2::theme(legend.position = "none",
                                                                                                                                                                                                                                            text=element_text(family="serif",size = 11),
                                                                                                                                                                                                                                            # axis.text.y = ggplot2::element_text(size = 8), 
                                                                                                                                                                                                                                            axis.ticks.y = ggplot2::element_line(size = 0), 
                                                                                                                                                                                                                                            # axis.title.x = ggplot2::element_text(size = 8),
                                                                                                                                                                                                                                            panel.grid.major =element_blank(),
                                                                                                                                                                                                                                            panel.grid.minor = element_blank(),
                                                                                                                                                                                                                                            panel.background = element_blank(),
                                                                                                                                                                                                                                            axis.line = element_line(colour = "grey")) +
                                                # ggplot2::labs(y = "", x = paste0("MR leave-one-out sensitivity analysis for\n",d$exposure[1]," on ", d$outcome[1]))
                                                ggplot2::labs(y = "", x = expression("MR leave-one-out sensitivity analysis for"~italic(TYK2)~"on CAD"))
                                            })
  res
}
TYK2_leaveoneout_CAD<-leaveoneout_plotanddata(TYK2_har_CAD,exposure_name = "TYK2",outcome_name = "CAD")
p <- cowplot::plot_grid(JAK1_leaveoneout$`Eduu0v.ukb-b-8714`,JAK2_leaveoneout$`4EhbuB.ukb-b-8714`,JAK3_leaveoneout$`q05Z7K.ukb-b-8714`,TYK2_leaveoneout$`uqup3e.ukb-b-8714`,JAK1_leaveoneout_CAD$`Eduu0v.ukb-b-15748`,JAK2_leaveoneout_CAD$`4EhbuB.ukb-b-15748`,JAK3_leaveoneout_CAD$`q05Z7K.ukb-b-15748`,TYK2_leaveoneout_CAD$`uqup3e.ukb-b-15748`,labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),label_fontfamily = "serif",nrow = 2,ncol = 4)
tiff(filename = paste0("tidy_reasult/JAKcombine_leaveoneout.tiff"),
     width = 17.3,height = 12.3,units = "in",
     pointsize = 5,res = 100)
par(cex = 2);
par(mar = c(4,4,4,4))
p
dev.off()

# Tidy table
JAK1_stroke <- read.csv("JAK1_strokeleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
JAK2_stroke <- read.csv("JAK2_strokeleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
JAK3_stroke <- read.csv("JAK3_strokeleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
TYK2_stroke <- read.csv("TYK2_strokeleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
JAK1_CAD <- read.csv("JAK1_CADleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
JAK2_CAD <- read.csv("JAK2_CADleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
JAK3_CAD <- read.csv("JAK3_CADleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
TYK2_CAD <- read.csv("TYK2_CADleaveoneout.csv") %>% .[,c("exposure","outcome","SNP","b","se","p")]
leaveoneout_result <- rbind(JAK1_stroke,JAK2_stroke,JAK3_stroke,TYK2_stroke,JAK1_CAD,JAK2_CAD,JAK3_CAD,TYK2_CAD)
names(leaveoneout_result)<-c("exposure","outcome","Removed SNP","Beta","SE","P")
write.csv(leaveoneout_result,"tidy_reasult/leaveoneout_result.csv",row.names = F,quote = F)



# Tidy table top 30% SNP
origin_JAK1_stroke_SNP <- read.csv("JAK1_SNP_asIV.csv")
origin_JAK2_stroke_SNP <- read.csv("JAK2_SNP_asIV.csv")
origin_JAK3_stroke_SNP <- read.csv("JAK3_SNP_asIV.csv")
origin_TYK2_stroke_SNP <- read.csv("TYK2_SNP_asIV.csv")
origin_JAK1_CAD_SNP <- read.csv("JAK1_SNP_asIV_CAD.csv")
origin_JAK2_CAD_SNP <- read.csv("JAK2_SNP_asIV_CAD.csv")
origin_JAK3_CAD_SNP <- read.csv("JAK3_SNP_asIV_CAD.csv")
origin_TYK2_CAD_SNP <- read.csv("TYK2_SNP_asIV_CAD.csv")

JAK1_stroke_SNP <- read.csv("30%subset_leaveoneout_JAK1_stroke.csv")
JAK2_stroke_SNP <- read.csv("30%subset_leaveoneout_JAK2_stroke.csv")
JAK3_stroke_SNP <- read.csv("30%subset_leaveoneout_JAK3_stroke.csv")
TYK2_stroke_SNP <- read.csv("30%subset_leaveoneout_TYK2_stroke.csv")
JAK1_CAD_SNP <- read.csv("30%subset_leaveoneout_JAK1_CAD.csv")
JAK2_CAD_SNP <- read.csv("30%subset_leaveoneout_JAK2_CAD.csv")
JAK3_CAD_SNP <- read.csv("30%subset_leaveoneout_JAK3_CAD.csv")
TYK2_CAD_SNP <- read.csv("30%subset_leaveoneout_TYK2_CAD.csv")

exclude_SNP_fun <- function(dat1,dat2,name){
  dat = dat1[which(!(dat1$SNP %in% dat2$SNP)),]
  dat = dat[,c("SNP","effect_allele.exposure","other_allele.exposure","beta.exposure","beta.outcome")]
  dat = dplyr::mutate(dat,ratio = beta.outcome/beta.exposure)
  dat = dplyr::mutate(dat,name_col = rep(name,nrow(dat)))
  write.csv(dat,paste0("30%exclude_list_",name,".csv"),row.names = F,quote = F)
  return(dat)
}
exclude_JAK1_SNP_stroke <- exclude_SNP_fun(origin_JAK1_stroke_SNP,JAK1_stroke_SNP,"JAK1_stroke")
exclude_JAK2_SNP_stroke <- exclude_SNP_fun(origin_JAK2_stroke_SNP,JAK2_stroke_SNP,"JAK2_stroke")
exclude_JAK3_SNP_stroke <- exclude_SNP_fun(origin_JAK3_stroke_SNP,JAK3_stroke_SNP,"JAK3_stroke")
exclude_TYK2_SNP_stroke <- exclude_SNP_fun(origin_TYK2_stroke_SNP,TYK2_stroke_SNP,"TYK2_stroke")
exclude_JAK1_SNP_CAD <- exclude_SNP_fun(origin_JAK1_CAD_SNP,JAK1_CAD_SNP,"JAK1_CAD")
exclude_JAK2_SNP_CAD <- exclude_SNP_fun(origin_JAK2_CAD_SNP,JAK2_CAD_SNP,"JAK2_CAD")
exclude_JAK3_SNP_CAD <- exclude_SNP_fun(origin_JAK3_CAD_SNP,JAK3_CAD_SNP,"JAK3_CAD")
exclude_TYK2_SNP_CAD <- exclude_SNP_fun(origin_TYK2_CAD_SNP,TYK2_CAD_SNP,"TYK2_CAD")

Tidy_exclude_30SNP <- rbind(exclude_JAK1_SNP_stroke,exclude_JAK2_SNP_stroke,
                            exclude_JAK3_SNP_stroke,
                            exclude_TYK2_SNP_stroke,
                            exclude_JAK1_SNP_CAD,
                            exclude_JAK2_SNP_CAD,
                            exclude_JAK3_SNP_CAD,
                            exclude_TYK2_SNP_CAD)
write.csv(Tidy_exclude_30SNP,"tidy_reasult/leavegroupout_30exclude_SNP.csv",row.names = F,quote = F)









