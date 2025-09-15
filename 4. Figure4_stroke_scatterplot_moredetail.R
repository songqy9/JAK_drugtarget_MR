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
result_smr_jak1_CAD <- read.delim("JAK1_CAD_Msmr.msmr")
result_smr_jak2_CAD <- read.delim("JAK2_CAD_Msmr.msmr")
result_smr_jak3_CAD <- read.delim("JAK3_CAD_Msmr.msmr")
result_smr_TYK2_CAD <- read.delim("TYK2_CAD_Msmr.msmr")
# result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3)
# MR_exposure
exposure_read <- function(name){
  JAK = read_exposure_data(filename = paste0(name,"_5e-0x_freq_Fstat.csv"),sep = ",",snp_col = "SNP",beta_col = "beta",se_col = "se",effect_allele_col = "effect_allele",other_allele_col = "other_allele",eaf_col = "eaf",samplesize_col = "samplesize",clump = FALSE)
  JAK = clump_data(JAK,clump_kb = 100,
                   clump_r2 = 0.3,
                   pop = "EUR")
  return(JAK)
}
JAK1_exposure <- exposure_read("JAK1")
JAK2_exposure <- exposure_read("JAK2")
JAK3_exposure <- exposure_read("JAK3")
TYK2_exposure <- exposure_read("TYK2")
# Scatter plot add smr and raps
harmoize_dat2 <- function(JAK,outcome_dataset){
  outcome = extract_outcome_data(snps = JAK$SNP,
                                 outcomes = outcome_dataset)
  dat <- harmonise_data(JAK,outcome)
  # dat2 <- dplyr::filter(dat,mr_keep=TRUE)
  # write.csv(dat2,paste0(name,"_SNP_asIV.csv"),row.names = F,quote = F)
  return(dat)
}
# scatter plot function mr_scatter_plot with legend
myfunction_mr_scatter_plot_withlegend <- function(mr_results, dat, plotnumber)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#e31a1c", "#fb9a99", "#33a02c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      # ggplot2::ggtitle(paste(plotnumber,"Scatterplot for",d$exposure[1],"on",d$outcome[1]))+
      ggplot2::labs(colour="", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::guides(colour=ggplot2::guide_legend(ncol=2,direction = "vertical")) +
      ggplot2::theme(text=element_text(family="serif"),legend.position="top", legend.direction="vertical",
                     panel.grid.major =element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black")) +
      scale_y_continuous(breaks = scales::breaks_pretty(4),labels = scales::scientific) #y轴采用科学计数法;规定tick的数量 并使得它们整齐
      # ggplot2::theme_classic()
  })
  mrres
}
# scatter plot function mr_scatter_plot without legend
myfunction_mr_scatter_plot_withoutlegend <- function(mr_results, dat, plotnumber)
{
  # dat <- subset(dat, paste(id.outcome, id.exposure) %in% paste(mr_results$id.outcome, mr_results$id.exposure))
  mrres <- plyr::dlply(dat, c("id.exposure", "id.outcome"), function(d)
  {
    d <- plyr::mutate(d)
    if(nrow(d) < 2 | sum(d$mr_keep) == 0)
    {
      return(blank_plot("Insufficient number of SNPs"))
    }
    d <- subset(d, mr_keep)
    index <- d$beta.exposure < 0
    d$beta.exposure[index] <- d$beta.exposure[index] * -1
    d$beta.outcome[index] <- d$beta.outcome[index] * -1
    mrres <- subset(mr_results, id.exposure == d$id.exposure[1] & id.outcome == d$id.outcome[1])
    mrres$a <- 0
    if("MR Egger" %in% mrres$method)
    {
      temp <- mr_egger_regression(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger"] <- temp$b_i
    }
    
    if("MR Egger (bootstrap)" %in% mrres$method)
    {
      temp <- mr_egger_regression_bootstrap(d$beta.exposure, d$beta.outcome, d$se.exposure, d$se.outcome, default_parameters())
      mrres$a[mrres$method == "MR Egger (bootstrap)"] <- temp$b_i
    }
    
    ggplot2::ggplot(data=d, ggplot2::aes(x=beta.exposure, y=beta.outcome)) +
      ggplot2::geom_errorbar(ggplot2::aes(ymin=beta.outcome-se.outcome, ymax=beta.outcome+se.outcome), colour="grey", width=0) +
      ggplot2::geom_errorbarh(ggplot2::aes(xmin=beta.exposure-se.exposure, xmax=beta.exposure+se.exposure), colour="grey", height=0) +
      ggplot2::geom_point() +
      ggplot2::geom_abline(data=mrres, ggplot2::aes(intercept=a, slope=b, colour=method), show.legend=TRUE) +
      ggplot2::scale_colour_manual(values=c("#a6cee3", "#1f78b4", "#b2df8a", "#e31a1c", "#fb9a99", "#33a02c", "#fdbf6f", "#ff7f00", "#cab2d6", "#6a3d9a", "#ffff99", "#b15928")) +
      ggplot2::labs(colour="", x=paste("SNP effect on", d$exposure[1]), y=paste("SNP effect on", d$outcome[1])) +
      ggplot2::theme(text=element_text(family="serif",size = 15),legend.position= "none",
                     panel.grid.major =element_blank(),
                     panel.grid.minor = element_blank(),
                     panel.background = element_blank(),
                     axis.line = element_line(colour = "black")) + #去掉legend
      # ggplot2::guides(colour=ggplot2::guide_legend(ncol=2)) +
      ggplot2::guides(colour="none") +
      scale_y_continuous(breaks = scales::breaks_pretty(4),labels = scales::scientific) #y轴采用科学计数法;规定tick的数量 并使得它们整齐
  })
  mrres
}
blank_plot <- function(message)
{
  ggplot2::ggplot(data.frame(a=0,b=0,n=message)) + 
    ggplot2::geom_text(ggplot2::aes(x=a,y=b,label=n)) + 
    ggplot2::labs(x=NULL,y=NULL) + 
    ggplot2::theme(axis.text=ggplot2::element_blank(), axis.ticks=ggplot2::element_blank())
}
# scatter plot add smr
mr_scatterplot_add_smr_raps <- function(JAK_dat,JAK_smr,name,plotnumber,withlegend,ylab_outcome){
  JAK_dat = dplyr::mutate(JAK_dat,outcome=rep(ylab_outcome,nrow(JAK_dat)))
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
  if(withlegend==TRUE){
    p = myfunction_mr_scatter_plot_withlegend(JAK,JAK_dat,plotnumber)
  }
  else{
    p = myfunction_mr_scatter_plot_withoutlegend(JAK,JAK_dat,plotnumber)
  }
  return(p)
}
#stroke plot
JAK1_dat <- harmoize_dat2(JAK1_exposure,"ukb-b-8714")
JAK2_dat <- harmoize_dat2(JAK2_exposure,"ukb-b-8714")
JAK3_dat <- harmoize_dat2(JAK3_exposure,"ukb-b-8714")
TYK2_dat <- harmoize_dat2(TYK2_exposure,"ukb-b-8714")
#CAD plot
JAK1_dat_CAD <- harmoize_dat2(JAK1_exposure,"ukb-b-15748")
JAK2_dat_CAD <- harmoize_dat2(JAK2_exposure,"ukb-b-15748")
JAK3_dat_CAD <- harmoize_dat2(JAK3_exposure,"ukb-b-15748")
TYK2_dat_CAD <- harmoize_dat2(TYK2_exposure,"ukb-b-15748")
# 1
JAK1_plot<-mr_scatterplot_add_smr_raps(JAK1_dat,result_smr_jak1,name = "JAK1",plotnumber = "(a)",withlegend = FALSE,ylab_outcome = "stroke")
tiff(filename = paste0("tidy_reasult/JAK1_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300 )
par(cex = 2);
par(mar = c(4,4,4,4))
JAK1_plot
dev.off()
# 2
JAK2_plot<-mr_scatterplot_add_smr_raps(JAK2_dat,result_smr_jak2,name = "JAK2",plotnumber = "(b)",withlegend = FALSE,ylab_outcome = "stroke")
tiff(filename = paste0("tidy_reasult/JAK2_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
JAK2_plot
dev.off()
# 3
JAK3_plot<-mr_scatterplot_add_smr_raps(JAK3_dat,result_smr_jak3,"JAK3",plotnumber = "(c)",withlegend = FALSE,ylab_outcome = "stroke")
tiff(filename = paste0("tidy_reasult/JAK3_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
JAK3_plot
dev.off()
# 4
TYK2_plot<-mr_scatterplot_add_smr_raps(TYK2_dat,result_smr_TYK2,"TYK2",plotnumber = "(d)",withlegend = FALSE,ylab_outcome = "stroke")
tiff(filename = paste0("tidy_reasult/TYK2_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
TYK2_plot
dev.off()
# CAD
# 1
JAK1_plot_CAD<-mr_scatterplot_add_smr_raps(JAK1_dat_CAD,result_smr_jak1_CAD,name = "JAK1",plotnumber = "(e)",withlegend = FALSE,ylab_outcome = "CAD")
tiff(filename = paste0("tidy_reasult/JAK1_CAD_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300 )
par(cex = 2);
par(mar = c(4,4,4,4))
JAK1_plot
dev.off()
# 2
JAK2_plot_CAD<-mr_scatterplot_add_smr_raps(JAK2_dat_CAD,result_smr_jak2_CAD,name = "JAK2",plotnumber = "(f)",withlegend = FALSE,ylab_outcome = "CAD")
tiff(filename = paste0("tidy_reasult/JAK2_CAD_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300 )
par(cex = 2);
par(mar = c(4,4,4,4))
JAK1_plot
dev.off()
# 3
JAK3_plot_CAD<-mr_scatterplot_add_smr_raps(JAK3_dat_CAD,result_smr_jak3_CAD,name = "JAK3",plotnumber = "(g)",withlegend = FALSE,ylab_outcome = "CAD")
tiff(filename = paste0("tidy_reasult/JAK3_CAD_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300 )
par(cex = 2);
par(mar = c(4,4,4,4))
JAK1_plot
dev.off()
# 4 
TYK2_plot_CAD<-mr_scatterplot_add_smr_raps(TYK2_dat_CAD,result_smr_TYK2_CAD,name = "TYK2",plotnumber = "(h)",withlegend = FALSE,ylab_outcome = "CAD")
tiff(filename = paste0("tidy_reasult/TYK2_CAD_eQTL.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300 )
par(cex = 2);
par(mar = c(4,4,4,4))
JAK1_plot
dev.off()
# 合并 
# 注意这个图形的名称每次都不一样 根据$的情况选
p <- cowplot::plot_grid(JAK1_plot$`bYgehS.ukb-b-8714`,JAK2_plot$`bVfECB.ukb-b-8714`,JAK3_plot$`tcrPSN.ukb-b-8714`,TYK2_plot$`SCk6Ut.ukb-b-8714`,JAK1_plot_CAD$`bYgehS.ukb-b-15748`,JAK2_plot_CAD$`bVfECB.ukb-b-15748`,JAK3_plot_CAD$`tcrPSN.ukb-b-15748`,TYK2_plot_CAD$`SCk6Ut.ukb-b-15748`,labels = c("(a)","(b)","(c)","(d)","(e)","(f)","(g)","(h)"),label_fontfamily = "serif",nrow = 2,ncol = 4)
tiff(filename = paste0("tidy_reasult/JAKcombine_eQTL.tiff"),
     width = 12,height = 6,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
p
dev.off()

# 保存legend
JAK2_plot<-mr_scatterplot_add_smr_raps(JAK2_dat,result_smr_jak2,name = "JAK2",plotnumber = "(b)",withlegend = TRUE)
tiff(filename = paste0("tidy_reasult/legend.tiff"),
     width = 5.5,height = 5.5,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
JAK2_plot
dev.off()




