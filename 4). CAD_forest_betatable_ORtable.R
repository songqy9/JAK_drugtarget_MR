setwd("D:/GitHub/MR/")
library(forestploter)
# library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
result <- read.csv("result_ukb-b-8714.csv",header = T)
result_smr_jak1 <- read.delim("JAK1_Msmr.msmr")
result_smr_jak2 <- read.delim("JAK2_Msmr.msmr")
result_smr_jak3 <- read.delim("JAK3_Msmr.msmr")
result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3)
# tidy for OR
tidy_function <- function(name){
  JAK = dplyr::filter(result,exposure==name)
  JAK_OR = dplyr::select(JAK,method,OR,OR_lci95,OR_uci95)
  JAK_OR = dplyr::mutate(JAK_OR," "=rep("                                                                                                              ",6))# 这个空格就是留着森林图的
  JAK_OR = dplyr::mutate(JAK_OR,estimate = paste0(round(JAK_OR$OR,digits = 2),"(",round(JAK_OR$OR_lci95,digits = 2),"-",round(JAK_OR$OR_uci95,digits = 2),")"))
  JAK_OR[nrow(JAK_OR)+1,] <- " "
  JAK_OR2 = JAK_OR[c(7,6,2,1,3,4,5),]
  nsnp = JAK[1,"nsnp"]
  title = paste0(name,"(nSNP=",nsnp,")")
  name_order = names(JAK_OR2)
  JAK_OR2 = mutate(JAK_OR2,exposure = c(title,rep(" ",6)))
  JAK_OR2 = dplyr::select(JAK_OR2,exposure,name_order)
  return(JAK_OR2)
}
JAK1_tidy <- tidy_function("JAK1")
JAK2_tidy <- tidy_function("JAK2")
JAK3_tidy <- tidy_function("JAK3")
# # add smr
# add_smr <- function(JAK_tidy,name){
#   smr = dplyr::filter(result_smr, Gene==name)
#   smr = dplyr::mutate(smr, beta_lci=b_SMR-se_SMR)
#   smr = dplyr::mutate(smr, beta_uci=b_SMR+se_SMR)
#   smr_row = c(" ","SMR",smr$b_SMR,smr$beta_lci,smr$beta_uci,"                                          ",paste0(round(smr$b_SMR,digits = 2),"(",round(smr$beta_lci,digits = 2),"-",round(smr$beta_uci,digits = 2),")"))
#   tidy = rbind(JAK_tidy,smr_row)
#   tidy = tidy[c(1,7,2,3,4,5,6),]# 把smr排在前面
#   return(tidy)
# }
# JAK1_tidy_smr <- add_smr(JAK1_tidy,"JAK1")
# JAK2_tidy_smr <- add_smr(JAK2_tidy,"JAK2")
# JAK3_tidy_smr <- add_smr(JAK3_tidy,"JAK3")
# dt <- read.csv("dt准备1.csv",check.names = F,
#                strip.white =F,
#                # 全部读取为字符变量
#                colClasses ="character"  )
# forest
dt <- rbind(JAK1_tidy,JAK2_tidy,JAK3_tidy)
dt <- dplyr::rename(dt,`OR(confident interval)`=estimate)
numvar=c("OR" ,"OR_lci95","OR_uci95")
dt[,c(numvar)]<-lapply(dt[,numvar],as.numeric)
{# 创建空白并移到
# dt$` ` <- paste(rep(" ", 20), # 这个空格就是留着森林图的
#                 collapse = " ")
# names(dt)
# # 移动下空白行的位置
# dt <- dt[ ,c(1:5,7,6)]
}

# Define theme
tm <- forest_theme(base_size = 8, # 更改整体字体大小
                   refline_col = "red",
                   # 树颜色
                   # ci_col="#228B22",
                   footnote_col = "#636363",
                   footnote_fontface = "italic")

p <- forest(
  #  筛选用到的列,8为第4列在excel为空
  dt[,c(1:2, 6:7)], 
  est = dt$OR,
  lower = dt$OR_lci95,
  upper = dt$OR_uci95,
  #sizes = dt$se,
  # 这边可以更改森林图的位置
  ci_column = 3, # 第四列为空，显示森林图
  # 数值虚线，一般OR/HR=1的位置
  ref_line = 1,
  # clip = c(-0.05,0.03),
  boxsize = 100, # 设置box 大小
  # arrow_lab = c("Placebo Better", "Treatment Better"),# 底下两个文本
  xlim = c(0.95,1.05),# X轴刻度范围
  ticks_at = c(0.96,0.98,1,1.02,1.04),# X轴刻度间隔
  # footnote = "forest plot of JAK1, JAK2, and JAK3",# 左下角的脚注
  theme = tm)# 自定义的主题

# eoffice::topptx(filename = "plot导出的.pptx") # Plot作图
tiff(filename = "tidy_reasult/CAD_OR_forestplot.tiff",
     width = 10,height = 10,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
plot(p)
dev.off()

# table
# table for β
table_function <- function(name){
  JAK = dplyr::filter(result,exposure==name)
  JAK_OR = dplyr::select(JAK,method,nsnp,beta,beta_lci,beta_uci,p)
  JAK_OR = dplyr::mutate(JAK_OR,CI = paste0(round(JAK_OR$beta_lci,digits = 3),":",round(JAK_OR$beta_uci,digits = 3)))
  JAK_OR = dplyr::mutate(JAK_OR,beta=round(as.numeric(beta),digits = 3))
  JAK_OR = dplyr::mutate(JAK_OR,p=round(as.numeric(p),digits = 3))
  name_order = names(JAK_OR)
  JAK_OR2 = mutate(JAK_OR,exposure = c(name,rep(" ",5)))
  JAK_OR2 = dplyr::select(JAK_OR2,exposure,name_order)
  JAK_OR2 = dplyr::rename(JAK_OR2,`Exposure traits`=exposure,NSNP=nsnp,`MR method`=method,Beta=beta,`p Val`=p)
  JAK_OR2 = dplyr::select(JAK_OR2,`Exposure traits`,NSNP,`MR method`,Beta,CI,`p Val`)
  # smr 
  # nsnp = JAK[1,"nsnp"]
  # smr = dplyr::filter(result_smr,Gene==name)
  # smr = dplyr::filter(result_smr, Gene==name)
  # smr = dplyr::mutate(smr, beta_lci=b_SMR-se_SMR)
  # smr = dplyr::mutate(smr, beta_uci=b_SMR+se_SMR)
  # smr = dplyr::mutate(smr, b_SMR=round(as.numeric(b_SMR),digits = 3))
  # smr = dplyr::mutate(smr, p_SMR=round(as.numeric(p_SMR),digits = 3))
  # smr_row = c(" ",nsnp,"SMR",round(as.numeric(smr$b_SMR),digits = 3),paste0(round(smr$beta_lci,digits = 3),"-",round(smr$beta_uci,digits = 3)),round(as.numeric(smr$p_SMR),digits = 3))
  # tidy = rbind(JAK_OR2,smr_row)
  tidy = JAK_OR2
  tidy = tidy[c(6,2,1,3,4,5),]
  tidy[1,1] = tidy[3,1] #调整“JAKx”的位置
  tidy[3,1] = tidy[2,1]
  tidy = dplyr::mutate(tidy,CI=as.character(CI))
  return(tidy)
}
JAK1_table <- table_function("JAK1")
JAK2_table <- table_function("JAK2")
JAK3_table <- table_function("JAK3")
JAK_table = rbind(JAK1_table,JAK2_table,JAK3_table)
write.csv(JAK_table,"tidy_reasult/CAD_scatterplot_ukb-b-8714_table.csv",row.names = F,quote = F)

# table for OR
table_function <- function(name){
  JAK = dplyr::filter(result,exposure==name)
  JAK_OR = dplyr::select(JAK,method,nsnp,OR,OR_lci95,OR_uci95,p)
  JAK_OR = dplyr::mutate(JAK_OR,CI = paste0(round(JAK_OR$OR_lci95,digits = 3),"-",round(JAK_OR$OR_uci95,digits = 3)))
  JAK_OR = dplyr::mutate(JAK_OR,OR=round(as.numeric(OR),digits = 3))
  JAK_OR = dplyr::mutate(JAK_OR,p=round(as.numeric(p),digits = 3))
  name_order = names(JAK_OR)
  JAK_OR2 = mutate(JAK_OR,exposure = c(name,rep(" ",5)))
  JAK_OR2 = dplyr::select(JAK_OR2,exposure,name_order)
  JAK_OR2 = dplyr::rename(JAK_OR2,`Exposure traits`=exposure,NSNP=nsnp,`MR method`=method,`Odd ratio`=OR,`p Val`=p)
  JAK_OR2 = dplyr::select(JAK_OR2,`Exposure traits`,NSNP,`MR method`,`Odd ratio`,CI,`p Val`)
  # smr 
  # nsnp = JAK[1,"nsnp"]
  # smr = dplyr::filter(result_smr,Gene==name)
  # smr = dplyr::filter(result_smr, Gene==name)
  # smr = dplyr::mutate(smr, beta_lci=b_SMR-se_SMR)
  # smr = dplyr::mutate(smr, beta_uci=b_SMR+se_SMR)
  # smr = dplyr::mutate(smr, b_SMR=round(as.numeric(b_SMR),digits = 3))
  # smr = dplyr::mutate(smr, p_SMR=round(as.numeric(p_SMR),digits = 3))
  # smr_row = c(" ",nsnp,"SMR",round(as.numeric(smr$b_SMR),digits = 3),paste0(round(smr$beta_lci,digits = 3),"-",round(smr$beta_uci,digits = 3)),round(as.numeric(smr$p_SMR),digits = 3))
  # tidy = rbind(JAK_OR2,smr_row)
  tidy = JAK_OR2
  tidy = tidy[c(6,2,1,3,4,5),]
  tidy[1,1] = tidy[3,1] #调整“JAKx”的位置
  tidy[3,1] = tidy[2,1]
  tidy = dplyr::mutate(tidy,CI=as.character(CI))
  return(tidy)
}
JAK1_table <- table_function("JAK1")
JAK2_table <- table_function("JAK2")
JAK3_table <- table_function("JAK3")
JAK_table = rbind(JAK1_table,JAK2_table,JAK3_table)
write.csv(JAK_table,"tidy_reasult/CAD_forest_ukb-b-8714_table.csv",row.names = F,quote = F)
