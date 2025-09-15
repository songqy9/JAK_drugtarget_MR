setwd("D:/GitHub/MR/")
library(forestploter)
library(forestplot)
# library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
library(cowplot)
# tidy for OR
tidy_function <- function(name,main_title){
  JAK = dplyr::filter(result,exposure==name)
  JAK_OR = dplyr::select(JAK,method,OR,OR_lci95,OR_uci95, p)
  # OR 添加
  # 自动保留三位并补齐
  # JAK_OR = dplyr::mutate(JAK_OR,`oddratio` = dplyr::if_else(round(JAK_OR$OR,digits = 3) == 1, paste0(round(JAK_OR$OR,digits = 3),".","0","0","0"),as.character(round(JAK_OR$OR,digits = 3))))
  # JAK_OR = dplyr::mutate(JAK_OR,`oddratio` = dplyr::if_else(nchar(round(JAK_OR$OR,digits = 3))!=5, dplyr::if_else(round(JAK_OR$OR,digits = 3) == 1,paste0(round(JAK_OR$OR,digits = 3),".","0","0","0"),paste0(round(JAK_OR$OR,digits = 3),"0")),as.character(round(JAK_OR$OR,digits = 3)))) # OR
  # JAK_OR = dplyr::mutate(JAK_OR,`oddratiolci` = dplyr::if_else(nchar(round(JAK_OR$OR_lci95,digits = 3))!=5, dplyr::if_else(round(JAK_OR$OR_lci95,digits = 3) == 1,paste0(round(JAK_OR$OR_lci95,digits = 3),".","0","0","0"),paste0(round(JAK_OR$OR_lci95,digits = 3),"0")),as.character(round(JAK_OR$OR_lci95,digits = 3)))) # OR lower CI
  # JAK_OR = dplyr::mutate(JAK_OR,`oddratiouci` = dplyr::if_else(nchar(round(JAK_OR$OR_uci95,digits = 3))!=5, dplyr::if_else(round(JAK_OR$OR_uci95,digits = 3) == 1,paste0(round(JAK_OR$OR_uci95,digits = 3),".","0","0","0"),paste0(round(JAK_OR$OR_uci95,digits = 3),"0")),as.character(round(JAK_OR$OR_uci95,digits = 3)))) # OR upper CI
  JAK_OR = dplyr::mutate(JAK_OR, `oddratio` = stringr::str_sub(OR,1,5))
  JAK_OR = dplyr::mutate(JAK_OR, `oddratiolci` = stringr::str_sub(OR_lci95,1,5))
  JAK_OR = dplyr::mutate(JAK_OR, `oddratiouci` = stringr::str_sub(OR_uci95,1,5))
  JAK_OR = dplyr::mutate(JAK_OR,`     OR(95% CI)` = paste0(oddratio,"(",oddratiolci,"~",oddratiouci,")"))
  # JAK_OR = dplyr::mutate(JAK_OR,`OR(CI)` = paste0(round(JAK_OR$OR,digits = 3),"(",round(JAK_OR$OR_lci95,digits = 3),"~",round(JAK_OR$OR_uci95,digits = 3),")"))
  JAK_OR = dplyr::mutate(JAK_OR," "=rep("                                                                                                                      ",6))# 这个空格就是留着森林图的
  JAK_OR = dplyr::mutate(JAK_OR,p = round(p,digits = 3))
  JAK_OR = dplyr::mutate(JAK_OR,p = as.character(p))
  JAK_OR = dplyr::mutate(JAK_OR,p = dplyr::if_else(p == "0", "<0.001",paste0("  ",p))) # 给小于0.001的p值赋大于号
  JAK_OR = dplyr::mutate(JAK_OR,p = dplyr::if_else(nchar(p) < 7 , dplyr::if_else(p == "<0.001",p, dplyr::if_else(p == "1", "1.000",paste0(p,"0"))),p)) # 对于不足3位的p值进行补齐 要注意上一步在前面加上的空格
  JAK_OR = dplyr::mutate(JAK_OR,p = dplyr::if_else(p == "  10","  1.000",p))# 修改p值为1的成为1.000
  JAK_OR = dplyr::relocate(JAK_OR,p,.after = last_col())
  # 修改MR method简写
  JAK_OR = dplyr::mutate(JAK_OR,method = dplyr::if_else(method=="Summary-data-based Mendelian Randomization","SMR", method))
  JAK_OR = dplyr::mutate(JAK_OR,method = dplyr::if_else(method=="Inverse variance weighted", "IVW", method))
  JAK_OR = dplyr::mutate(JAK_OR,method = dplyr::if_else(method=="Robust adjusted profile score", "RAPS", method))
  if(is.na(main_title) != TRUE){
    JAK_OR[nrow(JAK_OR)+1,] <- " "
    JAK_OR[nrow(JAK_OR)+1,] <- " "
    JAK_OR2 = JAK_OR[c(8,7,6,2,1,3,4,5),]
    nsnp = JAK[1,"nsnp"]
    title = paste0(name,"(n SNPs=",nsnp,")")
    name_order = names(JAK_OR2)
    JAK_OR2 = mutate(JAK_OR2,exposure = c(main_title,title,rep(" ",6)))
    JAK_OR2 = dplyr::select(JAK_OR2,exposure,name_order)
  }else{
    JAK_OR[nrow(JAK_OR)+1,] <- " "
    JAK_OR2 = JAK_OR[c(7,6,2,1,3,4,5),]
    nsnp = JAK[1,"nsnp"]
    title = paste0(name,"(n SNPs=",nsnp,")")
    name_order = names(JAK_OR2)
    JAK_OR2 = mutate(JAK_OR2,exposure = c(title,rep(" ",6)))
    JAK_OR2 = dplyr::select(JAK_OR2,exposure,name_order)
  }
  return(JAK_OR2)
}
# ---------------------------------------------------------------
## for stroke
result <- read.csv("result_ukb-b-8714.csv",header = T)
result_smr_jak1 <- read.delim("JAK1_Msmr.msmr")
result_smr_jak2 <- read.delim("JAK2_Msmr.msmr")
result_smr_jak3 <- read.delim("JAK3_Msmr.msmr")
result_smr_TYK2 <- read.delim("TYK2_Stroke_Msmr.msmr")
result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3)
JAK1_tidy <- tidy_function("JAK1","Stroke")
JAK2_tidy <- tidy_function("JAK2",NA)
JAK3_tidy <- tidy_function("JAK3",NA)
TYK2_tidy <- tidy_function("TYK2",NA)
# fix(JAK2_tidy) # 部分进行手动填充
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
dt <- rbind(JAK1_tidy,JAK2_tidy,JAK3_tidy,TYK2_tidy)
dt <- dt[,c(1:5,10,9,11)]
# dt <- dplyr::rename(dt,`OR(confident interval)`=estimate)
dt <- dplyr::rename(dt,`p Value`=p) %>% dplyr::rename(`MR method` = method) %>% dplyr::rename(Exposure = exposure)
# dt <- dplyr::rename(dt, "                                         stroke outcome" = " ")
numvar=c("OR" ,"OR_lci95","OR_uci95")
dt[,c(numvar)]<-lapply(dt[,numvar],as.numeric)
{# 创建空白并移到
  # dt$` ` <- paste(rep(" ", 20), # 这个空格就是留着森林图的
  #                 collapse = " ")
  # names(dt)
  # # 移动下空白行的位置
  # dt <- dt[ ,c(1:5,7,6)]
}
# ---------------------------------------------------------------

## for CAD
result <- read.csv("result_CADukb-b-15748.csv",header = T)
result_smr_jak1 <- read.delim("JAK1_CAD_Msmr.msmr")
result_smr_jak2 <- read.delim("JAK2_CAD_Msmr.msmr")
result_smr_jak3 <- read.delim("JAK3_CAD_Msmr.msmr")
result_smr_TYK2 <- read.delim("TYK2_CAD_Msmr.msmr")
result_smr <- rbind(result_smr_jak1,result_smr_jak2,result_smr_jak3,result_smr_TYK2)
JAK1_tidy <- tidy_function("JAK1","CAD")
JAK2_tidy <- tidy_function("JAK2",NA)
JAK3_tidy <- tidy_function("JAK3",NA)
TYK2_tidy <- tidy_function("TYK2",NA)
# forest
dt2 <- rbind(JAK1_tidy,JAK2_tidy,JAK3_tidy,TYK2_tidy)
dt2 <- dt2[,c(1:5,10,9,11)]
# dt <- dplyr::rename(dt,`OR(confident interval)`=estimate)
dt2 <- dplyr::rename(dt2,`p Value`=p) %>% dplyr::rename(`MR method` = method) %>% dplyr::rename(Exposure = exposure)
# dt2 <- dplyr::rename(dt2, "                                         CAD outcome" = " ")
numvar=c("OR" ,"OR_lci95","OR_uci95")
dt2[,c(numvar)]<-lapply(dt2[,numvar],as.numeric)
# fix(dt2) # 修改为“1”改为“1.000”
# ---------------------------------------------------------------
####dt combine cbind
# dt_combine <- cbind(dt,dt2)
# names(dt_combine)<-c("Exposure","MR method","OR_stroke","lci_stroke","uci_stroke","                  stroke outcome","stroke p Val","E","me","OR_CAD","lci_CAD","uci_CAD","                     CAD outcome","CAD p Val")
# dt_combine <- dplyr::select(dt_combine,2:7,10:14) %>% cbind(data.frame(Exposure = c("JAK1",rep(" ",6),"JAK2",rep(" ",6),"JAK3",rep(" ",6))),.)
# dt_combine <- dplyr::mutate(dt_combine,`stroke p Val`= paste0("   ",`stroke p Val`),`CAD p Val`=paste0("   ",`CAD p Val`))
####with forestplot
# dt2 <- rbind(colnames(dt),dt)
# p <- forestplot::forestplot(labeltext = as.matrix(dt2[,c(1,2,7)]),
#                             mean = dt2$OR,
#                             lower = dt2$OR_lci95,
#                             upper = dt2$OR_uci95,
#                             is.summary=c(T,T,F,F,F,F,F,F,T,F,F,F,F,F,F,F,T,F,F,F,F,F,F),
#                             graph.pos = 4,
#                             zero = 1,
#                             boxsize = 0.5,
#                             # xticks = c(0.96,0.98,1,1.02,1.04)
#                             )
####dt combine rbind
dt_combine <- rbind(dt,dt2)

####with forestploter
# Define theme
# 设定字体
windowsFonts(JP1 = windowsFont("Times New Roman"),
             JP2 = windowsFont("MS Gothic"),
             JP3 = windowsFont("Arial Unicode MS"),
             JP4 = windowsFont("宋体"),
             JP5 = windowsFont("微软雅黑"))
tm <- forest_theme(base_size = 15, # 更改整体字体大小
                   refline_col = "red",
                   base_family = "JP1",
                   ci_lwd = 1.5, # ci线宽
                   # 树颜色
                   # ci_col="#228B22",
                   # summary_col = "grey",
                   footnote_col = "#636363",
                   # footnote_col = "grey",
                   footnote_fontface = "italic")
# ----------------- verbal combine ---------------------------
p <- forest(
  #  筛选用到的列,8为第4列在excel为空
  dt_combine[,c(1:2, 6:8)],
  est = list(dt_combine$OR),
  lower = list(dt_combine$OR_lci95),
  upper = list(dt_combine$OR_uci95),
  sizes = 0.7,
  # 这边可以更改森林图的位置
  ci_column = 3, # 第四列为空，显示森林图
  # 数值虚线，一般OR/HR=1的位置
  ref_line = 1,
  # clip = c(-0.05,0.03),
  boxsize = 1, # 设置box 大小 不设置的话默认的会比较大
  # arrow_lab = c("Placebo Better", "Treatment Better"),# 底下两个文本
  xlim = c(0.95,1.05),# X轴刻度范围
  ticks_at = c(0.96,0.98,1,1.02,1.04),# X轴刻度间隔
  # footnote = "forest plot of JAK1, JAK2, and JAK3",# 左下角的脚注
  ci.vertices = TRUE,# 设置CI小竖线
  ci.vertices.height = 0.4,
  theme = tm)# 自定义的主题
# 对不同的组添加颜色区块
g <-  edit_plot(p, row = c(1:8), which = "background",
                gp = gpar(fill = "#f2f2f2")) #gpar是grid包中的函数 要library grid
# fill = "darkolivegreen1"

g2 <- edit_plot(g, row = c(9:15), which = "background",
                gp = gpar(fill = "white"))
g3 <- edit_plot(g2, row = c(16:22), which = "background",
                gp = gpar(fill = "#f2f2f2"))
g4 <- edit_plot(g3, row = c(23:29), which = "background",
                gp = gpar(fill = "white"))
g5 <- edit_plot(g4, row = c(30:37), which = "background",
                gp = gpar(fill = "#f2f2f2"))
g6 <- edit_plot(g5, row = c(38:44), which = "background",
                gp = gpar(fill = "white"))
g7 <- edit_plot(g6, row = c(45:51), which = "background",
                gp = gpar(fill = "#f2f2f2"))
g8 <- edit_plot(g7, row = c(52:58), which = "background",
                gp = gpar(fill = "white"))
g9 <- add_border(g8, part = "header", gp = gpar(lwd = 1.3))
g10 <- edit_plot(g9, row = c(1,30),gp = gpar(fontface = "bold")) #对应行加粗 
g11 <- edit_plot(g10, row = c(3,10,17,24,32,39,46,53),col = 3, which = "ci",gp = gpar(col = "#e31a1c")) %>% 
  edit_plot(row = c(4,11,18,25,33,40,47,54),col = 3, which = "ci",gp = gpar(col = "#a6cee3")) %>%
  edit_plot(row = c(5,12,19,26,34,41,48,55),col = 3, which = "ci",gp = gpar(col = "#1f78b4")) %>%
  edit_plot(row = c(6,13,20,27,35,42,49,56),col = 3, which = "ci",gp = gpar(col = "#fb9a99")) %>%
  edit_plot(row = c(7,14,21,28,36,43,50,57),col = 3, which = "ci",gp = gpar(col = "#33a02c")) %>%
  edit_plot(row = c(8,15,22,29,37,44,51,58),col = 3, which = "ci",gp = gpar(col = "#b2df8a"))
g12 <- edit_plot(g11,col = 5, row = c(10,11,12,13,14,15,40,41,42,44),gp = gpar(fontface = "bold"))
g13 <- edit_plot(g12,col = 1, row = c(2,9,16,23,31,38,45,52),gp = gpar(fontface = "italic"))
# ----------------------左右图----------------------------------
# # p <- forest(
# #   #  筛选用到的列,8为第4列在excel为空
# #   dt_combine[,c(1:2, 6:7,11:12)], 
# #   est = list(dt_combine$OR_stroke,
# #              dt_combine$OR_CAD),
# #   lower = list(dt_combine$lci_stroke,
# #                dt_combine$lci_CAD),
# #   upper = list(dt_combine$uci_stroke,
# #                dt_combine$uci_CAD),
# #   #sizes = dt$se,
# #   # 这边可以更改森林图的位置
# #   ci_column = c(3,5), # 第四列为空，显示森林图
# #   # 数值虚线，一般OR/HR=1的位置
# #   ref_line = 1,
# #   # clip = c(-0.05,0.03),
# #   boxsize = .5, # 设置box 大小 不设置的话默认的会比较大
# #   # arrow_lab = c("Placebo Better", "Treatment Better"),# 底下两个文本
# #   xlim = c(0.95,1.05),# X轴刻度范围
# #   ticks_at = c(0.96,0.98,1,1.02,1.04),# X轴刻度间隔
# #   # footnote = "forest plot of JAK1, JAK2, and JAK3",# 左下角的脚注
# #   ci.vertices = TRUE,# 设置CI小竖线
# #   ci.vertices.height = 0.4,
# #   theme = tm)# 自定义的主题
# # # 对不同的组添加颜色区块
# # g <-  edit_plot(p, row = c(1:7), which = "background",
# #                 gp = gpar(fill = "white")) #gpar是grid包中的函数 要library grid
# # # fill = "darkolivegreen1"
# # 
# # g2 <- edit_plot(g, row = c(8:14), which = "background",
# #                 gp = gpar(fill = "#e5e5e5"))
# # g3 <- edit_plot(g2, row = c(15:21), which = "background",
# #                 gp = gpar(fill = "white"))
# # g4 <- add_border(g3, part = "header", gp = gpar(lwd = 1.5))
# #--------------上下图---------------------------------------
# # for stroke
# p1 <- forest(
#   #  筛选用到的列,8为第4列在excel为空
#   dt[,c(1:2, 6:8)], 
#   est = list(dt$OR),
#   lower = list(dt$OR_lci95),
#   upper = list(dt$OR_uci95),
#   sizes = 0.7,
#   # 这边可以更改森林图的位置
#   ci_column = 3, # 第四列为空，显示森林图
#   # 数值虚线，一般OR/HR=1的位置
#   ref_line = 1,
#   # clip = c(-0.05,0.03),
#   boxsize = 1, # 设置box 大小 不设置的话默认的会比较大
#   # arrow_lab = c("Placebo Better", "Treatment Better"),# 底下两个文本
#   xlim = c(0.95,1.05),# X轴刻度范围
#   ticks_at = c(0.96,0.98,1,1.02,1.04),# X轴刻度间隔
#   # footnote = "forest plot of JAK1, JAK2, and JAK3",# 左下角的脚注
#   ci.vertices = TRUE,# 设置CI小竖线
#   ci.vertices.height = 0.4,
#   theme = tm)# 自定义的主题
# # 对不同的组添加颜色区块
# g_stroke <-  edit_plot(p1, row = c(1:7), which = "background",
#                 gp = gpar(fill = "white")) #gpar是grid包中的函数 要library grid
# # fill = "darkolivegreen1"
# 
# g2_stroke <- edit_plot(g_stroke, row = c(8:14), which = "background",
#                 gp = gpar(fill = "#e5e5e5"))
# g3_stroke <- edit_plot(g2_stroke, row = c(15:21), which = "background",
#                 gp = gpar(fill = "white"))
# g4_stroke <- add_border(g3_stroke, part = "header", gp = gpar(lwd = 1.5))
# # eoffice::topptx(filename = "plot导出的.pptx") # Plot作图
# # for CAD
# p2 <- forest(
#   #  筛选用到的列,8为第4列在excel为空
#   dt2[,c(1:2, 6:8)], 
#   est = list(dt2$OR),
#   lower = list(dt2$OR_lci95),
#   upper = list(dt2$OR_uci95),
#   sizes = 0.7,
#   # 这边可以更改森林图的位置
#   ci_column = 3, # 第四列为空，显示森林图
#   # 数值虚线，一般OR/HR=1的位置
#   ref_line = 1,
#   # clip = c(-0.05,0.03),
#   boxsize = 1, # 设置box 大小 不设置的话默认的会比较大
#   # arrow_lab = c("Placebo Better", "Treatment Better"),# 底下两个文本
#   xlim = c(0.95,1.05),# X轴刻度范围
#   ticks_at = c(0.96,0.98,1,1.02,1.04),# X轴刻度间隔
#   # footnote = "forest plot of JAK1, JAK2, and JAK3",# 左下角的脚注
#   ci.vertices = TRUE,# 设置CI小竖线
#   ci.vertices.height = 0.4,
#   theme = tm)# 自定义的主题
# # 对不同的组添加颜色区块
# g_CAD <-  edit_plot(p2, row = c(1:7), which = "background",
#                        gp = gpar(fill = "white")) #gpar是grid包中的函数 要library grid
# # fill = "darkolivegreen1"
# 
# g2_CAD <- edit_plot(g_CAD, row = c(8:14), which = "background",
#                        gp = gpar(fill = "#e5e5e5"))
# g3_CAD <- edit_plot(g2_CAD, row = c(15:21), which = "background",
#                        gp = gpar(fill = "white"))
# g4_CAD <- add_border(g3_CAD, part = "header", gp = gpar(lwd = 1.5))
# # eoffice::topptx(filename = "plot导出的.pptx") # Plot作图
# cow <- cowplot::plot_grid(g4_stroke,g4_CAD,labels = c("(a)","(b)"),label_fontfamily = "serif",nrow = 2,ncol = 1)
tiff(filename = "tidy_reasult/beta_forestplot.tiff",
     width = 12.7,height = 16,units = "in",
     # pointsize = 5,
     res = 300)
par(cex = 2);
par(mar = c(4,4,4,4),family = "JP1")
g13
dev.off()

