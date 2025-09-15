setwd("D:/GitHub/MR/")
# library(TwoSampleMR)
library(dplyr)
library(tidyr)
library(readr)
library(openxlsx)
library(stringr)
library(data.table)
library(ggrepel)
library(ggplot2)
# stroke
JAK1_IV_stroke <- read.csv("JAK1_SNP_asIV.csv",header = T)
JAK2_IV_stroke <- read.csv("JAK2_SNP_asIV.csv",header = T)
JAK3_IV_stroke <- read.csv("JAK3_SNP_asIV.csv",header = T)
TYK2_IV_stroke <- read.csv("TYK2_SNP_asIV.csv",header = T)
# CAD
JAK1_IV_CAD <- read.csv("JAK1_SNP_asIV_CAD.csv",header = T)
JAK2_IV_CAD <- read.csv("JAK2_SNP_asIV_CAD.csv",header = T)
JAK3_IV_CAD <- read.csv("JAK3_SNP_asIV_CAD.csv",header = T)
TYK2_IV_CAD <- read.csv("TYK2_SNP_asIV_CAD.csv",header = T)
# merge JAK2 ---------------------------------------------
JAK2_merge <- merge(JAK2_IV_stroke,JAK2_IV_CAD,by = "SNP") %>% 
  dplyr::select(SNP,
                chr.x,
                pos.x,
                effect_allele.exposure.x,
                other_allele.exposure.x,
                eaf.exposure.x,
                beta.exposure.x,
                pval.exposure.x,
                beta.outcome.x,
                pval.outcome.x,
                beta.outcome.y,
                pval.outcome.y)
names(JAK2_merge) <- c("SNP","Chr","position","EA","OA","Freq",
                       "Beta1","p1 value",
                       "Beta2","p2 value",
                       "Beta3","p3 value")
JAK2_merge_adjust <- dplyr::mutate(JAK2_merge,Chr = paste0("Chr",Chr))%>% 
  tidyr::unite(Chr,position,col = "Chr:Position") %>% 
  tidyr::unite(SNP,EA,col = "rsID_EA") %>% 
  dplyr::select(`Chr:Position`,`rsID_EA`,
                Beta1,`p1 value`,
                Beta2,`p2 value`,
                Beta3,`p3 value`)
# write.csv(JAK2_merge_adjust,"tidy_reasult/JAK2_commonSNP.csv",row.names = F,quote = F)
# digit adjust
JAK2_merge_adjust <- read.csv("tidy_reasult/JAK2_commonSNP.csv",header = T) %>% 
  dplyr::mutate(OR1 = exp(Beta1)) %>% 
  dplyr::mutate(OR2 = exp(Beta2)) %>% 
  dplyr::mutate(OR3 = exp(Beta3))
JAK2_adjust_num <- JAK2_adjust_num <- dplyr::mutate(JAK2_merge_adjust,Beta1 = sprintf("%0.2f",round(Beta1,digits = 2))) %>% 
     dplyr::mutate(p1.value = signif(p1.value,digits = 3)) %>% 
     dplyr::mutate(Beta2 = signif(Beta2,digits = 2)) %>% 
     dplyr::mutate(p2.value = sprintf("%0.3f",signif(p2.value,digits = 3))) %>% 
     dplyr::mutate(Beta3 = signif(Beta3,digits = 2)) %>% 
     dplyr::mutate(p3.value = sprintf("%0.3f",signif(p3.value,digits = 3)))
JAK2_adjust_num2 <- tidyr::separate(JAK2_adjust_num,col = `Chr.Position`,into = c("Chr","Position"),sep = "_")
# write.table(JAK2_adjust_num,"tidy_reasult/JAK2_commonSNP_adjustdigits.txt",row.names = F,quote = F) 
# JAK2_merge_inGene <- dplyr::filter(JAK2_adjust_num2, (Position <= 5129948 & Position >= 4984390))
# JAK2_merge_outGene <- dplyr::filter(JAK2_adjust_num2, (Position > 5129948 | Position < 4984390))
# merge JAK ----------------------------------------------------------
merge_JAK <- function(JAK_IV_stroke,JAK_IV_CAD,name){
  JAK_merge = merge(JAK_IV_stroke,JAK_IV_CAD,by = "SNP")
  JAK_merge = dplyr::select(JAK_merge,SNP,
                                      chr.x,
                                      pos.x,
                                      effect_allele.exposure.x,
                other_allele.exposure.x,
                eaf.exposure.x,
                beta.exposure.x,
                pval.exposure.x,
                beta.outcome.x,
                pval.outcome.x,
                beta.outcome.y,
                pval.outcome.y)
names(JAK_merge) = c("SNP","Chr","position","EA","OA","Freq",
                       "Beta1","p1 value",
                       "Beta2","p2 value",
                       "Beta3","p3 value")
JAK_merge_adjust = dplyr::mutate(JAK_merge,Chr = paste0("Chr",Chr))
JAK_merge_adjust =  tidyr::unite(JAK_merge_adjust,Chr,position,col = "Chr:Position")
JAK_merge_adjust =  tidyr::unite(JAK_merge_adjust,SNP,EA,col = "rsID_EA")
JAK_merge_adjust =  dplyr::select(JAK_merge_adjust,`Chr:Position`,`rsID_EA`,
                                  Beta1,`p1 value`,
                                  Beta2,`p2 value`,
                                  Beta3,`p3 value`)
write.csv(JAK_merge_adjust,paste0("tidy_reasult/",name,"_commonSNP.csv"),row.names = F,quote = F)
   return(JAK_merge_adjust)
}
JAK1_merge_adjust <- merge_JAK(JAK1_IV_stroke,JAK1_IV_CAD,"JAK1")
JAK3_merge_adjust <- merge_JAK(JAK3_IV_stroke,JAK3_IV_CAD,"JAK3")
TYK2_merge_adjust <- merge_JAK(TYK2_IV_stroke,TYK2_IV_CAD,"TYK2")
# digit adjust-----------------------------------------
JAK2_merge_adjust <- read.csv("tidy_reasult/JAK2_commonSNP.csv",header = T) %>% 
  dplyr::mutate(OR1 = exp(Beta1)) %>% 
  dplyr::mutate(OR2 = exp(Beta2)) %>% 
  dplyr::mutate(OR3 = exp(Beta3))
JAK2_adjust_num <- JAK2_adjust_num <- dplyr::mutate(JAK2_merge_adjust,Beta1 = sprintf("%0.2f",round(Beta1,digits = 2))) %>% 
  dplyr::mutate(p1.value = signif(p1.value,digits = 3)) %>% 
  dplyr::mutate(Beta2 = signif(Beta2,digits = 2)) %>% 
  dplyr::mutate(p2.value = sprintf("%0.3f",signif(p2.value,digits = 3))) %>% 
  dplyr::mutate(Beta3 = signif(Beta3,digits = 2)) %>% 
  dplyr::mutate(p3.value = sprintf("%0.3f",signif(p3.value,digits = 3)))
JAK2_adjust_num2 <- tidyr::separate(JAK2_adjust_num,col = `Chr.Position`,into = c("Chr","Position"),sep = "_")
write.table(JAK2_adjust_num,"tidy_reasult/JAK2_commonSNP_adjustdigits.txt",row.names = F,quote = F) 


# -------------------plot in the Chr9----------------------
JAK2_merge_data <- read.csv("tidy_reasult/JAK2_commonSNP.csv",header = T) %>% tidyr::separate(col = `rsID_EA`,into = c("rsID","EA"),sep = "_")
JAK2_eQTL <- data.table::fread("JAK2.txt",header = F) %>% .[,c(1,2,4)]
names(JAK2_eQTL) <- c("pval","rsID","Position")
JAK2_rsstroke <- read.csv("JAK2_SNP_asIV.csv",header = T)
JAK2_rsCAD <- read.csv("JAK2_SNP_asIV_CAD.csv",header = T)
JAK2_eQTL <- dplyr::mutate(JAK2_eQTL,SNP = dplyr::if_else((rsID %in% JAK2_merge_data$rsID),"Common IVs",dplyr::if_else((rsID %in% JAK2_rsstroke$SNP),"Unique IVs",dplyr::if_else((rsID %in% JAK2_rsCAD$SNP), "Unique IVs of CAD","Not IVs")))) %>% 
  dplyr::mutate(note = dplyr::if_else((Position <= 5129948 & Position >= 4984390),dplyr::if_else((rsID %in% JAK2_merge_data$rsID),"note","none"),"none")) #要标注的点
JAK2_eQTL <- dplyr::mutate(JAK2_eQTL,SNP = factor(JAK2_eQTL$SNP,levels = c("Common IVs","Unique IVs","Unique IVs of CAD","Not IVs")))
# set fontform
windowsFonts(JP1 = windowsFont("Times New Roman"),
             JP2 = windowsFont("MS Gothic"),
             JP3 = windowsFont("Arial Unicode MS"),
             JP4 = windowsFont("宋体"),
             JP5 = windowsFont("微软雅黑"))
p <- ggplot2::ggplot(data = JAK2_eQTL,aes(x = (Position/1000000),y = (-log10(pval)),color = SNP)) +
  ggplot2::geom_point(size = 2,alpha = 0.5)+
  scale_color_manual(values = c("red","green","#cecece"))+ # 调节点的颜
  scale_x_continuous(breaks = scales::breaks_pretty(5)) + # 调节x轴ticks
  # scale_shape_manual(values = c('Common SNP'=16,'Other SNP'=16))+ # 调节图例图标形状
  geom_vline(aes(xintercept = 4.984390),color = "#0000ff",linetype="solid",size = 0.25) +
  geom_vline(aes(xintercept = 5.129948),color = "blue",linetype="solid",size = 0.25) +
  geom_vline(aes(xintercept = 3.984390),color = "#0000ff",linetype="dashed") +
  geom_vline(aes(xintercept = 6.129948),color = "blue",linetype="dashed") +
  geom_segment(aes(x = 4.974390, y = 148,
               xend = 3.994390, yend = 148),
               size = 0.1,
               color = "#000062",
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_segment(aes(x = 5.139948, y = 148,
                   xend = 6.119948, yend = 148),
               size = 0.1,
               color = "#000062",
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_segment(aes(x = 4.994390, y = 148,
                   xend = 5.119948, yend = 148),
               # lineend = "round",
               size = 1,
               color = "#000062",
               arrow = arrow(length = unit(0.05, "inches"))) +
  annotate("text", x = 5.057169, y = 152, label = "JAK2",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  annotate("text", x =  4.73439, y = 152, label = "0.5Mb",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  annotate("text", x = 5.379948, y = 152, label = "0.5Mb",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  # geom_repel有点问题 会导致图例出错 要调整show.legend
  # geom_text_repel(data = dplyr::filter(JAK2_eQTL,note == "note"),aes(label = rsID),
  #                 box.padding = 0.5, #字到点的距离
  #                 point.padding = 0.5, #字到点的距离
  #                 segment.color = "red",
  #                 min.segment.length = 0, #加线段
  #                 size = 3.5,
  #                 # fontface = "bold",
  #                 show.legend = FALSE) + #去掉不是NA 是FALSE
  # geom_rect(aes(xmin=4.984390, xmax=5.129948, ymin=-Inf, ymax=Inf),fill='#d9d9d9',alpha = 0.002,color = NA)+
  xlab("Position on Chr9 (Mb)") + ylab("-log10(p-val)") +
  ggplot2::theme(text = element_text(family = "serif"),
                 # legend.position = c(0.9,0.8),
                 legend.direction = "horizontal",
                 legend.position = "top",
                 # legend.position= "none",
                 # legend.title = element_text("SNP"),
                 legend.box.background = element_blank(),
                 legend.key = element_blank(), #图例背景去掉
                 legend.title = element_blank(),
                 panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
p
# label手动调节吧
tiff(filename = paste0("tidy_reasult/JAKchrlocation_legend.tiff"),
     width = 6.4,height = 4.4,
     units = "in",
     # pointsize = 2,
     res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
p
dev.off()
p2 <- ggplot2::ggplot(data = JAK2_eQTL,aes(x = (Position/1000000),y = (-log10(pval)),color = SNP)) +
  ggplot2::geom_point(size = 2,alpha = 0.5)+
  scale_color_manual(values = c("red","green","#cecece"))+ # 调节点的颜色
  scale_x_continuous(breaks = scales::breaks_pretty(4)) + # 调节x轴ticks
  # scale_shape_manual(values = c('Common SNP'=16,'Other SNP'=16))+ # 调节图例图标形状
  geom_vline(aes(xintercept = 4.984390),color = "#0000ff",linetype="solid",size = 0.25) +
  geom_vline(aes(xintercept = 5.129948),color = "blue",linetype="solid",size = 0.25) +
  geom_vline(aes(xintercept = 4.054390),color = "#0000ff",linetype="dashed") +
  geom_vline(aes(xintercept = 6.059948),color = "blue",linetype="dashed") +
  geom_hline(aes(yintercept = -log10(5e-05)),color = "blue",linetype = "dashed")+
  geom_segment(aes(x = 4.974390, y = 148,
                   xend = 4.064390, yend = 148),
               size = 0.1,
               color = "#000062",
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_segment(aes(x = 5.139948, y = 148,
                   xend = 6.049948, yend = 148),
               size = 0.1,
               color = "#000062",
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_segment(aes(x = 4.994390, y = 148,
                   xend = 5.119948, yend = 148),
               # lineend = "round",
               size = 1,
               color = "#000062",
               arrow = arrow(length = unit(0.05, "inches"))) +
  annotate("text", x = 5.057169, y = 153, label = "JAK2",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  annotate("text", x = 4.51939, y = 153, label = "1Mb",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  annotate("text", x = 5.594948, y = 153, label = "1Mb",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  # geom_repel有点问题 会导致图例出错 要调整show.legend
  # geom_text_repel(data = dplyr::filter(JAK2_eQTL,note == "note"),aes(label = rsID),
  #                 box.padding = 0.5, #字到点的距离
  #                 point.padding = 0.5, #字到点的距离
  #                 segment.color = "red",
  #                 min.segment.length = 0, #加线段
  #                 size = 3.5,
  #                 # fontface = "bold",
  #                 show.legend = FALSE) + #去掉不是NA 是FALSE
  # geom_rect(aes(xmin=4.984390, xmax=5.129948, ymin=-Inf, ymax=Inf),fill='#d9d9d9',alpha = 0.002,color = NA)+
  xlab("Position on Chr9 (Mb)") + ylab("-log10(p-val)") +
  ggplot2::theme(text = element_text(family = "serif", size = 16),
                 # legend.position = c(0.9,0.8),
                 # legend.direction = "horizontal",
                 # legend.position = "top",
                 legend.position= "none",
                 # legend.title = element_text("SNP"),
                 # legend.box.background = element_blank(),
                 # legend.key = element_blank(), #图例背景去掉
                 # legend.title = element_blank(),
                 panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
p2 #without legend
# JAK plot --------------------------------------------------------
JAK_plot <- function(name,start_pos,end_pos){
JAK_merge_data = read.csv(paste0("tidy_reasult/",name,"_commonSNP.csv"),header = T)
JAK_merge_data = tidyr::separate(JAK_merge_data,col = `rsID_EA`,into = c("rsID","EA"),sep = "_")
JAK_eQTL = data.table::fread(paste0(name,".txt"),header = F)
JAK_eQTL = JAK_eQTL[,c(1,2,4)]
names(JAK_eQTL) = c("pval","rsID","Position")
JAK_rsstroke = read.csv(paste0(name,"_SNP_asIV.csv"),header = T)
JAK_rsCAD = read.csv(paste0(name,"_SNP_asIV_CAD.csv"),header = T)
JAK_eQTL = dplyr::mutate(JAK_eQTL,SNP = dplyr::if_else((rsID %in% JAK_merge_data$rsID),"Common IVs",dplyr::if_else((rsID %in% JAK_rsstroke$SNP),"Unique IVs",dplyr::if_else((rsID %in% JAK_rsCAD$SNP), "Unique IVs","Not IVs"))))
JAK_eQTL = dplyr::mutate(JAK_eQTL,note = dplyr::if_else((Position <= end_pos & Position >= start_pos),dplyr::if_else((rsID %in% JAK_merge_data$rsID),"note","none"),"none")) #要标注的点
JAK_eQTL <- dplyr::mutate(JAK_eQTL,SNP = factor(JAK_eQTL$SNP,levels = c("Common IVs","Unique IVs","Not IVs")))
return(JAK_eQTL)
}
JAK1_eQTL <- JAK_plot("JAK1",start_pos = 65298912,end_pos = 65533429) #转录顺序为反向的
JAK3_eQTL <- JAK_plot("JAK3",start_pos = 17935591, end_pos = 17958791) #转录顺序为反向的
TYK2_eQTL <- JAK_plot("TYK2",start_pos = 10461209,end_pos = 10491248) #转录顺序为反向的
#要标注的点
# set fontform
windowsFonts(JP1 = windowsFont("Times New Roman"),
             JP2 = windowsFont("MS Gothic"),
             JP3 = windowsFont("Arial Unicode MS"),
             JP4 = windowsFont("宋体"),
             JP5 = windowsFont("微软雅黑"))
#------------- JAK1 false plot function-------------------------------------
# JAK_eQTL = JAK1_eQTL 
# start_pos = 65533429/1000000 # 注意基因的转录顺序
# end_pos = 65298912/1000000
# name = "JAK1"
# Chr = "Chr1"
# y_top = 13
# y_text_top = 13.3
# 记得改p1 p3
plot_charplot <- function(JAK_eQTL,start_pos,end_pos,name,Chr,y_top,y_text_top){
pp=ggplot2::ggplot(data = JAK_eQTL,aes(x = (Position/1000000),y = (-log10(pval)),color = SNP)) +
  ggplot2::geom_point(size = 2,alpha = 0.5)+
  scale_color_manual(values = c("red","#cecece"))+ # 调节点的颜色
  scale_x_continuous(breaks = scales::breaks_pretty(4)) + # 调节x轴ticks
  # scale_shape_manual(values = c('Common SNP'=16,'Other SNP'=16))+ # 调节图例图标形状
  geom_vline(aes(xintercept = start_pos),color = "#0000ff",linetype="solid",size = 0.25) +
  geom_vline(aes(xintercept = end_pos),color = "blue",linetype="solid",size = 0.25) +
  geom_vline(aes(xintercept = (start_pos + 0.855)),color = "#0000ff",linetype="dashed") +
  geom_vline(aes(xintercept = (end_pos - 0.855)),color = "blue",linetype="dashed") +
  geom_hline(aes(yintercept = -log10(5e-05)),color = "blue",linetype = "dashed")+
  geom_segment(aes(x = (end_pos - 0.01), y = y_top,
                   xend = (end_pos - 0.855 + 0.01), yend = y_top),
               size = 0.1,
               color = "#000062",
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_segment(aes(x = (start_pos + 0.01), y = y_top,
                   xend = (start_pos + 0.855 - 0.01), yend = y_top),
               size = 0.1,
               color = "#000062",
               arrow = arrow(length = unit(0.1, "inches"))) +
  geom_segment(aes(x = (start_pos - 0.01), y = y_top,
                   xend = (end_pos + 0.01), yend = y_top),
               # lineend = "round",
               size = 1,
               color = "#000062",
               arrow = arrow(length = unit(0.05, "inches"))) +
  annotate("text", x = (((start_pos - 0.01)+(end_pos + 0.01))/2), y = y_text_top, label = name,
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  annotate("text", x =  (((end_pos - 0.01)+(end_pos - 0.855 + 0.01))/2), y = y_text_top, label = "1Mb",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  annotate("text", x = (((start_pos + 0.01)+(start_pos + 0.855 - 0.01))/2), y = y_text_top, label = "1Mb",
           family = "serif", fontface = "italic",
           colour = "#000062", size = 3.5) + 
  # geom_repel有点问题 会导致图例出错 要调整show.legend
  # geom_text_repel(data = dplyr::filter(JAK2_eQTL,note == "note"),aes(label = rsID),
  #                 box.padding = 0.5, #字到点的距离
  #                 point.padding = 0.5, #字到点的距离
  #                 segment.color = "red",
  #                 min.segment.length = 0, #加线段
  #                 size = 3.5,
  #                 # fontface = "bold",
  #                 show.legend = FALSE) + #去掉不是NA 是FALSE
  # geom_rect(aes(xmin=4.984390, xmax=5.129948, ymin=-Inf, ymax=Inf),fill='#d9d9d9',alpha = 0.002,color = NA)+
  xlab(paste0("Position on ",Chr," (Mb)")) + ylab("-log10(p-val)") +
  ggplot2::theme(text = element_text(family = "serif", size = 16),
                 # legend.position = c(0.9,0.8),
                 # legend.direction = "horizontal",
                 # legend.position = "top",
                 # legend.title = element_text("SNP"),
                 # legend.box.background = element_blank(),
                 # legend.key = element_blank(), #图例背景去掉
                 # legend.title = element_blank(),
                 legend.position= "none",
                 panel.grid.major =element_blank(),
                 panel.grid.minor = element_blank(),
                 panel.background = element_blank(),
                 axis.line = element_line(colour = "black"))
return(pp)}
p1 <- plot_charplot(JAK_eQTL = JAK1_eQTL, 
                    start_pos = 65533429/1000000 ,# 注意基因的转录顺序
                    end_pos = 65298912/1000000,
                    name = "JAK1",
                    Chr = "Chr1",
                    y_top = 13,
                    y_text_top = 13.5)
#-------------JAK3 false plot function-------------------------------------
# JAK_eQTL = JAK3_eQTL 
# start_pos = 17958791/1000000 # 注意基因的转录顺序
# end_pos = 17935591/1000000
# name = "JAK3"
# Chr = "Chr19"
# y_top = 9
# y_text_top = 9.3
p3 <- plot_charplot(JAK_eQTL = JAK3_eQTL, 
                    start_pos = 17958791/1000000, # 注意基因的转录顺序
                    end_pos = 17935591/1000000,
                    name = "JAK3",
                    Chr = "Chr19",
                    y_top = 9,
                    y_text_top = 9.3)
# 记得改p1 p3
#-------------TYK2 false plot function-------------------------------------
plot_charplot_TYK2 <- function(JAK_eQTL,start_pos,end_pos,name,Chr,y_top,y_text_top){
  pp=ggplot2::ggplot(data = JAK_eQTL,aes(x = (Position/1000000),y = (-log10(pval)),color = SNP)) +
    ggplot2::geom_point(size = 2,alpha = 0.5)+
    scale_color_manual(values = c("red","green","#cecece"))+ # 调节点的颜色
    scale_x_continuous(breaks = scales::breaks_pretty(4)) + # 调节x轴ticks
    # scale_shape_manual(values = c('Common SNP'=16,'Other SNP'=16))+ # 调节图例图标形状
    geom_vline(aes(xintercept = start_pos),color = "#0000ff",linetype="solid",size = 0.25) +
    geom_vline(aes(xintercept = end_pos),color = "blue",linetype="solid",size = 0.25) +
    geom_vline(aes(xintercept = (start_pos + 0.855)),color = "#0000ff",linetype="dashed") +
    geom_vline(aes(xintercept = (end_pos - 0.855)),color = "blue",linetype="dashed") +
    geom_hline(aes(yintercept = -log10(5e-05)),color = "blue",linetype = "dashed")+
    geom_segment(aes(x = (end_pos - 0.01), y = y_top,
                     xend = (end_pos - 0.855 + 0.01), yend = y_top),
                 size = 0.1,
                 color = "#000062",
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_segment(aes(x = (start_pos + 0.01), y = y_top,
                     xend = (start_pos + 0.855 - 0.01), yend = y_top),
                 size = 0.1,
                 color = "#000062",
                 arrow = arrow(length = unit(0.1, "inches"))) +
    geom_segment(aes(x = (start_pos - 0.01), y = y_top,
                     xend = (end_pos + 0.01), yend = y_top),
                 # lineend = "round",
                 size = 1,
                 color = "#000062",
                 arrow = arrow(length = unit(0.05, "inches"))) +
    annotate("text", x = (((start_pos - 0.01)+(end_pos + 0.01))/2), y = y_text_top, label = name,
             family = "serif", fontface = "italic",
             colour = "#000062", size = 3.5) + 
    annotate("text", x =  (((end_pos - 0.01)+(end_pos - 0.855 + 0.01))/2), y = y_text_top, label = "1Mb",
             family = "serif", fontface = "italic",
             colour = "#000062", size = 3.5) + 
    annotate("text", x = (((start_pos + 0.01)+(start_pos + 0.855 - 0.01))/2), y = y_text_top, label = "1Mb",
             family = "serif", fontface = "italic",
             colour = "#000062", size = 3.5) + 
    # geom_repel有点问题 会导致图例出错 要调整show.legend
    # geom_text_repel(data = dplyr::filter(JAK2_eQTL,note == "note"),aes(label = rsID),
    #                 box.padding = 0.5, #字到点的距离
    #                 point.padding = 0.5, #字到点的距离
    #                 segment.color = "red",
    #                 min.segment.length = 0, #加线段
    #                 size = 3.5,
    #                 # fontface = "bold",
    #                 show.legend = FALSE) + #去掉不是NA 是FALSE
    # geom_rect(aes(xmin=4.984390, xmax=5.129948, ymin=-Inf, ymax=Inf),fill='#d9d9d9',alpha = 0.002,color = NA)+
    xlab(paste0("Position on ",Chr," (Mb)")) + ylab("-log10(p-val)") +
    ggplot2::theme(text = element_text(family = "serif", size = 16),
                   # legend.position = c(0.9,0.8),
                   # legend.direction = "horizontal",
                   # legend.position = "top",
                   # legend.title = element_text("SNP"),
                   # legend.box.background = element_blank(),
                   # legend.key = element_blank(), #图例背景去掉
                   # legend.title = element_blank(),
                   legend.position= "none",
                   panel.grid.major =element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_blank(),
                   axis.line = element_line(colour = "black"))
  return(pp)}
# JAK_eQTL = TYK2_eQTL 
# start_pos = 10491248/1000000 # 
# end_pos = 10461209/1000000
# name = "TYK2"
# Chr = "Chr19"
# y_top = 9
# y_text_top = 9.3
p4 <- plot_charplot_TYK2(JAK_eQTL = TYK2_eQTL, 
                    start_pos = 10491248/1000000, # 注意基因的转录顺序
                    end_pos = 10461209/1000000,
                    name = "TYK2",
                    Chr = "Chr19",
                    y_top = 176,
                    y_text_top = 180.6)

# cowplot -------------------------------
pcow <- cowplot::plot_grid(p1,p2,p3,p4,labels = c("(a)","(b)","(c)","(d)"),label_fontfamily = "serif",nrow = 2,ncol = 2)
tiff(filename = paste0("tidy_reasult/JAKcombine_chrpostion_eQTL.tiff"),
     width = 8,height = 8,units = "in",
     pointsize = 5,res = 300)
par(cex = 2);
par(mar = c(4,4,4,4))
pcow
dev.off()
















