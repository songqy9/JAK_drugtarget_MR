setwd("D:/GitHub/MR/")
library(dplyr)

# JAK2 <- read.table("JAK2.txt",header = F)
# JAK_myquery <- function(name){
#   JAK = read.table(paste0(name,".txt"),header = F)
#   JAK = dplyr::select(JAK,V2,V3,V4,V5,V6,V8,V10,)
# }

ldsc_gwas <- read.csv("JAK2_snp.csv",header = T)
ldsc_gwas <- dplyr::filter(ldsc_gwas,pval.exposure <= 5e-08)
ldsc_gwas<-clump_data(ldsc_gwas,clump_kb = 100,
                      clump_r2 = 0.3,
                      pop = "EUR")
ldsc_gwas<-subset(ldsc_gwas,eaf.exposure>0.01)
ldsc_select <- dplyr::select(ldsc_gwas,SNP,chr.exposure,pos.exposure,effect_allele.exposure,other_allele.exposure,beta.exposure,pval.exposure) %>% 
  mutate(beta.exposure=exp(beta.exposure))
names(ldsc_select)<-c("SNPid","chr","bp","a1","a2","or","pval")
write.table(ldsc_select,"D:/Ubuntu/rootfs/home/songqy9/downloaddir/ldsc/data/trait1_JAKplasm.txt",row.names = F,quote = F)
