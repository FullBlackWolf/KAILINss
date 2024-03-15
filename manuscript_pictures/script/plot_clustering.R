library(ggplot2)
library(reshape2)
library(BuenColors)
wid = 6
heg = 5
setwd('F:\\CJP_cellfate\\result\\clustering\\')
simulation = read.csv('F:\\CJP_cellfate\\clustering\\clustering_simulation_no-c8.csv',skip = 1,row.names = 1)
simulation = melt(t(simulation))
simulation_ari = simulation[grep('ari',simulation$Var1),]
simulation_nmi = simulation[grep('nmi',simulation$Var1),]
pdf('simulation_ari.pdf',width = wid,height = heg)
ggplot(simulation_ari,aes(y=Var2,x=value,color=Var2))+geom_boxplot()+
  theme_classic()+xlab('ari')+ylab('methods')+
  scale_color_manual(values = jdb_palette("lawhoops")[c(1:5,10:20)])+
  theme(legend.position = "none")+theme(text = element_text(size=12))
dev.off()
pdf('simulation_nmi.pdf',width = wid,height = heg)
ggplot(simulation_nmi,aes(y=Var2,x=value,color=Var2))+geom_boxplot()+
  theme_classic()+xlab('nmi')+ylab('methods')+
  scale_color_manual(values = jdb_palette("lawhoops")[c(1:5,10:20)])+
  theme(legend.position = "none")+theme(text = element_text(size=12))
dev.off()

simulation = read.csv('F:\\CJP_cellfate\\clustering\\clustering_real_data_basic_noHB.csv',row.names = 1)
simulation = melt(t(simulation))
simulation_ari = simulation[grep('ari',simulation$Var1),]
simulation_nmi = simulation[grep('nmi',simulation$Var1),]
pdf('real_ari.pdf',width = wid,height = heg)
ggplot(simulation_ari,aes(y=Var2,x=value,color=Var2))+geom_boxplot()+
  theme_classic()+xlab('ari')+ylab('methods')+
  scale_color_manual(values = jdb_palette("lawhoops")[c(1:5,10:20)])+
  theme(legend.position = "none")+theme(text = element_text(size=12))
dev.off()
pdf('real_nmi.pdf',width = wid,height = heg)
ggplot(simulation_nmi,aes(y=Var2,x=value,color=Var2))+geom_boxplot()+
  theme_classic()+xlab('nmi')+ylab('methods')+
  scale_color_manual(values = jdb_palette("lawhoops")[c(1:5,10:20)])+
  theme(legend.position = "none")+theme(text = element_text(size=12))
dev.off()

library(Seurat)
obj = readRDS('F:\\CJP_cellfate\\clustering\\real_data_use_rds\\pbmc.rds')
obj = UpdateSeuratObject(obj)
pbmc <- read_excel("F:/CJP_cellfate/clustering/all_result/clustering/real/pbmc.xlsx")
pbmc = pbmc [order(pbmc $Var3),]
obj$pred = pbmc$Predict
UMAPPlot(obj,group.by='celltype')+theme_void()+ggtitle('')+
  scale_color_manual(values = jdb_palette("corona"))
UMAPPlot(obj,group.by='pred')+theme_void()+ggtitle('')+scale_color_manual(values = jdb_palette("corona"))
