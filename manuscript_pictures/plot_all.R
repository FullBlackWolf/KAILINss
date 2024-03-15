library(ggplot2)
library(pheatmap)
library(BuenColors)
setwd('F:\\CJP_cellfate\\result\\cell_fate_pred')
###########################
### larry
###########################
### pheatmap
larry_cluster <- read.csv("F:/CJP_cellfate/zoxueganxibao_zuizhong.csv",
                          header=FALSE)
larry_cluster[larry_cluster>0.0006]=0.0006
larry_cluster[larry_cluster<0.0002]=0.0002
colnames(larry_cluster) = paste0('Cluster',1:33)
rownames(larry_cluster) = paste0('Cluster',1:33)
pdf(file = 'larry_cluster.pdf',width = 6,height = 5)
pheatmap(larry_cluster,cluster_rows = F,cluster_cols = F,border_color = NA,
         color = colorRampPalette(c("#FAFAFC","#3F96B7"))(50))
dev.off()
### acc col plot
df_acc = data.frame('accuracy'=c(0.80,0.77,0.69,0.61,0.41,0.58),
                    'methods'=c('KAILIN',
                                'cospar_all_clone',
                                'cospar_state_info',
                                'Waddington-OT',
                                'K-means',
                                'hclust'))
pdf(file = 'larry_acc.pdf',width = 6,height = 5)
ggplot(df_acc,aes(y=accuracy,x=reorder(methods,accuracy,order = F),fill=methods))+geom_col(width = 0.9)+
  theme_classic()+scale_y_continuous(expand = c(0,0))+
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1:4,6:8)])+xlab('')+
  theme(axis.text.x = element_text(size=12,vjust = 0.5, hjust = 0.5,
                                   angle = 35))
dev.off()
### embedding
larry_emb <- read.csv("F:/CJP_cellfate/result/cell_fate_pred/larry_emb.csv",
                      header=T,row.names = 1)
pdf(file = 'larry_gt.pdf',width = 6,height = 5)
ggplot(larry_emb,aes(x=X0,y=X1,color=label))+geom_point() + theme_void() +
  scale_color_manual(values = c(jdb_palette("brewer_spectra")[c(1,8)],
                                'grey'))+theme(legend.position="none")
dev.off()
pdf(file = 'larry_pred.pdf',width = 6,height = 5)
ggplot(larry_emb,aes(x=X0,y=X1,color=pred))+geom_point() + theme_void() +
  scale_color_manual(values = c(jdb_palette("brewer_spectra")[c(1,8)],
                                'grey'))+theme(legend.position="none")
dev.off()
pdf(file = 'larry_descrip.pdf',width = 6,height = 5)
ggplot(larry_emb,aes(x=X0,y=X1,color=descrip))+geom_point() + theme_void() +
  scale_color_manual(values = c(jdb_palette("brewer_spectra")[c(1,8)],
                                jdb_palette("brewer_spectra")[6]))+
  theme(legend.position="none")
dev.off()
pdf(file = 'larry_pseudotime.pdf',width = 6,height = 5)
ggplot(larry_emb,aes(x=X0,y=X1,color=pseudotime))+geom_point() + theme_void()+
  scale_color_gradientn(colors = jdb_palette("Zissou"))
dev.off()
###########################
### repro
###########################
### pheatmap
larry_cluster <- read.csv("F:/CJP_cellfate/chongbiancheng_zuizhong.csv",
                          header=FALSE)
colnames(larry_cluster) = paste0('Cluster',1:9)
rownames(larry_cluster) = paste0('Cluster',1:9)
pdf(file = 'repro_cluster.pdf',width = 6,height = 5)
pheatmap(larry_cluster,cluster_rows = F,cluster_cols = F,border_color = NA,
         color = colorRampPalette(c("#FAFAFC","#3F96B7"))(50))
dev.off()
### acc col plot

df_acc = data.frame('accuracy'=c(0.88,0.64,0.57,0.31,0.52,0.52),
                    'methods'=c('KAILIN',
                                'cospar_all_clone',
                                'cospar_state_info',
                                'Waddington-OT',
                                'K-means',
                                'hclust'))
pdf(file = 'repro_acc.pdf',width = 6,height = 5)
ggplot(df_acc,aes(y=accuracy,x=reorder(methods,accuracy),fill=methods))+geom_col(width = 0.9)+
  theme_classic()+scale_y_continuous(expand = c(0,0))+
  theme(text = element_text(size = 15)) +
  scale_fill_manual(values = jdb_palette("brewer_spectra")[c(1:4,6:8)])+xlab('')+
  theme(axis.text.x = element_text(size=12,vjust = 0.5, hjust = 0.5,
                                   angle = 35))
dev.off()
### embedding
repro_emb <- read.table("F:/CJP_cellfate/result/cell_fate_pred/repro_embed.csv",
                      header=T,row.names = 1)

pdf(file = 'repro_gt.pdf',width = 6,height = 5)
ggplot(repro_emb,aes(x=X0,y=X1,color=label))+geom_point() + theme_void() +
  scale_color_manual(values = c(jdb_palette("brewer_spectra")[3],
                                'grey',jdb_palette("brewer_spectra")[7]))+
  theme(legend.position="none")
dev.off()
pdf(file = 'repro_pred.pdf',width = 6,height = 5)
ggplot(repro_emb,aes(x=X0,y=X1,color=pred))+geom_point() + theme_void() +
  scale_color_manual(values = c(jdb_palette("brewer_spectra")[3],
                                'grey',jdb_palette("brewer_spectra")[7]))+
  theme(legend.position="none")
dev.off()
pdf(file = 'repro_descrip.pdf',width = 6,height = 5)
ggplot(repro_emb,aes(x=X0,y=X1,color=descrip))+geom_point() + theme_void() +
  scale_color_manual(values = c(jdb_palette("brewer_spectra")[3],
                                jdb_palette("brewer_spectra")[6]
                                ,jdb_palette("brewer_spectra")[7]))+
  theme(legend.position="none")
dev.off()
pdf(file = 'repro_pseudotime.pdf',width = 6,height = 5)
ggplot(repro_emb,aes(x=X0,y=X1,color=pseudotime))+geom_point() + theme_void()+
  scale_color_gradientn(colors = jdb_palette("Zissou"))
dev.off()
