larry_clone_mt <- read.csv("F:/CJP_cellfate/data/larry_clone_mt.csv",row.names = 1)
larry_metadata <- read.csv("F:/CJP_cellfate/data/larry_metadata.csv",row.names = 1)

clone = c()
cell = c()
ct = c()
time = c()
for (i in 1:ncol(larry_clone_mt)) {
  clone_use = larry_clone_mt[,i]
  cell_use = which(clone_use=='True')
  ct_use = larry_metadata[which(clone_use=='True'),2]
  time_use = larry_metadata[which(clone_use=='True'),1]
  clone = c(clone,rep(i,length(cell_use)))
  cell = c(cell,cell_use)
  ct = c(ct,ct_use)
  time = c(time,time_use)
}
clone_df = data.frame(clone,ct,cell,time)
library(FateMapper)
colnames(clone_df)[1:2] = c('barcodes','cell_type')
fate_ct_heatmap(clone_df,cluster_row = T,show_rownames = F)
list1 = list()
for (i in unique(clone_df$barcodes)) {
  df_use = clone_df[clone_df$barcodes==i,]
  cell_name = df_use[,3]
  cell_type = df_use[,2]
  list1[[as.character(i)]] = c(paste0(cell_name,collapse = ','),paste0(cell_type,collapse = ','))
}
df_final = as.data.frame(t(as.data.frame(list1)))

write.table(df_final,'F:\\CJP_cellfate\\data/larry_barcode_ct.txt',quote = F,sep = '\t')
###superot task
a1 = larry_metadata[larry_metadata$NeuMon_fate_bias!=0.5,]
a1 = a1[a1$state_info %in% c('undiff','Monocyte','Neutrophil'),]
label = c()
for (i in 1:nrow(a1)) {
  if (a1[i,2]=='Neutrophil') {
    label = c(label,'Neutrophil')
  }else if(a1[i,2]=='Monocyte'){
    label = c(label,'Monocyte')
  }else{
    if (a1[i,3]>0.5) {
      label = c(label,'Neutrophil_prog')
    }else{
      label = c(label,'Monocyte_prog')
    }
  }
}





