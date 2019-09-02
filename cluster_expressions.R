library(limma)
library(Glimma)
library(edgeR)
library(data.table)

setwd("/Volumes/MackDisk/Wong lab- jointed/Samantha - Wound healing project/DE_davide_Alu")
      
data0 <- read.delim("kd2_vs_scr_0h")
data6 <- read.delim("kd2_vs_scr_6h")
data12 <- read.delim("kd2_vs_scr_12h")
data24 <- read.delim("kd2_vs_scr_24h")
data48 <- read.delim("kd2_vs_scr_48h")

target_gene <- rownames(data0[data0[5]<0.05,])

data0_filtered=subset(data0, rownames(data0) %in% target_gene)
data6_filtered=subset(data6, rownames(data6) %in% target_gene)
data12_filtered=subset(data12, rownames(data12) %in% target_gene)
data24_filtered=subset(data24, rownames(data24) %in% target_gene)
data48_filtered=subset(data48, rownames(data48) %in% target_gene)

data_merged_early <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(data0_filtered[2],data6_filtered[2],data12_filtered[2]))

######################### Early
for_clust_early=na.omit(data_merged_early) 
for_clust_with_gene_early=setDT(data.frame(for_clust_early), keep.rownames = TRUE)[]
names(for_clust_with_gene_early) <- c("gene","hr_00", "hr_06", "hr_12")
names(for_clust_early) <- c("hr_00", "hr_06", "hr_12")


for (n_clust in c(1:10)) {
  print(n_clust)
  ### kmeans
  max_itr <-  50
  #n_clust  <-  5  ## number of cluster 
  set.seed(123) ## reproduce the cluster 
  kmeans_out  <- kmeans(for_clust_early,n_clust,iter.max = max_itr)
  
  data_with_cust_info_early <- for_clust_with_gene_early %>% 
    mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))
  
  ## visualise  each cluster 
  temp_plot = data_with_cust_info_early %>% 
    gather(key = "variable" , value = "value", -c(1,5)) %>%  ### 1 is the index of column 'geneName' and 5 is the index of column 'clust'
    group_by(variable) %>%  
    mutate(row_num =  1:n()) %>% 
    ggplot(aes(x =  variable , y = value , group = row_num)) +   
    geom_point() +  
    geom_line(alpha = 1 , aes(col = as.character(clust))) + 
    theme_bw() +  
    theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
    facet_wrap(~clust)
  
  ggsave(temp_plot, file=paste('early_cluster_',n_clust,'.png', sep = ""))
  write.csv(data_with_cust_info_early,file=paste('early_cluster_',n_clust,'.csv', sep = ""))
}

######################### Late
data_merged_late <- Reduce(function(a,b){
  ans <- merge(a,b,by="row.names",all=T)
  row.names(ans) <- ans[,"Row.names"]
  ans[,!names(ans) %in% "Row.names"]
}, list(data0_filtered[2],data24_filtered[2],data48_filtered[2]))

for_clust_late=na.omit(data_merged_late) 
for_clust_with_gene_late=setDT(data.frame(for_clust_late), keep.rownames = TRUE)[]
names(for_clust_with_gene_late) <- c("gene","hr_00", "hr_24", "hr_48")
names(for_clust_late) <- c("hr_00", "hr_24", "hr_48")


for (n_clust in c(1:10)) {
  print(n_clust)
  ### kmeans
  max_itr <-  50
  #n_clust  <-  5  ## number of cluster 
  set.seed(123) ## reproduce the cluster 
  kmeans_out  <- kmeans(for_clust_late,n_clust,iter.max = max_itr)
  
  data_with_cust_info_late <- for_clust_with_gene_late %>% 
    mutate(clust = paste("clust_", kmeans_out$cluster,sep = ""))
  
  ## visualise  each cluster 
  temp_plot = data_with_cust_info_late %>% 
    gather(key = "variable" , value = "value", -c(1,5)) %>%  ### 1 is the index of column 'geneName' and 5 is the index of column 'clust'
    group_by(variable) %>%  
    mutate(row_num =  1:n()) %>% 
    ggplot(aes(x =  variable , y = value , group = row_num)) +   
    geom_point() +  
    geom_line(alpha = 1 , aes(col = as.character(clust))) + 
    theme_bw() +  
    theme(legend.position = "none" , axis.text.x = element_text(angle = 90 , vjust = 0.4)) +
    facet_wrap(~clust)
  
  ggsave(temp_plot, file=paste('late_cluster_',n_clust,'.png', sep = ""))
  write.csv(data_with_cust_info_late,file=paste('late_cluster_',n_clust,'.csv', sep = ""))
}
savehistory("late_early.R")

