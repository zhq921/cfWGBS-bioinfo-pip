library("glmnet")
library("randomForest")
library("pROC")
library(pheatmap)
library(ggplot2)
library(RColorBrewer)
library("caret") #for confusion matrix
library(e1071)
library("verification") # for roc p value
library(scatterplot3d)
library("ggsci")
library("Rtsne")


dmr_mat_total <- read.table("./met_mat/cfDNA_mat",sep = "\t",stringsAsFactors = F,header = T,row.names = 1)
dmr_mat_total <- dmr_mat_total[!grepl("^M",rownames(dmr_mat_total)),] #rm the rows of chrom MT 
dmr_mat_total <- dmr_mat_total[,!grepl("Z17",colnames(dmr_mat_total))] # remove the Z17 sample from the XH

#add chr to be a name
rownames(dmr_mat_total) = paste0("chr",rownames(dmr_mat_total))

{####filter by tissue
  #paste chr start end
  cpg.mat[,2] = as.character(cpg.mat[,2])
  cpg.mat[,3] = as.character(cpg.mat[,3])
  cpg.mat[,1] = apply(as.matrix(cpg.mat[,1:3]),1,function(x) paste(x,collapse = "_"))
  
  #rm these col  
  cpg.mat = cpg.mat[,-(4:6)]
  #rm the last col 
  cpg.mat = cpg.mat[,-dim(cpg.mat)[2]]
  
  #fill the na value
  nm.idx = grepl("N",colnames(cpg.mat))
  early.idx = grepl("T",colnames(cpg.mat))
  for(i in 1:dim(cpg.mat)[1])
  {
    na.idx = is.na(cpg.mat[i,])
    cpg.mat[i,nm.idx&na.idx] = median(unlist(cpg.mat[i,nm.idx]),na.rm = T)
    cpg.mat[i,early.idx&na.idx] = median(unlist(cpg.mat[i,early.idx]),na.rm = T)
  }
  
  #calculate the mean methylation for each region
  region = unique(cpg.mat[,1])
  region.mat = c()
  for(i in 1:length(region)) #collapse the positions into region
  {
    tmp = colMeans(cpg.mat[cpg.mat[,1]==region[i],4:dim(cpg.mat)[2]],na.rm = T)
    region.mat = rbind(region.mat,tmp)
    rownames(region.mat)[i] = region[i]
  }
  
  #Calculate the p value and mean differ
  nm.idx = which(grepl("N",colnames(region.mat))==T)
  early.idx = which(grepl("T",colnames(region.mat))==T)
  region.mat = cbind(region.mat,tissue_pvalue = 0,tissue_early_meth_mean = 0, tissue_normal_meth_mean= 0, tissue_meandiff = 0)
  for(i in 1:dim(region.mat)[1])
  {
    region.mat[i,"tissue_pvalue"] = t.test(x=region.mat[i,nm.idx],y = region.mat[i,early.idx])$p.value
    region.mat[i,"tissue_early_meth_mean"]  =  mean(region.mat[i,early.idx])
    region.mat[i,"tissue_normal_meth_mean"]  =  mean(region.mat[i,nm.idx])
    region.mat[i,"tissue_meandiff"] = mean(region.mat[i,early.idx])-mean(region.mat[i,nm.idx])
  }
  
  #read the training set's p and meandiffer
  train_dmr_differ_info = dmr_mat_total[,1:2]
  filtered_dmr_name = rownames(region.mat)[(abs(region.mat[,"tissue_meandiff"])>0.2)&
                                             (train_dmr_differ_info$meandiff*region.mat[,"tissue_meandiff"]>0)]
  filtered_dmr_name = paste0("chr",filtered_dmr_name)
  dmr_mat_total = dmr_mat_total[filtered_dmr_name,]
  
  ##output the 68 regions' detail info
  out.tmp <- t(apply(dmr_mat_total,1,function(x) 
  {
    mean_normal_in_cfdna = mean(x[grepl("Training_Normal",colnames(dmr_mat_total))])
    mean_early_in_cfdna = mean(x[grepl("Training_Early",colnames(dmr_mat_total))])
    return(c(mean_normal_in_cfdna,mean_early_in_cfdna))
  }))
  out.tmp <-cbind(out.tmp, dmr_mat_total[,1:2])
  colnames(out.tmp)[1:2] <-c("cfdna_normal_meth_mean","cfdna_early_meth_mean")
  rownames(region.mat) <-paste0("chr",rownames(region.mat))
  out.tmp <-cbind(out.tmp,region.mat[filtered_dmr_name,c(10,11,9,12)])
  write.table(out.tmp,"model/DMR68_info.txt",row.names = T,col.names = T,quote = F,sep="\t")
}
