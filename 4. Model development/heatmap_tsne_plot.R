library(ggplot2)
library(pheatmap)
library(Rtsne)
library(RColorBrewer)

cpg.mat = read.table("./met_mat/cfDNA_mat",header=T,sep="\t",row.names = 1,stringsAsFactors = F,na.strings = "-")
cpg.mat <- cpg.mat[!grepl("^M",rownames(cpg.mat)),] #rm the cpg from chrom MT 

#length distribution
out_DMR = rownames(cpg.mat)
tmp <- unlist(strsplit(out_DMR,split = "_"))
tmp <- as.numeric(tmp[-seq(1,length(tmp),by = 3)])
out_DMR <- matrix(tmp,ncol = 2,byrow = T)

dmr_info <- data.frame(length = out_DMR[,2]-out_DMR[,1])
ggplot(dmr_info, aes(x=length)) + 
  geom_density(alpha=.5, fill="#FF6666")+  
  xlab("Length")+ylab("Density")+
  geom_vline(aes(xintercept=mean(length, na.rm=T)),   # Ignore NA values for mean
             color="red", linetype="dashed", size=1)+
  theme(panel.background=element_rect(fill='transparent', color='black'),
        legend.position = "none",
        panel.grid=element_blank(),
        axis.ticks = element_line(color='black'),
        axis.text=element_text(color='black')
  )+geom_text(aes(x = 250,y = 0.008,label = "Mean = 82.07\n95% CI 77.57-86.56"))


#prepare data
#trainging set
region.mat <- cpg.mat[,grepl("Training",colnames(cpg.mat))]

#pheatmap
annotation = data.frame(Type = factor(ifelse(grepl("Normal",colnames(region.mat)),"Normal","Early")))
rownames(annotation) = colnames(region.mat)
annotation_colors = list(Type = c(Normal = rgb(24,53,103,maxColorValue = 255),
                                  Early = rgb(248,186,70,maxColorValue = 255)))

pheatmap(region.mat,annotation_col = annotation,show_rownames = F,
         show_colnames = F,annotation_colors = annotation_colors,
         color = colorRampPalette(rev(brewer.pal(n = 7, name =
                                                   "RdBu")))(100),
         clustering_method = "ward.D2")


#tsne
set.seed(888)
tsne_out = Rtsne(t(region.mat),perplexity = 5)# Run TSNE
# Show the objects in the 2D tsne representation
pdf("DMR_tsne.pdf",height = 6,width = 6)
plot(tsne_out$Y,col="black",
     xlab = "t-SNE_1",ylab = "t-SNE_2",pch=21,
     bg = ifelse(grepl("Normal",colnames(region.mat)),rgb(24,53,103,maxColorValue = 255),rgb(248,186,70,maxColorValue = 255)),
     cex=1.5)
dev.off()


#test set 1
region.mat <- cpg.mat[,grepl("Test_Normal|Test_Early",colnames(cpg.mat))]
set.seed(888)
tsne_out = Rtsne(t(region.mat),perplexity = 2)# Run TSNE
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col="black",
     xlab = "t-SNE_1",ylab = "t-SNE_2",pch=21,
     bg = ifelse(grepl("Normal",colnames(region.mat)),rgb(24,53,103,maxColorValue = 255),rgb(248,186,70,maxColorValue = 255)),
     cex=1.5)

#test set 2
region.mat <- cpg.mat[,grepl("Test_Normal|Test_Late",colnames(cpg.mat))]
set.seed(888)
tsne_out = Rtsne(t(region.mat),perplexity = 6)# Run TSNE
# Show the objects in the 2D tsne representation

plot(tsne_out$Y,col="black",
     xlab = "t-SNE_1",ylab = "t-SNE_2",pch=21,
     bg = ifelse(grepl("Normal",colnames(region.mat)),rgb(24,53,103,maxColorValue = 255),rgb(248,186,70,maxColorValue = 255)),
     cex=1.5)


#all samples
region.mat <- cpg.mat[,grepl("Training|Test",colnames(cpg.mat))]
col.point <- pal_nejm("default")(3)
col.points <- rep(col.point[3],163)
col.points[grepl("Normal",colnames(region.mat))] <- col.point[1]
col.points[grepl("Early",colnames(region.mat))] <- col.point[2]
set.seed(888)
tsne_out = Rtsne(t(region.mat),perplexity = 2)# Run TSNE
# Show the objects in the 2D tsne representation
plot(tsne_out$Y,col="black",
     xlab = "t-SNE_1",ylab = "t-SNE_2",pch=21,
     bg = col.points,
     cex=1.5)
legend(0,0,legend = c("Normal","Early","Late"),pt.bg = col.point,pch = 21)

