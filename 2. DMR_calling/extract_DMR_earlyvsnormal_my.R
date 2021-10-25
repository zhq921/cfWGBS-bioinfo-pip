date()
setwd("nmvsearly_multicenter/")

#extract sam names
sam.name = as.character(read.table("cpg_in_seg_all_sam",sep="\t",stringsAsFactors = F,nrows = 1))[7:170]
#cli
cli = read.table("../../SMART/data_preparation/sample_info_0724.txt",stringsAsFactors = F,header=T,sep="\t",row.names = 1)
all(sam.name%in%rownames(cli))
cli = cli[sam.name,]

#get sam name and the label
disease_label <- ifelse(cli$case1.control0==0,"Normal","Late")
disease_label[which(cli$advance.early0late1.==0)] <- "Early"
traintest_label <- rep("Test",164)
traintest_label[which((grepl("TZ",sam.name)&cli$advance.early0late1.==0)==T)] <- "Training" #training :TZ early
traintest_label[(grepl("TZ",sam.name)&cli$case1.control0==0)] <-"Training" #training :TZ normal
traintest_label[1:19] <- "Training" #training: normal with no source
sam.name_with_label <- paste(traintest_label,disease_label,sam.name,sep = "_")

#write header for the out put
write.table(matrix(c("region","pvalue","meandiff",sam.name_with_label),nrow = 1),"seg_mat_earlyvsNormal_p0.05_diff0.2",col.names = F,row.names = F,sep="\t",quote = F)

#sam's idx for na filling
train_normal_idx <- grepl("Training_Normal",sam.name_with_label)
train_early_idx <- grepl("Training_Early",sam.name_with_label)
test_normal_idx <- grepl("Test_Normal",sam.name_with_label)
test_early_idx <- grepl("Test_Early",sam.name_with_label)
test_late_idx <- grepl("Test_Late",sam.name_with_label)


con = file("cpg_in_seg_all_sam",open = "r")
line=readLines(con,n=1)
line=readLines(con,n=1) # from the second line
region.cpg.mat = rep(0,164)
region.pos = c("chrom","start","end")
while( length(line) != 0 ) {
    line.ele = unlist(strsplit(line,split = "\t"))
    line.region.name = line.ele[1:3]
    line.value = as.numeric(line.ele[7:170])
    if(all(line.region.name==region.pos))
    {
        region.cpg.mat = rbind(region.cpg.mat,line.value)
    }else
    {
        ###summary the last region
        #colmeans
        if(NCOL(region.cpg.mat)==1)
        {
            out.value = region.cpg.mat
        }else
        {
            out.value = colMeans(region.cpg.mat[,1:164],na.rm = T)
        }
        #fill the NA in region's value
        if((sum(is.na(out.value)&train_normal_idx) < (0.5*sum(train_normal_idx)))&  ##only fill the NA fraction <0.5
           (sum(is.na(out.value)&train_early_idx) < (0.5*sum(train_early_idx)))){
            out.value[is.na(out.value)&train_normal_idx] = median(out.value[train_normal_idx],na.rm = T)
            out.value[is.na(out.value)&train_early_idx] = median(out.value[train_early_idx],na.rm = T)
            pvalue = t.test(out.value[train_early_idx],out.value[train_normal_idx])$p.value
            meandiff = mean(out.value[train_early_idx])-mean(out.value[train_normal_idx])
            if(pvalue<0.05&abs(meandiff)>0.2)
            {
                #fill the test and late data
                out.value[is.na(out.value)&test_normal_idx] = median(out.value[test_normal_idx],na.rm = T)
                out.value[is.na(out.value)&test_early_idx] = median(out.value[test_early_idx],na.rm = T)
                out.value[is.na(out.value)&test_late_idx] = median(out.value[test_late_idx],na.rm = T)
                out.region.name = paste(region.pos,collapse = "_")
                write.table(rbind(c(),c(out.region.name,pvalue,meandiff,out.value)),"seg_mat_earlyvsNormal_p0.05_diff0.2",append = T,col.names = F,row.names = F,sep="\t",quote = F)
            }
        }


        ###the new region
        region.pos = line.region.name
        region.cpg.mat = line.value
    }
    line=readLines(con,n=1)
}
close(con)
date()
