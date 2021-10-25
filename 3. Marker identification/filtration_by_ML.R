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

#dataset extraction
train = dmr_mat_total[,grepl("Train",colnames(dmr_mat_total))]
train = t(train)
train = data.frame(train,Type = ifelse(grepl("Normal",rownames(train)),"Normal","Early"))
test = dmr_mat_total[,grepl("Test_Normal|Test_Early",colnames(dmr_mat_total))]
test = t(test)
test = data.frame(test,Type = ifelse(grepl("Normal",rownames(test)),"Normal","Early"))
late = dmr_mat_total[,grepl("Test_Late",colnames(dmr_mat_total))]
late = t(late)
late = data.frame(late,Type = "Late")

iter.n <- 53
tree.n <- 500
fold.n <- 10
lambda.n <- 500
min1se <-"lambda.1se"

#randomforest filtration
iter = iter.n
train.tmp = train

acc_sen_spe_f1 <-c()
acc_sen_spe_f1.test <-c()

sed <- seed.n
for(i in 1:iter)
{
  set.seed(sed)
  rf = randomForest(Type~.,data = train.tmp,importance=T,ntree = tree.n)
  train_metric <- caret::confusionMatrix(data = rf$predicted,reference = train.tmp$Type)
  test_metric <- caret::confusionMatrix(data = predict(rf,test),reference = test$Type)
  acc_sen_spe_f1 <- rbind(acc_sen_spe_f1,c(train_metric$overall["Accuracy"],
                                           train_metric$byClass["Sensitivity"],
                                           train_metric$byClass["Specificity"],
                                           train_metric$byClass["F1"]))
  acc_sen_spe_f1.test <- rbind(acc_sen_spe_f1.test,c(test_metric$overall["Accuracy"],
                                                     test_metric$byClass["Sensitivity"],
                                                     test_metric$byClass["Specificity"],
                                                     test_metric$byClass["F1"]))
  
  #print(rf$confusion)
  #print(table(predict(rf,test),test$Type))
  cutoff = sort(rf$importance[,4])[1] #remove one marker each time
  region_out = rownames(rf$importance)[rf$importance[,4]<=cutoff]
  train.tmp = train.tmp[,-match(region_out,colnames(train.tmp))]
  print(dim(train.tmp))
}
