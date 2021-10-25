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

########### lasso model #############
x = as.matrix(train.tmp[,-dim(train.tmp)[2]])
y = train.tmp[,dim(train.tmp)[2]]

# tuning the parameters(lambda)
set.seed(sed)
cvfit = cv.glmnet(x = x, y = y,standardize = F,family = "binomial",nfolds = fold.n,nlambda = lambda.n,
                  type.measure = "mse")
plot(cvfit)
coef(cvfit,s = cvfit$lambda.1se)
best_lambda = cvfit[[min1se]]
fit = glmnet(x= x,y =y,standardize = F,family = "binomial",lambda = best_lambda)
coef(fit)

# performance in training set
res = predict(fit, newx=x, type="response", s=best_lambda)
print(table(res[,1]>0.6,grepl("Normal",y)))#change threshold
y.pred <- factor(ifelse(res[,1]>0.6,"Normal","Early"))
a <-confusionMatrix(data = y.pred, reference = y)
normal_early_col <- c(rgb(24,53,103,maxColorValue = 255),
                      rgb(248,186,70,maxColorValue = 255))

# boxplot for probability's distribution
plot(res[,1],pch = 21,bg = c(rep(normal_early_col[1],30),rep(normal_early_col[2],38)),
     col = "black",xlab = "Samples_idx",ylab = "Probability")
abline(h = 0.6)
box.y <- res[,1]
box.x <- factor(c(rep("Normal",30),rep("Early stage",38)),levels = c("Normal","Early stage"))
boxplot(box.y~box.x,col = normal_early_col)
t.test(box.y~box.x)


# roc curve
roc.tmp=roc(response = y, predictor = res[,1],auc=T,ci =T,plot = T)
text(0.2,0.1,paste0("AUC: ",round(roc.tmp$auc,3)))
roc.area(obs = as.numeric(y)-1,pred = res[,1])



#predict test
x.test = as.matrix(test[,colnames(x)])
y.test = test[,dim(test)[2]]
res = predict(fit, newx=x.test, type="response", s=best_lambda)
print(table(res[,1]>0.6,grepl("Normal",y.test)))

#probability distribution
plot(res[,1],pch = 21,bg = c(rep(normal_early_col[1],10),rep(normal_early_col[2],15)),
     col = "black",xlab = "Samples_idx",ylab = "Probability")
box.y <- res[,1]
box.x <- factor(c(rep("Normal",10),rep("Early stage",15)),levels = c("Normal","Early stage"))
boxplot(box.y~box.x,col = normal_early_col)
t.test(box.y~box.x)

# roc curve for test set 1
roc.tmp=roc(response = y.test, predictor = res[,1],auc=T,ci =T,plot = T)
text(0.2,0.1,paste0("AUC: ",round(roc.tmp$auc,3)))
print(roc.tmp$auc)
roc.area(obs = as.numeric(y.test)-1,pred = res[,1])




#probability distribution
prob_distribution <- c()
prob_distribution.test <- c()
for(i in 1:100){
  set.seed(i)
  cvfit = cv.glmnet(x = x, y = y,standardize = F,family = "binomial",nfolds = 10,
                    type.measure = "mse")
  best_lambda = cvfit$lambda.1se
  fit = glmnet(x= x,y =y,standardize = F,family = "binomial",lambda = best_lambda)
  res = predict(fit, newx=x, type="response", s=best_lambda)
  prob_distribution <- rbind(prob_distribution,res[,1])
  res = predict(fit, newx=x.test, type="response", s=best_lambda)
  prob_distribution.test <- rbind(prob_distribution.test,res[,1])
}
boxplot(prob_distribution,col = c(rep(normal_early_col[1],30),rep(normal_early_col[2],38)))
boxplot(prob_distribution.test,col = c(rep(normal_early_col[1],10),rep(normal_early_col[2],15)))


#predict test_normal and late
x.test = rbind(as.matrix(test[grepl("Normal",rownames(test)),colnames(x)]),
               as.matrix(late[,colnames(x)]))
y.test = c(as.character(test[grepl("Normal",rownames(test)),dim(test)[2]]),
           as.character(late[,dim(late)[2]]))

res = predict(fit, newx=x.test, type="response", s=best_lambda)
table(res[,1]>0.6,grepl("Normal",y.test))

#probability_distribution
plot(res[,1],pch = 21,bg = c(rep(normal_early_col[1],10),rep(normal_early_col[2],70)),
     col = "black",xlab = "Samples_idx",ylab = "Probability")
box.y <- res[,1]
box.x <- factor(c(rep("Normal",10),rep("Late stage",70)),levels = c("Normal","Late stage"))
boxplot(box.y~box.x,col = normal_early_col)
t.test(box.y~box.x)

# roc curve for test set 2
roc.tmp=roc(response = y.test, predictor = res[,1],auc=T,ci =T,plot = T)
text(0.2,0.1,paste0("AUC: ",round(roc.tmp$auc,3)))
roc.area(obs = ifelse(grepl("Normal",y.test),1,0),pred = res[,1])

