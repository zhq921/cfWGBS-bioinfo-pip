# cfWGBS-bioinfo-pip
A comprehensive and rigid computational framework to identify reliable cfDNA methylation biomarkers from cfWGBS data

> This repository provides the code for the paper "Whole-genome circulating tumor DNA methylomes of minimal plasma for sensitive detection and molecular classification of cancer".

<img src="https://github.com/zhq921/cfWGBS-bioinfo-pip/blob/master/imgs/computational_workflow.png" width = "70%" />

##### 1) Identification of cfDNA recurrent regions in populations of two cohorts of normal and breast cancer samples. 
To reduce the impact of missing values of cfDNA fragment in the cohorts of samples, the cfDNA recurrent regions with population of samples were identified based on Poisson test in normal cfDNA samples and breast cancer cfDNA samples respectively. The recurrence ratio of each site was calculated using the percentage of samples, which covered by at least 1 read at the site. High confidence cfDNA recurrent regions were selected by the stringent threshold with p value < 0.01 and the recurrence ratio at each site >70%. The overlapping regions between recurrent regions of breast cancer cfDNA samples and those of normal cfDNA samples were extracted as reference recurrent regions for further analysis.

##### 2) Identification of de novo differentially methylated regions calling from cfDNA reference recurrent regions. 
To discover the reliable and precise cfDNA differentially methylated regions (cfDMRs) between case and control in training set, we performed de novo cfDMRs calling in the cfDNA reference recurrent regions that could directly reduce the effect of missing values owing to the increasing number of samples. We scanned to identify de novo cfDMRs within cfDNA reference recurrent regions based on the adjacent CpGs methylated patterns changes by our previous developed method [CpG_MPs][1] with the rigid threshold of absolute mean methylation difference of each region > 0.2 and p value < 0.01 (two sided t-test).

##### 3) Identification of optimal breast cancer-specific cfDNA methylation markers. 
To reduce the effect of methylation noise from other tissues in plasma, we used WGBS samples from the primary tumor tissue to extract the breast cancer-specific cfDMRs. The consistency of methylation pattern changes of cfDMRs in cfDNA and tissue were assessed by mean methylation difference. The cfDMRs were remained as cancer-specific cfDNA methylation markers by the consistent absolute mean methylation difference (> 0.2).
Furthermore, the backward stepwise strategy was implemented to find the optimal cfDMRs as cancer-specific markers. The method started with all cancer-specific cfDMRs ranked by the importance score, and iteratively thrown out the least important feature by one-at-a-time. The importance scores were evaluated using random forest algorithm. Finally, the optimal 15 cancer-specific cfDMRs were remained to construct predictive model.

##### 4) Model construction and validation. 
Least Absolute Shrinkage and Selection Operator (LASSO)-penalized logistic regression was used to construct model. To obtain the most regularized model, the largest lambda (known as “lambda.1se”) at which the error was within one standard error of the minimum was adopted using 10-fold cross-validation. The ultimate model derived from training set was applied into test set for independent validation.

[1]:https://academic.oup.com/nar/article/41/1/e4/1167588
