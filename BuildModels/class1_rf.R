library(doParallel)
library(caret)
library(randomForest)
library(MLmetrics)
library(pROC)
library(dplyr)
library(data.table)

cl <- makePSOCKcluster(20)
registerDoParallel(cl)

set.seed(12345)
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)

control <- trainControl(method="cv", number=10, classProbs= TRUE, summaryFunction = twoClassSummary)
metric <- "ROC"

input.csv1 <- args[1] #"featureMatrix_arboviruses_500bp_12.csv.gz"
input.csv2 <- args[2] #"featureMatrix_mosquito_500bp_12.csv.gz"
input.csv3 <- args[3] #"featureMatrix_otherviruses_500bp_12.csv.gz"
init <- args[5] 
final <- args[6]

class1 <- fread(input.csv1)
class1 <- subset(class1, select = -class) 
class2 <- fread(input.csv2)
class2 <- subset(class2, select = -class) 
class3 <- fread(input.csv3)
class3 <- subset(class3, select = -class) 

n_sample <- 500  #nrow(class2)
n_sample2 <- 2 * n_sample

for (i in init:final){

    class1.subset <- sample_n(class1, n_sample)
    class1.subset$class <- replicate(nrow(class1.subset), "mosquito.associated")
    class2.subset <- sample_n(class2, n_sample)
    class2.subset$class <- replicate(nrow(class2.subset), "mosquito.associated")
    class3.subset <- sample_n(class3, n_sample2)
    class3.subset$class <- replicate(nrow(class3.subset), "otherviruses")
    
    mn <- rbind(class1.subset, class2.subset)
    mn <- rbind(mn, class3.subset)
    mn <- mn[,-1]

    rf <- train(class~., data=mn, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))

    save (rf, file=paste0("./",args[4],"/model_rf_",i,".rda"))
    
   
    print (paste0("Class1. Iterat:", i, " input: ", input.csv1))
}

print ("All done. Input: ", input.csv1)

