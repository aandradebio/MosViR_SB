library(doParallel)
library(caret)
library(randomForest)
library(MLmetrics)
library(pROC)
library(dplyr)
library(data.table)

cl <- makePSOCKcluster(25)
registerDoParallel(cl)

set.seed(12345)
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)

control <- trainControl(method="cv", number=10, classProbs= TRUE, summaryFunction = twoClassSummary)
metric <- "ROC"

input.csv1 <- args[1] #"featureMatrix_arboviruses_500bp_12.csv.gz"
input.csv2 <- args[2] #"featureMatrix_mosquito_500bp_12.csv.gz"
init <- args[4] # first iteration
final <- args[5] # final iteration
#input.csv3 <- args[3] #"featureMatrix_otherviruses_500bp_12.csv.gz"

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)

n_sample <- 500  #nrow(class2)

for (i in init:final){

    print("starting")

    class1.subset <- sample_n(class1, n_sample)
    class2.subset <- sample_n(class2, n_sample)
    
    mn <- rbind(class1.subset, class2.subset)
    mn <- mn[,-1]

    rf <- train(class~., data=mn, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))

    save (rf, file=paste0("./",args[3],"/model_rf_",i,".rda"))
    
    print (paste0("Class2. Iterat:", i, " input: ", input.csv1))
}

print ("All done. Input: ", input.csv1)
