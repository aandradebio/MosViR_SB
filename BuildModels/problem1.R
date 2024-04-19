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

n_sample <- 1500  #nrow(class2)
n_sample2 <- 2 * n_sample

### init the variables

overall_rf <- data.frame(matrix(NA, nrow = 7, ncol = 100))
byClass_rf <- data.frame(matrix(NA, nrow = 11, ncol = 100))  # Update the dimensions to match byClass_tmp
AUC <- data.frame(matrix(NA, nrow = 4, ncol = 100))

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
    
    mn$class <- as.factor(mn$class)
    train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
    train <- mn[train_index, ]
    test <- mn[-train_index, ]
    train$class <- as.factor(train$class)

    rf <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
    save (rf, file=paste0("./",args[4],"/model_rf_",i,".rda"))
    
    p <- predict(rf, test, preProcess = c("center", "scale"),type = "prob")
    AUC_tmp <- auc(multiclass.roc(test$class, p))
    AUC <- rbind(AUC, AUC_tmp)
    print(paste0(AUC_tmp))

    p <- predict(rf, newdata = test)
    cf <- confusionMatrix(p, as.factor(test$class))
    byClass_rf_tmp <- cf$byClass
    overall_rf_tmp <- cf$overall
    print(paste0(overall_rf_tmp))
    print(paste0(byClass_rf_tmp))
    byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)
    overall_rf <- cbind(overall_rf, overall_rf_tmp)
    
    print (paste0("Class1. Iterat:", i, " input: ", input.csv1))
}

#print ("All done. Input: ", input.csv1)

write.csv(AUC, file=paste0("./",args[4],"/AUC.txt"), row.names=F)
write.csv(overall_rf, file=paste0("./",args[4],"/Accuracy.txt"), row.names=F)
write.csv(byClass_rf, file=paste0("./",args[4],"/byClass.txt"), row.names=F)
