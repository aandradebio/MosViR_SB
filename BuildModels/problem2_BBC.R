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

input.csv1 <- args[1] # path for the feature matrix
input.csv2 <- args[2] # path for the feature matrix

n_sample <- 1500  #nrow(class2)

### init the variables

overall_rf <- data.frame(matrix(NA, nrow = 7, ncol = 100))
byClass_rf <- data.frame(matrix(NA, nrow = 11, ncol = 100))  # Update the dimensions to match byClass_tmp
AUC <- data.frame(matrix(NA, nrow = 4, ncol = 100))

for (i in 1:100){

    class1 <- fread(paste0(input.csv1,"arboviruses_500bp_",i,"_BBC_k3.csv"))
 #   class1 <- subset(class1, select = -class) 
    class2 <- fread(paste0(input.csv2,"mosquito_subset_500bp_",i,"_BBC_k3.csv"))
 #   class2 <- subset(class2, select = -class) 

    class1.subset <- sample_n(class1, n_sample)
    class1.subset$class <- replicate(nrow(class1.subset), "arboviruses")
    class2.subset <- sample_n(class2, n_sample)
    class2.subset$class <- replicate(nrow(class2.subset), "mosquito")
    
    mn <- rbind(class1.subset, class2.subset)
    mn <- mn[,-1]
    
    mn$class <- as.factor(mn$class)
    train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
    train <- mn[train_index, ]
    test <- mn[-train_index, ]
    train$class <- as.factor(train$class)

    rf <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
    save (rf, file=paste0("./",args[3],"/model_rf_",i,"BBC_k3.rda"))
    
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
    
    print (paste0("Class2. BBC k3 Iterat:", i, " input: ", input.csv1))
    
    rm(train,test,rf,p,cf)
}

write.csv(AUC, file=paste0("./",args[3],"/AUC_BBC_k3.txt"), row.names=F)
write.csv(overall_rf, file=paste0("./",args[3],"/Accuracy_BBC_k3.txt"), row.names=F)
write.csv(byClass_rf, file=paste0("./",args[3],"/byClass_BBC_k3.txt"), row.names=F)

rm(byClass_rf,overall_rf,AUC)

### init the variables

overall_rf <- data.frame(matrix(NA, nrow = 7, ncol = 100))
byClass_rf <- data.frame(matrix(NA, nrow = 11, ncol = 100))  # Update the dimensions to match byClass_tmp
AUC <- data.frame(matrix(NA, nrow = 4, ncol = 100))

for (i in 1:100){

    class1 <- fread(paste0(input.csv1,"arboviruses_500bp_",i,"_BBC_k2.csv"))
#    class1 <- subset(class1, select = -class) 
    class2 <- fread(paste0(input.csv2,"mosquito_subset_500bp_",i,"_BBC_k2.csv"))
#    class2 <- subset(class2, select = -class) 

    class1.subset <- sample_n(class1, n_sample)
    class1.subset$class <- replicate(nrow(class1.subset), "arboviruses")
    class2.subset <- sample_n(class2, n_sample)
    class2.subset$class <- replicate(nrow(class2.subset), "mosquito")
    
    mn <- rbind(class1.subset, class2.subset)
    mn <- mn[,-1]
    
    mn$class <- as.factor(mn$class)
    train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
    train <- mn[train_index, ]
    test <- mn[-train_index, ]
    train$class <- as.factor(train$class)

    rf <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
    save (rf, file=paste0("./",args[3],"/model_rf_",i,"BBC_k2.rda"))
    
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
    
    print (paste0("Class2. BBC k2 Iterat:", i, " input: ", input.csv1))
    
    rm(train,test,rf,p,cf)
}

write.csv(AUC, file=paste0("./",args[3],"/AUC_BBC_k2.txt"), row.names=F)
write.csv(overall_rf, file=paste0("./",args[3],"/Accuracy_BBC_k2.txt"), row.names=F)
write.csv(byClass_rf, file=paste0("./",args[3],"/byClass_BBC_k2.txt"), row.names=F)

rm(byClass_rf,overall_rf,AUC)

### init the variables

overall_rf <- data.frame(matrix(NA, nrow = 7, ncol = 100))
byClass_rf <- data.frame(matrix(NA, nrow = 11, ncol = 100))  # Update the dimensions to match byClass_tmp
AUC <- data.frame(matrix(NA, nrow = 4, ncol = 100))

for (i in 1:100){

    class1 <- fread(paste0(input.csv1,"arboviruses_500bp_",i,"_BBC_k4.csv"))
  #  class1 <- subset(class1, select = -class) 
    class2 <- fread(paste0(input.csv2,"mosquito_subset_500bp_",i,"_BBC_k4.csv"))
  #  class2 <- subset(class2, select = -class) 

    class1.subset <- sample_n(class1, n_sample)
    class1.subset$class <- replicate(nrow(class1.subset), "arboviruses")
    class2.subset <- sample_n(class2, n_sample)
    class2.subset$class <- replicate(nrow(class2.subset), "mosquito")
    
    mn <- rbind(class1.subset, class2.subset)
    mn <- mn[,-1]
    
    mn$class <- as.factor(mn$class)
    train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
    train <- mn[train_index, ]
    test <- mn[-train_index, ]
    train$class <- as.factor(train$class)

    rf <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
    save (rf, file=paste0("./",args[3],"/model_rf_",i,"BBC_k4.rda"))
    
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
    
    print (paste0("Class2. BBC k4 Iterat:", i, " input: ", input.csv1))
    
    rm(train,test,rf,p,cf)
}

write.csv(AUC, file=paste0("./",args[3],"/AUC_BBC_k4.txt"), row.names=F)
write.csv(overall_rf, file=paste0("./",args[3],"/Accuracy_BBC_k4.txt"), row.names=F)
write.csv(byClass_rf, file=paste0("./",args[3],"/byClass_BBC_k4.txt"), row.names=F)

rm(byClass_rf,overall_rf,AUC)
