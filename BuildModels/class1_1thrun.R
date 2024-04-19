library(doParallel)
library(caret)
library(randomForest)
library(xgboost)
library(glmnet)
library(kknn)
library(e1071)
library(Biostrings)
library(mnmer)
library(MLmetrics)
library(pROC)
library(dplyr)
library(data.table)
cl <- makePSOCKcluster(50)
registerDoParallel(cl)
set.seed(12345)
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)

fragment <- args[1]
output_rdata <- args[2]

AUC_list <- list()
byClass_list <- list()
overall_list <- list()

control <- trainControl(method="cv", number=10, classProbs= TRUE, summaryFunction = multiClassSummary)
metric <- "AUC"

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_31.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_31.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_31.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf31 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf31, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_31 <- data.frame(mer = "31", AUC = AUC_rf)
byClass_31 <- data.frame(mer = "31", t(byClass_rf))
overall_31 <- data.frame(mer = "31", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[1]] <- AUC_31
byClass_list[[1]] <- byClass_31
overall_list[[1]] <- overall_31

print(paste0("31-mer: ", output_rdata))

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_21.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_21.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_21.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf21 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf21, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_21 <- data.frame(mer = "21", AUC = AUC_rf)
byClass_21 <- data.frame(mer = "21", t(byClass_rf))
overall_21 <- data.frame(mer = "21", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[2]] <- AUC_21
byClass_list[[2]] <- byClass_21
overall_list[[2]] <- overall_21

print(paste0("21-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_12.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_12.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_12.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf12 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf12, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_12 <- data.frame(mer = "12", AUC = AUC_rf)
byClass_12 <- data.frame(mer = "12", t(byClass_rf))
overall_12 <- data.frame(mer = "12", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[3]] <- AUC_12
byClass_list[[3]] <- byClass_12
overall_list[[3]] <- overall_12

print(paste0("12-mer: ", output_rdata))

###

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_13.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_13.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_13.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf13 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf13, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_13 <- data.frame(mer = "13", AUC = AUC_rf)
byClass_13 <- data.frame(mer = "13", t(byClass_rf))
overall_13 <- data.frame(mer = "13", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[4]] <- AUC_13
byClass_list[[4]] <- byClass_13
overall_list[[4]] <- overall_13

print(paste0("13-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_14.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_14.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_14.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf14 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf14, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_14 <- data.frame(mer = "14", AUC = AUC_rf)
byClass_14 <- data.frame(mer = "14", t(byClass_rf))
overall_14 <- data.frame(mer = "14", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[5]] <- AUC_14
byClass_list[[5]] <- byClass_14
overall_list[[5]] <- overall_14

print(paste0("14-mer: ", output_rdata))

####
  
  input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_20.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_20.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_20.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf20 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf20, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_20 <- data.frame(mer = "20", AUC = AUC_rf)
byClass_20 <- data.frame(mer = "20", t(byClass_rf))
overall_20 <- data.frame(mer = "20", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[6]] <- AUC_20
byClass_list[[6]] <- byClass_20
overall_list[[6]] <- overall_20

print(paste0("20-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_22.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_22.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_22.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf22 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf22, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_22 <- data.frame(mer = "22", AUC = AUC_rf)
byClass_22 <- data.frame(mer = "22", t(byClass_rf))
overall_22 <- data.frame(mer = "22", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[7]] <- AUC_22
byClass_list[[7]] <- byClass_22
overall_list[[7]] <- overall_22

print(paste0("22-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_23.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_23.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_23.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf23 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf23, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_23 <- data.frame(mer = "23", AUC = AUC_rf)
byClass_23 <- data.frame(mer = "23", t(byClass_rf))
overall_23 <- data.frame(mer = "23", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[8]] <- AUC_23
byClass_list[[8]] <- byClass_23
overall_list[[8]] <- overall_23

print(paste0("23-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_32.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_32.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_32.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf32 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf32, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_32 <- data.frame(mer = "32", AUC = AUC_rf)
byClass_32 <- data.frame(mer = "32", t(byClass_rf))
overall_32 <- data.frame(mer = "32", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[9]] <- AUC_32
byClass_list[[9]] <- byClass_32
overall_list[[9]] <- overall_32

print(paste0("32-mer: ", output_rdata))

###

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_41.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_41.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_41.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf41 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf41, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_41 <- data.frame(mer = "41", AUC = AUC_rf)
byClass_41 <- data.frame(mer = "41", t(byClass_rf))
overall_41 <- data.frame(mer = "41", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[10]] <- AUC_41
byClass_list[[10]] <- byClass_41
overall_list[[10]] <- overall_41

print(paste0("41-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_40.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_40.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_40.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf40 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf40, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_40 <- data.frame(mer = "40", AUC = AUC_rf)
byClass_40 <- data.frame(mer = "40", t(byClass_rf))
overall_40 <- data.frame(mer = "40", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[11]] <- AUC_40
byClass_list[[11]] <- byClass_40
overall_list[[11]] <- overall_40

print(paste0("40-mer: ", output_rdata))

###

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_50.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_50.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_50.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf50 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf50, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_50 <- data.frame(mer = "50", AUC = AUC_rf)
byClass_50 <- data.frame(mer = "50", t(byClass_rf))
overall_50 <- data.frame(mer = "50", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[12]] <- AUC_50
byClass_list[[12]] <- byClass_50
overall_list[[12]] <- overall_50

print(paste0("50-mer: ", output_rdata))

####

input.csv1 <- paste0("featureMatrix_arboviruses_", fragment, "_30.csv.gz")
input.csv2 <- paste0("featureMatrix_mosquito_", fragment, "_30.csv.gz")
input.csv3 <- paste0("featureMatrix_otherviruses_", fragment, "_30.csv.gz")

class1 <- fread(input.csv1)
class2 <- fread(input.csv2)
class3 <- fread(input.csv3)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)
class3.subset <- sample_n(class3, n_sample2)

mn <- rbind(class1.subset, class2.subset)
mn <- rbind(mn, class3.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

rf30 <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf30, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_30 <- data.frame(mer = "30", AUC = AUC_rf)
byClass_30 <- data.frame(mer = "30", t(byClass_rf))
overall_30 <- data.frame(mer = "30", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[13]] <- AUC_30
byClass_list[[13]] <- byClass_30
overall_list[[13]] <- overall_30

print(paste0("30-mer: ", output_rdata))
save.image(output_rdata)

#####

mn_prefix <- sub("\\.RData$", "", output_rdata)  

# Save AUC dataframe to CSV
auc.csv.gz <- paste0(mn_prefix, "_AUC.csv")
write.csv.gz(AUC_df, file = auc.csv.gz, row.names = FALSE)

# Save byClass dataframe to CSV
byclass.csv.gz <- paste0(mn_prefix, "_byClass.csv")
write.csv.gz(byClass_df, file = byclass.csv.gz, row.names = FALSE)

# Save overall dataframe to CSV
overall.csv <- paste0(mn_prefix, "_overall.csv")
write.csv.gz(overall_df, file = overall.csv, row.names = FALSE)

####
rm (test, train, class1, class2, class3, class3.subset, class1.subset, class2.subset, mn, AUC_df, AUC_svm, AUC_knn, AUC_rf, AUC_nb, p1, p2, p3, p4, p5, p6, cf, knn, svm, rf, glm, xg, nb)
print(paste0("Done", output_rdata))
save.image(output_rdata)

