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
cl <- makePSOCKcluster(10)
registerDoParallel(cl)
set.seed(12345)
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)

input_csv1 <- args[1]
input_csv2 <- args[2]
output_rdata <- args[3]

AUC_list <- list()
byClass_list <- list()
overall_list <- list()

control <- trainControl(method="cv", number=10, classProbs= TRUE, summaryFunction = multiClassSummary)
metric <- "AUC"

class1 <- fread(input_csv1)
class2 <- fread(input_csv2)
n_sample1 <- nrow(class2)
n_sample2 <- 2 * nrow(class2)

print(paste0("Matrices OK ", output_rdata))

class1.subset <- sample_n(class1, n_sample1)
class2.subset <- sample_n(class2, n_sample1)

mn <- rbind(class1.subset, class2.subset)
mn <- mn[,-1]
mn$class <- as.factor(mn$class)
train_index <- createDataPartition(mn$class, p=0.8, list=FALSE)
train <- mn[train_index, ]
test <- mn[-train_index, ]
train$class <- as.factor(train$class)

knn <- train(class~., data=train, method="kknn", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(knn, test, preProcess = c("center", "scale"), type = "prob")
AUC_knn <- multiclass.roc(test$class, p1)
AUC_knn <- AUC_knn$auc
p1 <- predict(knn, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_knn <- cf$byClass
overall_knn <- cf$overall

# Create dataframes for the metrics
AUC_df_knn <- data.frame(algorithm = "knn", AUC = AUC_knn)
byClass_df_knn <- data.frame(algorithm = "knn", t(byClass_knn))
overall_df_knn <- data.frame(algorithm = "knn", t(overall_knn))

# Add dataframes to the respective lists
AUC_list[[1]] <- AUC_df_knn
byClass_list[[1]] <- byClass_df_knn
overall_list[[1]] <- overall_df_knn

print(paste0("KNN: ", output_rdata))
  
svm <- train(class~., data=train, method="svmRadial", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(svm, test, preProcess = c("center", "scale"), type = "prob")
AUC_svm <- multiclass.roc(test$class, p1)
AUC_svm <- AUC_svm$auc
p1 <- predict(svm, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_svm <- cf$byClass
overall_svm <- cf$overall

# Create dataframes for the metrics
AUC_df_svm <- data.frame(algorithm = "svm", AUC = AUC_svm)
byClass_df_svm <- data.frame(algorithm = "svm", t(byClass_svm))
overall_df_svm <- data.frame(algorithm = "svm", t(overall_svm))

# Add dataframes to the respective lists
AUC_list[[2]] <- AUC_df_svm
byClass_list[[2]] <- byClass_df_svm
overall_list[[2]] <- overall_df_svm

print(paste0("svm: ", output_rdata))
  
nb <- train(class~., data=train, method="nb", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(nb, test, preProcess = c("center", "scale"), type = "prob")
AUC_nb <- multiclass.roc(test$class, p1)
AUC_nb <- AUC_nb$auc
p1 <- predict(nb, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_nb <- cf$byClass
overall_nb <- cf$overall

# Create dataframes for the metrics
AUC_df_nb <- data.frame(algorithm = "nb", AUC = AUC_nb)
byClass_df_nb <- data.frame(algorithm = "nb", t(byClass_nb))
overall_df_nb <- data.frame(algorithm = "nb", t(overall_nb))

# Add dataframes to the respective lists
AUC_list[[3]] <- AUC_df_nb
byClass_list[[3]] <- byClass_df_nb
overall_list[[3]] <- overall_df_nb

print(paste0("nb: ", output_rdata))
  
rf <- train(class~., data=train, method="rf", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(rf, test, preProcess = c("center", "scale"), type = "prob")
AUC_rf <- multiclass.roc(test$class, p1)
AUC_rf <- AUC_rf$auc
p1 <- predict(rf, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_rf <- cf$byClass
overall_rf <- cf$overall

# Create dataframes for the metrics
AUC_df_rf <- data.frame(algorithm = "rf", AUC = AUC_rf)
byClass_df_rf <- data.frame(algorithm = "rf", t(byClass_rf))
overall_df_rf <- data.frame(algorithm = "rf", t(overall_rf))

# Add dataframes to the respective lists
AUC_list[[4]] <- AUC_df_rf
byClass_list[[4]] <- byClass_df_rf
overall_list[[4]] <- overall_df_rf

print(paste0("rf: ", output_rdata))
  
glm <- train(class~., data=train, method="glmnet", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(glm, test, preProcess = c("center", "scale"), type = "prob")
AUC_glm <- multiclass.roc(test$class, p1)
AUC_glm <- AUC_glm$auc
p1 <- predict(glm, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_glm <- cf$byClass
overall_glm <- cf$overall

# Create dataframes for the metrics
AUC_df_glm <- data.frame(algorithm = "glm", AUC = AUC_glm)
byClass_df_glm <- data.frame(algorithm = "glm", t(byClass_glm))
overall_df_glm <- data.frame(algorithm = "glm", t(overall_glm))

# Add dataframes to the respective lists
AUC_list[[5]] <- AUC_df_glm
byClass_list[[5]] <- byClass_df_glm
overall_list[[5]] <- overall_df_glm

print(paste0("glm: ", output_rdata))
  
xg <- train(class~., data=train, method="xgbTree", metric=metric, trControl=control, preProcess = c("center", "scale"))
p1 <- predict(xg, test, preProcess = c("center", "scale"), type = "prob")
AUC_xg <- multiclass.roc(test$class, p1)
AUC_xg <- AUC_xg$auc
p1 <- predict(xg, newdata = test)
cf <- confusionMatrix(p1, as.factor(test$class))
byClass_xg <- cf$byClass
overall_xg <- cf$overall

# Create dataframes for the metrics
AUC_df_xg <- data.frame(algorithm = "xg", AUC = AUC_xg)
byClass_df_xg <- data.frame(algorithm = "xg", t(byClass_xg))
overall_df_xg <- data.frame(algorithm = "xg", t(overall_xg))

# Add dataframes to the respective lists
AUC_list[[6]] <- AUC_df_xg
byClass_list[[6]] <- byClass_df_xg
overall_list[[6]] <- overall_df_xg

print(paste0("xg: ", output_rdata))

AUC_df <- do.call(rbind, AUC_list)
byClass_df <- do.call(rbind, byClass_list)
overall_df <- do.call(rbind, overall_list)

mn_prefix <- sub("\\.RData$", "", output_rdata)  

# Save AUC dataframe to CSV
auc_csv <- paste0(mn_prefix, "_AUC.csv")
write.csv(AUC_df, file = auc_csv, row.names = FALSE)

# Save byClass dataframe to CSV
byclass_csv <- paste0(mn_prefix, "_byClass.csv")
write.csv(byClass_df, file = byclass_csv, row.names = FALSE)

# Save overall dataframe to CSV
overall_csv <- paste0(mn_prefix, "_overall.csv")
write.csv(overall_df, file = overall_csv, row.names = FALSE)

mn_names <- c("knn", "svm", "nb", "rf", "glm", "xg")
mn_names <- paste(mn_prefix, mn_names, sep = "_")  # Combine with prefix

assign(mn_names[1], knn, envir = .GlobalEnv)
assign(mn_names[2], svm, envir = .GlobalEnv)
assign(mn_names[3], nb, envir = .GlobalEnv)
assign(mn_names[4], rf, envir = .GlobalEnv)
assign(mn_names[5], glm, envir = .GlobalEnv)
#assign(mn_names[6], xg, envir = .GlobalEnv)
  
rm (test, train, class1, class2, class3, class3.subset, class1.subset, class2.subset, mn, AUC_df, AUC_svm, AUC_knn, AUC_rf, AUC_nb, p1, p2, p3, p4, p5, p6, cf, knn, svm, rf, glm, xg, nb)
print(paste0("Done", output_rdata))
save.image(output_rdata)
