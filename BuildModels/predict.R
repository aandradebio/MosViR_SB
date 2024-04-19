library(doParallel)
library(caret)
library(randomForest)
library(MLmetrics)
library(pROC)
library(mnmer)
library(genocut)
library(dplyr)
library(data.table)

cl <- makePSOCKcluster(20)
registerDoParallel(cl)

set.seed(12345)
options(scipen=999)

args <- commandArgs(trailingOnly = TRUE)

control <- trainControl(method="cv", number=10, classProbs= TRUE, summaryFunction = twoClassSummary)
metric <- "ROC"

input.csv1 <- args[1] #feature matriz de test
input.csv2 <- args[2] #feature matriz de test

class1 <- fread(input.csv1)
#lass1 <- subset(class1, select = -class) 
class2 <- fread(input.csv2)
#lass2 <- subset(class2, select = -class) 

test <- rbind(class1,class2)

### init the variables

overall_rf <- data.frame(matrix(NA, nrow = 7, ncol = 100))
byClass_rf <- data.frame(matrix(NA, nrow = 11, ncol = 100))  # Update the dimensions to match byClass_tmp
AUC <- data.frame(matrix(NA, nrow = 4, ncol = 100))

for (i in 1:100){

    load(paste0(args[3],"/model_rf_",i,".rda"))
    p <- predict(rf, test, preProcess = c("center", "scale"),type = "prob")
    AUC_tmp <- auc(multiclass.roc(test$class, p))
    AUC <- rbind(AUC, AUC_tmp)

    p <- predict(rf, newdata = test)
    cf <- confusionMatrix(p, as.factor(test$class))
    byClass_rf_tmp <- cf$byClass
    overall_rf_tmp <- cf$overall
    byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)
    overall_rf <- cbind(overall_rf, overall_rf_tmp)
    
    print (paste0("Class2. Iterat:", i, " input: ", input.csv1))
}

write.csv(AUC, file=paste0(args[3],"/AUC.txt"), row.names=F)
write.csv(overall_rf, file=paste0(args[3],"/Accuracy.txt"), row.names=F)
write.csv(byClass_rf, file=paste0(args[3],"/byClass.txt"), row.names=F)
