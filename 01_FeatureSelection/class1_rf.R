# Load required libraries
library(doParallel)   # For parallel processing
library(caret)        # For machine learning modeling
library(randomForest) # For random forest modeling
library(MLmetrics)    # For evaluating machine learning models
library(pROC)         # For ROC curve analysis
library(dplyr)        # For data manipulation
library(data.table)   # For efficient data manipulation

# Create a parallel cluster with 20 workers
cl <- makePSOCKcluster(20)
registerDoParallel(cl)

# Set the seed for reproducibility
set.seed(12345)

# Set options to suppress scientific notation
options(scipen = 999)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define train control parameters for cross-validation
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)

# Define evaluation metric
metric <- "ROC"

### init the variables

overall_rf <- data.frame(matrix(NA, nrow = 7, ncol = 100))
byClass_rf <- data.frame(matrix(NA, nrow = 11, ncol = 100))  # Update the dimensions to match byClass_tmp
AUC <- data.frame(matrix(NA, nrow = 4, ncol = 100))


# Get input file paths from command-line arguments
input.csv1 <- args[1] #"featureMatrix_arboviruses_500bp_12_train.csv.gz"
input.csv2 <- args[2] #"featureMatrix_mosquito_500bp_12_train.csv.gz"
input.csv3 <- args[3] #"featureMatrix_otherviruses_500bp_12_train.csv.gz"
input.csv4 <- args[4] #"featureMatrix_arboviruses_500bp_12_test.csv.gz"
input.csv5 <- args[5] #"featureMatrix_mosquito_500bp_12_test.csv.gz"
input.csv6 <- args[6] #"featureMatrix_otherviruses_500bp_12_test.csv.gz"

# Read and subset data for each class
class1_train <- fread(input.csv1)
class2_train <- fread(input.csv2)
class3_train <- fread(input.csv3)

class1_test <- fread(input.csv4)
class2_test <- fread(input.csv5)
class3_test <- fread(input.csv6)
class1_test$class <- replicate(nrow(class1_test), "mosquito.associated")
class2_test$class <- replicate(nrow(class2_test), "mosquito.associated")
class3_test$class <- replicate(nrow(class3_test), "otherviruses")
test <- rbind(class1_test,class2_test)
test <- rbind(test,class3_test)

# Define number of samples to select from each class
n_sample <- 500  # Number of samples to select from class 2
n_sample2 <- 2 * n_sample

# Loop through iterations
for (i in 1:100){

    # Sample from each class and assign class labels
    class1.subset <- sample_n(class1_train, n_sample)
    class1.subset$class <- replicate(nrow(class1.subset), "mosquito.associated")
    class2.subset <- sample_n(class2_train, n_sample)
    class2.subset$class <- replicate(nrow(class2.subset), "mosquito.associated")
    class3.subset <- sample_n(class3_train, n_sample2)
    class3.subset$class <- replicate(nrow(class3.subset), "otherviruses")
    
    # Combine subsets into a single data frame
    mn <- rbind(class1.subset, class2.subset)
    mn <- rbind(mn, class3.subset)
    mn <- mn[,-1]

    # Train random forest model
    rf <- train(class ~ ., data = mn, method = "rf", metric = metric, trControl = control, preProcess = c("center", "scale"))

    # Save trained model
    save(rf, file = paste0("./", args[7], "/model_rf_", i, ".rda"))
    
    # Predict probabilities for the testing set using the trained random forest model
    p <- predict(rf, test, preProcess = c("center", "scale"), type = "prob")

    # Calculate the Area Under the ROC Curve (AUC) for the multiclass classification
    AUC_tmp <- auc(multiclass.roc(test$class, p))

    # Append the calculated AUC to the existing AUC matrix
    AUC <- rbind(AUC, AUC_tmp)

    # Make predictions on the testing set using the trained random forest model
    p <- predict(rf, newdata = test)

    # Calculate confusion matrix and related metrics
    cf <- confusionMatrix(p, as.factor(test$class))

    # Extract class-specific metrics (e.g., sensitivity, specificity) from the confusion matrix
    byClass_rf_tmp <- cf$byClass

    # Extract overall metrics (e.g., accuracy, kappa) from the confusion matrix
    overall_rf_tmp <- cf$overall

    # Combine the class-specific metrics with the existing metrics matrix
    byClass_rf <- cbind(byClass_rf, byClass_rf_tmp)

    # Combine the overall metrics with the existing metrics matrix
    overall_rf <- cbind(overall_rf, overall_rf_tmp)

    print (paste0("Classification Step 1. Iterat:", i, " input: ", input.csv1))
}

write.csv(AUC, file=paste0(args[7],"/AUC.txt"), row.names=F)
write.csv(overall_rf, file=paste0(args[7],"/Accuracy.txt"), row.names=F)
write.csv(byClass_rf, file=paste0(args[7],"/byClass.txt"), row.names=F)


# Print completion message
print("All done. Input: ", input.csv1)

