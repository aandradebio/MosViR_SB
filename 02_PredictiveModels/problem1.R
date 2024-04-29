# Load libraries
library(doParallel)   # For parallel processing
library(caret)        # For machine learning modeling
library(dplyr)        # For data manipulation
library(data.table)   # For efficient data manipulation

# Create a parallel cluster with 10 workers
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

# Set the seed for reproducibility
set.seed(12345)

# Set options to suppress scientific notation
options(scipen = 999)

# Get command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Define train control parameters for cross-validation
control <- trainControl(method = "cv", number = 10, classProbs = TRUE, summaryFunction = twoClassSummary)
metric <- "ROC"

# Define input file paths
csv1 <- args[1]  # Path to the first CSV file containing mosquito-associated viruses data
csv2 <- args[2]  # Path to the second CSV file containing arboviruses data
csv3 <- args[3]  # Path to the third CSV file containing other viruses data

# Extract other arguments
num <- as.numeric(args[4])   # Number of iterations
algo <- args[5]               # Algorithm used by caret (e.g., "rf")
out <- args[6]                # Output filename prefix

# Read and preprocess data for each class
class1 <- fread(csv1)
class1$class <- "mosquito.associated"

class2 <- fread(csv2)
class2$class <- "mosquito.associated"
class2 <- class2[, -1]

class3 <- fread(csv3)
class3$class <- "otherviruses"
class3 <- class3[, -1]

n_sample <- nrow(class1)

# Initialize list to store models
mdlist <- list()

# Loop through iterations
for (i in 1:num) {
  # Sample from each class
  class1.subset <- sample_n(class1, n_sample)
  class2.subset <- sample_n(class2, n_sample)
  class3.subset <- sample_n(class3, 2 * n_sample)
  
  # Combine subsets into a single data frame
  mn <- rbind(class1.subset, class2.subset, class3.subset)
  
  # Convert class variable to factor
  mn$class <- as.factor(mn$class)
  
  # Train model
  fit <- train(class ~ ., data = mn, method = algo, metric = metric, trControl = control, preProcess = c("center", "scale"))
  
  # Store trained model in the list
  mdlist[[paste0("model_", algo, "_", i)]] <- fit
  
  # Print iteration details
  ROC_tmp <- max(fit$results$ROC)
  print(paste0("Class1. Iterat: ", i, " Algorithm: ", algo, " ROC: ", ROC_tmp))
}

# Save models to file
save(mdlist, file = paste0("./models/models_", out, "_", algo, ".rda"))

# Print completion message
print("OK")

