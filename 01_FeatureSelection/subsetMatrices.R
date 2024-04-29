# Load required libraries
library(data.table)  # For efficient data manipulation
library(doParallel)   # For parallel processing
library(caret)       # For data partitioning

# Create a parallel cluster with 50 workers
cl <- makePSOCKcluster(50)
registerDoParallel(cl)

# Set the seed for reproducibility
set.seed(12345)

# Set options to suppress scientific notation
options(scipen = 999)

# Define a function to subset a matrix and write to CSV files
subsetMatrix <- function(matrix_info, output_name) {
  
  # Read the matrix from the input file
  matrix <- fread(matrix_info, header = TRUE, stringsAsFactors = TRUE)
  
  # Create a training set and a test set using stratified sampling
  train_index <- createDataPartition(matrix$class, p = 0.7, list = FALSE)
  train <- matrix[train_index, ]
  test <- matrix[-train_index, ]
  
  # Write the training set and test set to CSV files
  write.csv(train, file = paste0(output_name, "_train.csv"), row.names = FALSE)
  write.csv(test, file = paste0(output_name, "_test.csv"), row.names = FALSE)
}

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

# Call the subsetMatrix function with the specified input and output names
subsetMatrix(args[1], args[2])

