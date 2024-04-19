library(data.table)
library(doParallel)
library(caret)

cl <- makePSOCKcluster(50)
registerDoParallel(cl)
set.seed(12345)
options(scipen=999)

subsetMatrix <- function(matrix_info, output_name) {
  matrix <- fread(matrix_info, header = TRUE, stringsAsFactors = TRUE)
  
  train_index <- createDataPartition(matrix$class, p = 0.7, list = FALSE)
  train <- matrix[train_index, ]
  test <- matrix[-train_index, ]
  
  write.csv(train, file = paste0(output_name, "_train.csv"), row.names = FALSE)
  write.csv(test, file = paste0(output_name, "_test.csv"), row.names = FALSE)
}

# Get the command-line arguments
args <- commandArgs(trailingOnly = TRUE)

subsetMatrix(args[1], args[2])

