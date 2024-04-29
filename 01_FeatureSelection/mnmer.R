# Load the doParallel package to enable parallel processing
library(doParallel)

# Create a parallel cluster with 10 workers
cl <- makePSOCKcluster(10)
registerDoParallel(cl)

# Set the seed for reproducibility
set.seed(12345)

# Load required libraries
library(Biostrings)
library(mnmer)

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Extract input file path, m, and n values from command line arguments
input_file <- args[1]
m <- as.numeric(args[2])  # Convert m to numeric type
n <- as.numeric(args[3])  # Convert n to numeric type

# Read the DNA sequences from the input FASTA file
fas <- readDNAStringSet(input_file)

# Extract the prefix from the input file name
prefix <- gsub("\\.fasta$", "", basename(input_file))

# Print status message indicating the start of processing for the current file and parameters
print(paste0("Processing ", prefix, "_", m, n))

# Generate the feature matrix using the mnmer function
mntab <- mnmer(fas, m, n)

# Write the feature matrix to a CSV file
write.csv(mntab, file = paste0("featureMatrix_", prefix, "_", m, n, ".csv"), row.names = FALSE, col.names = FALSE)

# Print status message indicating the completion of processing for the current file and parameters
print(paste0("Done: ", prefix, "_", m, n))

# Close the parallel cluster
stopCluster(cl)

