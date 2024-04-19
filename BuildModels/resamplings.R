# Install and load Bioconductor's Biostrings package
# if not installed already
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
library(Biostrings)
library(doParallel)
cl <- makePSOCKcluster(50)
registerDoParallel(cl)
set.seed(12345)

# Function to perform random resampling
perform_resampling <- function(input_fasta, output_folder, num_resamples, num_sequences) {
  # Read the original FASTA file
  sequences <- readDNAStringSet(input_fasta)

  # Create the output folder if it doesn't exist
  if (!dir.exists(output_folder)) dir.create(output_folder)
  prefix <- gsub("\\.fasta$", "", basename(input_fasta))

  # Perform resampling and save to separate files
  for (i in 1:num_resamples) {
    # Randomly sample a fixed number of sequences
    sampled_indices <- sample(length(sequences), size = num_sequences)
    sampled_sequences <- sequences[sampled_indices]

    # Create the output file name
    output_file <- file.path(output_folder, paste0(prefix,"_", i, ".fasta"))

    # Write the sampled sequences to a new FASTA file
    writeXStringSet(sampled_sequences, output_file, format = "fasta")
  }
}

# Retrieve command line arguments
args <- commandArgs(trailingOnly = TRUE)

# Check if there are enough arguments
if (length(args) != 4) {
  cat("Usage: Rscript script.R input_fasta output_folder num_resamples num_sequences\n")
  quit(status = 1)
}

# Extract command line arguments
input_fasta_file <- args[1]
output_folder <- args[2]
num_resamples <- as.integer(args[3])
num_sequences <- as.integer(args[4])

# Perform random resampling
perform_resampling(input_fasta_file, output_folder, num_resamples, num_sequences)

