library(rda)
library(data.table)
library(foreign)
library(dplyr)

process_rda_files <- function(label, fragment_length, mer, input_dir) {
  mntab_list <- list()
  files <- list.files(paste0(input_dir), "*.rda", full.names = TRUE)
  
  for (file in files) {
    load(file)
    mntab_list[[file]] <- mntab
  }
  
  combined_df <- do.call(rbind, mntab_list)
  combined_df$class <- label
  cat("Feature Matrix OK")
  output_filename <- paste0("./featureMatrix_", label, "_", fragment_length, "_", mer, ".csv")
  write.csv(combined_df, output_filename, row.names = FALSE)
  cat("Processing complete. CSV saved to:", output_filename, "\n")
} 

args <- commandArgs(trailingOnly = TRUE)

label <- args[1]
fragment_length <- args[2]
mer <- args[3]
input_dir <- args[4]

process_rda_files(label, fragment_length, mer, input_dir)
