#!/usr/bin/bash

# Set the path to the dataset directory provided as the first argument
dataset_path=$1

# Set the base path where the R script mnmer.R is located
base_path="/labinfo/home/aandrade/data/datasets/featureMatrices/"

# Change directory to the dataset path
cd $dataset_path

# Loop through each FASTA file in the dataset directory
for file in *fasta; do

	# Run mnmer.R script with different parameters in parallel for faster processing
	# It seems that mnmer.R takes three arguments: the input FASTA file, k value, and l value
	# You are running mnmer.R with multiple combinations of k and l values
	# '&' at the end of each command runs the commands in parallel
	Rscript "${base_path}mnmer.R" "$file" 1 1 &
	Rscript "${base_path}mnmer.R" "$file" 1 2 &
	Rscript "${base_path}mnmer.R" "$file" 2 0 &
	Rscript "${base_path}mnmer.R" "$file" 1 3 &
	Rscript "${base_path}mnmer.R" "$file" 3 0 &
	Rscript "${base_path}mnmer.R" "$file" 2 1 &
	Rscript "${base_path}mnmer.R" "$file" 2 2 &
	Rscript "${base_path}mnmer.R" "$file" 3 1 &
	Rscript "${base_path}mnmer.R" "$file" 4 0 &

done

