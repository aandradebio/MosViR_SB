#!/usr/bin/bash

dataset_path=$1
base_path="/labinfo/home/aandrade/data/datasets/featureMatrices/"
cd $dataset_path

for file in *fasta; do

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
