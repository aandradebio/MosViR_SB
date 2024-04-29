#!/usr/bin/bash

# Accept the dataset path as the first argument
dataset_path=$1

# Change directory to the dataset path
cd $dataset_path

# Create a directory for storing feature matrices
mkdir featureMatrices
feature_path="featureMatrices/"

# Loop through each FASTA file in the dataset directory
for file in *fasta; do

    # Extract the filename without extension
    fasta=$(echo "$file" | cut -d '.' -f 1) 

    # Calculate various features for the current FASTA file using Python scripts and output to CSV files
    calc_bbc.py --fasta "${fasta}.fasta" --molecule dna --k 2 --out "${feature_path}${fasta}_BBC_k2.csv"; echo "Done: $fasta BBC K2" &
    calc_fcgr.py --fasta "${fasta}.fasta" --word_size 2 --out "${feature_path}${fasta}_FCGR_k2.csv"; echo "Done: $fasta FCGR K2" &
    calc_word_ffp.py --fasta "${fasta}.fasta" --molecule dna --word_size 2 --out "${feature_path}${fasta}_FFP_k2.csv"; echo "Done: $fasta FFP K2" &
    calc_word_rtd.py --fasta "${fasta}.fasta" --word_size 2 --out "${feature_path}${fasta}_RTD_k2.csv"; echo "Done: $fasta RTD K2" &
    calc_word_cv.py --fasta "${fasta}.fasta" --word_size 3 --out "${feature_path}${fasta}_CV_k3.csv"; echo "Done: $fasta CV K3" &
    calc_bbc.py --fasta "${fasta}.fasta" --molecule dna --k 3 --out "${feature_path}${fasta}_BBC_k3.csv"; echo "Done: $fasta BBC K3" &
    calc_fcgr.py --fasta "${fasta}.fasta" --word_size 3 --out "${feature_path}${fasta}_FCGR_k3.csv"; echo "Done: $fasta FCGR K3" & 
    calc_word_ffp.py --fasta "${fasta}.fasta" --molecule dna --word_size 3 --out "${feature_path}${fasta}_FFP_k3.csv"; echo "Done: $fasta FFP K3" &
    calc_word_rtd.py --fasta "${fasta}.fasta" --word_size 3 --out "${feature_path}${fasta}_RTD_k3.csv"; echo "Done: $fasta RTD K3" &
    calc_word_cv.py --fasta "${fasta}.fasta" --word_size 4 --out "${feature_path}${fasta}_CV_k4.csv"; echo "Done: $fasta CV K4" & 
    calc_bbc.py --fasta "${fasta}.fasta" --molecule dna --k 4 --out "${feature_path}${fasta}_BBC_k4.csv"; echo "Done: $fasta BBC K4" & 
    calc_fcgr.py --fasta "${fasta}.fasta" --word_size 4 --out "${feature_path}${fasta}_FCGR_k4.csv"; echo "Done: $fasta FCGR K4" &
    calc_word_ffp.py --fasta "${fasta}.fasta" --molecule dna --word_size 4 --out "${feature_path}${fasta}_FFP_k4.csv"; echo "Done: $fasta FFP K4" & 
    calc_word_rtd.py --fasta "${fasta}.fasta" --word_size 4 --out "${feature_path}${fasta}_RTD_k4.csv"; echo "Done: $fasta RTD K4" &
    
    # Remove the header and perform some text replacements in the generated CSV files
    sed -i '1d' "${feature_path}${fasta}*.csv"
    sed -i 's/>//g' "${feature_path}${fasta}*.csv"
    sed -i 's/ /,/g' "${feature_path}${fasta}*.csv"
    
done

