#!/bin/bash

# Set the parent directory containing the folders with Bowtie log files
parent_directory="/labinfo/home/aandrade/contigs_test/00_Alignment"

# Get a list of folders containing Bowtie log files
folders=($(find "$parent_directory" -type d -name "PRJ*"))

# Loop through each folder
for folder in "${folders[@]}"; do
    # Set the path to Bowtie log files in the current folder
    bowtie_logs="$folder/*_reads.mapping.bowtie.log"

    # Initialize a variable to store the total unaligned reads
    total_unaligned_reads=0

    # Loop through each Bowtie log file in the current folder
    for log_file in $bowtie_logs; do
        # Extract the number of reads aligned concordantly 0 times
        unaligned_reads=$(awk '/aligned concordantly 0 times/ {print $1}' "$log_file")

        # Add the unaligned reads to the total
        total_unaligned_reads=$((total_unaligned_reads + unaligned_reads))
    done

    # Extract the folder name from the path
    folder_name=$(basename "$folder")

    # Print the folder name and total unaligned reads
    echo "$folder_name total unaligned reads: $total_unaligned_reads"
done

