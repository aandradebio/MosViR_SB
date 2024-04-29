#!/bin/bash

# Check if an accession number is provided as a command line argument
if [ -z "$1" ]; then
    echo "Error: Please provide an accession number as a command line argument."
    exit 1
fi

# Accession number from command line
accession_number="$1"

# Output CSV file
output_file="taxonomic_info2.csv"

# Use esearch to retrieve the taxonomic ID
taxid=$(esearch -db protein -query "$accession_number" | efetch -format docsum | xtract -pattern DocumentSummary -element TaxId)

# Use efetch to retrieve taxonomic information
lineage=$(efetch -db taxonomy -id "$taxid" -format xml | xtract -pattern Taxon -element Lineage)
name=$(efetch -db taxonomy -id "$taxid" -format xml | xtract -pattern Taxon -element ScientificName)

# Output to CSV file
#echo "Accession number,Taxid,Lineage,Scientific_Name" > "$output_file"
echo "$accession_number,$taxid,$lineage,$name" >> "$output_file"

echo "Processed: $accession_number"
echo "Taxonomic information retrieval completed. Output saved to $output_file."

