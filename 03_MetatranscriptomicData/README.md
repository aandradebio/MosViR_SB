## Retrieval and Pre-processing of metatranscriptomic mosquito data

A crucial step in testing the MosViR pipeline involved retrieving real-world data from mosquito metatranscriptomics. Figure 1 presents a comprehensive diagram describing the analysis conducted to prepare the metatranscriptomic data for classification by the MosViR pipeline. 

**Figure 1**. A step-by-step diagram of the comprehensive similarity-based methodology employed for preparing previously published metatranscriptomic data to test the MosViR pipeline. 


**Dependencies**: To run this analysis, you should have a recent instalation of the following: 

- BBDUK (BBmap suite)
- Bowtie2
- CD-HIT
- SPAdes
- HMMER (version 3.4)
- Diamond (version BLASTX for similarity-based search)
- Prodigal
- Prodigal-gv
- pFam database 
- NCBI Non-Redundant database 
- NCBI viral RefSeq database 

We initially pre-selected 24 projects publicly available at the NCBI Sequence Read Archive (SRA). These selections were based on thorough literature surveys.

1. **Retrieval procedure**

Prepare BioProject IDs: We wrote the PRJNA.txt file containing the list of BioProject IDs.

Generate Script: We ran the following command in the terminal to generate the script for downloading SRA data for each BioProject ID:

```
   while read line; do echo "./getSRA.sh ${line}" >> run.sh; done < PRJNA.txt
```

Set Permissions and run the `getSRA.sh` script: 

```
   chmod 777 *sh
   ./run.sh
``` 

This script automates the download of SRA data for the provided BioProject IDs. It calls the getSRA.sh script for each BioProject ID, passing it as an argument. The `getSRA.sh` script utilizes NCBI Entrez Direct utilities (esearch and efetch) to fetch the list of SRR accessions associated with each BioProject ID. It then downloads the data using prefetch and fasterq-dump, storing the fastq files in separate folders named after each BioProject.

2. **Read Quality Control and removal of host reads**

For each set of reads retrieved from the SRA, we labeled them according to their corresponding BioProject IDs. Post-retrieval, we subjected the reads to quality control using the bbduk software. This encompassed adapter trimming, filtering out low-quality reads (with a Phred quality score > 20 and a minimum length of 100), and removing contaminants. 

To exclude host sequences and endogenous viral elements (EVEs), we employed Bowtie2 to align quality reads to mosquito genomes from the genera Culex sp. and Aedes sp. Unmapped reads underwent de novo assembly into contigs using SPAdes v3.16 with default settings. 


```
projects_file="${contigs_path}PRJNA.txt"
projects=$(cat "$projects_file")
index_basename="index_genomes"

# Loop over projects
for project in $projects; do
    echo "Processing project: $project"
    
    cd "$fastq_path$project"

    total_reads=0
    aligned_reads=0

    for a in *1.fastq.gz; do
        echo "Processing file: $a"
        
        bbduk.sh in="${a}" in2=${a%_1.fastq.gz}_2.fastq.gz \
            out="bbduk_${a}" out2="bbduk_${a%_1.fastq.gz}_2.fastq.gz" \
            outs="bbduk_S_${a%_1.fastq.gz}.fastq" stats="bbduk.stats" statscolumns=5 \
            ordered=f threads=48 qtrim=rl minlength=50 -Xmx30g \
            ref=adapter_phix.fa minavgquality=20 ktrim=l mink=8 \
            qout=auto k=31 overwrite=true;
        
        # Align the preprocessed reads
        bowtie2 -x "$contigs_path$index_basename" -1 "bbduk_${a}" -2 "bbduk_${a%_1.fastq.gz}_2.fastq.gz" \
            -D 20 -R 3 -N 0 -L 20 -i S,1,0.50 -S /dev/null -p 48 --un-conc "${alignment_path}unmapped_${a%_1.fastq.gz}.fastq" \
            2> "reads.mapping.bowtie.log"
        
        # Update statistics
       total_reads=$((total_reads + $(grep "reads processed" "bbduk.stats" | cut -f2)))
        aligned_reads=$((aligned_reads + $(grep "aligned" "reads.mapping.bowtie.log" | cut -f6)))
                
        # Perform assembly with SPAdes for each set of unmapped reads
        spades.py -o "$assembly_path$project/${a%_1.fastq.gz}" -t 40 \
            -s "${alignment_path}unmapped_${a%_1.fastq.gz}.fastq" \
            --meta -k 21,33,55,77
        
        echo "Assembly of $a done"
    done

```

`sum_bowtie.sh` was used to collect alignment metrics from the log files. 

Subsequently, we concatenated the contigs for all reads into a single file named: contigs_ID. For each BioProjects ID. The contigs exceeding 500 bp and containing less than 2% non-ACTG base content were clustered to remove redundancy at 90% nucleotide identity using CD-HIT-EST. These processed contigs served as the basis for subsequent steps.

```
cat "$assembly_path$project/contigs_${project}.fasta" | awk '$0 ~ ">" {if (NR > 1) {print c;} c=0;printf substr($0,2,100) "\t"; } $0 !~ ">" {c+=length($0);} END { print c; }' > sizes.csv
cat sizes.csv | awk '{print $NF,$0}'| sort -nr | cut -f2- -d' ' > sort.csv
cat sort.csv | awk -F'\t' '{ if ($2 >= 500 && $2 <=100000) { print $1} }' > names1.tab
filterbyname.sh in="$assembly_path$project/contigs_${project}.fasta" out=${project}.filt.fasta names=names1.tab include=t
cd-hit-est -i ${project}.filt.fasta -o "${project}.nr.fasta" -c 0.999 -n 10 -M 10000 -G 1 -T 20 
```

3. **Similarity-based search**

To identify viral contigs, we conducted a similarity-based search against the NCBI viral RefSeq database, which was downloaded in February 2024. We utilized Diamond BLASTX with an e-value threshold set at 1e-10 for this purpose. Based on the BLASTX outcomes, we categorized our contigs into three distinct groups:

**Known Viruses**: Contigs demonstrating at least 90% amino acid identity to previously described viruses were classified as Known viruses.

**Putative Novel Viruses**: Contigs exhibiting less than 90% identity to known viruses were categorized as Putative novel viruses.

**Unknown Contigs**: Contigs that lacked significant matches to any sequences in the database were labeled as Unknown contigs.

```
diamond blastx -d viral_refseq_01152024_protein.faa.dmnd -q "${project}.nr.fasta" -o "${project}_blastx.m8" --sensitive --threads 20
cat "${project}_blastx.m8" | sort -k1,1 -k11,11n | awk '$11 < 1e-10 && !a[$1]++' > "${project}_blastx_sorted.m8
mkdir 01_noMatch 02_knownViruses 03_novelViruses

cut -f 1 "${project}_blastx_sorted.m8" > "${project}_names.tab"
filterbyname.sh in="${project}.nr.fasta" out="./01_noMatch/${project}_nomatch.fasta" names="${project}_names.tab" include=f

awk '$3 > 90' "${project}_blastx_sorted.m8" | cut -f 1 > "${project}_knownviruses.txt"
filterbyname.sh in="${project}.nr.fasta" out="./02_knownViruses/${project}_known.fasta" names="${project}_knownviruses.txt" include=t
```

The contigs labeled as known viruses were used to assess the generalization of our pipeline with real-world data. We confirmed the viral origin and aminoacid identity values of all contigs by conducting a similarity-based search against the NCBI non-redundant database, downloaded in February 2024, using Diamond BLASTX with an e-value threshold of 1e-10. 

```
diamond blastx -d nr_07022024.dmnd -q "./02_knownViruses/${project}_known.fasta" -o "${project}_blastx_nr.m8" --sensitive --threads 50
cat "${project}_blastx_nr.m8" | sort -k1,1 -k11,11n | awk '$11 < 1e-10 && !a[$1]++' > "${project}_blastx_sorted_nr.m8"

awk '$3 > 90' "${project}_blastx_sorted_nr.m8" | cut -f 1 > "${project}_knownviruses_nr.txt"
filterbyname.sh in=${project}.nr.fasta" out="./02_knownViruses/${project}_known_nr.fasta" names="${project}_knownviruses_nr.txt" include=t

awk '$3 < 90' "${project}_blastx_sorted.m8" | cut -f 1 > "${project}_putativenovelviruses.tab"
filterbyname.sh in="${project}.nr.fasta" out="./03_novelViruses/${project}_novel.fasta" names="${project}_putativenovelviruses.tab" include=t

```

We assigned classes to each contig based on the BLASTX output table, prioritizing hits with lower e-values and higher bit scores for greater statistical significance and sequence similarity. Taxonomic information from these hits was cross-referenced with previously selected viral species to ensure accurate class assignment. These labeled contigs formed feature matrices for testing the MosViR pipeline. 

The script `mnmer.R` demonstrates the feature extraction process for the known contigs, while `computeMetrics_class1.R` and `computeMetrics_class2.R` were employed to evaluate the contigs within the MosViR pipeline. All findings can be accessed in the publication by Andrade et al., 2024. 

3. Functional annotation of metatranscriptomic dark matter

To assess the capability of our pipeline in predicting novel mosquito-associated viruses, we detected divergent contigs bearing the RNA-dependent RNA polymerase (RdRp) protein from both the Putative novel viruses and the Unknown contigs groups.

```
# Navigate to the directory containing the noMatch contigs
cd "./01_noMatch"

# Run Prodigal for gene prediction
prodigal -i "${project}_nomatch.fasta" -o my.genes -a "${project}_nomatch_proteins.faa" -p meta

# Run Prodigal-gv for gene prediction in anonymous mode
prodigal-gv -p meta -i "${project}_nomatch.fasta" -a "${project}_nomatch_proteins_gv.fasta" -o my_genes_gv 

# Concatenate the protein sequences from Prodigal and Prodigal-gv
cat "${project}_nomatch_proteins.faa" "${project}_nomatch_proteins_gv.fasta" > tmp

# Replace the original protein file with the concatenated one
mv tmp "${project}_nomatch_proteins.faa"

# Navigate to the directory containing the novelViruses contigs
cd "./03_novelViruses"

# Run Prodigal for gene prediction
prodigal -i "${project}_novel.fasta" -o my.genes -a "${project}_novel_proteins.faa" -p meta

# Run Prodigal-gv for gene prediction in anonymous mode
prodigal-gv -p meta -i "${project}_novel.fasta" -a "${project}_novel_proteins_gv.fasta" -o my_genes_gv 

# Concatenate the protein sequences from Prodigal and Prodigal-gv
cat "${project}_novel_proteins.faa" "${project}_novel_proteins_gv.fasta" > tmp

# Replace the original protein file with the concatenated one
mv tmp "${project}_novel_proteins.faa"

# Concatenate the protein files from noMatch and novelViruses directories
cat "./01_noMatch/${project}_nomatch_proteins.faa" "./03_novelViruses/${project}_novel_proteins.faa" > "${project}_proteins.faa"

# Cluster the protein sequences using CD-HIT-EST
cd-hit-est -i "${project}_proteins.faa" -o "${project}_proteins_clstr.faa" -c 1 -T 50 -M 200000

# Run HMMER for functional annotation using Pfam
hmmsearch --tblout "${hmmer_path}${project}_pFam_hmmer.out" -E 1e-10 "Pfam-A.hmm" "${project}_proteins_clstr.faa" --cpu 40

# Run Diamond BLASTP against the NCBI nr database
diamond blastp -d nr_07022024.dmnd -q "${project}_proteins_clstr.faa" -o "${project}_blastp_nr.m8" --very-sensitive --threads 40

# Sort and filter the Diamond BLASTP results
cat "${project}_blastp_nr.m8" | sort -k1,1 -k11,11n | awk '$11 < 1e-10 && !a[$1]++' > "${project}_blastp_sorted_nr.m8"

```

The next folder contains the phylogeny analysis. Due to the heavy size of most generated files, they are not included in this repository. However, we are happy to provide them upon request (atrv@lncc.br or aandradebio@gmail.com).

For further details, refer to the publication by Andrade et al., 2024.

