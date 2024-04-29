## Phylogeny analysis of unknown RdRp proteins 

To assess the capability of our pipeline in predicting novel mosquito-associated viruses, we detected divergent contigs bearing the RNA-dependent RNA polymerase (RdRp) protein from both the Putative novel viruses and the Unknown contigs groups. 

Building on the previous functional annotation, we subsetted the novel RdRp contigs according to their homology with specific RdRp domains: RdRP_1 (PF00680), RdRP_2 (PF00978), RdRP_3 (PF00998), RdRP_4 (PF02123), RdRP_5 (PF07925), Birna_RdRp (PF04197), Flavi_NS5 (PF00972), Mitovir_RNA_pol (PF05919), Bunya_RdRp (PF04196), Arena_RNA_pol (PF06317), Mononeg_RNA_pol (PF00946), and Flu_PB1 (PF00602). We subsetted the novel RdRp contigs into 12 viral groups based on their homology to these RdRp domains. 

All novel RdRp sequences are stored in the `Novel_sequences.zip.`

1. **Retrieve representative sequences for each RdRp domain**

We identified the representative RdRp sequences for each domain based on the RdRp-scan database. The accession numbers were used to retrieve the sequences and their taxonomy from the NCBI database.

1. Prepare Accession Numbers: The file named `accession_numbers.txt` containing the list of accession numbers to process.

2. Generate Script: Run the following command in the terminal to generate the script to retrieve taxonomic information for each accession number:
   
   ```bash
   while read line; do echo "./getTax.sh ${line}" >> run.sh; done < accession_numbers.txt

3. Set Permissions and run

   ```bash
   chmod 777 *sh
   ./run.sh

This script automates the retrieval of taxonomic information for the provided accession numbers. It calls the getTax.sh script for each accession number, passing it as an argument. The getTax.sh script utilizes the Entrez Direct utilities (esearch and efetch) to fetch taxonomic information from NCBI databases.
The getTax.sh script generates a CSV file named taxonomic_info.csv, which contains the following information for each accession number:

    Accession number
    Taxonomic ID (Taxid)
    Lineage
    Scientific Name

The taxonomy information was used to annotate the Itool Trees. All representative sequences are available at `Representative_sequences.zip`

The novel RdRp contigs, together with their corresponding homologs and representative species from each order or family retrieved from the RdRp-Scan sequence database, formed the following groups: Picornavirales-like and Nidovirales-like, Tymovirales-like and Hepe-Virga-like, Tombusviridae-like and Nodaviridae-like, Toti-, Luteo-, and Sobemoviridae-like, Reoviridae-like, Birnaviridae-like, Flaviviridae-like, Narnaviridae-like, Bunyavirales-like, Arenaviridae-like, Mononega- and Chuviridae-like, Orthomyxoviridae-like.

2. **Alignment**

For each RdRp domain, we performed amino acid sequence alignment using MAFFT v.7.407 (https://mafft.cbrc.jp/alignment/software/) and the E-INS-I algorithm. We removed ambiguously aligned regions using TrimAl v. 2.0 (https://github.com/inab/trimal), employing an automated trimming heuristic. 

```
mafft --auto --globalpair --thread 50 FlaviNS5.faa > FlaviNS5.aln
trimal -in FlaviNS5.aln -out FlaviNS5_trim.aln -automated1
```

3. **Tree**

All multi-sequence alignments were manually checked to ensure the presence of at least two conserved motifs in the RdRp domain. Each multi-sequence alignment was subjected to maximum likelihood phylogenetic analysis as implemented inIQ-TREE v. 2 (https://www.iqtree.org/), with substitution model selection was carried out by the ModelFinder algorithm and branch support assessed by 1,000 bootstrap replicates. The resulting phylogenies were annotated using iTol (https://itol.embl.de). All trees are available at `Trees.zip`

```
iqtree2 -s FlaviNS5_trim.aln -B 10000 -m MFP -nt AUTO --seqtype AA -pre FlaviNS5_trim
```

The text files used to annotate the tree are available in the file `itol.txt`

