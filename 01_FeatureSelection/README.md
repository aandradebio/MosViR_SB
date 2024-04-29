## Feature Selection and Single-model classification

The FASTA files with the sequences used for training and testing our predictive models can be found on Zenodo (https://zenodo.org/records/10950999). These files cover fragment lengths of 500 bp, 1000 bp, 3000 bp, 5000 bp, and 10000 bp for Other viruses, Mosquito-associated viruses, and Arboviruses (https://zenodo.org/records/10975789). Additionally, we have hosted the complete feature matrices in the Zenodo repository (https://zenodo.org/records/10975801).

We used five alignment-free feature extraction methods from the FASTA files: 1) K-mers, 2) (m,n)-mers, 3) Frequency Chaos Game Representation, 4) Feature Frequency Profiles, and 5) Composition vector. More details on this analysis can be found in the manuscript and in the Supplementary Material.

1. **Extract K-mers and (m,n)-mers using the mnmer R package**

```
	./run_mnmer.sh arbo_1000bp &
	./run_mnmer.sh arbo_10000bp &
	./run_mnmer.sh arbo_3000bp &
	./run_mnmer.sh arbo_5000bp &
	./run_mnmer.sh arbo_500bp &
	./run_mnmer.sh mos_1000bp &
	./run_mnmer.sh mos_10000bp &
	./run_mnmer.sh mos_3000bp &
	./run_mnmer.sh mos_5000bp &
	./run_mnmer.sh mos_500bp &
	./run_mnmer.sh other_1000bp &
	./run_mnmer.sh other_10000bp &
	./run_mnmer.sh other_3000bp &
	./run_mnmer.sh other_5000bp &
	./run_mnmer.sh other_500bp &
```

The `run_mnmer.sh` script serves as a wrapper to execute the `mnmer.R` script with different input files and parameter combinations, facilitating parallel processing of feature extraction for various biological sequence datasets. By providing different input files and parameters, the run_mnmer.sh script orchestrates the generation of feature matrices for distinct sequence datasets, enabling efficient analysis and comparison of biological sequences.
   
2. **Alfree software to extract alignment-free features**

```./alfpy.sh arbo_1000bp &
  	./alfpy.sh arbo_10000bp &
  	./alfpy.sh arbo_3000bp &
  	./alfpy.sh arbo_5000bp &
  	./alfpy.sh arbo_500bp &
  	./alfpy.sh mos_1000bp &
  	./alfpy.sh mos_10000bp &
  	./alfpy.sh mos_3000bp &
  	./alfpy.sh mos_5000bp &
  	./alfpy.sh mos_500bp &
  	./alfpy.sh other_1000bp &
  	./alfpy.sh other_10000bp &
  	./alfpy.sh other_3000bp &
  	./alfpy.sh other_5000bp &
  	./alfpy.sh other_500bp &
```

3. After obtaining the feature matrices, we used the caret package in R to subset it into train and testing sets with the `subsetMatrices.R` script. 

```cd /labinfo/home/aandrade/data/datasets/featureMatrices/mosquito
	while read line; do Rscript ../subsetMatrices.R "${line}.csv.gz" ${line};  done < samples.tab &
	cd /labinfo/home/aandrade/data/datasets/featureMatrices/arboviruses
	while read line; do Rscript ../subsetMatrices.R "${line}.csv.gz" ${line} ; done < samples.tab &
	cd /labinfo/home/aandrade/data/datasets/featureMatrices/otherviruses
	while read line; do Rscript ../subsetMatrices.R "${line}.csv.gz" ${line} ; done < samples.tab &
	wait
```

4. The classification is run by the `class1_rf.R` and `class2_rf.R` scripts that automates the following steps:

**Class Assignment**: Assigns classes to the data based on predefined criteria.
   
**Removal of Contiguous IDS**: Filters out contiguous IDS to improve data quality.
   
**Random Resampling (50 times)**: Performs random resampling of the data to generate diverse training and testing datasets.

**Subsetting the Matrices into Train and Test Sets**: Divide the datasets into separate subsets for training and testing purposes.

**Classification using the Random Forest Algorithm**: Utilizes the Random Forest algorithm to classify the biological samples based on the input feature matrices.

**Saving the Models**: Saves the trained Random Forest models for future use and reference.

To run the classification for the 3000bp fragment length and the (2,2)-mer, execute the following command:
   
```
   mkdir class1_3000bp_22
   Rscript class1_rf.R "featureMatrix_arboviruses_3000bp_22_filtered_train.csv.gz" "featureMatrix_mosquito_3000bp_22_train.csv.gz" "featureMatrix_otherviruses_3000bp_22_train.csv.gz" "featureMatrix_arboviruses_3000bp_22_filtered_test.csv.gz" "featureMatrix_mosquito_3000bp_22_test.csv.gz" "featureMatrix_otherviruses_3000bp_22_test.csv.gz" "class1_3000bp_22" &
   mkdir class2_3000bp_22
   Rscript class2_rf.R "featureMatrix_arboviruses_3000bp_22_filtered_train.csv.gz" "featureMatrix_mosquito_3000bp_22_train.csv.gz"  "featureMatrix_arboviruses_3000bp_22_filtered_test.csv.gz" "featureMatrix_mosquito_3000bp_22_test.csv.gz" "class2_3000bp_22" &
```

This sequence of commands creates directories for storing classification results and executes the Random Forest classification scripts (`class1_rf.R` and `class2_rf.R`) with the corresponding input and output file paths. These results were used by the `featureSelection_plots.R` script to generate figures 1 and 2 above (Supplementary Material 1 for the publication Andrade et al., 2024). 

![](https://github.com/aandradebio/MosViR_SB/blob/main/01_FeatureSelection/FeatureSelection_Fig1_SupplementaryMaterial1.png)
**Figure 1**. Boxplot for the AUC values retrieved through 50 random resampling for each feature extraction method. 

![](https://github.com/aandradebio/MosViR_SB/blob/main/01_FeatureSelection/SupplementaryMaterial1_AUCs.png)
**Figure 2**. A comparison of the classification performances of (m,n)-mer and k-mer after 50 random resamplings. A) Performance for step 1 (Mosquito-associated vs Other viruses), and B) Performance for step 2 (Arboviruses vs Mosquito-specific viruses). The figure shows the mean and confidence interval for AUC values. In most cases, the confidence intervals are barely visible on the current scale.
