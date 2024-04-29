## Building predictive models

This repository contains scripts used to construct predictive models for the MosViR pipeline. All predictive models are available at https://zenodo.org/records/10950999. 

1. Evaluate multiple algorithms

We began by selecting algorithms based on the best-performing (m,n)-mers identified during feature selection. The following algorithms were evaluated:

- Random Forest (rf)
- K-Nearest Neighbors (knn)
- Support Vector Machine with Radial Kernel (svmRadial)
- Linear Discriminant Analysis (lda)
- Quadratic Discriminant Analysis (qda)
- Boosted Logistic Regression (LogitBoost)

We selected the (1,2)-mer and the (2,1)-mer based on their performance in previous testing (01_FeatureSelection folder). These particular (m,n)-mers demonstrated superior predictive capabilities compared to other feature extraction methods evaluated during our feature selection process. he scripts `problem1.R` and `problem2.R` execute classification tasks for distinguishing between Mosquito-associated viruses and Other viruses, as well as Arboviruses and Mosquito-specific viruses. These scripts construct predictive models and store the ROC metrics, facilitating model selection in subsequent stages.

```
mn=(12 21)
alg=("rf" "svmRadial" "knn" "LogitBoost" "lda" "qda")

for aa in ${alg[@]}
do
	for m in ${mn[@]}
	do
		Rscript problem1.R featureMatrix_mosquito_500bp_${m}_train.csv.gz featureMatrix_arboviruses_500bp_${m}_train.csv.gz featureMatrix_otherviruses_500bp_${m}_train.csv.gz 100 $aa 500bp_10k_${m} &
		Rscript problem2.R featureMatrix_mosquito_500bp_${m}_train.csv.gz featureMatrix_arboviruses_500bp_${m}_train.csv.gz 100 $aa 500bp_10k_${m} &
	done

	wait
done
```

All models were preserved for each (m,n)-mer. From these, the top 50 models with the highest ROC metrics were singled out for use in the pipeline through a soft voting mechanism. In total, 50 predictive models were selected for each fragment length (500bp, 1000bp, 3000bp, 5000bp, and 10000bp) for both classification steps, resulting in a total of 500 predictive models in our pipeline! The code for the ensemble approach is contained in the package `appmodels.R` function (https://github.com/aandradebio/MosViR/). 

