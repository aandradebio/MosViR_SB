library(data.table)
library(pROC)
library(caret)
library(ggplot2)
library(doParallel)

args <- commandArgs(trailingOnly = TRUE)

num_cores <- 80  # Set the number of cores you want to use
cl <- makeCluster(num_cores)
registerDoParallel(cl)

# Define the fragment lengths and corresponding directories
length <- args[1]

# Print current fragment length and step
cat("Running for fragment length:", args[1], "\n")

# Load data for the current fragment length
cat("Loading data...\n")

tab <- fread(paste("featureMatrices_knownContigs_", args[1], "_12_step1.csv", sep = ""), header = TRUE)
tab <- as.data.frame(tab)
tab <- tab[, !names(tab) %in% "seqid", drop = FALSE]
tab$classes <- factor(tab$classes, levels = c("mosquito.associated", "otherviruses"))

cat("Loading model...\n")
load(paste0("/labinfo/apps/labinfo/R/4.3.0_bio/lib64/R/library/MosVir/models/models_class1_",args[1],"_12.rda"))

cat("Predicting...\n")
pma <- predict(models[[1]], tab, preProcess = c("center", "scale"), type = "prob")
  for (j in 2:length(models)) {
    p <- predict(models[[j]], tab, preProcess = c("center", "scale"), type = "prob")
    pma <- pma + p
  }
pma <- pma / length(models)

cat("Computing metrics...\n")
roc_obj <- roc(tab$classes, pma[, "mosquito.associated"])
auc_score <- auc(roc_obj)
predictions <- ifelse(pma[, "mosquito.associated"] > 0.5, "mosquito.associated", "otherviruses")
cf <- confusionMatrix(data = factor(predictions), reference = factor(tab$classes))
byClass_metrics <- as.numeric(cf$byClass)
overall_metrics <- as.numeric(cf$overall)
precision <- byClass_metrics[6]
accuracy <- overall_metrics[1]

rm(models)
save.image(paste0(args[1],"_metrics_Class1_knownContigs.RData"))

results <- data.table(fragment_length = as.numeric(args[1]), AUC = as.numeric(auc_score), Sensitivity = as.numeric(byClass_metrics[1]),
                            Specificity = as.numeric(byClass_metrics[2]), F1 = as.numeric(byClass_metrics[5]),
                            Precision = as.numeric(precision), Accuracy = as.numeric(accuracy))

write.csv(results,"metrics_Class1_knownContigs.csv", row.names=F, append=T)
#m(pma,models,cf,tab,tab_mosquito,tab_otherviruses)

cat("Fragment length:", length, "completed!\n\n")
save.image(paste0(length,"_metrics_Class1_knownContigs.RData"))

stopCluster(cl)
