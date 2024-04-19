library(doParallel)
cl <- makePSOCKcluster(10)
registerDoParallel(cl)
set.seed(12345)
library (Biostrings)
library (mnmer)

args <- commandArgs(trailingOnly = TRUE)
input_file <- args[1]
m <- args[2]
n <- args[3]

fas <- readDNAStringSet (input_file)
prefix <- gsub("\\.fasta$", "", basename(input_file))
print(paste0("Doing", prefix, "_", m, n))
mntab <- mnmer (fas, m, n)
#mntab$class <- replicate(nrow(mntab), "mosquito")
write.csv(mntab, file = paste0("featureMatrix_",prefix,"_",m,n,".csv"), row.names = FALSE, col.names = FALSE)
#print(paste0("Done:", prefix, "_", m, n))
#save.image("teste.RData")
