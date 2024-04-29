library(Biostrings)
library(mnmer)
library(data.table)

sequences <- readDNAStringSet("knownContigs.fasta")
sequence_lengths <- width(sequences)
min_length <- 500
max_length <- 1000
filtered_sequences <- sequences[sequence_lengths >= min_length]
filtered_sequences <- filtered_sequences[filtered_sequences =< max_length]
writeXStringSet(filtered_sequences, format = "fasta", file = "filtered_sequences.fasta")

mn <- mnmer(filtered_sequences,1,2)
classes <- fread("knownContigs_classes_step1.csv")
merged_df <- merge(mn, classes, by.x = "seqid", by.y = "contigs", all.x = TRUE)
cleaned_df <- na.omit(merged_df)
write.csv(cleaned_df, file="featureMatrices_knownContigs_10000bp_12_step1.csv", row.names=F)

