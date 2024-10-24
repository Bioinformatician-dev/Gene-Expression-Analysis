# Load necessary libraries
library(DESeq2)

# Load the count data and sample information
count_data <- read.csv("gene_counts.csv", row.names = 1)
sample_info <- read.csv("sample_info.csv")

# Create a DESeq2 dataset object
dds <- DESeqDataSetFromMatrix(countData = count_data,
                              colData = sample_info,
                              design = ~ condition)

# Pre-filtering: Remove rows with low counts
dds <- dds[rowSums(counts(dds)) > 1, ]

# Run the DESeq2 analysis
dds <- DESeq(dds)

# Extract results
results <- results(dds)

# View the results
head(results)

# Save results to a CSV file
write.csv(as.data.frame(results), file = "differential_expression_results.csv")
