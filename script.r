# Load required libraries
install.packages("DESeq2")  # Uncomment if DESeq2 is not installed
library(DESeq2)

# Step 1: Load your RNA-seq count data
# Assuming you have a counts matrix where rows are genes and columns are samples
# Replace with your actual data path
counts <- read.csv("path/to/counts_matrix.csv", row.names = 1)

# Step 2: Create a sample information data frame
# This should include the condition for each sample
# Replace with your actual sample conditions
coldata <- data.frame(
  row.names = colnames(counts),
  condition = c("healthy", "healthy", "diseased", "diseased")  # Adjust according to your data
)

# Step 3: Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = coldata,
                              design = ~ condition)

# Step 4: Pre-filtering (optional)
# Filter out low count genes
dds <- dds[rowSums(counts(dds)) > 1, ]

# Step 5: Run DESeq2
dds <- DESeq(dds)

# Step 6: Results extraction
res <- results(dds)

# Step 7: Summary of results
summary(res)

# Step 8: Adjust p-values for multiple testing
res <- lfcShrink(dds, coef="condition_diseased_vs_healthy", type="apeglm")

# Step 9: Identify significant genes
# Adjust the threshold as needed
sig_genes <- res[which(res$padj < 0.05), ]
head(sig_genes)

# Step 10: Visualization
# Plot a MA plot
plotMA(res, main="MA Plot", ylim=c(-2, 2))

# Plot a volcano plot
library(EnhancedVolcano)  # Uncomment to install if needed
# install.packages("EnhancedVolcano")
EnhancedVolcano(res,
                lab = rownames(res),
                x = 'log2FoldChange',
                y = 'padj',
                xlim = c(-5, 5),
                ylim = c(0, 10),
                title = "Volcano Plot")

# Step 11: Save results
write.csv(as.data.frame(res), file = "differential_expression_results.csv")
write.csv(as.data.frame(sig_genes), file = "significant_genes.csv")
