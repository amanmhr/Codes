library(DESeq2)

# Example count matrix (replace this with your own data)
counts <- matrix(c(100, 150, 120, 90, 30, 50,
                   80, 110, 130, 100, 40, 60),
                 nrow = 2,
                 ncol = 6,
                 byrow = TRUE)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = c("Sample1", "Sample2", "Sample3", "Sample4", "Sample5", "Sample6"),
  Group = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")
)

# Set row names of the count matrix to be the gene names
rownames(counts) <- c("GeneA", "GeneB")



# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = counts,
                              colData = sample_info,
                              design = ~ Group)



# Estimate size factors and dispersions
dds <- DESeq(dds)

# Test for differential expression
results <- results(dds)

# Display the top significant differentially expressed genes
top_genes <- head(results[order(results$padj), ], 10)
print(top_genes)




# Assuming you already have your DESeqDataSet object 'dds', which you created earlier

dds <- nbinomLRT(dds, reduced = ~ 1)

# Estimate size factors and gene-wise dispersions
dds <- estimateSizeFactors(dds)
dds <- estimateDispersionsGeneEst(dds)

# Set the gene-wise estimates as final dispersion estimates
dispersions(dds) <- mcols(dds)$dispGeneEst

# Continue with differential expression testing using nbinomWaldTest or nbinomLRT
dds <- nbinomWaldTest(dds)
# OR
dds <- nbinomLRT(dds)

# Further downstream analysis and visualization
# For example, getting the results table
results_table <- as.data.frame(results(dds))
print(results_table)







# Assuming you have already performed the differential expression analysis and have the results table 'results_table'

# Load required packages
install.packages("ggplot2")
library(ggplot2)

# Create the volcano plot
ggplot(results_table, aes(x = log2FoldChange, y = -log10(padj))) +
  geom_point(color = ifelse(results_table$padj < 0.05, "red", "black")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  labs(x = "Log2 Fold Change", y = "-log10(Adjusted p-value)") +
  theme_minimal()











# Load the required packages
library(DESeq2)

# Replace 'metabolomics_data.csv' with the actual file path of your CSV data file
file_path <- "C:/desktop/Yashwant_sir/final_data_ner_multitissue_duplicate removed.csv"
# Read the CSV data into a data frame
metabolomics_data <- read.csv(file_path, header = TRUE, row.names = 1)

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = colnames(metabolomics_data),
  Group = rep(1:4, each = 11)
)

# Create a DESeqDataSet object
dds <- DESeqDataSetFromMatrix(countData = metabolomics_data,
                              colData = sample_info,
                              design = ~ Group)

# Perform differential expression analysis
dds <- DESeq(dds)

# Test for differential expression using likelihood ratio test (LRT)
dds <- nbinomLRT(dds, reduced = ~ 1)

# Get the results table
results_table <- as.data.frame(results(dds))
print(results_table)











# Load the MetaboAnalystR package
install.packages("MetaboAnalystR")
library(MetaboAnalystR)

# Replace 'metabolomics_data.csv' with the actual file path of your CSV data file
file_path <- "C:/desktop/Yashwant_sir/final_data_ner_multitissue_duplicate removed.csv"

# Read the CSV data into a data frame
metabolomics_data <- read.csv(file_path, header = TRUE, row.names = 1)

# Assuming you have already read the metabolomics data into the 'metabolomics_data' data frame

# Create a data frame with sample information
sample_info <- data.frame(
  Sample = colnames(metabolomics_data),
  Group = rep(1:4, each = 11)
)


# Load the data into MetaboAnalystR
inputData <- createMSMSSet(inputMatrix = metabolomics_data, sampleInfo = sample_info, dataScale = TRUE)

# Perform differential expression analysis (e.g., ANOVA)
results <- anovaMS(inputData)

# Display the top significant metabolites
top_metabolites <- head(results, 10)
print(top_metabolites)







# Load the "limma" package
library(limma)

# Example continuous gene expression data
# Replace this with your actual data
data <- matrix(c(10.5, 12.3, 8.7, 9.8, 11.2, 10.1,
                 5.6, 7.2, 6.8, 5.9, 7.5, 8.1,
                 4.3, 3.9, 5.2, 6.1, 4.7, 5.9),
               nrow = 3,
               ncol = 6,
               byrow = TRUE)

# Create a design matrix (replace this with your actual design)
design <- model.matrix(~ Group, data=data.frame(Group = c("Control", "Control", "Control", "Treatment", "Treatment", "Treatment")))

# Fit linear models
fit <- lmFit(data, design)

# Perform empirical Bayes moderated t-statistics
fit <- eBayes(fit)

#Extract differentially expressed genes with FDR < 0.05
DEGs <- topTable(fit, coef = "GroupTreatment", adjust = "fdr", number = Inf)

# View the results
print(DEGs)










install.packages("readxl")
# Load necessary libraries

library(limma)
library(readxl)

# Load the CSV data file
data <- read.csv("C:/desktop/Yashwant_sir/final_data_ner_multitissue_duplicate removed.csv", header = TRUE, row.names = 1)

# Extract sample names from the column names of the data matrix
sample_names <- colnames(data)

# Determine the number of samples
num_samples <- ncol(data)

# Set up the experimental design matrix
# Assuming the first 6 columns are control and the rest are case in each group
group <- rep(1:4, each = 11)
condition <- rep(c(rep("Control", 6), rep("Case", 5)), times = 4)

# Print lengths of group and condition variables, and num_samples
print(paste("Length of group:", length(group)))
print(paste("Length of condition:", length(condition)))
print(paste("Number of samples:", num_samples))

# Verify that the lengths of 'group' and 'condition' match the number of samples
stopifnot(length(group) == num_samples)
stopifnot(length(condition) == num_samples)

# Combine group and condition into the design matrix
design <- model.matrix(~0 + factor(group) + condition)

# Fit linear models
# Fit linear models
# Fit linear models
# Fit linear models
# Fit linear models
fit <- lmFit(data, design)

# Define valid contrast names for each group's condition
valid_contrast_names <- c("Control1_vs_Case1", "Control2_vs_Case2", "Control3_vs_Case3", "Control4_vs_Case4")

# Create a list of contrasts
contrast_list <- list()
for (i in 1:length(valid_contrast_names)) {
  contrast_name <- valid_contrast_names[i]
  control_col <- paste0("conditionControl", i)
  case_col <- paste0("conditionCase", i)
  contrast_list[[contrast_name]] <- c(0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, -1)
}

# Create the contrast matrix
contrast_matrix <- do.call("cbind", contrast_list)

# Perform empirical Bayes moderated t-statistics
fit <- eBayes(fit)

# Extract differentially expressed metabolites with FDR < 0.05 for the first contrast
DEMs <- topTable(fit, coef = 1, adjust = "fdr", number = Inf)

output_file <- "differential_expression_results.csv"
write.csv(DEMs, file = output_file, row.names = TRUE)

# Print a message indicating where the file was saved
cat("Results saved to:", output_file, "\n")








library(ggplot2)

# Load the differential expression results
DEMs <- read.csv("differential_expression_results.csv", header = TRUE, row.names = 1)

# Create a volcano plot
volcano_plot <- ggplot(DEMs, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(abs(logFC) > 1 & P.Value < 0.05, "red", "black")), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(x = "log2 Fold Change", y = "-log10(P-value)", title = "Volcano Plot")

# Display the volcano plot
print(volcano_plot)


output_file <- "volcano_plot.png"
ggsave(output_file, plot = volcano_plot, width = 8, height = 6)

# Print a message indicating where the file was saved
cat("Volcano plot saved to:", output_file, "\n")






#network analysis

library(igraph)
library(ggplot2)

DEMs <- read.csv("differential_expression_results.csv", header = TRUE, row.names = 1)

significant_DEGs <- DEMs$adj.P.Val < 0.05
correlation_matrix <- cor(data[significant_DEGs, ], method = "spearman")

cor_threshold <- 0.8
adjacency_matrix <- ifelse(correlation_matrix > cor_threshold, 1, 0)
network <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")

# Provide gene labels to the graph
gene_labels <- rownames(correlation_matrix)
V(network)$label <- gene_labels

# Layout for better visualization
layout <- layout_with_fr(network)

# Plot the network using the layout
plot(network, layout = layout, 
     vertex.size = 5,            # Adjust vertex size
     vertex.label.cex = 0.7,     # Adjust vertex label size
     vertex.label.dist = 0.8,    # Adjust vertex label distance from vertex
     edge.color = "gray",
     edge.width = 0.5,           # Adjust edge width
     main = "Co-Expression Network")

# Save the plot as an image file
output_file <- "coexpression_network.png"
dev.copy(png, file = output_file, width = 800, height = 800)
dev.off()

# Print a message indicating where the file was saved
cat("Network image saved to:", output_file, "\n")




#compare each group 

library(limma)
library(readxl)

# Load the CSV data file
data <- read.csv("C:/desktop/Yashwant_sir/final_data_ner_multitissue_duplicate removed.csv", header = TRUE, row.names = 1)

# Extract sample names from the column names of the data matrix
sample_names <- colnames(data)

# Determine the number of samples
num_samples <- ncol(data)

# Set up the experimental design matrix
# Assuming the first 6 columns are control and the rest are case in each group
group <- rep(1:4, each = 11)
condition <- rep(c(rep("Control", 6), rep("Case", 5)), times = 4)

# Print lengths of group and condition variables, and num_samples
print(paste("Length of group:", length(group)))
print(paste("Length of condition:", length(condition)))
print(paste("Number of samples:", num_samples))

# Verify that the lengths of 'group' and 'condition' match the number of samples
stopifnot(length(group) == num_samples)
stopifnot(length(condition) == num_samples)

# Combine group and condition into the design matrix
design <- model.matrix(~0 + factor(group) + condition)

# Fit linear models
fit <- lmFit(data, design)

# Perform empirical Bayes moderated t-statistics
fit <- eBayes(fit)

# Define valid contrast names for each group's condition
valid_group_ids <- 1:4
valid_contrast_names <- paste("Group", valid_group_ids, "_Case_vs_Control", sep = "")

# Create a list of contrasts
contrast_list <- list()
for (group_id in valid_group_ids) {
  contrast_name <- valid_contrast_names[group_id]
  control_col <- paste0("conditionControl", group_id)
  case_col <- paste0("conditionCase", group_id)
  contrast_vector <- rep(0, length(sample_names))
  contrast_vector[condition == "Case" & group == group_id] <- 1
  contrast_vector[condition == "Control" & group == group_id] <- -1
  contrast_list[[contrast_name]] <- contrast_vector
}

# Create the contrast matrix
contrast_matrix <- do.call("cbind", contrast_list)

# Extract differentially expressed metabolites with FDR < 0.05 for each contrast
DEMs_list <- list()
for (i in 1:length(valid_contrast_names)) {
  DEMs <- topTable(fit, coef = i, adjust = "fdr", number = Inf)
  DEMs_list[[valid_contrast_names[i]]] <- DEMs
}

# Save differentially expressed metabolites for each contrast
output_dir <- "separate_differential_expression_results"
dir.create(output_dir, showWarnings = FALSE)
for (i in 1:length(valid_contrast_names)) {
  contrast_name <- valid_contrast_names[i]
  output_file <- file.path(output_dir, paste0(contrast_name, "_results.csv"))
  write.csv(DEMs_list[[contrast_name]], file = output_file, row.names = TRUE)
  cat("Results saved to:", output_file, "\n")
}







#gg plot 


library(ggplot2)

# Load the differential expression results
DEMs <- read.csv("Group4_Case_vs_Control_results.csv", header = TRUE, row.names = 1)

# Create a volcano plot
volcano_plot <- ggplot(DEMs, aes(x = logFC, y = -log10(P.Value))) +
  geom_point(aes(color = ifelse(abs(logFC) > 1 & P.Value < 0.05, "red", "black")), size = 2) +
  scale_color_manual(values = c("black", "red")) +
  theme_minimal() +
  labs(x = "log2 Fold Change", y = "-log10(P-value)", title = "Volcano Plot")

# Display the volcano plot
print(volcano_plot)


output_file <- "group_4volcano_plot.png"
ggsave(output_file, plot = volcano_plot, width = 8, height = 6)

# Print a message indicating where the file was saved
cat("Volcano plot saved to:", output_file, "\n")







#network analysis

library(igraph)
library(ggplot2)

DEMs <- read.csv("Group4_Case_vs_Control_results.csv", header = TRUE, row.names = 1)

significant_DEGs <- DEMs$adj.P.Val < 0.05
correlation_matrix <- cor(data[significant_DEGs, ], method = "spearman")

cor_threshold <- 0.8
adjacency_matrix <- ifelse(correlation_matrix > cor_threshold, 1, 0)
network <- graph_from_adjacency_matrix(adjacency_matrix, mode = "undirected")

# Provide gene labels to the graph
gene_labels <- rownames(correlation_matrix)
V(network)$label <- gene_labels

# Layout for better visualization
layout <- layout_with_fr(network)

# Plot the network using the layout
plot(network, layout = layout, 
     vertex.size = 5,            # Adjust vertex size
     vertex.label.cex = 0.7,     # Adjust vertex label size
     vertex.label.dist = 0.8,    # Adjust vertex label distance from vertex
     edge.color = "gray",
     edge.width = 0.5,           # Adjust edge width
     main = "Co-Expression Network")

# Save the plot as an image file
output_file <- "grp_4coexpression_network.png"
dev.copy(png, file = output_file, width = 800, height = 800)
dev.off()

# Print a message indicating where the file was saved
cat("Network image saved to:", output_file, "\n")



#validate results

# Example data for Condition A and Condition B
condition_a <- c(17.56607592,42.15436665,545.2188824,96.22080573,93.42128015,76.04241234)
condition_b <- c(7.785475552,6.772316217,397.4998612,10.73966038,8.231525898)

# Calculate log-fold change
logFC <- log2(condition_a / condition_b)

# Print the log-fold change values
print(logFC)

