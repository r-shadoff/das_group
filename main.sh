# ----- Load Libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)
library(ggrepel)
library(DESeq2)
library(tidyr)
library(pheatmap)
library(tibble)

# ----- Import Data
customized_read_tsv <- function(file){
  read_tsv(file) %>%
    mutate(fileName = file)}

cancer_df <- list.files(path = "/Users/rachel/Desktop/binf5503/data", full.names = TRUE) %>%
  lapply(customized_read_tsv) %>% # Reads in all tsv files
  reduce(bind_rows) %>% # Combines files into a single dataframe
  select(gene_id, gene_name, read_count, fileName)


# ----- Data Cleaning

# Extracting patient ID and day number from filename string
cancer_df[c('path', 'day')] <- str_split_fixed(cancer_df$fileName, '_', 2)
cancer_df[c('day_number', 'extra')] <- str_split_fixed(cancer_df$day, '_', 2)
cancer_df[c('extra2', 'sample_id')] <- str_split_fixed(cancer_df$path,'\\/(?!.*/)', 2)

# Drop unnecessary columns created during splitting
cancer_df_updated <- subset(cancer_df, select = (-c(extra2, extra, day, path, fileName)))

# Check for NA values
sum(is.na(cancer_df_updated)) # 40 NA values - These correspond to 5 lines of summary metrics included in every file
cancer_df_clean <- na.omit(cancer_df_updated) # Removal of summary metric rows from dataframe

# Split df into D0 and D6 matrices 
cancer_df_clean$day_number <- as.factor(cancer_df_clean$day_number)
cancer_df_test <- split(cancer_df_clean, cancer_df_clean$day_number)
day0 <- as.data.frame(cancer_df_test$D0)
day6 <- as.data.frame(cancer_df_test$D6)

## Exploration of Data

# Number of genes in dataset
sample_list = c("A778", "A899", "A820", "A870")

# Check that records are complete/D0 and D6 have data for the same number of genes
for (x in sample_list) {
    print(dim(subset(cancer_df_clean, sample_id == x & day_number == "D0"))) }

for (x in sample_list) {
  print(dim(subset(cancer_df_clean, sample_id == x & day_number == "D6"))) }

# Genes with 0 reads on D0
for (x in sample_list) {
  print(x)
  print(nrow(subset(cancer_df_clean, sample_id == x & day_number == "D0" & read_count == 0))) }

# Genes with 0 reads on D6
for (x in sample_list) {
  print(x)
  print(nrow(subset(cancer_df_clean, sample_id == x & day_number == "D6" & read_count == 0))) }

# Genes with 0 reads on D0 and D6
for (x in sample_list){
  # Subset for sample x
  sample_data = cancer_df_clean[cancer_df_clean$sample_id == x, ]
  # Filter subset for genes with 0 reads on D0 and D6
  gene_list <- sample_data %>%
    group_by(gene_name) %>%
    filter(any(day_number == "D0" & read_count == 0) & any(day_number == "D6" & read_count == 0)) %>%
    distinct(gene_name) %>%
    pull(gene_name)
  list_name <- paste0("gene_list_",x)
  assign(list_name, gene_list)
}

gene_lists <- list(
  "A780" = gene_list_A870,
  "A778" = gene_list_A778,
  "A820" = gene_list_A820,
  "A899" = gene_list_A899
)

common_0read_genes = Reduce(intersect, gene_lists) # Compares lists and takes only genes that appear in all lists
print(length(common_0read_genes)) # 25,476 genes with 0 reads

cancer_df_filtered <- cancer_df_clean %>%
  filter(!gene_name %in% common_0read_genes) # Filters out data associated with gene names that appear on the common_genes list

dim(cancer_df_filtered) # 253,520 rows remaining - this number is so low because 25,476 genes were removed, but each gene had 8 observations (2 per sample x 4 samples)

# Filtering out low variance genes 
variance_threshold <- 0.99 # Selects 25,035 genes at 0.99 threshold

high_variance_genes <- cancer_df_filtered %>%
  group_by(gene_name) %>%
  summarise(var_read_count = var(read_count)) %>% # calculate read count variance
  filter(var_read_count > variance_threshold) %>% # compares variance to threshold
  pull(gene_name) # pulls gene names that correspond with variances > threshold

# Filter dataframe to only contain genes with high variance
genes_df_filtered <- cancer_df_filtered %>% 
  filter(gene_name %in% high_variance_genes)

genes_df_filtered <- genes_df_filtered %>%
  group_by(gene_name, sample_id, day_number) %>% 
  summarise(read_count=sum(read_count), .groups = "drop")# Ensures that for each gene, sample, and day number there is only one read count value

# Create a count matrix
count_matrix <- genes_df_filtered %>%
  pivot_wider(names_from = c(sample_id, day_number), values_from= read_count) %>%
  column_to_rownames("gene_name")

# Create metadata object 
metadata <- data.frame(
  row.names = colnames(count_matrix), # Format is sample_day
  sample_id = sub("_.*", "", colnames(count_matrix)),  # Extract sample id
  day_number = as.factor(sub(".*_", "", colnames(count_matrix)))  # Extract day 
)

# Differential gene analysis

# create Dseq2 object 
dds <- DESeqDataSetFromMatrix(
  countData = count_matrix,
  colData = metadata,
  design = ~day_number
)

# Run DESeq2
dds <- DESeq(dds)

# Store and sort results
DESeq_results <- results(dds, contrast = c("day_number", "D0", "D6"))
DESeq_results_sorted <- DESeq_results[order(DESeq_results$padj, na.last = TRUE), ] # Creates table with most sigificant genes at the top

# Create Volcano Plot 
DESeq_results_df <- as.data.frame(DESeq_results)
DESeq_results_df <- subset(DESeq_results_df, !is.na(padj)) # Drops p values that are NA

# Add significance column for colour coding plot
DESeq_results_df <- DESeq_results_df %>%  mutate(significance = ifelse(padj < 0.05 & abs(log2FoldChange) > 1, "Significant", "Not Significant"))

# Add expression labels for colour coding
DESeq_results_df <- DESeq_results_df %>%   mutate(expression = case_when(
  log2FoldChange > 1  ~ "Overexpressed",
  log2FoldChange < -1 ~ "Underexpressed",
  TRUE ~ "Neither" # TRUE = when log2fold is between -1 and 1
))

# Take the top 20 genes and label them on the volcano plot 
most_sig_genes <- DESeq_results_df %>%
  filter(padj<0.01 & abs(log2FoldChange) > 3) %>%
  arrange(padj) %>%
  slice_head(n =10) # Takes top 20 most significant genes 

# Plot 1: Volcano Plot

ggplot(DESeq_results_df, aes(x = log2FoldChange, y = -log10(padj), color = padj)) +
  geom_point(alpha = 0.75, size = 2) +
  scale_color_viridis_c(option = "mako", direction = 1, limits = c(0,0.05)) + 
  theme_minimal() +
  labs(
    x = "Log2 Fold Change",
    y = "-Log10 Adjusted p-value",
    colour = "Adjusted p-value"
  ) +
  geom_vline(xintercept = c(-1, 1), linetype = "dashed", color = "black", size =0.3) +  # FC threshold lines
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", size = 0.3) +
  geom_text_repel(data = most_sig_genes, aes(label = rownames(most_sig_genes)), size = 2, colour = 'black')

# Plot 2 : PCA Analysis
# Perform variance stabilizing transformation (VST) on the count data
vsd <- vst(dds, blind = FALSE)

# Get the transformed data for PCA
pca_data <- plotPCA(vsd, intgroup = c("day_number"), returnData = TRUE)

# Create the PCA plot
ggplot(pca_data, aes(PC1, PC2, color = day_number, shape = day_number, label = rownames(pca_data))) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(size = 3, color = "black") +
  scale_color_manual(values = c("D0" = "blue", "D6" = "red")) +  # keeping with the color scheme
  labs(
    title = "PCA of Gene Expression by Day Number", 
    x = "Principal Component 1", 
    y = "Principal Component 2", 
    color = "Day", shape = "Day",
    caption = "Figure 1: PCA plot of gene expression profiles colored by day number (D0 and D6). Points represent individual samples, with day_number used to distinguish between groups."
    ) +
  theme_minimal() +
  theme(
    legend.position = "top",
    plot.title = element_text(hjust = 0.5, face = "bold"),
    legend.text = element_text(size = 12),
    plot.caption = element_text(size = 10, face = "italic", hjust = 0.5)
    )
  # annotate("text", x = -50, y = 100, size = 4, hjust = 0, label = "Figure 1: PCA plot of gene expression profiles colored by day number (D0 and D6). Points represent individual samples, with day_number used to distinguish between groups. The separation along PC1 and PC2 indicates variation in gene expression between the time points.")
  # annotate("text", x = 0, y = -80, size = 4, hjust = 0.5, 
  #           label = "Figure 1: PCA plot of gene expression profiles colored by day number (D0 and D6).\nPoints represent individual samples, with day_number used to distinguish between groups. \nThe separation along PC1 and PC2 indicates variation in gene expression between the time points.")
  # 

# Plot 3: Heatmap 

# aggregate data to get mean read count per gene per day
deg_data <- cancer_df_filtered %>%
  group_by(gene_name, day_number) %>%
  summarise(mean_read_count = mean(read_count), .groups = 'drop') %>%
  spread(key = day_number, value = mean_read_count, fill = 0)

# ensure D0 and D6 exist
if (!all(c("D0", "D6") %in% colnames(deg_data))) {
  stop("Error: D0 or D6 column missing. Check your data.")
}

# calculate log fold change (D6 vs D0)
deg_data <- deg_data %>%
  mutate(logFC = log2(D6 + 1) - log2(D0 + 1)) %>%
  arrange(desc(abs(logFC)))  # Sort by highest absolute logFC

# select top 20 DEGs
top_20_genes <- head(deg_data, 20)

# prepare matrix for heatmap: Aggregate read_count for each gene_name and sample_id
heatmap_data <- cancer_df_filtered %>%
  filter(gene_name %in% top_20_genes$gene_name) %>%
  group_by(gene_name, sample_id) %>%
  summarise(total_read_count = sum(read_count), .groups = 'drop') %>%
  spread(key = sample_id, value = total_read_count, fill = 0) 

# Ccnvert gene names to rownames
heatmap_matrix <- column_to_rownames(heatmap_data, var = "gene_name")

# convert to matrix
heatmap_matrix <- as.matrix(heatmap_matrix)

# log transform counts
heatmap_matrix <- log2(heatmap_matrix + 1)

# generate heatmap with custom blue and brown color palette and axis labels
pheatmap(heatmap_matrix, scale = "row", clustering_distance_rows = "euclidean",
         clustering_distance_cols = "euclidean", clustering_method = "complete",
         color = colorRampPalette(c("brown", "white", "blue"))(50),  # blue and brown color palette, colour-blind friendly 
         main = "Top 20 Differentially Expressed Genes",
         labels_row = top_20_genes$gene_name,  # label y-axis with gene names
         labels_col = colnames(heatmap_matrix), # label x-axis with sample ids
         fontsize_row = 8,  # adjust row label size
         fontsize_col = 8)  # adjust column label size

# Table of GSEA results 
# Filter for upregulated and downregulated genes
upregulated_genes <- DESeq_results_df %>%
  filter(padj < 0.05 & log2FoldChange > 1) %>%
  arrange(padj)

downregulated_genes <- DESeq_results_df %>%
  filter(padj < 0.05 & log2FoldChange < -1) %>%
  arrange(padj)

# Combine upregulated and downregulated genes into one gsea table
gsea_results <- bind_rows(
  upregulated_genes %>%
    mutate(regulation = "Upregulated"),
  downregulated_genes %>%
    mutate(regulation = "Downregulated")
)
# str(gsea_results) # to check column names

# Select relevant columns for the GSEA table
gsea_table <- gsea_results %>%
  rownames_to_column("gene_name") %>%
  select(gene_name, log2FoldChange, padj, regulation)

# View the gsea table
head(gsea_table)
