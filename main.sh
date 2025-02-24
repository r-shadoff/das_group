# ----- Load Libraries
library(tidyverse)
library(dplyr)
library(stringr)
library(ggplot2)

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
# ----- Plot 1
# ----- Plot 2

# ----- Table

# ----- Plot 3

print(length(high_variance_genes))
