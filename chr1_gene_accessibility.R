# Library Packages
library(dplyr)
library(data.table)
library(tidyr)
library(ggplot2)
library(tidyverse)
library(bigreadr)
library(GenomicRanges)
library(IRanges)
library(plotly)
library(ggpubr)
library(gtools)

# Datasets
fpath <- "/home/bdsi2024/Genomics/Data/scMultiome/Trevino_dataset/"
fpath2 <- "/home/bdsi2024/Genomics/Data/"
fpath3 <- "/home/bdsi2024/Genomics/Data/scMultiome/Trevino_dataset/scATAC_pcw16/"
fpath4 <- "/home/bdsi2024/Genomics/Data/scMultiome/Trevino_dataset/scATAC_pcw20/"
fpath5 <- "/home/bdsi2024/Genomics/Data/scMultiome/Trevino_dataset/scATAC_pcw21/"
fpath6 <- "/home/bdsi2024/Genomics/Data/scMultiome/Trevino_dataset/scATAC_pcw24/"
atac_metadata <- fread(paste0(fpath,"/GSE162170_atac_cell_metadata.txt"))
atac_cp <- fread(paste0(fpath,"/GSE162170_atac_consensus_peaks.bed.gz"))
atac_sqc <- fread(paste0(fpath,"/GSE162170_atac_sample_qc.tsv.gz"))
atac_counts_chr1_pcw16 <- fread2(paste0(fpath3,"chr1_scATAC_pcw16.tsv.gz"))
atac_counts_chr1_pcw20 <- fread2(paste0(fpath4,"chr1_scATAC_pcw20.tsv.gz"))
atac_counts_chr1_pcw21 <- fread2(paste0(fpath5,"chr1_scATAC_pcw21.tsv.gz"))
atac_counts_chr1_pcw24 <- fread2(paste0(fpath6,"chr1_scATAC_pcw24.tsv.gz"))
multiome_cell_metadata <- fread(paste0(fpath,"/GSE162170_multiome_cell_metadata.txt.gz"))
multiome_cluster_names <- fread(paste0(fpath,"/GSE162170_multiome_cluster_names.txt.gz"))
multiome_consensus_peaks <- fread(paste0(fpath,"/GSE162170_multiome_atac_consensus_peaks.txt.gz"))
gene_coordinates <- fread(paste0(fpath2,"gene_coordinates_GRCh38.p14.txt"))

# Filter for chromosome 1
chr1_peaks_atac <- atac_cp %>%
  filter(V1 == "chr1")
setnames(chr1_peaks_atac, old = c("V1", "V2", "V3"), new = c("Chromosome", "PeakStart", "PeakEnd"))

multiome_chr1_peaks <- multiome_consensus_peaks %>%
  filter(seqnames == "chr1")

gene_coordinates_chr1 <- gene_coordinates %>%
  filter(chromosome_name == 1)

# Preliminary Gene Accessibility Score

# Filtering for peaks using GenomicRanges package; same as what is done in loop
# Set basepair ranges for genes
gene_coordinates_chr1 <- gene_coordinates_chr1 %>%
  mutate(start = start_position - 100000,
         end = end_position + 100000,
         seq_names = paste0('chr', chromosome_name)) %>%
  makeGRangesFromDataFrame(start.field = "start",
                           end.field = "end",
                           seqnames.field = "seq_names",
                           keep.extra.columns = TRUE)

#chr1_peaks_atac <- cbind(chr1_peaks_atac, w16_counts_binary, w20_counts_binary, w21_counts_binary, w24_counts_binary)
# Set basepair ranges for peaks
chr1_peaks <- chr1_peaks_atac %>%
  select(Chromosome, PeakStart, PeakEnd) %>%
  makeGRangesFromDataFrame(start.field = "PeakStart",
                           end.field = "PeakEnd",
                           seqnames.field = "Chromosome",
                           keep.extra.columns = TRUE)

# Find peaks within gene range
x = findOverlaps(gene_coordinates_chr1, chr1_peaks, select="all")
x

# Display genes and the peaks that are associated with them
gr1.matched = gene_coordinates_chr1[queryHits(x)]
mcols(gr1.matched) = cbind.data.frame(mcols(gr1.matched),
                                      mcols(chr1_peaks[subjectHits(x)]))
gr1.matched = as.data.frame(gr1.matched)

# Indices of peaks; aligns with atac_counts_chr1_pcw16 file's number of observations
gene_matched <- gr1.matched %>%
  mutate(index = subjectHits(x))

# WEEK 16
# Sample unique values for splitting (replace this with actual data)
atac_metadata_w16 <- atac_metadata %>%
  filter(Age == "pcw16")

unique_values <- unique(atac_metadata_w16$Iterative.LSI.Clusters)

# Split the indices based on unique values
indices_list <- split(seq_along(unique_values), unique_values)

# Initialize an empty matrix to store the results
result_matrix <- matrix(NA, nrow = nrow(atac_counts_chr1_pcw16), ncol = length(indices_list))

# Calculate the row means
for (i in seq_along(indices_list)) {
  idx <- indices_list[[i]]
  result_matrix[, i] <- rowMeans(atac_counts_chr1_pcw16[, idx, drop = FALSE])
}

# Convert to dataframe and set column names for clarity
w16_pac <- as.data.frame(result_matrix)
indices_df <- as.data.frame(indices_list)
colnames(w16_pac) <- (colnames(indices_df))

# Initialize an empty dataframe to store the results
gene_accessibility_scores_w16 <- data.frame()

# Define the index variable (assuming you have it)
index <- 1:nrow(w16_pac) # Modify this according to your actual index variable

# Group by gene id and calculate the column means for each gene
gene_accessibility_scores_w16 <- gene_matched %>%
  group_by(ensembl_gene_id) %>%
  summarize(
    start = min(index),
    end = max(index)
  ) %>%
  rowwise() %>%
  mutate(
    gac = list(colMeans(w16_pac[start:end, ]))
  ) %>%
  select(ensembl_gene_id, gac) %>%
  unnest(cols = c(gac))

# Add a column for the repetition index within each group of 20
gene_accessibility_scores_w16 <- gene_accessibility_scores_w16 %>%
  group_by(ensembl_gene_id) %>%
  mutate(index = row_number()) %>%
  ungroup()

# Use pivot_wider to transform the dataframe
gac_wide_w16 <- gene_accessibility_scores_w16 %>%
  pivot_wider(
    names_from = ensembl_gene_id,
    values_from = gac
  ) %>%
  select(-index)

# Gives gene accessibility scores in gene by cluster format (3144x20)
gac_wide_w16 <- t(gac_wide_w16)
colnames(gac_wide_w16) <- (colnames(indices_df))
gac_wide_w16 <- as.data.frame(gac_wide_w16)
gac_wide_w16 <- gac_wide_w16 %>%
  mutate(timepoint = "Week 16")

# END OF WEEK 16

# BEGIN WEEK 20

# Sample unique values for splitting (replace this with actual data)
atac_metadata_w20 <- atac_metadata %>%
  filter(Age == "pcw20")

unique_values <- unique(atac_metadata_w20$Iterative.LSI.Clusters)

# Split the indices based on unique values
indices_list_ <- split(seq_along(unique_values), unique_values)

# Initialize an empty matrix to store the results
result_matrix <- matrix(NA, nrow = nrow(atac_counts_chr1_pcw20), ncol = length(indices_list))

# Calculate the row means
for (i in seq_along(indices_list)) {
  idx <- indices_list[[i]]
  result_matrix[, i] <- rowMeans(atac_counts_chr1_pcw20[, idx, drop = FALSE])
}

# Convert to dataframe and set column names for clarity
w20_pac <- as.data.frame(result_matrix)
indices_df <- as.data.frame(indices_list)
colnames(w20_pac) <- (colnames(indices_df))

# Initialize an empty dataframe to store the results
gene_accessibility_scores_w20 <- data.frame()

# Define the index variable (assuming you have it)
index <- 1:nrow(w20_pac) # Modify this according to your actual index variable

# Group by gene id and calculate the column means for each gene
gene_accessibility_scores_w20 <- gene_matched %>%
  group_by(ensembl_gene_id) %>%
  summarize(
    start = min(index),
    end = max(index)
  ) %>%
  rowwise() %>%
  mutate(
    gac = list(colMeans(w20_pac[start:end, ]))
  ) %>%
  select(ensembl_gene_id, gac) %>%
  unnest(cols = c(gac))

# Add a column for the repetition index within each group of 20
gene_accessibility_scores_w20 <- gene_accessibility_scores_w20 %>%
  group_by(ensembl_gene_id) %>%
  mutate(index = row_number()) %>%
  ungroup()

# Use pivot_wider to transform the dataframe
gac_wide_w20 <- gene_accessibility_scores_w20 %>%
  pivot_wider(
    names_from = ensembl_gene_id,
    values_from = gac
  ) %>%
  select(-index)

# Gives gene accessibility scores in gene by cluster format (3144x20)
gac_wide_w20 <- t(gac_wide_w20)
colnames(gac_wide_w20) <- (colnames(indices_df))
gac_wide_w20 <- as.data.frame(gac_wide_w20)
gac_wide_w20 <- gac_wide_w20 %>%
  mutate(timepoint = "Week 20")

# END OF WEEK 20

# BEGIN WEEK 21

# Sample unique values for splitting (replace this with actual data)
atac_metadata_w21 <- atac_metadata %>%
  filter(Age == "pcw21")

unique_values <- unique(atac_metadata_w21$Iterative.LSI.Clusters)

# Split the indices based on unique values
indices_list_ <- split(seq_along(unique_values), unique_values)

# Initialize an empty matrix to store the results
result_matrix <- matrix(NA, nrow = nrow(atac_counts_chr1_pcw21), ncol = length(indices_list))

# Calculate the row means
for (i in seq_along(indices_list)) {
  idx <- indices_list[[i]]
  result_matrix[, i] <- rowMeans(atac_counts_chr1_pcw21[, idx, drop = FALSE])
}

# Convert to dataframe and set column names for clarity
w21_pac <- as.data.frame(result_matrix)
indices_df <- as.data.frame(indices_list)
colnames(w21_pac) <- (colnames(indices_df))

# Initialize an empty dataframe to store the results
gene_accessibility_scores_w21 <- data.frame()

# Define the index variable (assuming you have it)
index <- 1:nrow(w21_pac) # Modify this according to your actual index variable

# Group by gene id and calculate the column means for each gene
gene_accessibility_scores_w21 <- gene_matched %>%
  group_by(ensembl_gene_id) %>%
  summarize(
    start = min(index),
    end = max(index)
  ) %>%
  rowwise() %>%
  mutate(
    gac = list(colMeans(w21_pac[start:end, ]))
  ) %>%
  select(ensembl_gene_id, gac) %>%
  unnest(cols = c(gac))

# Add a column for the repetition index within each group of 20
gene_accessibility_scores_w21 <- gene_accessibility_scores_w21 %>%
  group_by(ensembl_gene_id) %>%
  mutate(index = row_number()) %>%
  ungroup()

# Use pivot_wider to transform the dataframe
gac_wide_w21 <- gene_accessibility_scores_w21 %>%
  pivot_wider(
    names_from = ensembl_gene_id,
    values_from = gac
  ) %>%
  select(-index)

# Gives gene accessibility scores in gene by cluster format (3144x20)
gac_wide_w21 <- t(gac_wide_w21)
colnames(gac_wide_w21) <- (colnames(indices_df))
gac_wide_w21 <- as.data.frame(gac_wide_w21)
gac_wide_w21 <- gac_wide_w21 %>%
  mutate(timepoint = "Week 21")

# END OF WEEK 21

# BEGIN WEEK 24

# Sample unique values for splitting (replace this with actual data)
atac_metadata_w24 <- atac_metadata %>%
  filter(Age == "pcw24")

unique_values <- unique(atac_metadata_w24$Iterative.LSI.Clusters)

# Split the indices based on unique values
indices_list_ <- split(seq_along(unique_values), unique_values)

# Initialize an empty matrix to store the results
result_matrix <- matrix(NA, nrow = nrow(atac_counts_chr1_pcw24), ncol = length(indices_list))

# Calculate the row means
for (i in seq_along(indices_list)) {
  idx <- indices_list[[i]]
  result_matrix[, i] <- rowMeans(atac_counts_chr1_pcw24[, idx, drop = FALSE])
}

# Convert to dataframe and set column names for clarity
w24_pac <- as.data.frame(result_matrix)
indices_df <- as.data.frame(indices_list)
colnames(w24_pac) <- (colnames(indices_df))

# Initialize an empty dataframe to store the results
gene_accessibility_scores_w24 <- data.frame()

# Define the index variable (assuming you have it)
index <- 1:nrow(w24_pac) # Modify this according to your actual index variable

# Group by gene id and calculate the column means for each gene
gene_accessibility_scores_w24 <- gene_matched %>%
  group_by(ensembl_gene_id) %>%
  summarize(
    start = min(index),
    end = max(index)
  ) %>%
  rowwise() %>%
  mutate(
    gac = list(colMeans(w24_pac[start:end, ]))
  ) %>%
  select(ensembl_gene_id, gac) %>%
  unnest(cols = c(gac))

# Add a column for the repetition index within each group of 20
gene_accessibility_scores_w24 <- gene_accessibility_scores_w24 %>%
  group_by(ensembl_gene_id) %>%
  mutate(index = row_number()) %>%
  ungroup()

# Use pivot_wider to transform the dataframe
gac_wide_w24 <- gene_accessibility_scores_w24 %>%
  pivot_wider(
    names_from = ensembl_gene_id,
    values_from = gac
  ) %>%
  select(-index)

# Gives gene accessibility scores in gene by cluster format (3144x20)
gac_wide_w24 <- t(gac_wide_w24)
colnames(gac_wide_w24) <- (colnames(indices_df))
gac_wide_w24 <- as.data.frame(gac_wide_w24)
gac_wide_w24 <- gac_wide_w24 %>%
  mutate(timepoint = "Week 24")

# END OF WEEK 24

# Combine all four time points
all_gacs_chr1_wide <- rbindlist(list(gac_wide_w16, gac_wide_w20,
                                gac_wide_w21, gac_wide_w24), idcol = "ID")
all_gacs_chr1_wide$ID = rep(rownames(gac_wide_w16), 4)

all_gacs_chr1_long <- all_gacs_chr1_wide %>%
  pivot_longer(cols = starts_with("c"),
               names_to = "Cluster_Type",
               values_to = "gac")

# Violin Plot  
all_gacs_chr1_long %>%
  group_by(ID) %>%
  ggplot(aes(x = timepoint, y = gac, fill = timepoint)) + 
  geom_violin(trim = FALSE,
              linewidth = 0) + 
  #geom_boxplot(width = 0.1,
  #             outlier.size = 0.5) +
  theme_minimal() +
  labs(title = "Gene Accessibility Scores by Cluster Type in Chromosome 1",
       x = "Timepoint",
       y = "Gene Accessibility Score") +
  scale_fill_discrete(name = "Timepoint") +
  facet_wrap(mixedsort(vars(Cluster_Type)),
             nrow = 4,
             ncol = 5) +
  theme(axis.text.x = element_text(size = 7))


# Filter out genes that exhibit consistent zero accessibility across all four timepoints
# Identify IDs to be filtered out
ids_to_remove <- all_gacs_chr1_long %>%
  group_by(ID) %>%
  summarise(all_zero = all(gac == 0)) %>%
  filter(all_zero) %>%
  pull(ID)

# Filter out those IDs
nonzero_gacs_chr1_long <- all_gacs_chr1_long %>%
  filter(!ID %in% ids_to_remove)
length(unique(nonzero_gacs_chr1_long$ID)) # drops from 3144 to 3126

# Violin Plot  
nonzero_gacs_chr1_long %>%
  group_by(ID) %>%
  ggplot(aes(x = timepoint, y = gac, fill = timepoint)) + 
  geom_violin(trim = FALSE,
              linewidth = 0) + 
  #geom_boxplot(width = 0.1,
  #             outlier.size = 0.5) +
  theme_minimal() +
  labs(title = "Gene Accessibility Scores by Cluster Type in Chromosome 1",
       x = "Timepoint",
       y = "Gene Accessibility Score") +
  scale_fill_discrete(name = "Timepoint") +
  facet_wrap(mixedsort(vars(Cluster_Type)),
             nrow = 4,
             ncol = 5) +
  theme(axis.text.x = element_text(size = 7))
# This shows virtually no difference than all genes but still worth exploring