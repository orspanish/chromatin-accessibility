library(tidyverse)
library(dtwclust)
library(dplyr)
library(ggplot2)

# Read in combined file
combined_data <- readRDS("allClusterByGene.rds")

# Run this for sensitivity testing dataset; add other ensembl IDs for further 
# sensitivity analyses
chr19_gacs_no_outs <- all_gacs_chr19_wide %>%
  filter(!(ID == "ENSG00000167785"))

# Compute mean gene accessibility scores at each time point for each gene
all_genes_avg_gacs <- combined_data %>%
  rowwise() %>%
  mutate(avg_gac = mean(c_across(1:20), na.rm = TRUE)) %>%
  ungroup() %>%
  select(gene_id, avg_gac, timepoint)

# Convert data frame to wide format
all_avg_gacs_wide <- all_genes_avg_gacs %>%
  pivot_wider(names_from = timepoint, values_from = avg_gac)

# Convert to a matrix for dtwclust compatibility
ts_data <- as.matrix(all_avg_gacs_wide[, -1]) # Exclude the ID column
rownames(ts_data) <- all_avg_gacs_wide$gene_id

# Time Series Clustering
# USE THIS ONE
pc <- tsclust(ts_data, k = 2L:4L, type = "partitional",  
              distance = "sbd", centroid = "pam", 
              seed = 3247L, trace = TRUE,
              args = tsclust_args(dist = list(window.size = 20L))) 
names(pc) <- paste0("k_", 2L:4L)
sapply(pc, cvi, type = "internal")

# Extract cluster assignments
cluster_assignments <- pc@cluster

# Extract centroids (for partitional clustering)
centroids <- pc@centroids

# Convert the time series data into a long format
ts_data_long <- data.frame(
  ID = rep(1:length(ts_data), each = length(ts_data[[1]])),
  Time = rep(1:length(ts_data[[1]]), length(ts_data)),
  Value = unlist(ts_data),
  Cluster = rep(cluster_assignments, each = length(ts_data[[1]]))
)

ts_data_long <- ts_data_long %>%
  pivot_longer(
    cols = starts_with("Value."),
    names_to = "Timepoint",
    values_to = "Value"
  ) %>%
  mutate(Timepoint = sub("Value.", "", Timepoint))

ts_data_long$Timepoint <- gsub("\\.", " ", ts_data_long$Timepoint)

# Convert the centroids into a long format
centroids_long <- data.frame(
  Cluster = rep(1:length(centroids), each = length(centroids[[1]])),
  Time = rep(1:length(centroids[[1]]), length(centroids)),
  Value = unlist(centroids)
)

# Rewrite t as Week Timepoints
centroids_long <- centroids_long %>%
  mutate(Week = case_when(
    Time == 1 ~ "Week 16",
    Time == 2 ~ "Week 20",
    Time == 3 ~ "Week 21",
    Time == 4 ~ "Week 24",
    TRUE ~ as.character(Time) # In case there are values other than 1, 2, 3, 4
  )) %>%
  select(Cluster, Time, Value, Week)

# Plot time series data colored by cluster assignments
ggplot() +
  geom_line(data = ts_data_long, aes(x = Timepoint, y = Value, group = ID, color = ID), alpha = 0.5) +
  geom_line(data = centroids_long, aes(x = Week, y = Value, group = Cluster), color = "black", size = 1) +
  facet_wrap(~Cluster, scales = "free_y") +
  coord_cartesian(ylim = c(0, 0.30)) +
  theme_minimal() +
  labs(title = "Time Series Clustering with Centroids",
       x = "Time Point",
       y = "Average Gene Accessibility Score",
       color = "Cluster")