# File: nano_analysis/scripts/plot_coverage.R

# Load required libraries
library(tidyverse)
library(ggplot2)

# Get parameters and inputs
bedgraph_files <- snakemake@input
results_dir <- snakemake@params$results_dir

# Function to read bedgraph file and add sample information
read_bedgraph <- function(file) {
  sample_name <- basename(file) %>% 
    str_replace(".bedgraph", "")
  
  read_tsv(file, col_names = c("chromosome", "start", "end", "coverage"), 
           col_types = "ciid") %>%
    mutate(sample = sample_name)
}

# Read all bedgraph files
coverage_data <- map_dfr(bedgraph_files, read_bedgraph)

# Create coverage plot for all samples
p <- ggplot(coverage_data, aes(x = start, y = coverage, color = sample)) +
  geom_line(alpha = 0.7) +
  facet_wrap(~chromosome, scales = "free_x") +
  labs(title = "Coverage Across All Samples",
       x = "Position (bp)",
       y = "Coverage Depth") +
  theme_minimal() +
  theme(legend.position = "bottom")

# Save plot
ggsave(snakemake@output[[1]], plot = p, width = 12, height = 8)

# Create per-sample coverage plots
samples <- unique(coverage_data$sample)

sample_plots <- ggplot(coverage_data, aes(x = start, y = coverage)) +
  geom_line(alpha = 0.7) +
  facet_grid(sample~chromosome, scales = "free") +
  labs(title = "Per-Sample Coverage",
       x = "Position (bp)",
       y = "Coverage Depth") +
  theme_minimal()

# Save per-sample plots
per_sample_plot_path <- file.path(dirname(snakemake@output[[1]]), "per_sample_coverage.pdf")
ggsave(per_sample_plot_path, plot = sample_plots, width = 14, height = 10)

# Calculate summary statistics
coverage_stats <- coverage_data %>%
  group_by(sample, chromosome) %>%
  summarize(
    mean_coverage = mean(coverage),
    median_coverage = median(coverage),
    max_coverage = max(coverage),
    min_coverage = min(coverage),
    regions_below_10x = sum(coverage < 10),
    .groups = "drop"
  )

# Save summary statistics
stats_path <- file.path(dirname(snakemake@output[[1]]), "coverage_statistics.csv")
write_csv(coverage_stats, stats_path)