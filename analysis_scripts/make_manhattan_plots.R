# Load libraries
library(dplyr)
library(ggplot2)


# # Discovery vss replication
# # Read the discovery file
# discovery_file <- "/Users/user/Desktop/beta_beta/merged_pit_volume_left_irn_primary.txt"
# disc <- read.table(discovery_file, header = TRUE)

# # Read the replication file
# replication_file <- "/Users/user/Desktop/beta_beta/merged_pit_volume_left_irn_replication.txt"
# repl <- read.table(replication_file, header = TRUE)


# # Merge by SNP ID (you can change this to CHR + POS if needed)
# merged <- inner_join(disc, repl, by = "ID", suffix = c("_disc", "_repl"))
# merged <- merged %>% select(ID, BETA_disc, BETA_repl)

# # Convert -log10(P) to P
# disc <- disc %>%
#   mutate(P_true = 10^(-LOG10P))  # assuming column P is -log10(P)

# # Filter for genome-wide significance
# disc_sig <- disc %>% filter(P_true < 5e-8)


# # Merge on ID (or CHR + POS if needed)
# merged_sig <- inner_join(disc_sig, repl, by = "ID", suffix = c("_disc", "_repl"))

# # Optional sanity check
# cat("Number of significant SNPs in discovery found in replication:", nrow(merged_sig), "\n")

# # Plot
# ggplot(merged_sig, aes(x = BETA_disc, y = BETA_repl)) +
#   geom_point(color = "blue", alpha = 0.7) +
#   geom_smooth(method = "lm", se = FALSE, color = "red") +
#   labs(
#     title = "Beta-Beta Plot for Discovery Significant SNPs",
#     x = "Discovery BETA",
#     y = "Replication BETA"
#   ) +
#   theme_minimal()

# # Correlation
# cor_test <- cor.test(merged_sig$BETA_disc, merged_sig$BETA_repl)
# print(cor_test)




# PLOT
# Load required libraries
library(data.table)
library(qqman)
library(ggplot2)
library(dplyr)
library(viridisLite)
library(data.table)
library(ggplot2)
library(viridis)
library(data.table)
library(CMplot)
library(topr)
library(stringr)

# List of traits to process
traits <- c("rim_height", "pit_volume", "cft", "pit_depth", "mean_slope", "rim_radius")

for (trait_name in traits) {
  
  # Load the Regenie summary statistics file
  input_path <- sprintf("merged_%s_left.regenie_filtered", trait_name)
  df <- fread(input_path, sep = " ", header = TRUE)
  
  # Ensure required columns exist and are renamed for topr
  df <- df %>%
    mutate(P = 10^(-LOG10P)) %>%
    rename(CHR = CHROM,  # no need to rename if keeping CHROM (topr accepts both CHR or CHROM)
         BP = GENPOS,
         SNP = ID) %>%
    select(CHR, BP, SNP, P)
  
  # Check chromosome formatting (topr requires numeric CHROM)
  df$CHR <- as.numeric(df$CHR)
  
  # Remove rows with missing or malformed values
  df_clean <- df %>%
    filter(!is.na(CHR), !is.na(P), !is.na(BP), !is.na(SNP))
  
  # Create Manhattan plot
  p <- manhattan(df_clean,
                 build = "37",
                 sign_thresh = 5e-8,
                 sign_thresh_label_size = 0,
                 annotate=5e-09)
  
  # Save plot
  ggsave(filename = sprintf("topr_manhattan_%s.png", trait_name),
         plot = p, width = 12, height = 6, dpi = 300)




  
  # Print SNP stats
  total_snps <- nrow(df_clean)
  sig_snps <- sum(df_clean$P < 5e-8)
  
  cat(sprintf("\nSummary for %s:\n", trait_name))
  cat(sprintf("Total SNPs: %s\n", format(total_snps, big.mark = ",")))
  cat(sprintf("Genome-wide significant SNPs (P < 5e-8): %s\n", format(sig_snps, big.mark = ",")))
}







####### Make the miami plot for the main text
# Load required libraries
library(data.table)
library(dplyr)
library(ggplot2)
library(topr)

setwd("~/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/manhattan")

# Traits to compare
traits <- c("pit_volume", "rim_radius")

# Function to load and format a trait's data
load_trait_data <- function(trait_name) {
  input_path <- sprintf("merged_%s_left.regenie_filtered", trait_name)
  df <- fread(input_path, sep = " ", header = TRUE)
  
  df %>%
    mutate(P = 10^(-LOG10P)) %>%
    rename(CHROM = CHROM,
           POS = GENPOS,
           SNP = ID) %>%
    select(CHROM, POS, SNP, P) %>%
    mutate(CHROM = case_when(
  CHROM %in% c("X", "23", 23) ~ "23",
  TRUE ~ as.character(CHROM)
)) %>%
    filter(!is.na(CHROM), !is.na(P), !is.na(POS), !is.na(SNP))
}

highlight_snps <- c("rs1042602")

# Load and prepare both datasets
df_list <- lapply(traits, load_trait_data)

# Plot both using TopR
p <- manhattan(
  df_list,
  sign_thresh = 5e-8,
  legend_labels = c("Pit volume", "Rim radius"),
  ntop = 1,
  sign_thresh_label_size = 0,
  build = "37"
)

# Save plot
  ggsave(filename = sprintf("topr_manhattan_main.png"),
         plot = p, width = 12, height = 6, dpi = 300)

ggsave("topr_manhattan_pit_vs_rim.png", plot = p, width = 12, height = 7, dpi = 300)




