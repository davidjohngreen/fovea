# -------------------------------------------------------------
# This script performs a broad set of downstream phenotypic and
# epidemiological analyses on OCT-derived foveal parameters.
#
# It serves as a working analysis script for integrating foveal
# morphology data with participant-level phenotypic information,
# generating descriptive statistics and figures, and running a
# range of exploratory and multivariable models.
#
# The script includes the following major components:
#
# 1. Data loading and integration
#    - imports cleaned left-eye foveal parameter data
#    - merges in ancestry assignments, spherical equivalent,
#      visual acuity (logMAR), sex, retinal pigmentation score,
#      scan dates, and other participant-level variables
#    - excludes participants with relevant retinal disorders
#      where required
#    - restricts analyses to unrelated individuals
#
# 2. Outlier handling
#    - trims extreme values for selected foveal traits by setting
#      the upper and lower tails to missing values on a
#      per-variable basis
#    - generates random samples of retained individuals for
#      manual review or QC checking
#
# 3. Derived variable construction
#    - calculates age at scan from scan timestamps and birth data
#    - prepares summary metrics such as fold-range statistics
#      across foveal traits
#
# 4. Descriptive visualisation
#    - produces histograms of key foveal parameters
#    - generates violin/boxplots stratified by sex and ancestry
#    - exports publication-style figures as PDF files
#
# 5. Summary statistics and group comparisons
#    - calculates trait-level means, medians, and standard
#      deviations
#    - compares parameter distributions between sexes
#    - summarises trait distributions across ancestry groups
#
# 6. Biological association analyses
#    - tests the relationship between foveal morphology and:
#         * ancestry
#         * retinal pigmentation score (RPS)
#         * eye size / spherical equivalent
#         * sex
#         * visual acuity
#    - fits multivariable linear regression models for the main
#      foveal traits
#    - extracts effect sizes, confidence intervals, p-values,
#      and model fit statistics
#    - generates figure-ready outputs from regression results
#
# 7. Averaged-eye analysis
#    - loads right-eye trait data
#    - applies equivalent filtering and trimming steps
#    - merges left and right eye data
#    - creates “smart-averaged” trait values by averaging both
#      eyes where available, while retaining single-eye values
#      where only one measurement is present
#    - writes an averaged analysis dataset for downstream use
#
# Overall, this script acts as a flexible working file for
# exploring the phenotypic architecture of foveal morphology and
# generating figures, tables, and model outputs used in the
# associated project.
#
# Note:
# This script is intentionally broad and contains a mixture of
# data preparation, exploratory analyses, figure generation, and
# model-fitting code developed iteratively during the project.
# Some sections are standalone, some are optional, and some are
# retained as commented reference code for alternative analyses.
# -------------------------------------------------------------



# Load necessary libraries
library(tidyverse)
library(ggplot2)
library(reshape2)
library(GGally)
library(broom)
library(tibble)
library(dplyr)
library(showtext)
library(tidyr)
library(ggpubr)
library(viridisLite)
library(effectsize)
library(purrr)
library(lubridate)
library(tidyverse)
library(showtext)
library(sysfonts)
library(ggplot2)
library(dplyr)
library(tidyr)
library(showtext)
library(sysfonts)





dates <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/dates_of_scan.csv')
dates <- dates %>% janitor::clean_names()



# ================================
#         Data Preparation
# ================================
## This refers to the left-sided data
# Set your directory containing the _fixed.csv files
data_dir <- "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/retimat_data"  # <-- replace this with your actual path

# Read and merge all *_fixed.csv files
all_files <- list.files(data_dir, pattern = "_fixed\\.csv$", full.names = TRUE)

df_list <- lapply(all_files, read.csv)

merged_df <- bind_rows(df_list)


### Read in the info about non_european_popultations
non_euro_pops <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/non_eur_pops.csv')


# Read in the file with genetic ethnic grouping
euro_pops <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/genetic_ethnic_group.csv') %>%
  rename(IID='Participant.ID')

merged_df <- left_join(merged_df, non_euro_pops, by='IID')
merged_df <- left_join(merged_df, euro_pops, by='IID')

merged_df <- merged_df %>%
  mutate(Ancestry = case_when(
    Genetic.ethnic.grouping == 1 ~ "EUR",
    Assigned == "AFR" ~ "AFR",
    Assigned == "SAS" ~ "SAS",
    Assigned == "EAS" ~ "EAS", 
    TRUE ~ NA_character_
  ))

merged_df <- merged_df %>%
  filter(!is.na(Ancestry))

# Select relevant columns
selected_df <- merged_df %>%
  select(IID, cft, pit_volume, pit_depth, mean_slope, rim_radius, Ancestry, rim_height)

# Convert to numeric just in case
selected_df <- selected_df %>%
  mutate(across(c(cft, pit_volume, pit_depth, mean_slope, rim_radius, rim_height), as.numeric))

# Import spherical error and do some work to produce one value of SER per individual
SPD <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/3__gwas/3__african/spherical_power_both_instances.csv')

# Join spherical equivalent data to the main dataframe
selected_df <- left_join(selected_df, SPD, by = 'IID')

# Import the logmar values
logmar <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/logmar_data.csv') %>%
  janitor::clean_names()

logmar <- logmar %>%
  mutate(
    logmar_diff = abs(log_mar_final_left_instance_1 - log_mar_final_left_instance_0),
    
    log_mar_final_left = case_when(
      !is.na(log_mar_final_left_instance_0) & !is.na(log_mar_final_left_instance_1) & logmar_diff < 0.1 ~
        (log_mar_final_left_instance_0 + log_mar_final_left_instance_1) / 2,
      
      !is.na(log_mar_final_left_instance_1) ~ log_mar_final_left_instance_1,
      !is.na(log_mar_final_left_instance_0) ~ log_mar_final_left_instance_0,
      
      TRUE ~ NA_real_
    )
  )

logmar <- logmar %>% rename(IID='iid') %>% select(IID, log_mar_final_left)

selected_df <- left_join(selected_df, logmar, by='IID')

sex <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/3__gwas/1__discovery/covariate_file.txt', sep='\t') %>% 
  select(IID, sex)

selected_df <- left_join(selected_df, sex, by='IID')

RPS <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/retinal_pigmentation_score.csv') %>%
  mutate(IID = str_extract(Name, "(?<=/)[0-9]{7}(?=_)")) %>%
  select(IID, pigmentation)

RPS$IID <- as.numeric(RPS$IID)

selected_df <- left_join(selected_df, RPS, by = 'IID')



# Read in list of people with H35
# This line is not used when preparing the cox data, as we need to retain those with disease
H35 <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/retinal_disorders_H35.csv')

selected_df <- selected_df %>% anti_join(H35, by = "IID")





# Read in the kinship data
oct_data <- read.csv("/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/1__image_QC/imaging_quality_metrics_both_instances.csv") %>%
  janitor::clean_names() %>% rename(IID = 'participant_id')


removed_kinship <- oct_data %>%
  filter(is.na(genetic_kinship_to_other_participants) |
         genetic_kinship_to_other_participants != 0)

oct_unrelated <- oct_data %>%
  filter(genetic_kinship_to_other_participants == 0)


selected_df$IID <- as.character(selected_df$IID)
oct_unrelated$IID <- as.character(oct_unrelated$IID)

selected_df_unrelated <- selected_df %>%
  semi_join(oct_unrelated, by = "IID")




# ================================
#        Outlier Filtering
# ================================
# Variables to trim
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius",
            "rim_height")

# Function to trim outliers to NA per variable
trim_outliers_to_na <- function(df, vars, trim = 0.01) {
  for (var in vars) {
    lower <- quantile(df[[var]], trim, na.rm = TRUE)
    upper <- quantile(df[[var]], 1 - trim, na.rm = TRUE)
    
    df[[var]] <- ifelse(!is.na(df[[var]]) & (df[[var]] < lower | df[[var]] > upper), NA, df[[var]])
  }
  return(df)
}

# Apply the trimming function
trimmed_df <- trim_outliers_to_na(selected_df_unrelated, params)




# Random sample included
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", 
            "rim_radius", "rim_height")

set.seed(123)  # optional for reproducibility

samples_list <- lapply(params, function(trait) {
  
  # rows where the trait exists
  available <- trimmed_df[!is.na(trimmed_df[[trait]]), ]
  
  # sample up to 100 (or fewer if not enough rows)
  n_take <- min(100, nrow(available))
  
  picked <- available[sample(nrow(available), n_take), "IID", drop = FALSE]
  
  picked$trait <- trait
  picked
})

final_samples <- do.call(rbind, samples_list)

write.csv(final_samples, "/Users/user/Desktop/appeal_work/outliers/sampled_IIDs_per_trait.csv", row.names = FALSE)


# write.csv(trimmed_df, file = '/Users/user/Desktop/appeal_work/cox_and_logistic/new_cox_without_related/for_cox_left.csv', quote = F)




#Time calculations
# Step 1: Convert timestamps to POSIXct
dates <- dates %>%
  mutate(
    scan_date_0 = ymd_hms(eye_measures_sign_off_timestamp_instance_0, quiet = TRUE),
    scan_date_1 = ymd_hms(eye_measures_sign_off_timestamp_instance_1, quiet = TRUE)
  )

# Step 2: Pick the earliest available scan date
dates <- dates %>%
  mutate(
    scan_date_used = pmin(scan_date_0, scan_date_1, na.rm = TRUE),
    scan_year = year(scan_date_used),
    scan_month = month(scan_date_used)
  )

# Step 3: Calculate age at scan
dates <- dates %>%
  mutate(
    age_at_scan = scan_year + (scan_month / 12) - (year_of_birth + (month_of_birth / 12))
  )

# Optional: Keep just participant ID and final age at scan
age_at_scan_df <- dates %>%
  select(participant_id, age_at_scan) %>%
  rename(IID = participant_id)

age_at_scan_df <- age_at_scan_df %>%
  mutate(IID = as.character(IID))

# Optional: Join to your trimmed_df
trimmed_df<- left_join(trimmed_df, age_at_scan_df, by = "IID")








# ================================
#     Calculate the fold range
# ================================
fold_range_summary <- sapply(params, function(trait) {
  vals <- trimmed_df[[trait]]
  max_val <- max(vals, na.rm = TRUE)
  min_val <- min(vals, na.rm = TRUE)
  q1 <- quantile(vals, 0.25, na.rm = TRUE)
  q3 <- quantile(vals, 0.75, na.rm = TRUE)
  
  c(
    max_min = max_val / min_val,
    iqr_ratio = q3 / q1
  )
})

# Transpose for better viewing
fold_range_df <- as.data.frame(t(fold_range_summary))
fold_range_df <- round(fold_range_df, 2)

fold_range_df





# ================================
#            Plotting
# ================================
# ================================
#           Histograms
# ================================

# Parameters in desired order and labels
params <- c("pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft")
param_labels <- c(
  "pit_volume" = "Pit volume (mm³)",
  "rim_radius" = "Rim radius (mm)",
  "rim_height" = "Rim height (µm)",
  "mean_slope" = "Mean slope (°)",
  "pit_depth" = "Pit depth (µm)",
  "cft" = "CFT (µm)"
)

# Prepare data
long_df <- trimmed_df |>
  select(all_of(params)) |>
  pivot_longer(cols = everything(), names_to = "Parameter", values_to = "Value") |>
  mutate(Parameter = factor(Parameter, levels = params, labels = param_labels))

# Plot
p_hist <- ggplot(long_df, aes(x = Value, fill = Parameter)) +
  geom_histogram(bins = 30, color = "white", size = 0.2, na.rm = TRUE) +
  facet_wrap(~Parameter, scales = "free", ncol = 3) +
  scale_x_continuous(expand = expansion(mult = c(0.05, 0.05))) +  # Add space on x-axis
  scale_fill_viridis_d(option = "D", begin = 0.2, end = 0.9, direction = -1) +
  coord_cartesian(clip = "off") +  # Allow tick labels to overflow
  theme_minimal(base_family = "Arial", base_size = 12) +
  labs(
    x = NULL,
    y = "Count"
  ) +
  theme(
    legend.position = "none",
    strip.text = element_text(family = "Arial", face = "bold", size = 12),
    plot.title = element_text(size = 15, face = "bold", hjust = 0.5),
    axis.title.y = element_text(margin = margin(r = 10)),
    axis.text = element_text(size = 12),
    panel.grid.major = element_line(size = 0.2),
    panel.grid.minor = element_blank(),
    plot.margin = unit(c(1, 2, 1, 1), "cm")  # top, right, bottom, left
  )

# Load and configure font (macOS)
font_add("Arial", regular = "/System/Library/Fonts/Supplemental/Arial.ttf")
showtext_auto()

# Save to PDF
ggsave(
  "/Users/user/Desktop/appeal_work/histograms_new.pdf",
  plot = p_hist,
  width = 10,
  height = 6,
  device = "pdf"
)


# # ================================
# #        Raw plots scatter
# # ================================
# # Variables to plot against logmar
# params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")

# plot_df_ser <- trimmed_df %>%
#   select(spherical_equivalent_left_final, all_of(params)) %>%
#   pivot_longer(cols = all_of(params), names_to = "parameter", values_to = "value") %>%
#   filter(!is.na(spherical_equivalent_left_final) & !is.na(value)) %>%
#   group_by(parameter) %>%
#   filter(n() > 0) %>%
#   ungroup()


# # Plot
# ggplot(plot_df_ser, aes(x = spherical_equivalent_left_final, y = value)) +
#   geom_point(alpha = 0.05, size = 0.8) +
#   geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
#   facet_wrap(~ parameter, scales = "free_y") +
#   labs(
#     x = "Spherical equivalent",
#     y = "Foveal parameter",
#   ) +
#   theme_minimal(base_size = 12)



# # Make the raw plot of plotting foveal parameters against logmar
# plot_df <- trimmed_df %>%
#   select(log_mar_final_left, all_of(params)) %>%
#   pivot_longer(cols = all_of(params), names_to = "parameter", values_to = "value") %>%
#   filter(!is.na(log_mar_final_left) & !is.na(value))

# ggplot(plot_df, aes(x = log_mar_final_left, y = value)) +
#   geom_point(alpha = 0.05, size = 0.8) +
#   geom_smooth(method = "lm", se = TRUE, color = "firebrick") +
#   facet_wrap(~ parameter, scales = "free_y") +
#   labs(
#     x = "logMAR (VA)",
#     y = "Foveal Parameter Value",
#     title = "Foveal Parameters vs logMAR"
#   ) +
#   theme_minimal(base_size = 12)




# Correlation between pit volume and logMAR
cor_pit <- cor(trimmed_df$pit_volume, trimmed_df$log_mar_final_left, use = "complete.obs")

# Correlation between rim radius and logMAR
cor_radius <- cor(trimmed_df$rim_radius, trimmed_df$log_mar_final_left, use = "complete.obs")

# Print results
cat("Correlation (pit_volume ~ logMAR):", cor_pit, "\n")
cat("Correlation (rim_radius ~ logMAR):", cor_radius, "\n")


# ================================
#        Sex violin plot
# ================================
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")
vcols <- viridis(n = 10, option = "viridis")
sex_colors <- c("Female" = vcols[2], "Male" = vcols[7])

# === Prepare data ===
sex_df <- trimmed_df %>%
  mutate(Sex = factor(sex, levels = c(0, 1), labels = c("Female", "Male"))) %>%
  pivot_longer(cols = all_of(params), names_to = "Parameter", values_to = "Value") %>%
  filter(!is.na(Value), !is.na(Sex))

# Set desired facet order
sex_df$Parameter <- factor(sex_df$Parameter, levels = c(
  "pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft"
))

# Define pretty names for facets
param_labels <- c(
  "pit_volume" = "Pit volume (mm³)",
  "rim_radius" = "Rim radius (mm)",
  "rim_height" = "Rim height (µm)",
  "mean_slope" = "Mean slope (°)",
  "pit_depth" = "Pit depth (µm)",
  "cft" = "CFT (µm)"
)

# === Plot ===
p_sex <- ggplot(sex_df, aes(x = Sex, y = Value, fill = Sex)) +
  geom_violin(trim = FALSE, alpha = 0.9) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  facet_wrap(~Parameter, scales = "free_y", labeller = labeller(Parameter = param_labels)) +
  scale_fill_manual(values = sex_colors) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "plain", size = 14),
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )



ggsave(
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/fig3.pdf",
  plot = p_sex,
  width = 10,
  height = 6,
  device = "pdf"
)



# ================================
#             Sex stats
# ================================
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")

# Compute summary stats and p-values
summary_stats <- lapply(params, function(p) {
  female_vals <- trimmed_df %>% filter(sex == 0) %>% pull(p)
  male_vals   <- trimmed_df %>% filter(sex == 1) %>% pull(p)

  test <- wilcox.test(female_vals, male_vals)

  tibble(
    Parameter = p,
    Female = sprintf("%.2f (%.2f)", mean(female_vals, na.rm = TRUE), median(female_vals, na.rm = TRUE)),
    Male = sprintf("%.2f (%.2f)", mean(male_vals, na.rm = TRUE), median(male_vals, na.rm = TRUE)),
    P_value = signif(test$p.value, 3)
  )
}) %>%
  bind_rows()

summary_stats


# ================================
#         Sex stats (split)
# ================================
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")

# Compute summary stats and p-values
summary_stats <- lapply(params, function(p) {
  female_vals <- trimmed_df %>% filter(sex == 0) %>% pull(p)
  male_vals   <- trimmed_df %>% filter(sex == 1) %>% pull(p)

  test <- wilcox.test(female_vals, male_vals)

  tibble(
    Parameter = p,
    Female_Mean = round(mean(female_vals, na.rm = TRUE), 3),
    Female_Median = round(median(female_vals, na.rm = TRUE), 3),
    Male_Mean = round(mean(male_vals, na.rm = TRUE), 3),
    Male_Median = round(median(male_vals, na.rm = TRUE), 3),
    P_value = signif(test$p.value, 3)
  )
}) %>%
  bind_rows()

# Print to console
print(summary_stats)
# Write to CSV
write.csv(summary_stats, "/Users/user/Desktop/appeal_work/new figures and tables/sex_difference_stats.csv", row.names = FALSE)





# ==========================
# 🌍    Ancestry Plots
# ==========================
# ==========================
# 🌍         FULL
# ==========================

library(dplyr)
library(ggplot2)
library(tidyr)

# Define trait order and pretty labels
params <- c("pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft")
param_labels <- c(
  "pit_volume" = "Pit volume (mm³)",
  "rim_radius" = "Rim radius (mm)",
  "rim_height" = "Rim height (µm)",
  "mean_slope" = "Mean slope (°)",
  "pit_depth" = "Pit depth (µm)",
  "cft" = "CFT (µm)"
)

# Define ancestry colours
ancestry_colors <- c("EUR" = vcols[2], "AFR" = vcols[5], "SAS" = vcols[9], "EAS" = vcols[7])

# Prepare data: reshape and apply ordered factors
trimmed_df_long <- trimmed_df %>%
  pivot_longer(cols = all_of(params), names_to = "Parameter", values_to = "Value") %>%
  filter(!is.na(Value) & !is.na(Ancestry)) %>%
  mutate(
    Ancestry = factor(Ancestry, levels = c("EUR", "AFR", "SAS", "EAS")),
    Parameter = factor(Parameter, levels = params)
  )

# Plot
# Plot
p_ancestry <- ggplot(trimmed_df_long, aes(x = Ancestry, y = Value, fill = Ancestry)) +
  geom_violin(trim = FALSE, alpha = 0.7) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  facet_wrap(~Parameter, scales = "free_y", labeller = labeller(Parameter = param_labels)) +
  scale_fill_manual(values = ancestry_colors) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1, size = 14),
    axis.text.y = element_text(size = 14),
    axis.title = element_text(face = "bold", size = 14),
    strip.text = element_text(face = "plain", size = 14),  # <-- make facet titles plain
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
  )

# Save the plot
ggsave(
  "/Users/user/Desktop/appeal_work/new figures and tables/suppfig1.pdf",
  plot = p_ancestry,
  width = 10,
  height = 6,
  device = "pdf"
)



# ================================
#     Pivoted ancestry stats
# ================================
params <- c("pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft")
ancestry_levels <- c("EUR", "AFR", "SAS", "EAS")

# Collect mean/median per ancestry × parameter
ancestry_stats <- lapply(params, function(p) {
  lapply(ancestry_levels, function(anc) {
    vals <- trimmed_df %>%
      filter(Ancestry == anc) %>%
      pull(p)

    tibble(
      Parameter = p,
      Ancestry = anc,
      Mean = round(mean(vals, na.rm = TRUE), 3),
      Median = round(median(vals, na.rm = TRUE), 3)
    )
  }) %>% bind_rows()
}) %>% bind_rows()

# Pivot wider
pivoted_stats <- ancestry_stats %>%
  pivot_wider(
    names_from = Ancestry,
    values_from = c(Mean, Median),
    names_glue = "{Ancestry}_{.value}"
  ) %>%
  select(Parameter, everything())  # Keep order clean


# Write to CSV
write.csv(pivoted_stats, "/Users/user/Desktop/appeal_work/new figures and tables/ancestry_trait_summary_stats.csv", row.names = FALSE)






# Trait order and parameter setup
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")
trait_order <- c(
  "Pit volume (mm³)", 
  "Rim radius (mm)", 
  "Rim height (µm)", 
  "Mean slope (°)", 
  "Pit depth (µm)", 
  "CFT (µm)"
)
vcols <- viridis(n = 10, option = "viridis")
ancestry_colors_eur_afr <- c("EUR" = vcols[4], "AFR" = vcols[9])

# === Prepare data ===
eur_afr_df <- trimmed_df %>%
  filter(Ancestry %in% c("EUR", "AFR")) %>%
  pivot_longer(cols = all_of(params), names_to = "Parameter", values_to = "Value") %>%
  filter(!is.na(Value)) %>%
  mutate(
    Ancestry = factor(Ancestry, levels = c("EUR", "AFR")),
    Parameter = recode(Parameter,
                        pit_volume = "Pit volume (mm³)",
                        rim_radius = "Rim radius (mm)",
                        rim_height = "Rim height (µm)",
                        mean_slope = "Mean slope (°)",
                        pit_depth = "Pit depth (µm)",
                        cft = "CFT (µm)"),
    Parameter = factor(Parameter, levels = trait_order)
  )



# Define facet labels
param_labels <- setNames(trait_order, trait_order)

# === Plot ===
p_ancestry_eur_afr <- ggplot(eur_afr_df, aes(x = Ancestry, y = Value, fill = Ancestry)) +
  geom_violin(trim = FALSE, alpha = 0.8) +
  geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
  facet_wrap(~Parameter, scales = "free_y", labeller = labeller(Parameter = param_labels)) +
  scale_fill_manual(values = ancestry_colors_eur_afr) +
  labs(
    x = NULL,
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
  legend.position = "none",
  axis.text.x = element_text(angle = 0, hjust = 0.5, size = 14),
  axis.text.y = element_text(size = 14),
  axis.title = element_text(face = "bold", size = 14),
  strip.text = element_text(face = "plain", size = 14),
  plot.title = element_text(face = "bold", hjust = 0.5, size = 14)
)



ggsave(
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/fig4.pdf",
  plot = p_ancestry_eur_afr,
  width = 10,
  height = 6,
  device = "pdf"
)



# # ---------------------------------
# # Statistical analysis
# # ---------------------------------
# # Parameters and predictors
# # Produce the plot showing effect sizes for visual acuity, SER, and RPS
# foveal_params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")
# main_predictors <- list(
#   SER = "spherical_equivalent_left_final",
#   RPS = "pigmentation",
#   VA  = "log_mar_final_left"
# )

# # Collect model results
# model_results <- list()

# for (param in foveal_params) {
#   for (label in names(main_predictors)) {
#     predictor <- main_predictors[[label]]
#     formula_str <- paste(param, "~", paste(c(predictor, "spherical_equivalent_left_final", "sex", "Ancestry"), collapse = " + "))
#     if (label == "SER") {
#       formula_str <- paste(param, "~", paste(c(predictor, "sex", "Ancestry"), collapse = " + "))
#     }

#     model <- lm(as.formula(formula_str), data = trimmed_df)
#     coefs <- summary(model)$coefficients
#     r2 <- summary(model)$adj.r.squared

#     if (predictor %in% rownames(coefs)) {
#       beta <- coefs[predictor, "Estimate"]
#       se   <- coefs[predictor, "Std. Error"]
#       pval <- coefs[predictor, "Pr(>|t|)"]
#     } else {
#       beta <- NA
#       se   <- NA
#       pval <- NA
#     }

#     model_results[[length(model_results) + 1]] <- tibble(
#       Parameter   = param,
#       Predictor   = label,
#       Beta        = beta,
#       SE          = se,
#       Lower_CI    = beta - 1.96 * se,
#       Upper_CI    = beta + 1.96 * se,
#       P_value     = signif(pval, 3),
#       R_squared   = signif(r2, 3)
#     )
#   }
# }

# summary_df <- bind_rows(model_results)


# # Do a plot for this
# summary_df <- summary_df %>%
#   mutate(Parameter = factor(Parameter, levels = rev(foveal_params)))

# ggplot(summary_df, aes(x = Beta, y = Parameter)) +
#   geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
#   geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI), size = 0.3) +
#   facet_wrap(~ Predictor, scales = "free_x") +
#   labs(
#     title = "Effect of Eye Size, Pigmentation, and Visual Acuity on Foveal Parameters",
#     x = "Effect Size (Beta ± 95% CI)",
#     y = "Foveal Parameter"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(strip.text = element_text(face = "bold"))




# # RPS only, not using ancestry asd covariate
# foveal_params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")

# rps_only_results <- lapply(foveal_params, function(p) {
#   model <- lm(reformulate(c("pigmentation", "spherical_equivalent_left_final", "sex"), p), data = trimmed_df)
#   coefs <- summary(model)$coefficients
#   r2 <- summary(model)$adj.r.squared

#   tibble(
#     Parameter   = p,
#     Beta_RPS    = coefs["pigmentation", "Estimate"],
#     P_RPS       = signif(coefs["pigmentation", "Pr(>|t|)"], 3),
#     R_squared   = signif(r2, 3)
#   )
# }) %>%
#   bind_rows()



### get simple summary stats for the traits
library(dplyr)

traits <- c("cft", "pit_volume", "rim_radius", 
            "rim_height", "mean_slope", "pit_depth")

summary_stats <- trimmed_df %>%
  summarise(across(all_of(traits),
                   list(mean = ~mean(.x, na.rm = TRUE),
                        sd   = ~sd(.x, na.rm = TRUE)),
                   .names = "{.col}_{.fn}"))

summary_stats





# -----------------------------------------------------------------------------------
# BIOLOGICAL QUESTION #1: what aspects of foveal morphology are associated with ancestry 
# and to what extent this signal is driven by pigmentation (RPS)
# -----------------------------------------------------------------------------------
# Relevel ancestry to set EUR as reference
trimmed_df$Ancestry <- relevel(factor(trimmed_df$Ancestry), ref = "EUR")

# Define traits and initialize storage
foveal_traits <- c("pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft")
results_list <- list()
effectsize_list <- list()

# Loop through traits
for (trait in foveal_traits) {
  # Model 1: RPS only
  model_rps <- lm(reformulate(c("pigmentation", "spherical_equivalent_left_final", "sex", "age_at_scan"), trait), data = trimmed_df)
  r2_rps <- summary(model_rps)$adj.r.squared

  # Model 2: Ancestry only
  model_ancestry <- lm(reformulate(c("Ancestry", "spherical_equivalent_left_final", "sex", "age_at_scan"), trait), data = trimmed_df)
  r2_ancestry <- summary(model_ancestry)$adj.r.squared

  # Model 3: Combined
  model_combined <- lm(reformulate(c("pigmentation", "Ancestry", "spherical_equivalent_left_final", "sex", "age_at_scan"), trait), data = trimmed_df)
  r2_combined <- summary(model_combined)$adj.r.squared

  # Extract standardized betas
  std_betas <- standardize_parameters(model_combined) %>%
    filter(grepl("^pigmentation$|^Ancestry", Parameter)) %>%
    mutate(
      Trait = trait,
      Predictor = case_when(
        Parameter == "pigmentation" ~ "RPS",
        grepl("^Ancestry", Parameter) ~ paste0("Ancestry: ", gsub("Ancestry", "", Parameter))
      )
    ) %>%
    select(Trait, Predictor, Std_Coefficient)

  # Store results
  results_list[[trait]] <- tibble(
    Parameter = trait,
    R2_RPS_Only = r2_rps,
    R2_Ancestry_Only = r2_ancestry,
    R2_RPS_and_Ancestry = r2_combined,
    R2_RPS_Unique = r2_combined - r2_ancestry,
    R2_Ancestry_Unique = r2_combined - r2_rps
  )

  effectsize_list[[trait]] <- std_betas
}

# Combine R2 and beta outputs
final_r2_table <- bind_rows(results_list)
std_beta_df <- bind_rows(effectsize_list)

# Plot 1: Unique R2
r2_stacked <- final_r2_table %>%
  select(Parameter, R2_RPS_Unique, R2_Ancestry_Unique) %>%
  pivot_longer(cols = -Parameter, names_to = "Source", values_to = "R2") %>%
  mutate(Source = recode(Source,
                         R2_RPS_Unique = "RPS (unique)",
                         R2_Ancestry_Unique = "Ancestry (unique)"))

pretty_names <- c(
  cft = "CFT",
  pit_volume = "Pit volume",
  pit_depth = "Pit depth",
  mean_slope = "Mean slope",
  rim_radius = "Rim radius",
  rim_height = "Rim height"
)

r2_stacked$Parameter <- factor(r2_stacked$Parameter, levels = c("pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft"))


p_r2 <- ggplot(r2_stacked, aes(x = Parameter, y = R2, fill = Source)) +
  geom_col(width = 0.7, color = "white", linewidth = 0.2) +
  # geom_text(aes(label = round(R2, 3)), position = position_stack(vjust = 0.5), size = 3.5, color = "black") +
  scale_x_discrete(labels = pretty_names) +  # <--- HERE
  scale_fill_viridis_d(option = "viridis", begin = 0.5, end = 0.9) +
  coord_cartesian(ylim = c(-0.025, 0.05)) + # adjust upper limit as needed 
  labs(
    y = "Unique adjusted R²",
    x = "Foveal parameter",
    fill = "Predictor"
  ) +
  theme_minimal(base_size = 13) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top"
  )



ggsave(
  "/Users/user/Desktop/appeal_work/new figures and tables/fig5.pdf",
  plot = p_r2,
  width = 10,
  height = 6,
  device = "pdf"
)



#GRAB the table for suppl table s5
# Function to extract overall model p-value from lm summary
get_model_p <- function(model) {
  s <- summary(model)
  pf(s$fstatistic[1], s$fstatistic[2], s$fstatistic[3], lower.tail = FALSE)
}

# Initialize lists
results_list <- list()

for (trait in foveal_traits) {
  # Model 1: RPS only
  model_rps <- lm(reformulate(c("pigmentation", "spherical_equivalent_left_final", "sex", "age_at_scan"), trait), data = trimmed_df)
  r2_rps <- summary(model_rps)$adj.r.squared
  p_rps <- get_model_p(model_rps)

  # Model 2: Ancestry only
  model_ancestry <- lm(reformulate(c("Ancestry", "spherical_equivalent_left_final", "sex", "age_at_scan"), trait), data = trimmed_df)
  r2_ancestry <- summary(model_ancestry)$adj.r.squared
  p_ancestry <- get_model_p(model_ancestry)

  # Model 3: Combined
  model_combined <- lm(reformulate(c("pigmentation", "Ancestry", "spherical_equivalent_left_final", "sex", "age_at_scan"), trait), data = trimmed_df)
  r2_combined <- summary(model_combined)$adj.r.squared
  p_combined <- get_model_p(model_combined)

  # Store results
  results_list[[trait]] <- tibble(
    Parameter = trait,
    R2_RPS_Only = r2_rps,
    P_RPS_Only = p_rps,
    R2_Ancestry_Only = r2_ancestry,
    P_Ancestry_Only = p_ancestry,
    R2_RPS_and_Ancestry = r2_combined,
    P_RPS_and_Ancestry = p_combined,
    R2_RPS_Unique = r2_combined - r2_ancestry,
    R2_Ancestry_Unique = r2_combined - r2_rps
  )
}

results_df_final <- bind_rows(results_list)




# -----------------------------------------------------------------------------------
# BIOLOGICAL QUESTION #2: what aspects of foveal morphology are associated with eye size
# and sex
# -----------------------------------------------------------------------------------
# Define traits to analyse
foveal_traits <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")

# Initialise result holder
sex_ser_results <- list()

# Loop through each trait
for (trait in foveal_traits) {
  # Linear model with just sex and eye size (SER)
  model <- lm(reformulate(c("sex", "spherical_equivalent_left_final"), trait), data = trimmed_df)
  s <- summary(model)
  
  # Extract key stats
  sex_ser_results[[trait]] <- tibble(
    Parameter = trait,
    Beta_Sex = coef(s)["sex", "Estimate"],
    P_Sex = coef(s)["sex", "Pr(>|t|)"],
    Beta_SER = coef(s)["spherical_equivalent_left_final", "Estimate"],
    P_SER = coef(s)["spherical_equivalent_left_final", "Pr(>|t|)"],
    R2_Adjusted = s$adj.r.squared
  )
}

# Combine
sex_ser_table <- bind_rows(sex_ser_results)
print(sex_ser_table)







# -----------------------------------------------------------------------------------
# BIOLOGICAL QUESTION #3: can visual acuity be predicted from pit morphology?
# -----------------------------------------------------------------------------------

library(dplyr)
library(ggplot2)
library(readr)
library(tibble)

# Define foveal parameters and pretty names
foveal_traits <- c("pit_volume", "rim_radius", "rim_height", "mean_slope", "pit_depth", "cft")
pretty_names <- c(
  cft = "CFT",
  pit_volume = "Pit volume",
  pit_depth = "Pit depth",
  mean_slope = "Mean slope",
  rim_radius = "Rim radius",
  rim_height = "Rim height"
)

# Create a list to store results
va_model_results <- list()

# Loop over each foveal trait as predictor of visual acuity
for (trait in foveal_traits) {
  formula <- as.formula(paste("log_mar_final_left ~", paste(c(trait, "spherical_equivalent_left_final", "sex", "Ancestry", 'age_at_scan'), collapse = " + ")))
  model <- lm(formula, data = trimmed_df)
  summary_model <- summary(model)
  
  coef_data <- summary_model$coefficients
  r2 <- summary_model$adj.r.squared
  
  if (trait %in% rownames(coef_data)) {
    beta <- coef_data[trait, "Estimate"]
    se <- coef_data[trait, "Std. Error"]
    pval <- coef_data[trait, "Pr(>|t|)"]
  } else {
    beta <- NA
    se <- NA
    pval <- NA
  }
  
  va_model_results[[trait]] <- tibble(
    Parameter = trait,
    Beta = beta,
    SE = se,
    Lower_CI = beta - 1.96 * se,
    Upper_CI = beta + 1.96 * se,
    P_value = signif(pval, 3),
    R_squared = signif(r2, 3)
  )
}

# Combine into one dataframe
va_summary_df <- bind_rows(va_model_results)

# Apply correct plotting order and pretty names
va_summary_df <- va_summary_df %>%
  mutate(Parameter = factor(Parameter, levels = rev(foveal_traits), labels = pretty_names[rev(foveal_traits)]))

# Set axis limits
xmax <- max(abs(c(va_summary_df$Lower_CI, va_summary_df$Upper_CI)), na.rm = TRUE) * 1.2

# Plot
p <- ggplot(va_summary_df, aes(x = Beta, y = Parameter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI), size = 0.4) +
  scale_x_continuous(limits = c(-xmax, xmax)) +
  labs(
    x = "Effect Size (Beta ± 95% CI)",
    y = NULL
  ) +
  theme_minimal(base_size = 13) +
  theme(
    strip.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(face = "bold", size = 15),
    plot.subtitle = element_text(size = 12)
  )

# Save plot
ggsave(
  "/Users/user/Desktop/appeal_work/new figures and tables/suppfig2.pdf",
  plot = p,
  width = 8,
  height = 8,
  device = "pdf"
)






# Plot
ggplot(va_summary_df, aes(x = Beta, y = Parameter)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
  geom_pointrange(aes(xmin = Lower_CI, xmax = Upper_CI), size = 0.4) +
  labs(
    title = "Foveal Morphology Predicting Visual Acuity",
    subtitle = "Adjusted for sex, SER, and ancestry",
    x = "Effect Size (Beta ± 95% CI)",
    y = "Foveal Parameter"
  ) +
  theme_minimal(base_size = 13) +
  theme(strip.text = element_text(face = "bold"))


# # AMD and Glaucoma
# # Set up the parameters and outcomes
# foveal_traits <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")
# binary_outcomes <- c("amd", "glaucoma", "erm", "macular_hole", "retinal_detachment")

# # Collect results
# logistic_results <- list()

# for (outcome in binary_outcomes) {
#   for (trait in foveal_traits) {
    
#     formula_str <- paste0(outcome, " ~ ", trait, " + sex + spherical_equivalent_left_final + Ancestry")
#     model <- glm(as.formula(formula_str), data = trimmed_df, family = binomial)
#     coefs <- summary(model)$coefficients

#     # Check that the trait was retained in the model
#     if (trait %in% rownames(coefs)) {
#       beta <- coefs[trait, "Estimate"]
#       se <- coefs[trait, "Std. Error"]
#       p_val <- coefs[trait, "Pr(>|z|)"]
      
#       logistic_results[[length(logistic_results) + 1]] <- tibble(
#         Outcome = outcome,
#         Parameter = trait,
#         OR = exp(beta),
#         Lower_CI = exp(beta - 1.96 * se),
#         Upper_CI = exp(beta + 1.96 * se),
#         P_value = signif(p_val, 3)
#       )
#     }
#   }
# }

# # Combine all results
# logistic_summary <- bind_rows(logistic_results)

# # View the results
# print(logistic_summary)

# library(ggplot2)
# library(dplyr)

# logistic_summary %>%
#   mutate(Parameter = gsub("_", " ", Parameter),
#          Outcome = tools::toTitleCase(Outcome)) %>%
#   ggplot(aes(x = OR, y = reorder(Parameter, OR))) +
#   geom_vline(xintercept = 1, linetype = "dashed", color = "grey60") +
#   geom_point(color = "darkgreen", size = 2) +
#   geom_errorbarh(aes(xmin = Lower_CI, xmax = Upper_CI), height = 0.2, color = "darkgreen") +
#   facet_wrap(~ Outcome, scales = "free_x") +
#   scale_x_log10() +
#   labs(
#     x = "Odds Ratio (log scale)",
#     y = "Foveal Parameter",
#     title = "Association Between Foveal Parameters and Ocular Diseases"
#   ) +
#   theme_minimal(base_size = 13) +
#   theme(strip.text = element_text(face = "bold"))





#### Multivariable linear regression for the six foveal traits
library(dplyr)
library(broom)
library(tibble)

run_multivariable_models <- function(
  df,
  output_csv,
  traits = c("cft", "pit_volume", "pit_depth",
             "mean_slope", "rim_radius", "rim_height"),
  predictors = "age_at_scan + spherical_equivalent_left_final + sex + Ancestry"
) {
  # ---- 1. Prepare categorical variables ----
  df$sex <- factor(df$sex)
  df$Ancestry <- factor(df$Ancestry)

  # Set EUR as reference
  df$Ancestry <- relevel(df$Ancestry, ref = "EUR")

  # ---- 2. Run models ----
  all_results <- list()

  for (trait in traits) {

    # Z-scale outcome
    df$y_std <- scale(df[[trait]])

    # Build formula dynamically
    f <- as.formula(paste("y_std ~", predictors))

    # Fit model
    model <- lm(f, data = df)

    # Tidy
    tb <- broom::tidy(model)

    # Confidence intervals
    ci <- confint(model)
    ci_tb <- as_tibble(ci, rownames = "term") |>
      rename(conf.low = `2.5 %`, conf.high = `97.5 %`)

    merged <- tb |>
      left_join(ci_tb, by = "term") |>
      mutate(
        outcome = trait,
        n_used = model$df.residual + length(model$coefficients)
      )

    all_results[[trait]] <- merged
  }

  # ---- 3. Bind and save ----
  final_table <- bind_rows(all_results) |>
    select(outcome, term, estimate, std.error, statistic, p.value,
           conf.low, conf.high, n_used)

  write.csv(final_table, output_csv, row.names = FALSE)

  return(final_table)
}

df_to_use <- trimmed_df
output_file <- "/Users/user/Desktop/appeal_work/new figures and tables/multivariable_models_foveal_traits_left.csv"

results <- run_multivariable_models(
  df = df_to_use,
  output_csv = output_file
)




library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(dplyr)
library(ggplot2)
library(readr)
library(patchwork)
library(dplyr)
library(readr)

df <- read_csv("/Users/user/Desktop/multivariable_models_foveal_traits_averaged.csv")

# Rows your supervisor wants
terms_needed <- c(
  "age_at_scan",
  "spherical_equivalent_left_final",
  "sex1",          # will be recoded to female vs male
  "AncestryAFR"
)

out <- df %>%
  filter(term %in% terms_needed) %>%
  mutate(
    term_clean = recode(term,
      "age_at_scan" = "Age (1 year)",
      "spherical_equivalent_left_final" = "Refraction (SE, 1 dioptre)",
      
      # flip the sex effect
      "sex1" = "Sex (female)",
      
      "AncestryAFR" = "Ancestry (AFR v EUR)"
    ),
    
    # invert male→female contrast
    estimate  = ifelse(term == "sex1", -estimate, estimate),
    conf.low  = ifelse(term == "sex1", -conf.low, conf.low),
    conf.high = ifelse(term == "sex1", -conf.high, conf.high)
  ) %>%
  select(outcome, term_clean, estimate, conf.low, conf.high, p.value)

out




# ---- 1. Clean outcome labels ----
plot_df <- out %>%
  mutate(
    outcome = recode(outcome,
      "cft"        = "CFT",
      "pit_volume" = "Pit volume",
      "pit_depth"  = "Pit depth",
      "mean_slope" = "Mean slope",
      "rim_radius" = "Rim radius",
      "rim_height" = "Rim height"
    )
  )

# ---- 2. Order predictors ----
plot_df$term_clean <- factor(
  plot_df$term_clean,
  levels = c(
    "Age (1 year)",
    "Refraction (SE, 1 dioptre)",
    "Sex (female)",
    "Ancestry (AFR v EUR)"
  )
)

# ---- 3. Order traits ----
plot_df$outcome <- factor(
  plot_df$outcome,
  levels = c(
    "CFT",
    "Pit volume",
    "Pit depth",
    "Mean slope",
    "Rim radius",
    "Rim height"
  )
)

# ---- 4. Panel function ----
make_panel <- function(trait_name) {
  ggplot(
    filter(plot_df, outcome == trait_name),
    aes(x = term_clean, y = estimate)
  ) +
    geom_point(size = 2, colour = "black") +
    geom_errorbar(
      aes(ymin = conf.low, ymax = conf.high),
      width = 0, colour = "black"
    ) +
    geom_hline(
      yintercept = 0, linetype = 2, colour = "grey40"
    ) +
    coord_flip() +
    labs(
      title = trait_name,
      x = "",
      y = "Effect (SD units)"
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      panel.grid.minor = element_blank()
    )
}

# ---- 5. Build panels ----
p1 <- make_panel("CFT")
p2 <- make_panel("Pit volume")
p3 <- make_panel("Pit depth")
p4 <- make_panel("Mean slope")
p5 <- make_panel("Rim radius")
p6 <- make_panel("Rim height")


# Combine using patchwork (2 rows × 3 columns)
final_fig <- (p2 | p5 | p6) /
             (p4 | p3 | p1)

final_fig


# ---- 6. Save to PDF ----
ggsave(
  filename = "/Users/user/Desktop/multivariable_models_foveal_traits_averaged.pdf",
  plot     = final_fig,
  width    = 12,      # adjust as needed
  height   = 8,       # adjust as needed
  units    = "in"
)





##########################################################
##########################################################
##########################################################

# New averaged analysis

##########################################################
##########################################################
##########################################################
path <- "/Users/user/Desktop/appeal_work/right_vs_left/retimat_output"

right_eye_traits <- list.files(
  path,
  pattern = "_fixed\\.csv$",
  full.names = TRUE
) %>%
  lapply(read.csv, header = TRUE) %>%
  bind_rows()

right_eye_traits

right_eye_traits$IID <- sub("_.*", "", right_eye_traits$filename)

traits <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")
right_sub <- right_eye_traits %>%
  select(IID, all_of(traits))




# Convert to numeric just in case
right_sub <- right_sub %>%
  mutate(across(c(IID), as.numeric))



right_sub <- left_join(right_sub, non_euro_pops, by='IID')
right_sub <- left_join(right_sub, euro_pops, by='IID')

right_sub <- right_sub %>%
  mutate(Ancestry = case_when(
    Genetic.ethnic.grouping == 1 ~ "EUR",
    Assigned == "AFR" ~ "AFR",
    Assigned == "SAS" ~ "SAS",
    Assigned == "EAS" ~ "EAS", 
    TRUE ~ NA_character_
  ))

right_sub <- right_sub %>%
  filter(!is.na(Ancestry))



# Join spherical equivalent data to the main dataframe
right_sub <- left_join(right_sub, SPD, by = 'IID')


right_sub <- left_join(right_sub, logmar, by='IID')

# SEX
right_sub <- left_join(right_sub, sex, by='IID')

# RPS
right_sub <- left_join(right_sub, RPS, by = 'IID')




age_at_scan_df <- age_at_scan_df %>%
  mutate(across(c(IID), as.numeric))
right_sub <- left_join(right_sub, age_at_scan_df, by = "IID")



# Read in the kinship data
oct_data <- read.csv("/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/1__image_QC/imaging_quality_metrics_both_instances.csv") %>%
  janitor::clean_names() %>% rename(IID = 'participant_id')


removed_kinship <- oct_data %>%
  filter(is.na(genetic_kinship_to_other_participants) |
         genetic_kinship_to_other_participants != 0)



oct_unrelated <- oct_data %>%
  filter(genetic_kinship_to_other_participants == 0)


right_sub$IID <- as.character(right_sub$IID)
oct_unrelated$IID <- as.character(oct_unrelated$IID)

right_sub <- right_sub %>%
  semi_join(oct_unrelated, by = "IID")


right_sub <- right_sub %>% select(-Genetic.ethnic.grouping, -Assigned)


right_sub$IID <- as.integer(right_sub$IID)
#Remove the retinal disorders
right_sub <- right_sub %>% anti_join(H35, by = "IID")




# need to add logmar, pigmentation, SER, sex, ancestry







# ================================
#        Outlier Filtering for Right Eye
# ================================
# Variables to trim
params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius",
            "rim_height")

# Function to trim outliers to NA per variable
trim_outliers_to_na <- function(df, vars, trim = 0.01) {
  for (var in vars) {
    lower <- quantile(df[[var]], trim, na.rm = TRUE)
    upper <- quantile(df[[var]], 1 - trim, na.rm = TRUE)
    
    df[[var]] <- ifelse(!is.na(df[[var]]) & (df[[var]] < lower | df[[var]] > upper), NA, df[[var]])
  }
  return(df)
}

# Apply the trimming function
trimmed_right_df <- trim_outliers_to_na(right_sub, params)



left_df <- trimmed_df %>%
  rename_with(~ paste0(.x, "_L"),
              .cols = !c(IID, Ancestry, spherical_equivalent_left_final, log_mar_final_left, sex, pigmentation, age_at_scan))

right_df <- trimmed_right_df %>%
  rename_with(~ paste0(.x, "_R"),
              .cols = !c(IID, Ancestry, spherical_equivalent_left_final, log_mar_final_left, sex, pigmentation, age_at_scan))


left_df$IID  <- as.character(left_df$IID)
right_df$IID <- as.character(right_df$IID)






#### SMART AVERAGING
merged <- full_join(left_df, right_df, by = "IID")



merged <- merged %>%
  mutate(
    Ancestry = coalesce(Ancestry.x, Ancestry.y),
    spherical_equivalent_left_final = coalesce(spherical_equivalent_left_final.x,
                                               spherical_equivalent_left_final.y),
    log_mar_final_left = coalesce(log_mar_final_left.x, log_mar_final_left.y),
    sex = coalesce(sex.x, sex.y),
    pigmentation = coalesce(pigmentation.x, pigmentation.y),
    age_at_scan = coalesce(age_at_scan.x, age_at_scan.y)
  )

merged <- merged %>%
  select(
    -ends_with(".x"),
    -ends_with(".y")
  )


averaged <- merged %>%
  mutate(across(
    ends_with("_L"),
    ~ {
      left_val  <- .
      right_val <- merged[[sub("_L$", "_R", cur_column())]]
      # fallback logic:
      ifelse(
        !is.na(left_val) & !is.na(right_val),  (left_val + right_val) / 2,  # both present
        ifelse(!is.na(left_val), left_val, right_val)                        # one present
      )
    },
    .names = "{sub('_L$', '', .col)}"
  ))



write.csv(averaged, file = '/Users/user/Desktop/appeal_work/cox_and_logistic/new_cox_without_related/for_cox_smart_avg.csv', quote = F)













# # --------------------------------------------------------
# # Start with the new averaged dataset
# # --------------------------------------------------------
# df <- final_df   # rename for convenience
# df$IID <- as.character(df$IID)   # ensure IID is character everywhere
# SPD <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/3__gwas/3__african/spherical_power_both_instances.csv')

# SPD$IID <- as.character(SPD$IID)

# df <- left_join(df, SPD, by = "IID")
# logmar <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/logmar_data.csv') %>%
#   janitor::clean_names()

# logmar <- logmar %>%
#   mutate(
#     logmar_diff = abs(log_mar_final_left_instance_1 - log_mar_final_left_instance_0),
    
#     log_mar_final_left = case_when(
#       !is.na(log_mar_final_left_instance_0) & 
#       !is.na(log_mar_final_left_instance_1) & 
#       logmar_diff < 0.1 ~ (log_mar_final_left_instance_0 + log_mar_final_left_instance_1) / 2,
      
#       !is.na(log_mar_final_left_instance_1) ~ log_mar_final_left_instance_1,
#       !is.na(log_mar_final_left_instance_0) ~ log_mar_final_left_instance_0,
      
#       TRUE ~ NA_real_
#     )
#   ) %>%
#   rename(IID = iid) %>%
#   mutate(IID = as.character(IID)) %>%
#   select(IID, log_mar_final_left)

# df <- left_join(df, logmar, by = "IID")

# sex <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/3__gwas/1__discovery/covariate_file.txt',
#                 sep = '\t') %>%
#   select(IID, sex)

# sex$IID <- as.character(sex$IID)

# df <- left_join(df, sex, by = "IID")

# RPS <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/retinal_pigmentation_score.csv') %>%
#   mutate(IID = str_extract(Name, "(?<=/)[0-9]{7}(?=_)")) %>%
#   mutate(IID = as.character(IID)) %>%
#   select(IID, pigmentation)

# df <- left_join(df, RPS, by = "IID")

# H35 <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/retinal_disorders_H35.csv')

# H35$IID <- as.character(H35$IID)

# df <- df %>% anti_join(H35, by = "IID")

# params <- c("cft", "pit_volume", "pit_depth", "mean_slope", "rim_radius", "rim_height")

# trim_outliers_to_na <- function(df, vars, trim = 0.01) {
#   for (var in vars) {
#     lower <- quantile(df[[var]], trim, na.rm = TRUE)
#     upper <- quantile(df[[var]], 1 - trim, na.rm = TRUE)
    
#     df[[var]] <- ifelse(!is.na(df[[var]]) & 
#                         (df[[var]] < lower | df[[var]] > upper),
#                         NA, df[[var]])
#   }
#   df
# }

# trimmed_df <- trim_outliers_to_na(df, params)

