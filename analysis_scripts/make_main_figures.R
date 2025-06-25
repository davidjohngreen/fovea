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




dates <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/dates_of_scan.csv')
dates <- dates %>% janitor::clean_names()



# ================================
#         Data Preparation
# ================================

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

sex <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/3__gwas/1__primary/covariate_file.txt', sep='\t') %>% 
  select(IID, sex)

selected_df <- left_join(selected_df, sex, by='IID')

RPS <- read.csv('/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/retinal_pigmentation_score.csv') %>%
  mutate(IID = str_extract(Name, "(?<=/)[0-9]{7}(?=_)")) %>%
  select(IID, pigmentation)

RPS$IID <- as.numeric(RPS$IID)

selected_df <- left_join(selected_df, RPS, by = 'IID')


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
trimmed_df <- trim_outliers_to_na(selected_df, params)



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
# Load necessary packages
library(tidyverse)
library(showtext)
library(sysfonts)

library(ggplot2)
library(dplyr)
library(tidyr)
library(showtext)
library(sysfonts)

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
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/histograms.pdf",
  plot = p_hist,
  width = 10,
  height = 6,
  device = "pdf"
)


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
    y = "Foveal Parameter"
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
    y = "Foveal Parameter"
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
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/suppfig1.pdf",
  plot = p_ancestry,
  width = 10,
  height = 6,
  device = "pdf"
)



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
    y = "Foveal parameter"
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
  model_rps <- lm(reformulate(c("pigmentation", "spherical_equivalent_left_final", "sex"), trait), data = trimmed_df)
  r2_rps <- summary(model_rps)$adj.r.squared

  # Model 2: Ancestry only
  model_ancestry <- lm(reformulate(c("Ancestry", "spherical_equivalent_left_final", "sex"), trait), data = trimmed_df)
  r2_ancestry <- summary(model_ancestry)$adj.r.squared

  # Model 3: Combined
  model_combined <- lm(reformulate(c("pigmentation", "Ancestry", "spherical_equivalent_left_final", "sex"), trait), data = trimmed_df)
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
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/fig5.pdf",
  plot = p_r2,
  width = 10,
  height = 6,
  device = "pdf"
)



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
  formula <- as.formula(paste("log_mar_final_left ~", paste(c(trait, "spherical_equivalent_left_final", "sex", "Ancestry"), collapse = " + ")))
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
    title = "Foveal Morphology Predicting Visual Acuity",
    subtitle = "Adjusted for sex, SER, and ancestry",
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
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/params_vs_acuity.pdf",
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



