library(ggplot2)
library(dplyr)
library(forcats)
library(readr)



# Load your cleaned data
df <- read_csv("/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/4__epidemiological_studies/david_cox_results_fixed.csv")

# Define significance colouring and reverse trait order
df <- df %>%
  mutate(
    Significance = ifelse(reject_fdr, "FDR < 0.05", "Not significant"),
    Trait = factor(Trait, levels = rev(c(
      "Pit volume", "Rim radius", "Rim height",
      "Mean slope", "Pit depth", "CFT"
    )))
  )

# Plot
p <- ggplot(df, aes(x = HR, y = Trait, xmin = LCI_95, xmax = UCI_95, color = Significance)) +
  geom_errorbarh(height = 0.3, linewidth = 0.9) +
  geom_point(size = 3) +
  geom_vline(xintercept = 1, linetype = "dotted", color = "gray40", linewidth = 0.7) +
  scale_x_log10(
    breaks = c(0.1, 0.5, 1, 2, 10, 100),
    labels = scales::label_number(accuracy = 0.1)
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = "#B22222", "Not significant" = "black")
  ) +
  facet_wrap(~ Disease, ncol = 1, scales = "free_y") +
  labs(
    x = "Hazard Ratio (log scale)",
    y = NULL,
    color = NULL,
    title = "Associations Between Retinal Traits and Eye Disease Risk"
  ) +
  theme_minimal(base_size = 14) +
  theme(
    legend.position = "top",
    strip.text = element_text(size = 14, face = "bold"),
    axis.text.y = element_text(face = "bold"),
    plot.title = element_text(size = 16, face = "bold"),
    panel.grid.minor = element_blank()
  )

# Save as PDF
# Save to PDF
ggsave(
  "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/cox.pdf",
  plot = p,
  width = 6,
  height = 10,
  device = "pdf"
)





