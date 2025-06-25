# ---- Parameters ----
trait <- "cft"  # Change this to your trait of interest
base_path <- "/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/7__LDSC/input_files"

# ---- File paths ----
discovery_file <- file.path(base_path, "discovery", paste0("merged_", trait, "_left.ldsc.tsv"))
replication_file <- file.path(base_path, "replication", paste0("merged_", trait, "_left.ldsc.tsv"))

# ---- Read in data ----
disco <- read_tsv(discovery_file, col_types = cols())
repli <- read_tsv(replication_file, col_types = cols())

# ---- Join on SNP ----
merged <- inner_join(disco, repli, by = "SNP", suffix = c("_disc", "_repl"))

sig_merged <- merged %>% filter(P_disc < 5e-8)

cor_test <- cor.test(sig_merged$BETA_disc, sig_merged$BETA_repl)
cor_test

# Extract R² and p-value
r2 <- round(cor_test$estimate^2, 3)
pval <- signif(cor_test$p.value, 3)

print(r2)
print(pval)

# ---- Plot ----
p <- ggplot(sig_merged, aes(x = BETA_disc, y = BETA_repl)) +
  geom_point(alpha = 0.7) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray50") +
  geom_smooth(method = "lm", se = FALSE, color = "firebrick", linewidth = 1) +
  labs(
    x = "Effect size (Discovery)",
    y = "Effect size (Replication)",
  ) +
  theme_minimal(base_size = 14)

print(p)

# ---- Save plot ----
ggsave(
  filename = file.path(paste0("/Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/8__figures/beta/beta_beta_", trait, ".png")),
  plot = p,
  width = 10, height = 5, dpi = 600
)