#!/bin/bash
set -euo pipefail

echo "📁 Listing downloaded inputs:"
ls -lh /home/dnanexus/in/

# Define output directory (must match output name in dxapp.json later)
output_dir="/home/dnanexus/out/LAVA_RESULTS"
mkdir -p "$output_dir"

echo "📁 Creating working directory..."
mkdir -p lava_env && cd lava_env

# 🔍 Identify and assign trait files
trait1_path=$(find /home/dnanexus/in/trait1/ -name "*.sumstats.gz" | head -n 1)
trait2_path=$(find /home/dnanexus/in/trait2/ -name "*.sumstats.gz" | head -n 1)

if [[ -z "$trait1_path" || -z "$trait2_path" ]]; then
  echo "❌ ERROR: One or both trait files not found!" >&2
  exit 1
fi

trait1_name=$(basename "$trait1_path" .sumstats.gz)
trait2_name=$(basename "$trait2_path" .sumstats.gz)

# 🔍 Identify loci chunk
loci_path=$(find /home/dnanexus/in/loci_chunk/ -name "loci_chunk_*.txt" | head -n 1)
if [[ -z "$loci_path" ]]; then
  echo "❌ ERROR: Loci chunk file not found!" >&2
  exit 1
fi

echo "✅ Trait 1: $trait1_name"
echo "✅ Trait 2: $trait2_name"
echo "✅ Loci file: $loci_path"



# ---------------------------------------------
# 🔽 Step 1: Download and unzip 1000G reference files
# ---------------------------------------------
echo "🔽 Downloading 1000G reference files (LAVA-compatible)..."
wget -q -O ref_part1.zip https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download
wget -q -O ref_part2.zip https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download
wget -q -O ref_part3.zip https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download

echo "📦 Unzipping reference files..."
unzip -Bq ref_part1.zip
unzip -Bq ref_part2.zip
unzip -Bq ref_part3.zip

dx download file-J1G0fp0JgK1x0pXJqgx45Bp8
dx download file-J1G09x8JgK1p22gq9vqvvzVX


# ---------------------------------------------
# 🧰 Step 2: Install LAVA and dependencies
# ---------------------------------------------
echo "🧰 Installing LAVA and R packages..."
Rscript - <<'EOF'
packages <- c("data.table", "matrixsampling", "R.utils")
install.packages(packages, repos = "https://cloud.r-project.org")
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install("snpStats", ask = FALSE)
if (!require("LAVA", character.only = TRUE)) {
  if (!dir.exists("LAVA")) system("git clone https://github.com/josefin-werme/LAVA.git")
  system("R CMD INSTALL ./LAVA")
}
EOF


# Export environmental variables to use inside the R script
export TRAIT1="$trait1_name"
export TRAIT2="$trait2_name"
export LOCI_CHUNK="$loci_path"
export OUTPUT_DIR="$output_dir"

# ---------------------------------------------
# 🧠 Step 3: Run lava
# ---------------------------------------------
echo "🧠 Running LAVA analysis in R..."

Rscript - <<EOF
trait1 <- Sys.getenv("TRAIT1")
trait2 <- Sys.getenv("TRAIT2")
traits <- c(trait1, trait2)

cat("🔬 Using traits: ", traits[1], "and", traits[2], "\n")

# Subset input info
input_info <- read.table("lava_input_info.txt", header=TRUE, sep="\t", stringsAsFactors=FALSE)
subset_info <- input_info[input_info\$phenotype %in% traits, ]
write.table(subset_info, "lava_input_info_subset.txt", sep="\t", quote=FALSE, row.names=FALSE)

# Subset sample overlap matrix
overlap_full <- read.table("lava_sample_overlap.csv", header=TRUE, row.names=1, check.names=FALSE)
overlap_subset <- overlap_full[traits, traits]
write.table(overlap_subset, "lava_sample_overlap_subset.txt", sep="\t", quote=FALSE, col.names=FALSE)

library(LAVA)

# Step 1: Process input
cat("⚙️ Processing input...\n")
input <- process.input(
  input.info.file = "lava_input_info_subset.txt",
  sample.overlap.file = "lava_sample_overlap_subset.txt",
  ref.prefix = "g1000_eur",
  phenos = traits
)

# Step 2: Load LD blocks
cat("📦 Reading LD blocks...\n")
loci <- read.loci(Sys.getenv("LOCI_CHUNK"))

# Step 3: Loop over loci
univ_out <- list()
bivar_out <- list()

for (i in 1:nrow(loci)) {
  cat(sprintf("🔄 Locus %d / %d...\n", i, nrow(loci)))
  locus <- process.locus(loci[i, ], input)
  if (!is.null(locus)) {
    res <- run.univ.bivar(locus)
    info <- data.frame(locus = locus\$id, chr = locus\$chr, start = locus\$start, stop = locus\$stop, n.snps = locus\$n.snps, n.pcs = locus\$K)
    univ_out[[i]] <- cbind(info, res\$univ)
    if (!is.null(res\$bivar)) {
      bivar_out[[i]] <- cbind(info, res\$bivar)
    }
  }
}

# Step 4: Save results
dir.create(Sys.getenv("OUTPUT_DIR"), showWarnings=FALSE)
write.table(do.call(rbind, univ_out), file.path(Sys.getenv("OUTPUT_DIR"), "univ.txt"), sep="\t", quote=FALSE, row.names=FALSE)
write.table(do.call(rbind, bivar_out), file.path(Sys.getenv("OUTPUT_DIR"), "bivar.txt"), sep="\t", quote=FALSE, row.names=FALSE)

cat("✅ Done.\n")
EOF


