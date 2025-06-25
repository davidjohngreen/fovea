#!/bin/bash
set -euo pipefail

# Install parallel
apt-get update && apt-get install -y parallel

# ---= Make an output folder for the final file
output_dir="/home/dnanexus/out/MAGMA_OUTPUTS"
mkdir -p "$output_dir"

# ---- Step 0: Download input ----
echo "📥 Downloading inputs..."
dx-download-all-inputs

echo "📁 Listing downloaded inputs:"
ls -lh /home/dnanexus/in/

# ---- Step 1: Find input file and extract trait name ----
INPUT_FILE=$(find /home/dnanexus/in/regenie_input/ -name "*.regenie_nofreq" | head -n 1)
if [[ ! -f "$INPUT_FILE" ]]; then
  echo "❌ ERROR: Regenie input file not found!" >&2
  exit 1
fi

BASENAME=$(basename "$INPUT_FILE")
TRAIT=$(echo "$BASENAME" | sed 's/^merged_//' | sed 's/\.regenie_nofreq$//')
Output_Prefix=${TRAIT}
echo "🔍 Detected trait: $TRAIT"

# ---- Step 2: Download MAGMA reference files and gene sets ----
echo "🔽 Downloading MAGMA reference files..."
wget -q -O ref_part1.zip https://vu.data.surfsara.nl/index.php/s/zkKbNeNOZAhFXZB/download
wget -q -O ref_part2.zip https://vu.data.surfsara.nl/index.php/s/Pj2orwuF2JYyKxq/download
wget -q -O ref_part3.zip https://vu.data.surfsara.nl/index.php/s/VZNByNwpD8qqINe/download

echo "📦 Unzipping reference files..."
unzip -Bq ref_part1.zip
unzip -Bq ref_part2.zip
unzip -Bq ref_part3.zip

echo "🧬 Downloading gene set annotations..."
dx download project-GyXx7GQJgK1Xv9Q1QPbG5jYq:/DaveGreen_temp/magma_gene_sets/c2.cp.biocarta.v2024.1.Hs.entrez.gmt.txt

# ---- Step 3: Prepare GWAS input ----
echo "📂 Preparing GWAS input files..."
mkdir -p ./GWAS_results/

awk -v OFS='\t' 'NR>1 {print $3, $1, $2}' "$INPUT_FILE" > ./GWAS_results/${Output_Prefix}.magma.input.snp.chr.pos.txt
awk -v OFS='\t' 'NR>1 && $13 != "NA" {p = 10^(-$13); if (p > 0 && p <= 1) print $3, p}' "$INPUT_FILE" > ./GWAS_results/${Output_Prefix}.magma.input.p.txt
N=$(awk 'NR==2 {print $8}' "$INPUT_FILE")

echo "🔍 P-value preview (head):"
head ./GWAS_results/${Output_Prefix}.magma.input.p.txt | cat -A

echo "🔍 SNP/CHR/POS preview (head):"
head ./GWAS_results/${Output_Prefix}.magma.input.snp.chr.pos.txt | cat -A

mkdir -p ./temp_annot

# ---- Step 4: Annotate SNPs to genes ----
echo "🧬 Annotating SNPs to genes..."
/home/dnanexus/magma \
  --annotate window=100,20 \
  --snp-loc ./GWAS_results/${Output_Prefix}.magma.input.snp.chr.pos.txt \
  --gene-loc /home/dnanexus/NCBI37.3.gene.loc \
  --out ./GWAS_results/${Output_Prefix}

# Step 5: Chunked MAGMA gene analysis with parallel
N_BATCHES=14
echo "🚀 Running MAGMA gene analysis with $N_BATCHES parallel batches..."

parallel /home/dnanexus/magma \
  --batch {} $N_BATCHES \
  --bfile /home/dnanexus/g1000_eur \
  --gene-annot ./GWAS_results/${Output_Prefix}.genes.annot \
  --gene-model snp-wise=mean \
  --pval ./GWAS_results/${Output_Prefix}.magma.input.p.txt N=$N \
  --out ./temp_annot/${Output_Prefix} \
  ::: $(seq 1 $N_BATCHES)

# Step 6: Merge results
/home/dnanexus/magma \
  --merge ./temp_annot/${Output_Prefix} \
  --out "$output_dir"/${Output_Prefix}


# ---- Step 7: Run gene-set enrichment ----
echo "🧪 Running gene-set enrichment..."
/home/dnanexus/magma \
  --gene-results "$output_dir"/${Output_Prefix}.genes.raw \
  --set-annot ./c2.cp.biocarta.v2024.1.Hs.entrez.gmt.txt \
  --out "$output_dir"/${Output_Prefix}

# ---- Step 8: Upload results ----
echo "📤 Uploading results to DNAnexus..."
# dx upload ./GWAS_results/* --destination /DaveGreen_temp/1__GWAS/1__DISCOVERY/MAGMA_OUTPUT/
dx-upload-all-outputs
