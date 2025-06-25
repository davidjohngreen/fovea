#!/bin/bash
set -e -x -o pipefail

# === Setup ===
work_dir="/home/dnanexus/in"
output_dir="/home/dnanexus/out/admixture_results"
mkdir -p "$output_dir"

# Download all input files
dx-download-all-inputs

# === Install ADMIXTURE ===
# Download and extract
dx download file-J0VX4KjJgK1QjzqXgQZFkj17 -o admixture_linux-1.3.0.tar.gz
tar -xzf admixture_linux-1.3.0.tar.gz

# Make the binary executable
chmod +x ./dist/admixture_linux-1.3.0/admixture

# Add it to PATH
export PATH="$PWD/dist/admixture_linux-1.3.0:$PATH"

# === Locate input files ===
bed_file=$(find "$work_dir" -name "*.bed" | head -n 1)
bim_file=$(find "$work_dir" -name "*.bim" | head -n 1)
fam_file=$(find "$work_dir" -name "*.fam" | head -n 1)
pop_file=$(find "$work_dir" -name "*.pop" | head -n 1)

if [[ -z "$bed_file" || -z "$bim_file" || -z "$fam_file" || -z "$pop_file" ]]; then
  echo "❌ One or more input files missing."
  exit 1
fi

prefix=$(basename "$bed_file" .bed)

# === Copy to working dir ===
cp "$bed_file" "${prefix}.bed"
cp "$bim_file" "${prefix}.bim"
cp "$fam_file" "${prefix}.fam"
cp "$pop_file" "${prefix}.pop"

# === Run supervised ADMIXTURE ===
admixture --supervised -j$(nproc) "${prefix}.bed" 4

# === Move outputs ===
cp "${prefix}.4.Q" "$output_dir/"
cp "${prefix}.4.P" "$output_dir/"

# Upload outputs
dx-upload-all-outputs
