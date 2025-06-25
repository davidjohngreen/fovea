#!/bin/bash
set -e -x -o pipefail

# === Setup ===
work_dir="/home/dnanexus/in"
output_dir="/home/dnanexus/out/cojo_results"
mkdir -p "$output_dir"

# Download all input files
dx-download-all-inputs

# === Install GCTA ===
dx download file-Gzq3p5jJgK1Z8QBqJ2g0kJ2F -o gcta.zip
unzip gcta.zip -d gcta
chmod +x gcta/gcta-1.94.4-linux-kernel-3-x86_64/gcta64
gcta_exec="gcta/gcta-1.94.4-linux-kernel-3-x86_64/gcta64"

# === Locate .ma file ===
ma_file=$(find "$work_dir" -name '*.ma' | head -n 1)
base_ma=$(basename "$ma_file" .ma)
echo "Found MA file: $ma_file"

# === Prepare PLINK file symlinks by chr ===
plink_dir="/home/dnanexus/plink_chr_inputs"
mkdir -p "$plink_dir"

for bed in $(find "$work_dir" -name "*.bed" | sort); do
  prefix=$(basename "$bed" .bed)
  chr=$(echo "$prefix" | grep -oE 'CHR[0-9XY]+' | sed 's/CHR//')

  bim=$(find "$work_dir" -name "${prefix}.bim" | head -n 1)
  fam=$(find "$work_dir" -name "${prefix}.fam" | head -n 1)

  if [[ -f "$bim" && -f "$fam" ]]; then
    ln -s "$bed" "${plink_dir}/chr${chr}.bed"
    ln -s "$bim" "${plink_dir}/chr${chr}.bim"
    ln -s "$fam" "${plink_dir}/chr${chr}.fam"
    echo "✅ Linked files for chr${chr}"
  else
    echo "⚠️ Missing .bim or .fam for $prefix, skipping..."
  fi
done

for chr in {1..22} X; do
    prefix="${plink_dir}/chr${chr}"
    if [[ -f "${prefix}.bed" && -f "${prefix}.bim" && -f "${prefix}.fam" ]]; then
        echo "🚀 Running COJO for chr$chr"
        out_prefix="/home/dnanexus/${base_ma}_chr${chr}_cojo"

        if [[ "$chr" == "X" ]]; then
            echo "🛠 Converting CHR column from X to 23 in .bim for GCTA"
            cp "${prefix}.bim" "${prefix}.bim.bak"
            awk '{$1 = ($1 == "X" ? 23 : $1); print}' OFS='\t' "${prefix}.bim.bak" > "${prefix}.bim"
        fi

        # ✅ Check for SNP overlap using temp files to reduce memory and terminal spam
        bim_tmp=$(mktemp)
        ma_tmp=$(mktemp)

        awk '{print $2}' "${prefix}.bim" | sort > "$bim_tmp"
        awk 'NR > 1 {print $1}' "$ma_file" | sort > "$ma_tmp"

        common_snps=$(comm -12 "$bim_tmp" "$ma_tmp" | wc -l)

        rm "$bim_tmp" "$ma_tmp"

        if [[ "$common_snps" -eq 0 ]]; then
            echo "⏭️  Skipping chr$chr — no overlapping SNPs with .ma file"
            continue
        fi

        if [[ "$chr" == "X" ]]; then
            "$gcta_exec" --bfile $prefix --maf 0 \
               --cojo-wind 1500 \
               --cojo-file $ma_file \
               --cojo-slct \
               --out $out_prefix
        else
            "$gcta_exec" --bfile $prefix --chr $chr --maf 0 \
               --cojo-wind 1500 \
               --cojo-file $ma_file \
               --cojo-slct \
               --out $out_prefix
        fi

        cp "${out_prefix}.jma.cojo" "$output_dir/" || true
        cp "${out_prefix}.log" "$output_dir/" || true
    else
        echo "⚠️  Missing PLINK files for chr$chr, skipping."
    fi
done


# Upload outputs
dx-upload-all-outputs
