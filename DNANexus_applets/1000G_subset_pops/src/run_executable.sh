#!/bin/bash
set -e -x -o pipefail

# === Download all inputs ===
dx-download-all-inputs

# === Setup ===
output_dir="/home/dnanexus/out/1000g_subsets"
mkdir -p "$output_dir"

# === Download static dependencies ===
dx download file-J0Jg308JgK1VJK9b4vqZYFgj  # bcftools
dx download file-J0Jg328JgK1gjP0Z7Jf0JfQP  # htslib (tabix)
dx download file-J0JgFQjJgK1f3x8Z89XkK9kQ -o panel.txt  # panel file

# Copy panel for use inside output dir
cp panel.txt "$output_dir/"

# === Install bcftools ===
tar -xjf bcftools-1.21.tar.bz2
cd bcftools-1.21
make
export PATH=$PWD:$PATH
cd ..

# === Install htslib (tabix) ===
tar -xjf htslib-1.21.tar.bz2
cd htslib-1.21
make
export PATH=$PWD:$PATH
cd ..

# === Go to output dir ===
cd "$output_dir"
panel_file="panel.txt"

# === Create sample lists for super-populations ===
POPS=("AFR" "EUR" "SAS" "EAS")
for POP in "${POPS[@]}"; do
    awk -v pop="$POP" '$3 == pop {print $1}' "$panel_file" > "${POP}_samples.txt"
    echo "✅ Sample list for $POP created: $(wc -l < ${POP}_samples.txt) samples"
done

# === Read input file ===
input_txt=$(find /home/dnanexus/in/ -name '*.txt' | head -n 1)

# === Download all listed files ===
while read -r file_id; do
  dx download "$file_id"
done < "$input_txt"

# === Build list of downloaded VCFs ===
vcf_ids=()
for chr in {1..22}; do
  vcf_name="ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  if [[ -f "$vcf_name" ]]; then
    vcf_ids+=("$vcf_name")
  else
    echo "⚠️ Missing expected VCF: $vcf_name"
    exit 1
  fi

done

# === Merge and filter all chromosomes at once ===
merged_vcf="allchr_merged.vcf.gz"
bcftools concat -Oz -o "$merged_vcf" "${vcf_ids[@]}"
tabix -p vcf "$merged_vcf"

# === Subset the merged VCF to just 4 populations ===
subset_vcf="subset_4pop_allchr.vcf.gz"



cat AFR_samples.txt EUR_samples.txt SAS_samples.txt EAS_samples.txt > all_samples_to_keep.txt
bcftools view -S all_samples_to_keep.txt -Oz -o "$subset_vcf" "$merged_vcf"


tabix -p vcf "$subset_vcf"

# Clean up before upload
rm -f "$panel_file"

# === Remove downloaded VCFs and TBIs ===
for chr in {1..22}; do
  rm -f "ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz"
  rm -f "ALL.chr${chr}.phase3_shapeit2_mvncall_integrated_v5b.20130502.genotypes.vcf.gz.tbi"
done


# === Upload outputs ===
dx-upload-all-outputs
