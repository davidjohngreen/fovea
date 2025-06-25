### Set a variable for which folder you want to operate in
run="1__DISCOVERY"


### REGENIE STEP 1 SETUP AND FILTERING OF BED ARRAY FILES ###
### CONVERT THE BED FILES INTO PGEN AND PERFORM FILTERING, EXCLUDE NON-INCLUDED INDIVIDUALS ETC ###
array_dir="/Bulk/Genotype Results/Genotype calls"
keep_file="/DaveGreen_temp/1__GWAS/${run}/input_files/samples_to_keep.txt"
output_dir="/DaveGreen_temp/1__GWAS/${run}/pgen_array_data"
base_name="ukb22418_c"

for chr in {1..22} X; do
  dx run swiss-army-knife \
    -iin="${array_dir}/${base_name}${chr}_b0_v2.bed" \
    -iin="${array_dir}/${base_name}${chr}_b0_v2.bim" \
    -iin="${array_dir}/${base_name}${chr}_b0_v2.fam" \
    -iin="${keep_file}" \
    --destination="${output_dir}" \
    --tag="chr${chr}_to_pgen" \
    --name="chr${chr}_to_pgen" \
    --instance-type="mem2_ssd2_x40" \
    --priority="high" \
    --brief \
    --yes \
    -icmd="plink2 \
      --bfile ${base_name}${chr}_b0_v2 \
      --keep samples_to_keep.txt \
      --maf 0.01 \
      --mac 100 \
      --geno 0.1 \
      --hwe 1e-15 \
      --mind 0.1 \
      --make-pgen \
      --out pgen_chr${chr}"
done


# NEXT YOU NEED TO MERGE THE ARRAY PGEN FILES IN ORDER FOR REGENIE STEP 1 TO WORK PROPERLY
# this requires making a merge list for the pgen files, rather than using multiple merge flags
# You need to make the merge list and then add it to the relevant folder in RAP
dx run swiss-army-knife \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/pgen_chr1.pgen" \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/pgen_chr1.pvar" \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/pgen_chr1.psam" \
  # List all other chromosome files as triplets here
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/pgen_merge_list.txt" \
  --destination="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data" \
  --tag="merged_array_pgen" \
  --name="merged_array_pgen" \
  --instance-type="mem2_ssd2_x40" \
  --priority="high" \
  --brief \
  --yes \
  -icmd="plink2 \
    --pfile pgen_chr1 \
    --pmerge-list pgen_merge_list.txt \
    --make-pgen \
    --out ukb_array_filtered_allchr"


# RUN REGENIE STEP 1 ON THE FILTERED ARRAY DATA
dx run swiss-army-knife \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/ukb_array_filtered_allchr.pgen" \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/ukb_array_filtered_allchr.pvar" \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_array_data/ukb_array_filtered_allchr.psam" \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/input_files/covariate_file.txt" \
  -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/input_files/phenotype_file.txt" \
  --destination="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output" \
  --tag="regenie_step1" \
  --name="regenie_step1" \
  --instance-type="mem2_ssd2_x40" \
  --priority="high" \
  --brief \
  --yes \
  -icmd="regenie \
    --step 1 \
    --pgen ukb_array_filtered_allchr \
    --phenoFile phenotype_file.txt \
    --covarFile covariate_file.txt \
    --bsize 1000 \
    --lowmem \
    --lowmem-prefix tmp_regenie_step1 \
    --out ukb_step1"






# ----------------------------------------------------------------------------
# CONVERT THE BGEN FILES INTO PGEN TO ALLOW AUTOMATIC HANDLING OF X CHROMOSOME
# INCORPORATE THE FILTERS INTO THIS STEP TO REDUCE REGENIE RUNTIME
# Also use this stage to extract the info scores (imputation metric) from the .mfi files
# ----------------------------------------------------------------------------
chroms=({1..22} X)


# Define shared paths
imputed_dir="/Bulk/Imputation/UKB imputation from genotype"
output_dir="/DaveGreen_temp/1__GWAS/${run}/pgen_imputed_data"
base_name="ukb22828"
keep_file="/DaveGreen_temp/1__GWAS/${run}/input_files/samples_to_keep.txt"

# Loop through each chromosome
for chr in "${chroms[@]}"; do
  dx run swiss-army-knife \
    -iin="${imputed_dir}/${base_name}_c${chr}_b0_v3.bgen" \
    -iin="${imputed_dir}/${base_name}_c${chr}_b0_v3.bgen.bgi" \
    -iin="${imputed_dir}/${base_name}_c${chr}_b0_v3.sample" \
    -iin="${imputed_dir}/${base_name}_c${chr}_b0_v3.mfi.txt" \
    -iin="${keep_file}" \
    --destination="${output_dir}" \
    --tag="filter_chr${chr}" \
    --name="filter_chr${chr}" \
    --instance-type="mem2_ssd2_x40" \
    --priority="high" \
    --brief \
    --yes \
    -icmd="awk '\$8 > 0.8 {print \$2}' ${base_name}_c${chr}_b0_v3.mfi.txt > snps_info_chr${chr}.txt && \
           plink2 \
             --bgen ${base_name}_c${chr}_b0_v3.bgen ref-first \
             --sample ${base_name}_c${chr}_b0_v3.sample \
             --keep samples_to_keep.txt \
             --maf 0.01 \
             --mac 20 \
             --hwe 1e-15 \
             --geno 0.1 \
             --mind 0.1 \
             --rm-dup force-first \
             --make-pgen \
             --out pgen_chr${chr}"
done






## RUN REGENIE STEP 2 on all chromosomes using the .loco files which are listed in the pred.list file
## remember to go inside the pred.list file and remove the part about dnanexus paths as it will mess up the next stage
for chr in {1..22} X; do
  dx run swiss-army-knife \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_imputed_data/pgen_chr${chr}.pgen" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_imputed_data/pgen_chr${chr}.pvar" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/pgen_imputed_data/pgen_chr${chr}.psam" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_pred.list" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_1.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_2.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_3.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_4.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_5.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_6.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_7.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_8.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_9.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_10.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_11.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_12.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_13.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_14.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_15.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_16.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_17.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_18.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step1_output/ukb_step1_19.loco" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/input_files/covariate_file.txt" \
    -iin="mvn:/DaveGreen_temp/1__GWAS/${run}/input_files/phenotype_file.txt" \
    --destination="mvn:/DaveGreen_temp/1__GWAS/${run}/regenie_step2_output" \
    --tag="regenie_step2_chr${chr}" \
    --name="regenie_chr${chr}" \
    --instance-type mem2_ssd2_x40 \
    --priority high \
    --brief \
    --yes \
    -icmd="regenie \
      --step 2 \
      --pgen pgen_chr${chr} \
      --phenoFile phenotype_file.txt \
      --covarFile covariate_file.txt \
      --pred ukb_step1_pred.list \
      --bsize 400 \
      --out ukb_step2_BT_chr${chr}"
done



# ---------------------------------
# TO FILTER REGENIE OUTPUTS YOU CAN USE TTYD, IT IS EASIER AND SIMPLER. JUST RUN THE REGENIE_FILTER_FINAL.sh script
# Then to upload files to a specific folder, make the folder first and then empty contents of output there recursively
#-----------------------------------

#!/bin/bash
set -euo pipefail

run="1__DISCOVERY"

regenie_dir="/mnt/project/DaveGreen_temp/1__GWAS/${run}/regenie_step2_output"
snps_info_dir="/mnt/project/DaveGreen_temp/1__GWAS/${run}/pgen_imputed_data/imputation_scores"
output_dir="/home/dnanexus/filtered_regenie_files"
mkdir -p "$output_dir"

# Step 1: Find all regenie files
mapfile -t regenie_files < <(find "$regenie_dir" -maxdepth 1 -name "*.regenie")

# Step 2: Group by trait name
declare -A trait_map
for file in "${regenie_files[@]}"; do
    filename=$(basename "$file")
    trait=$(echo "$filename" | sed -E 's/.*_chr[0-9XY]+_(.*)\.regenie/\1/')
    trait_map["$trait"]+="$file "
done

echo "🧠 Trait-to-files mapping:"
for trait in "${!trait_map[@]}"; do
    echo "Trait: $trait"
    for file in ${trait_map[$trait]}; do
        echo "  ↳ $file"
    done
    echo ""
done

# Parallel processing
max_jobs=5
current_jobs=0
pids=()

for trait in "${!trait_map[@]}"; do
    (
        echo "🔬 Processing trait: $trait"
        merged_ma="${output_dir}/merged_${trait}.ma"
        merged_regenie="${output_dir}/merged_${trait}.regenie_filtered"
        merged_ldsc="${output_dir}/merged_${trait}.ldsc.tsv"
        merged_regenie_nofreq="${output_dir}/merged_${trait}.regenie_nofreq"

        echo -e "SNP A1 A2 freq b se p N" > "$merged_ma"
        head -n 1 $(echo ${trait_map[$trait]} | awk '{print $1}') > "$merged_regenie"
        head -n 1 $(echo ${trait_map[$trait]} | awk '{print $1}') > "$merged_regenie_nofreq"
        echo -e "SNP\tA1\tA2\tFRQ\tBETA\tSE\tP\tN" > "$merged_ldsc"

        total_lines=0

        for file in ${trait_map[$trait]}; do
            chr=$(echo "$file" | sed -E 's/.*_chr([0-9XY]+)_.*/\1/')
            snp_list="${snps_info_dir}/snps_info_chr${chr}.txt"

            if [[ ! -f "$snp_list" ]]; then
                echo "⚠️  SNP list not found for chr${chr}, skipping: $file"
                continue
            fi

            lines_in_file=$(($(wc -l < "$file") - 1))
            total_lines=$((total_lines + lines_in_file))

            awk -v ma_out="$merged_ma" \
                -v regen_out="$merged_regenie" \
                -v ldsc_out="$merged_ldsc" \
                -v nofreq_out="$merged_regenie_nofreq" \
                -v snps="$snp_list" '
                BEGIN {
                    OFS = "\t"
                    while ((getline line < snps) > 0) keep[line] = 1
                }
                NR > 1 && keep[$3] && $8 >= 1000 {
                    if (!seen[$3]++) {
                        pval = 10^(-$13)
                        if ($6 > 0.05 && $6 < 0.95) {
                            printf "%s %s %s %.6f %.6f %.6f %.6g %d\n", $3, $5, $4, $6, $10, $11, pval, $8 >> ma_out
                            print >> regen_out
                            printf "%s\t%s\t%s\t%.6f\t%.6f\t%.6f\t%.6g\t%d\n", $3, $5, $4, $6, $10, $11, pval, $8 >> ldsc_out
                        }
                        print >> nofreq_out
                    }
                }
            ' "$file"
        done

        regen_count=$(($(wc -l < "$merged_regenie") - 1))
        regen_nofreq_count=$(($(wc -l < "$merged_regenie_nofreq") - 1))
        ldsc_count=$(($(wc -l < "$merged_ldsc") - 1))
        cojo_count=$(wc -l < "$merged_ma")

        echo "📊 Variant counts for $trait:"
        echo "   • Total SNPs before filtering: $total_lines"
        echo "   • After filtering (regenie_filtered):     $regen_count"
        echo "   • Without freq filter (regenie_nofreq):   $regen_nofreq_count"
        echo "   • After filtering (ldsc.tsv):             $ldsc_count"
        echo "   • COJO SNPs in .ma file:                  $cojo_count"
        echo "✅ Created:"
        echo "   • $merged_ma"
        echo "   • $merged_regenie"
        echo "   • $merged_ldsc"
        echo "   • $merged_regenie_nofreq"
        echo ""
    ) &

    pids+=($!)
    current_jobs=$((current_jobs + 1))

    if [[ "$current_jobs" -ge "$max_jobs" ]]; then
        wait -n
        current_jobs=$((current_jobs - 1))
    fi
done

wait


dx mkdir -p "/DaveGreen_temp/1__GWAS/${run}/regenie_step2_merged_output"
dx upload filtered_regenie_files/ -r  --destination "/DaveGreen_temp/1__GWAS/${run}/regenie_step2_merged_output/"


# ---------------------------------
# Make an LD reference dataset for each of the chromosomes
# using both the same filters as for step 2 AND an info-score filter
#-----------------------------------
#!/bin/bash
imputed_file_dir="/Bulk/Imputation/UKB imputation from genotype"
info_list_dir="/DaveGreen_temp/1__GWAS/${run}/pgen_imputed_data/imputation_scores"
output_dir="/DaveGreen_temp/1__GWAS/${run}/LD_references"
base_name="ukb22828"
keep_file="/DaveGreen_temp/1__GWAS/${run}/input_files/samples_to_keep.txt"
chroms=(5)

for chr in "${chroms[@]}"; do
  dx run swiss-army-knife \
    -iin="${imputed_file_dir}/${base_name}_c${chr}_b0_v3.bgen" \
    -iin="${imputed_file_dir}/${base_name}_c${chr}_b0_v3.bgen.bgi" \
    -iin="${imputed_file_dir}/${base_name}_c${chr}_b0_v3.sample" \
    -iin="${info_list_dir}/snps_info_chr${chr}.txt" \
    -iin="${keep_file}" \
    --destination="${output_dir}" \
    --tag="bed_chr${chr}" \
    --name="bed_chr${chr}" \
    --instance-type="mem2_ssd2_x40" \
    --priority="high" \
    --brief \
    --yes \
    -icmd="plink2 \
      --bgen ${base_name}_c${chr}_b0_v3.bgen ref-first \
      --sample ${base_name}_c${chr}_b0_v3.sample \
      --keep samples_to_keep.txt \
      --extract snps_info_chr${chr}.txt \
      --maf 0.01 \
      --mac 20 \
      --geno 0.1 \
      --mind 0.1 \
      --rm-dup force-first \
      --make-bed \
      --out LD_reference_CHR${chr}"
done


# -------------------------------------------------------------
# Run COJO on each chromosome-trait pair
#---------------------------------------------------------------
# Add all trait names here (no .ma extensions)
traits=(
    pit_volume_left_irn
    cft_left
    max_slope_disk_area_left
    max_slope_disk_area_left_irn
    max_slope_disk_perim_left
    max_slope_disk_perim_left_irn
    max_slope_height_left
    max_slope_left
    max_slope_radius_left
    mean_slope_left
    min_height_left
    pit_area_left
    pit_depth_left
    pit_volume_left
    rim_disk_area_left
    rim_disk_perim_left
    rim_disk_perim_left_irn
    rim_height_left
    rim_radius_left
)


chroms=(X)
ma_dir="/DaveGreen_temp/1__GWAS/${run}/regenie_step2_merged_output"
ld_dir="/DaveGreen_temp/1__GWAS/${run}/LD_references"
output_dir="/DaveGreen_temp/1__GWAS/${run}/COJO_output"

for trait in "${traits[@]}"; do
    ma_file="${ma_dir}/merged_${trait}.ma"
    echo "ma_file = ${ma_file}"

    for chr in "${chroms[@]}"; do
        dx run applet-J0vzYXQJgK1k371Xbp0X7g7j \
            --input in="${ma_file}" \
            --input in="${ld_dir}/LD_reference_CHR${chr}.bed" \
            --input in="${ld_dir}/LD_reference_CHR${chr}.bim" \
            --input in="${ld_dir}/LD_reference_CHR${chr}.fam" \
            --name "cojo_${trait}_chr${chr}" \
            --destination "$output_dir" \
            --tag "cojo" \
            --yes
    done
done


# -----------------------------------
# Merge all COJO outputs and produce vcf files for annotation using the merge_cojo_outputs.sh script
# -----------------------------------
run='1__DISCOVERY'

#!/bin/bash
# Check if input directory was provided
if [ -z "$1" ]; then
    echo "Usage: $0 <input_directory_with_cojo_files>"
    exit 1
fi

input_dir="$1"
input_dir="${input_dir%/}"  # Remove trailing slash if present

# Check directory exists
if [ ! -d "$input_dir" ]; then
    echo "Error: Directory not found: $input_dir"
    exit 1
fi

# Loop over unique trait prefixes in the input directory
for trait in $(ls "$input_dir"/*_cojo.jma.cojo 2>/dev/null | \
               xargs -n 1 basename | \
               sed 's/_chr[0-9]*_cojo\.jma\.cojo//' | \
               sort | uniq); do
    echo "Merging files for trait: $trait"

    output_file="${trait}_cojo_merged.jma.cojo"
    first=1

    for file in "$input_dir"/${trait}_chr*_cojo.jma.cojo; do
        if [[ $first -eq 1 ]]; then
            cat "$file" > "$output_file"
            first=0
        else
            tail -n +2 "$file" >> "$output_file"
        fi
    done

    echo " → Output written to: $output_file"
done


dx mkdir -p "/DaveGreen_temp/1__GWAS/${run}/COJO_output_merged"
dx upload merged* --destination "/DaveGreen_temp/1__GWAS/${run}/COJO_output_merged/"




# -----------------------------------
# Run SNPEff on the vcfs to get annotated variants for each trait
# -----------------------------------
#!/bin/bash

traits=(
    pit_volume_left_irn
    cft_left
    max_slope_disk_area_left
    max_slope_disk_area_left_irn
    max_slope_disk_perim_left
    max_slope_disk_perim_left_irn
    max_slope_height_left
    max_slope_left
    max_slope_radius_left
    mean_slope_left
    min_height_left
    pit_area_left
    pit_depth_left
    pit_volume_left
    rim_disk_area_left
    rim_disk_perim_left
    rim_disk_perim_left_irn
    rim_height_left
    rim_radius_left
)

input_dir="/DaveGreen_temp/1__GWAS/${run}/vep_input"
output_dir="${input_dir}/annotated"
genome="GRCh37.75 release from ENSEMBL"

for trait in "${traits[@]}"; do
    vcf_file="${input_dir}/merged_${trait}_vcf_input.vcf"

    echo "▶ Submitting SnpEff for: $vcf_file"

    dx run app-GyJjkqQ9kgKP61X8Zz8KVzZ7 \
        -ivariants_vcf="${vcf_file}" \
        -igenome="${genome}" \
        --destination="${output_dir}" \
        --priority high \
        --name="snpeff_${trait}" \
        --tag="snpeff" \
        --brief \
        --yes
done



# -----------------------------------
# Run LDSC heritability estimation with docker image
# This is quite light on resources and can be run locally
# So more convenient to downlaod the regenie .ldsc.tsv file produced earlier
# And run it on local machine. Can make plotting easier.
# -----------------------------------


cd /Users/user/Dropbox/work/projects/2024/241125_foveal_parameters/foveal-analysis_rap/7__LDSC
docker run -v "$PWD":/mnt -it zijingliu/ldsc

# mkdir -p /mnt/results

trait="pit_volume_left"

python /ldsc/munge_sumstats.py \
  --sumstats "/mnt/input_files/discovery/merged_${trait}.ldsc.tsv" \
  --N 29255 \
  --a1 A1 \
  --a2 A2 \
  --snp SNP \
  --p P \
  --out "/mnt/lava_files/${trait}"
  # --merge-alleles /mnt/w_hm3.snplist

python /ldsc/ldsc.py \
--h2 "/mnt/results/${trait}.sumstats.gz" \
--ref-ld-chr /mnt/eur_w_ld_chr/ \
--out "/mnt/h2_files/${trait}_h2" \
--w-ld-chr /mnt/eur_w_ld_chr/ 



# -----------------------------------
# Merge the h2 files to make a table
# -----------------------------------
# Change this to your folder path
LOG_DIR="/mnt/h2_files"

# Output file
OUTPUT="/mnt/ldsc_summary_table.tsv"
echo -e "Trait\th2\th2_SE\tLambda_GC\tMean_Chi2\tIntercept\tIntercept_SE\tRatio\tRatio_SE" > "$OUTPUT"

# Loop through each log file
for file in "$LOG_DIR"/*_h2.log; do
    trait=$(basename "$file" | sed 's/_h2.log//')
    h2=$(grep "Total Observed scale h2:" "$file" | awk '{print $5}')
    h2_se=$(grep "Total Observed scale h2:" "$file" | awk '{print $6}' | tr -d '()')
    lambda=$(grep "Lambda GC:" "$file" | awk '{print $3}')
    chi2=$(grep "Mean Chi^2:" "$file" | awk '{print $4}')
    intercept=$(grep "Intercept:" "$file" | awk '{print $2}')
    intercept_se=$(grep "Intercept:" "$file" | awk '{print $3}' | tr -d '()')
    ratio=$(grep "Ratio:" "$file" | awk '{print $2}')
    ratio_se=$(grep "Ratio:" "$file" | awk '{print $3}' | tr -d '()')

    echo -e "$trait\t$h2\t$h2_se\t$lambda\t$chi2\t$intercept\t$intercept_se\t$ratio\t$ratio_se" >> "$OUTPUT"
done



# -----------------------------------
# Run LDSC correlations with docker image
# -----------------------------------
# Define the path to LDSC script
LDSC="python ldsc/ldsc.py"

# Define LD reference path
LDREF="mnt/eur_w_ld_chr/"

# Create output directory for rg files
RG_DIR="mnt/lava/rg_files"
mkdir -p "$RG_DIR"

# Step 1: Gather sumstats files from mnt/results/
FILES=($(ls mnt/lava_files/*.sumstats.gz))
N=${#FILES[@]}

echo "🔍 Found $N sumstats files:"
for f in "${FILES[@]}"; do echo "  • $f"; done
echo ""

# Step 2: Run all pairwise correlations (including self)
for I in "${FILES[@]}"; do
    PHEN1=$(basename "$I" .sumstats.gz)
    for Z in "${FILES[@]}"; do
        PHEN2=$(basename "$Z" .sumstats.gz)
        OUT="${RG_DIR}/rg_${PHEN1}__${PHEN2}"
        echo "📊 Running LDSC genetic correlation: $PHEN1 vs $PHEN2"
        $LDSC \
            --rg "$I","$Z" \
            --ref-ld-chr "$LDREF" \
            --w-ld-chr "$LDREF" \
            --out "$OUT"
    done
done

echo ""
echo "📦 Collating LDSC rg results into $RG_DIR/FOVEA_all.rg"
rm -f "$RG_DIR/FOVEA_all.rg"
first=1

for LOGFILE in "$RG_DIR"/rg_*.log; do
    BASE=$(basename "$LOGFILE" .log)

    awk '
        /Summary of Genetic Correlation Results/ {found=1; next}
        found && NF==0 {exit}
        found {print}
    ' "$LOGFILE" > "$RG_DIR/${BASE}.rg"

    if [[ $first -eq 1 ]]; then
        cat "$RG_DIR/${BASE}.rg" > "$RG_DIR/FOVEA_all.rg"
        first=0
    else
        tail -n +2 "$RG_DIR/${BASE}.rg" >> "$RG_DIR/FOVEA_all.rg"
    fi
done

echo "✅ Done. Results saved to: $RG_DIR/FOVEA_all.rg"



# -----------------------------------
# Run RPS applet with a file containing a list of images
# -----------------------------------
for batch in {1..35}; do
  dx run applet-J0x2G98JgK1X6Y4xk3Zykk8B \
    --input image_list_file="/DaveGreen_temp/2__RPS/batches/batch_${batch}.txt" \
    --name "rps_${batch}" \
    --tag "rps" \
    --priority high \
    --yes
done


# -----------------------------------
# Run MAGMA applet on one trait regenie file
# -----------------------------------
traits=(
    pit_volume_left_irn
    cft_left
    mean_slope_left
    pit_depth_left
    rim_height_left
    rim_radius_left
)

run="1__DISCOVERY"
output_dir="/DaveGreen_temp/1__GWAS/${run}/MAGMA_OUTPUT"

for trait in "${traits[@]}"; do
  echo "🚀 Submitting MAGMA job for: $trait"
  dx run applet-J18g040JgK1yp2jyJqYgX4XP \
  --input regenie_input="/DaveGreen_temp/1__GWAS/1__DISCOVERY/regenie_step2_merged_output/merged_${trait}.regenie_nofreq" \
  --name "MAGMA" \
  --destination $output_dir \
  --tag "magma" \
  --priority high \
  --name "MAGMA_${trait}" \
  --yes
done




# -----------------------------------
# Create a VEP input from the .cojo files
# -----------------------------------
import pandas as pd
from pathlib import Path

# 🔁 Update this to your directory
input_dir = Path("/path/to/your/cojo_files")

for file in input_dir.glob("*.cojo"):
    df = pd.read_csv(file, sep="\t")
    vep_lines = []

    for _, row in df.iterrows():
        snp = str(row["SNP"])
        chrom = str(row["Chr"])
        pos = str(row["bp"])
        ref = row["refA"]
        alt = "."  # placeholder

        if snp.startswith("rs"):
            # Use . for REF and ALT (VEP will resolve rsID to actual alleles)
            vep_lines.append(f"{chrom}\t{pos}\t{snp}\t.\t.\t.\t.\t.")
        else:
            # e.g. "1:183018960_CA_C" → chrom=1, pos=183018960, ref=CA, alt=C
            parts = snp.replace(":", "_").split("_")
            if len(parts) == 4:
                chrom, pos, ref, alt = parts
            vep_lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\t.\t.")

    # Write .vep file
    output_file = file.with_suffix(".vep")
    with open(output_file, "w") as f:
        for line in vep_lines:
            f.write(line + "\n")

    print(f"✅ VEP input written: {output_file.name}")



# -----------------------------------
# Run LAVA applet on one trait regenie file
# -----------------------------------
# List of traits
traits=("cft_left" "pit_depth_left" "pit_volume_left" "mean_slope_left" "rim_height_left" "rim_radius_left" "RPS_left")

# Loop over all trait pairs
for (( i=0; i<${#traits[@]}; i++ )); do
  for (( j=i+1; j<${#traits[@]}; j++ )); do
    trait1=${traits[i]}
    trait2=${traits[j]}

    # Loop over 25 loci chunks
    for chunk in {1..25}; do

      dx run lava_locus_chunk \
        --name "${trait1}_vs_${trait2}_chunk${chunk}" \
        --input trait1="/DaveGreen_temp/lava/${trait1}.sumstats.gz" \
        --input trait2="/DaveGreen_temp/lava/${trait2}.sumstats.gz" \
        --input loci_chunk="/DaveGreen_temp/lava_inputs/loci_chunk_${chunk}.txt" \
        --destination "/DaveGreen_temp/results/${trait1}_vs_${trait2}" \
        --instance-type mem2_ssd1_v2_x16 \
        --yes

    done
  done
done



# 🧪 Test run for a single trait pair and one loci chunk
trait1="cft_left"
trait2="pit_depth_left"
chunk=1

dx run applet-J1GX4P0JgK1px643zPYbV9Jb \
  --name "${trait1}_vs_${trait2}_chunk${chunk}_TEST" \
  --input trait1="/DaveGreen_temp/lava/${trait1}.sumstats.gz" \
  --input trait2="/DaveGreen_temp/lava/${trait2}.sumstats.gz" \
  --input loci_chunk="/DaveGreen_temp/lava/loci_chunk_${chunk}.txt" \
  --destination "/DaveGreen_temp/results/${trait1}_vs_${trait2}" \
  --instance-type mem2_ssd1_v2_x16 \
  --priority high \
  --yes



