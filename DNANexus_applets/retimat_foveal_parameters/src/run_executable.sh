#!/bin/bash
set -e -x -o pipefail

# Define runtime and paths
runtime_dir="/home/dnanexus/matlab_runtime/installed_runtime"
installer_path="/home/dnanexus/matlab_runtime/install"
installer_input="/home/dnanexus/matlab_runtime/installer_input.txt"
# executable_dir="/home/dnanexus/in/files"
output_dir="/home/dnanexus/out/OCT_OUTPUTS_${DX_JOB_ID}"
final_output="${output_dir}/final_output.txt"  # Define final concatenated output file

# Download all inputs
dx-download-all-inputs

# Create necessary directories
mkdir -p "$runtime_dir" "$output_dir"
mkdir downloads

# Debugging: Show current structure
pwd
ls -l /
ls -l /home/dnanexus

# Download MATLAB Runtime
echo "Downloading MATLAB Runtime..."
dx download file-GxG4J2jJyQ6XB1BgV0PPyGfZ -o /home/dnanexus/matlab_runtime/matlab_runtime.zip


# Unzip MATLAB Runtime
unzip -q /home/dnanexus/matlab_runtime/matlab_runtime.zip -d /home/dnanexus/matlab_runtime

# Create installer input file
cat <<EOL > "$installer_input"
agreeToLicense=yes
destinationFolder=$runtime_dir
outputFile=/home/dnanexus/matlab_runtime/install_log.txt
EOL

# Install MATLAB Runtime
if [ -f "$installer_path" ]; then
  bash "$installer_path" -inputfile "$installer_input" | tee /home/dnanexus/matlab_runtime/install_log.txt
else
  echo "MATLAB Runtime installer not found."
  exit 1
fi

# Verify Runtime Installation
if [ ! -d "$runtime_dir" ]; then
  echo "MATLAB Runtime installation failed."
  exit 1
fi

ls $runtime_dir


# Set up MATLAB Runtime environment variables once
export LD_LIBRARY_PATH=.:${runtime_dir}/R2024b/runtime/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${runtime_dir}/R2024b/bin/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${runtime_dir}/R2024b/sys/os/glnxa64
export LD_LIBRARY_PATH=${LD_LIBRARY_PATH}:${runtime_dir}/R2024b/sys/opengl/lib/glnxa64

echo "MATLAB Runtime environment set up."
echo "LD_LIBRARY_PATH is ${LD_LIBRARY_PATH}"


# Get the path to the file list containing the relative paths
file_list_path=$(find /home/dnanexus/in/ -name "*.txt")

# Remove the relative path and replace with project
cat $file_list_path | sed 's|/mnt/project/|mvn:/|' > file_list.txt

# Show head for debugging
# head -n 10 file_list.txt

# Perform parallel downloads of files using tom's xargs method
cores=16  # Set max parallel downloads (DNANexus throttles above ~20)
cat file_list.txt | \
xargs -P $cores -I {} dx download "{}" --overwrite -o ./downloads/

# Confirm location of the downloads
# ls downloads

# Create file with list of samples in the downloads folder


# Create file with list of samples, ensuring full paths are included
ls downloads/*.fda > sample_file.txt


# Ensure scripts are executable
chmod +x /my_executable

# Run the MATLAB compiled executable ONCE and process all files in one go
echo "Running MATLAB executable for all samples..."

# Pass the sample file list to the executable
if /my_executable sample_file.txt "$output_dir"; then
  echo "Batch processing complete."
else
  echo "Error: Some files failed. Check logs."
fi

# Find all text files to concatenate
txt_files=($(find "$output_dir" -type f -name "*.csv" ! -name "final_output.txt" ! -name "failed_samples.log"))

# Check if there are any files to process
if [ ${#txt_files[@]} -eq 0 ]; then
    echo "No text files found for concatenation."
    exit 1
fi

# Extract the header from the first file and add 'filename' as the first column
first_file="${txt_files[0]}"
header=$(head -n 1 "$first_file")
echo -e "filename\t$header" > "$final_output"

# Process each file: Add the filename as a new column
for file in "${txt_files[@]}"; do
    # Extract filename without path
    filename=$(basename "$file")
    
    # Skip the first line (header) for all files except the first one
    tail -n +2 "$file" | awk -v fname="$filename" '{print fname "\t" $0}' >> "$final_output"
done

# Delete all files inside the output folder except for the concatenated one and the failed samples log file
find "$output_dir" -type f -name "*.csv" ! -name "final_output.txt" ! -name "failed_samples.log" -delete

# Upload results
dx-upload-all-outputs

echo "All processing complete."
