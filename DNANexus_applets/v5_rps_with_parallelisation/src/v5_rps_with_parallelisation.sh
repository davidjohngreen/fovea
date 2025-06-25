#!/bin/bash
main() {

trap 'echo "💥 Script killed by signal at $(date)"' SIGTERM SIGKILL

echo "📅 Script started at: $(date)"
ulimit -n 65536

# Download inputs
dx-download-all-inputs

# Install Miniconda
mkdir -p ~/miniconda3
wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh -O ~/miniconda3/miniconda.sh
bash ~/miniconda3/miniconda.sh -b -u -p ~/miniconda3 || { echo "❌ Miniconda install failed"; exit 1; }

# Activate Conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda create --name rps python=3.8 -y
conda activate rps

# Locate the input list
image_list_path=$(find /home/dnanexus/in/ -name "*.txt")

# Extract batch name from input filename
batch_name=$(basename "$image_list_path" .txt)  # e.g. batch_7

# Clean workspace for repo
rm -rf retinal-pigmentation-score
ls -l

# Clone the repo
git clone https://github.com/uw-biomedical-ml/retinal-pigmentation-score.git

# Create required directories
mkdir -p /home/dnanexus/retinal-pigmentation-score/src/test_images_UKBB
mkdir -p /home/dnanexus/retinal-pigmentation-score/src/test_outputs

# ---------------------------------------------------------
# Install dependencies with robust logging and error traps
# ---------------------------------------------------------

echo "📦 Step: Installing Python requirements.txt at $(date)"
if pip install -r retinal-pigmentation-score/requirements.txt; then
  echo "✅ requirements.txt installed successfully"
else
  echo "❌ requirements.txt install failed" >&2
  exit 1
fi

echo "📦 Step: Installing Torch stack at $(date)"
if pip install torch==1.9.0+cu111 torchvision==0.10.0+cu111 torchaudio==0.9.0 -f https://download.pytorch.org/whl/torch_stable.html; then
  echo "✅ Torch stack installed"
else
  echo "❌ Torch install failed" >&2
  exit 1
fi

# ---------------------------------------------------------
# Download config (log and trap)
# ---------------------------------------------------------
echo "📄 Downloading config.py at $(date)"
if dx download --overwrite mvn:/Users/TomJ/Retinal_images/config.py -o retinal-pigmentation-score/src/; then
  echo "✅ config.py downloaded"
else
  echo "❌ dx download config.py failed" >&2
  exit 1
fi

# Prepare file list (convert /mnt/project paths to mvn:/ for dx)
cat "$image_list_path" | sed 's|/mnt/project/|mvn:/|' > file_list.txt

cores=10
out_dir="/home/dnanexus/retinal-pigmentation-score/src/test_images_UKBB"

cat file_list.txt | xargs -P "$cores" -I {} dx download "{}" --overwrite -o "$out_dir/"
echo "📈 Files downloaded: $(find "$out_dir" -type f | wc -l)"

# List some files
echo "📂 First few downloaded files:"
find "$out_dir" -type f | head -n 5

# Optional: force CPU mode for PyTorch (if GPU is failing)
# echo ">>> Forcing CPU mode for debugging"
# sed -i 's/torch.device("cuda" if torch.cuda.is_available() else "cpu")/torch.device("cpu")/' retinal-pigmentation-score/src/main.py

# Reactivate Conda env (in case it was lost)
source ~/miniconda3/etc/profile.d/conda.sh
conda activate rps || { echo "❌ Conda activate failed at runtime"; exit 1; }

# Run the RPS model
echo "🚀 Starting RPS model at: $(date)"
python retinal-pigmentation-score/src/main.py
echo "✅ Finished RPS model at: $(date)"

echo "📦 Final output directory listing:"
find /home/dnanexus/out || echo "ℹ️ No output directory present"

# Upload images
echo "RPS script is completed, uploading images"
echo "Contents of /home/dnanexus/retinal-pigmentation-score/src/test_outputs is:"
ls /home/dnanexus/retinal-pigmentation-score/src/test_outputs

# Clean unintended outputs
rm -rf /home/dnanexus/out/*

# Rename the output file to include batch name
csv_orig="retinal-pigmentation-score/src/test_outputs/retinal_pigmentation_score.csv"
csv_renamed="retinal-pigmentation-score/src/test_outputs/retinal_pigmentation_score_${batch_name}.csv"
mv "$csv_orig" "$csv_renamed"

# Upload renamed file
dx upload "$csv_renamed" \
  -o mvn:/DaveGreen_temp/2__RPS/outputs/


}
main
