#!/bin/sh
# Script for execution of deployed MATLAB applications

exe_name=$0
exe_dir=$(dirname "$0")
echo "------------------------------------------"

if [ "x$1" = "x" ]; then
  echo "Usage:"
  echo "  $0 <input_file> <output_file>"
else
  echo "Running MATLAB executable..."

  input_file=$1  # First argument is input file
  output_file=$2 # Second argument is output file

  echo "Input file: $input_file"
  echo "Output file: $output_file"

  ls -lh
  
  # Run MATLAB function with passed arguments
  /my_executable "$input_file" "$output_file"
fi
exit