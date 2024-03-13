#!/bin/bash
#SBATCH -p short                     # short partition
#SBATCH --mem=10G
#SBATCH --job-name=multiplexing      # Job name\
#SBATCH -c 1                    # Run all processes on a single node\
#SBATCH --ntasks=1                   # Run a single task\
#SBATCH --output=multiplexing_output_test.txt  # Standard output and error log\

python multiplexing_VH_decoder.py
