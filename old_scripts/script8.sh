#!/bin/bash
#SBATCH -p compute
#SBATCH --time=15:00:00
#SBATCH --mem=10G
#SBATCH --job-name=QM_sc8      # Job name\
#SBATCH -c 1                    # Run all processes on a single node\
#SBATCH --ntasks=1                   # Run a single task\
#SBATCH --output=multiplexing_output_test.txt  # Standard output and error log\

python multiplexing_VH_decoder.py 8 1 0 10000 1.00 0.00 100 16 2 4 0