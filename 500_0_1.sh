#!/bin/bash
#SBATCH -p compute
#SBATCH --time=19:00:00
#SBATCH --mem=20G
#SBATCH --job-name=500_0_1  # Job name\
#SBATCH -c 1              # Run all processes on a single node\
#SBATCH --ntasks=1        # Run a single task\
#SBATCH --output=output_500_0_1.txt  # Standard output and error log\

python multiplexing_VH_decoder.py 1 10000 1.00 0.00 100 20 2 4 0