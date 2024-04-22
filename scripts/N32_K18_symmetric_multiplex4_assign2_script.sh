#!/bin/bash
#SBATCH -p compute
#SBATCH --time=99:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=4_100000_1_0_101_2_N32_K18_symmetric
#SBATCH --output=output_4_100000_1_0_101_2_N32_K18_symmetric.txt
python ../multiplexing_VH_decoder.py 4 100000 1 0 101 2 N32_K18_symmetric
