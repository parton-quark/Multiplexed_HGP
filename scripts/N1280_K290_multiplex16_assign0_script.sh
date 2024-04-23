#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=16_100000_1_0_101_0_N1280_K290
#SBATCH --output=output_16_100000_1_0_101_0_N1280_K290.txt
python ../multiplexing_VH_decoder.py 16 100000 1 0 101 0 N1280_K290
