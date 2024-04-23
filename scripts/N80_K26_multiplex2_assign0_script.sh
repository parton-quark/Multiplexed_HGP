#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=2_100000_1_0_101_0_N80_K26
#SBATCH --output=output_2_100000_1_0_101_0_N80_K26.txt
python ../multiplexing_VH_decoder.py 2 100000 1 0 101 0 N80_K26
