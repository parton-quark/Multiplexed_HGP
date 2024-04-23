#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=2_100000_1_0_101_1_N320_K82
#SBATCH --output=output_2_100000_1_0_101_1_N320_K82.txt
python ../multiplexing_VH_decoder.py 2 100000 1 0 101 1 N320_K82
