#!/bin/bash
#SBATCH -p compute
#SBATCH --time=90:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=4_100000_1_0_101_3_N320_K82
#SBATCH --output=output_4_100000_1_0_101_3_N320_K82.txt
python ../multiplexing_VH_decoder.py 4 100000 1 0 101 3 N320_K82
