#!/bin/bash
#SBATCH -p compute
#SBATCH --time=99:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=8_100000_1_0_101_0_N320_K82
#SBATCH --output=output_8_100000_1_0_101_0_N320_K82.txt
python ../multiplexing_VH_decoder.py 8 100000 1 0 101 0 N320_K82
