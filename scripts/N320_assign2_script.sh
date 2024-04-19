#!/bin/bash
#SBATCH -p compute
#SBATCH --time=99:00:00
#SBATCH --mem=20G
#SBATCH -c 1
#SBATCH --ntasks=1
#SBATCH --job-name=8_100_1_0_10_2_N320
#SBATCH --output=output_8_100_1_0_10_2_N320.txt
python ../multiplexing_VH_decoder.py 8 100 1 0 10 2 N320
