#!/bin/bash
#SBATCH -p compute
#SBATCH --time=19:00:00
#SBATCH --mem=20G
#SBATCH --job-name=QM_sc11      # Job name\
#SBATCH -c 1                    # Run all processes on a single node\
#SBATCH --ntasks=1                   # Run a single task\
#SBATCH --output=multiplexing_output_test.txt  # Standard output and error log\

<< COMMENTOUT
# Function Inputs
num_multiplexing=int(sys.argv[1])  # number of qubits per photon
num_trials=int(sys.argv[2])
max_erasure_rate=float(sys.argv[3])
min_erasure_rate=float(sys.argv[4])
num_steps=int(sys.argv[5])
assignment_type=int(sys.argv[6])   # 
input_matrices_filename=sys.argv[7]   # The input HGP matrix json file names
COMMENTOUT


python multiplexing_VH_decoder.py 16 1 0 10000 1.00 0.00 100 16 2 4 0