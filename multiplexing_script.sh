{\rtf1\ansi\ansicpg1252\cocoartf2759
\cocoatextscaling0\cocoaplatform0{\fonttbl\f0\fswiss\fcharset0 Helvetica;}
{\colortbl;\red255\green255\blue255;}
{\*\expandedcolortbl;;}
\paperw11900\paperh16840\margl1440\margr1440\vieww11520\viewh8400\viewkind0
\pard\tx720\tx1440\tx2160\tx2880\tx3600\tx4320\tx5040\tx5760\tx6480\tx7200\tx7920\tx8640\pardirnatural\partightenfactor0

\f0\fs24 \cf0 #!/bin/bash\
#SBATCH --job-name=multiplexing      # Job name\
#SBATCH --nodes=1                    # Run all processes on a single node\
#SBATCH --ntasks=1                   # Run a single task\
#SBATCH --cpus-per-task=1            # Number of CPU cores per task\
#SBATCH --time=100:00:00              # Time limit hrs:min:sec\
#SBATCH --output=multiplexing_output_test.txt  # Standard output and error log\
\
python multiplexing_VH_decoder.py}