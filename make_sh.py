# make shell script file for slurm
def make_shellscript(inputs):
    # inputs = [num_multiplexing, num_trials, max_erasure_rate, min_erasure_rate, num_steps, total_bits, bit_node_deg, check_node_deg, assignment_type]
    configurations = str(inputs[0]) + '_' + str(inputs[1]) + '_' + str(inputs[2]) + '_' + str(inputs[3]) + '_' + str(inputs[4]) + '_' + str(inputs[5]) + '_' + str(inputs[6]) + '_' + str(inputs[7]) + '_' + str(inputs[8])
    file_name = configurations + '.sh'
    
    f = open(file_name, 'w')
    datalist = ['#!/bin/bash\n', '#SBATCH -p compute\n', '#SBATCH --time=99:00:00\n', '#SBATCH --mem=20G\n', '#SBATCH -c 1\n', '#SBATCH --ntasks=1\n', ]
    # job_name
    job_name = '#SBATCH --job-name=' + str(configurations) + '\n'
    datalist.append(job_name)
    # output text file
    output_txt = '#SBATCH --output=' + 'output_' + str(configurations) + '.txt\n'
    datalist.append(output_txt)
    # run py
    run = 'python multiplexing_VH_decoder.py' + str(inputs[0]) + ' ' + str(inputs[1]) + ' ' + str(inputs[2]) + ' ' + str(inputs[3]) + ' ' + str(inputs[4]) + ' ' + str(inputs[5]) + ' ' + str(inputs[6]) + ' ' + str(inputs[7]) + ' ' + str(inputs[8]) + '\n'
    datalist.append(run)
    
    f.writelines(datalist)
    f.close()
    return 0    