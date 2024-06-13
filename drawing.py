import json
import numpy as np
import matplotlib.pyplot as plt

def agresti_coull_intetrval(pair):
    # num trial n
    success = pair[0]
    fail = pair[1]
    rate = success / (success + fail)
    n = success + fail
    z = 2
    n_tilda = n  + z ** 2
    p_tilda = (1 / n_tilda) * (success + (z**2/2) )
    dif = (z * (np.sqrt(  (p_tilda / n_tilda) * (1 - p_tilda) )))
    return dif

def read_multiple_jsons(json_file_names):
    list_of_res = []
    for i in json_file_names:
        curve = json_to_data(i)
        list_of_res.append(curve)
    return list_of_res

def json_to_data(json_file_name):
    data = open(json_file_name)
    data = json.load(data)
    erasure_rates = data['erasure_rates']
    failure_rates = data['failure_rates']
    failure_rates_error_bars = data['failure_rates_error_bars']  
    decoding_failure_rates = data['DFrate']
    non_df_le_rates = data['nonDFLErate']

    label = 'm=' + str(data['num_multiplexing']) + ' , strategy=' + str(data['assignment type'])
    return label, erasure_rates, failure_rates,failure_rates_error_bars,decoding_failure_rates, non_df_le_rates

def plot_multiple_data(list_of_res, phys_err, save = False, grid = True, log = False, with_le_and_df=False):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    for i in list_of_res:
        ax.errorbar(phys_err, i[2], i[3], linewidth= 2, label = i[0])
        if with_le_and_df == True:
            ax.plot(phys_err, i[4], linewidth = 2, label = i[0] + 'DF') # decoding failure rates
            ax.plot(phys_err, i[5], linewidth = 2, label = i[0] + 'non DF LE') # non decoding failure logical error rates 

    if grid == True:
        ax.grid(which="major", alpha=0.6)
        ax.grid(which="minor", alpha=0.3)
    if log == True:
        ax.set_yscale('log')

    ax.set_xlabel("photon loss probability", fontsize = 15)
    ax.set_ylabel(r"logical $Z$ error probability", fontsize = 15)
    # ax.set_xlim(0, 1)
    plt.style.use('tableau-colorblind10')
    plt.legend(loc='lower right', fontsize = 8)
    if save == True:
        filename = "multiplexed_toric"
        plt.savefig(filename + ".pdf")
    plt.show()
