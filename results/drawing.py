import json
import matplotlib.pyplot as plt
import numpy as np
from scipy import sparse
from scipy.sparse import csc_matrix

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

def read_jsons(json_file_names):
    data = []
    for i in json_file_names:
        curve = json_to_data(i)
        data.append(curve)
    return data

def json_to_data(json_file_name):
    assignment = ['deterministic','random','row col']
    
    data = open(json_file_name)
    data = json.load(data)
    res = data['res']
    rate = data['rate']
    error = data['error']
    label = 'm=' + str(data['num_multiplexing']) + ' , ' + assignment[data['assignment type']] 
    return res,rate,error,label

def plot_multiple_data(list_of_res, max_e, min_e, step,n,k,save = False):
    fig = plt.figure()
    ax = fig.add_subplot(1,1,1)
    
    phys_err = np.arange(min_e ,max_e, step)
    for i in list_of_res:
        ax.errorbar(phys_err, i[1], i[2], linewidth= 1, label = i[3])
    
    file_name = "[[" + str(n) + "," +str(k) + "]] LDPC HGP code"
    
    ax.set_yscale('log')
    ax.set_title(file_name, fontsize = 15)
    ax.set_xlabel("photon loss probability", fontsize = 15)
    ax.set_ylabel(r"decoder failure & logical $X$ error prob", fontsize = 15)
    ax.set_xlim(min_e, max_e)
    plt.style.use('tableau-colorblind10')
    plt.legend(loc='lower right', fontsize = 13)
    if save == True:
        file_name = file_name + '.pdf'
        plt.savefig(file_name)
    plt.show()
    return 0

# example 
# files = ['results_HGP_n100_m=1_rate0.001-0.997_time20240313141419_deterministic_assignment.json','results_HGP_n100_m=2_rate0.006-0.983_time20240313142542_random_assignment.json','results_HGP_n100_m=2_rate0.016-0.98_time20240313143125_row_col_assignment.json','results_HGP_n100_m=2_rate0.023-0.982_time20240313141958_deterministic_assignment.json']

# list_of_res = read_jsons(files)
# plot_multiple_data(list_of_res=list_of_res, max_e=0.55, min_e=0.2, step=0.01, n=100,k=4,save=True)