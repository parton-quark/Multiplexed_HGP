import numpy as np
import random
import datetime
import matplotlib.pyplot as plt
from itertools import chain, combinations
import json 
import sys
from bidict import bidict
import math

from Hypergraph_Product_Code_Construction_v3 import HGP_code, Toric3, compute_adjacency_list
from peeling_cluster_decoder import combined_peeling_and_cluster_decoder
from utilities import generate_erasure_pattern_index_set, generate_random_H_matrix, HGP_code, generate_random_error_index_set_with_erasure_support, index_to_biindex, biindex_to_index
from asssignment_strategies import *

def make_erasure_vec_from_ph(assignment, erasure_pattern):
    erased_qubits = []
    erasure_pattern = list(erasure_pattern)
    
    for index in erasure_pattern:
        erased_qubits  = erased_qubits + assignment[index]
    return erased_qubits

def get_success_prob(HGP, assignment, num_trials, num_photons, erasure_rate):
    num_success = 0
    num_failure = 0
    trials = 0
    while trials < num_trials:
        trials += 1
        try:
            errors_on_photons = generate_erasure_pattern_index_set(num_photons, erasure_rate)
            erasure_on_qubits = make_erasure_vec_from_ph(assignment = assignment, erasure_pattern = errors_on_photons)
            # print('erasure_on_qubits:')
            # print(erasure_on_qubits)
            random_pauli = generate_random_error_index_set_with_erasure_support(HGP.num_qubits, erasure_index_set = set(erasure_on_qubits),error_rate = 0.5)
            # print('random_pauli:')
            # print(random_pauli)
            syndrome = HGP.Hz_syn_index_set_for_X_err(random_pauli)
            result = combined_peeling_and_cluster_decoder(HGP_code=HGP,E_index_set_input=set(erasure_on_qubits),s_index_set_input=syndrome)
            
            errors_after_correction = random_pauli.symmetric_difference(result[0])
            
            if HGP.is_non_trivial_X_logical_error_index_set(errors_after_correction) == False:
                num_success += 1
            else:
                num_failure += 1
        except:
            # print("Decoder failed to decode")
            num_failure += 1
    res_trials = [num_success, num_failure]
    # print('results:' + str(res_trials))
    return res_trials


def get_DF_and_LE_and_failure_prob(HGP, assignment, num_trials, num_photons, erasure_rate):
    # get decoder failure & logical error probability
    num_success = 0
    num_failure = 0
    num_DF = 0
    num_non_DF_LE = 0
    
    trials = 0
    while trials < num_trials:
        trials += 1
        try:
            errors_on_photons = generate_erasure_pattern_index_set(num_photons, erasure_rate)
            erasure_on_qubits = make_erasure_vec_from_ph(assignment = assignment, erasure_pattern = errors_on_photons)
            # print('erasure_on_qubits:')
            # print(erasure_on_qubits)
            random_pauli = generate_random_error_index_set_with_erasure_support(HGP.num_qubits, erasure_index_set = set(erasure_on_qubits),error_rate = 0.5)
            # print('random_pauli:')
            # print(random_pauli)
            syndrome = HGP.Hz_syn_index_set_for_X_err(random_pauli)
            result = combined_peeling_and_cluster_decoder(HGP_code=HGP,E_index_set_input=set(erasure_on_qubits),s_index_set_input=syndrome)
            
            errors_after_correction = random_pauli.symmetric_difference(result[0])
            
            if HGP.is_non_trivial_X_logical_error_index_set(errors_after_correction) == False:
                num_success += 1
            else:
                # logical error
                num_failure += 1
                num_non_DF_LE += 1
        except:
            # decoder failure
            num_failure += 1
            num_DF += 1
    res_trials = [num_success, num_failure, num_DF, num_non_DF_LE]
    return res_trials


def agresti_coull_interval(pair):
    # num trial n
    success = pair[0]
    fail = pair[1]
    rate = success / (success + fail)
    n = success + fail
    z = 2
    n_tilda = n  + z ** 2
    p_tilda = (1 / n_tilda) * (success + (z**2/2) )
    # conf_int_min = p_tilda + (z * (np.sqrt(  (p_tilda / n_tilda) * (1 - p_tilda) )))
    # conf_int_max = p_tilda - (z * (np.sqrt(  (p_tilda / n_tilda) * (1 - p_tilda) )))
    dif = (z * (np.sqrt(  (p_tilda / n_tilda) * (1 - p_tilda) )))
    # return [conf_int_min, rate, conf_int_max]
    # dif = conf_int_max - conf_int_min 
    return dif


def rate_and_error(results, num):
    rates = []
    errors = []
    for i in range(num):
        num_success = results[i][0]
        num_fail = results[i][1]
        rate = num_fail / (num_success + num_fail)
        rates.append(rate)
        error = agresti_coull_interval([num_success, num_fail])
        errors.append(error)
    return rates, errors

def success_DF_nonDFLE_rate_and_error_bar(results, num):
    success_rates = []
    success_errors = []

    failure_rates = []
    failure_errors = []
    
    DFrates = []
    DFerrors = []
    
    nonDFLErates = []
    nonDFLEerrors = []
    
    for i in range(num):
        num_success = results[i][0]
        num_fail = results[i][1]
        num_DF = results[i][2]
        num_nonDFLE = results[i][3]
        
        
        success_rate = num_success / (num_success + num_fail)
        success_rates.append(success_rate)
        failure_rate = num_fail / (num_success + num_fail)
        failure_rates.append(failure_rate)
        DFrate = num_DF / (num_success + num_fail)
        DFrates.append(DFrate)
        nonDFLErate = num_nonDFLE / (num_success + num_fail)
        nonDFLErates.append(nonDFLErate)
        
        success_error = agresti_coull_interval([num_success, num_fail])
        success_errors.append(success_error)
        failure_error = agresti_coull_interval([num_fail, num_success])
        failure_errors.append(failure_error)
        DFerror = agresti_coull_interval([num_DF, (num_success + num_fail - num_DF)])
        DFerrors.append(DFerror)
        nonDFLEerror = agresti_coull_interval([num_nonDFLE, (num_success + num_fail - num_nonDFLE)])
        nonDFLEerrors.append(nonDFLEerror)
    return success_rates, success_errors, failure_rates, failure_errors, DFrates, DFerrors, nonDFLErates, nonDFLEerrors




def run_decoder_with_assignment_with_DFLE(
    code,num_multiplexing,assignment_type,num_trials,max_erasure_rate,min_erasure_rate,num_steps):
    num_photons = code.num_qubits//num_multiplexing
    
    if assignment_type == 0:
        # no multiplexing
        assignment = [[i] for i in range(code.num_qubits)]
    elif assignment_type == 1:
        # random
        assignment = randomly_assign_qubits_to_photons(num_multiplexing,num_photons)
    elif assignment_type == 2:
        # Stabilizer assignment
        assignment = photon_assigment_by_stabilizer_support(
            H=np.concatenate((code.Hz,code.Hx),axis=0),num_multiplexing=num_multiplexing)
    elif assignment_type == 3:
        # random with row col constraint
        # Sudoku assignment
        assignment = HGP_different_row_and_col_assign_qubits_to_photons(
        num_multiplexing=num_multiplexing,num_photons=num_photons,HGP_code=code)
    elif assignment_type == 4:
        # deterministic, it can also be used for the case without multiplexing
        # Row-Col assignment
        assignment = deterministically_assign_qubits_to_photons(num_multiplexing,num_photons)
    elif assignment_type == 5:
        assignment = HGP_diagonal_assign_qubits_to_photons(
            num_multiplexing=num_multiplexing,HGP_code=code)
    else:
        print("invalid assignment")
    
    res = []
    step_size = (max_erasure_rate-min_erasure_rate)/num_steps
    erasure_rates = [i*step_size+min_erasure_rate for i in range(num_steps)]
    
    for i in erasure_rates:
        # res_trials = get_success_prob(HGP=code,assignment=assignment,num_trials=num_trials,num_photons=num_photons,erasure_rate=i)
        res_trials = get_DF_and_LE_and_failure_prob(HGP=code,assignment=assignment,num_trials=num_trials,num_photons=num_photons,erasure_rate=i)
        res.append(res_trials)
        
    success_rates, success_errors, failure_rates, failure_errors, DFrates, DFerrors, nonDFLErates, nonDFLEerrors = success_DF_nonDFLE_rate_and_error_bar(res, num_steps) 
    return res, assignment, success_rates, success_errors, failure_rates, failure_errors, DFrates, DFerrors, nonDFLErates, nonDFLEerrors, erasure_rates


def save_results_with_DFLE(
        assignment_type,
        assignment,
        res,
        success_rates,
        success_errors,
        failure_rates,
        failure_errors,
        DFrates,
        DFerrors,
        nonDFLErates,
        nonDFLEerrors,
        code,
        num_multiplexing,
        max_erasure_rate,
        min_erasure_rate,
        num_steps,num_trials,
        erasure_rates,
        dt_start,
        dt_finished):
    
    num_photons = code.num_qubits//num_multiplexing
    dt_now = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    step_size = (max_erasure_rate-min_erasure_rate)/num_steps

    dt_diff = dt_finished - dt_start
    dt_duration = (divmod(round(dt_diff.total_seconds()),3600)[0], divmod(round(dt_diff.total_seconds()),3600)[1]//60)

    dt_start_str = dt_start.strftime("%Y%m%d%H%M%S")
    dt_finished_str = dt_finished.strftime("%Y%m%d%H%M%S")


    
    # Data to be written
    failure_reasonlist = []
    num_DFlist = []
    num_nonDFLElist = []
    
    for simulation in res:
        num_DFlist.append(simulation[2])   # the number of decoder failures
        num_nonDFLElist.append(simulation[3])   # the number of non-decoder failures giving a logical error
        failure_reasonlist.append((simulation[2],simulation[3]))
        
    results_dictionary = {
        "Code_length": code.num_qubits,
        "Code_dimension": code.dim,
        "assignment type": assignment_type,
        "max_erasure_rate": max_erasure_rate,
        "min_erasure_rate": min_erasure_rate,
        "num_steps": num_steps,
        "stepsize" : step_size, 
        "num_multiplexing": num_multiplexing,
        "num_photons": num_photons,
        "num_trials": num_trials,
        "dt_start": dt_start_str, 
        "dt_finished": dt_finished_str,
        "dt_now": dt_now,
        "dt_duration": dt_duration,
        "erasure_rates": erasure_rates,
        "results": res,
        "success_rates": success_rates,
        "success_rates_error_bars": success_errors,
        "failure_rates": failure_rates,
        "failure_rates_error_bars": failure_errors,
        "failure_reason": failure_reasonlist,
        "num_DF": num_DFlist,
        "num_nonDFLE": num_nonDFLElist,
        "DFrate": DFrates,
        "DFerrors": DFerrors,
        "nonDFLErate": nonDFLErates,
        "nonDFLEerrors": nonDFLEerrors,
        "assignment": assignment,
        "HGP.Hx": code.Hx.tolist(),
        "HGP.Hz": code.Hx.tolist(),
        "H1": code.H1.tolist(),
        "H2": code.H2.tolist()
    }

    # Serializing json
    json_object = json.dumps(results_dictionary, indent=4)
    
    folder_name = "../results/"
    file_name_base = folder_name + "results_HGP_n=" + str(code.num_qubits) + "_m="+str(num_multiplexing) + "_r=" + str(min_erasure_rate)+"-"+str(max_erasure_rate) # +"_time"+str(dt_now)
    file_name_base = file_name_base +"_as=" + str(assignment_type)
    # if assignment_type == 0:
        # deterministic, it can also be used for the case without multiplexing
        # file_name = file_name_base+"_deterministic_assignment"
    # elif assignment_type == 1:
        # random
        # file_name = file_name_base+"_random_assignment"
    # elif assignment_type == 2:
        # random with row col constraint
        # file_name = file_name_base+"_row_col_assignment"
    # else:
        # print("invalid assignment")
    # Writing to json file
    
    file_name = file_name_base + ".json"
    
    with open(file_name, "w") as outfile:
        outfile.write(json_object)
        
    return results_dictionary, file_name


def import_matrices(file_name):
    file_name = file_name + ".json"
    jsonfile = open(file_name, 'r')
    json_load = json.load(jsonfile)
    H1 = np.array(json_load['H1'])
    H2 = np.array(json_load['H2'])
    return [H1, H2]

def main_with_LE():
    # inputs
    
    print('inputs:')
    print(sys.argv)
    # sys.argv[0] is the file name
    
    num_multiplexing=int(sys.argv[1])
    num_trials=int(sys.argv[2])
    max_erasure_rate=float(sys.argv[3])
    min_erasure_rate=float(sys.argv[4])
    num_steps=int(sys.argv[5])
    assignment_type=int(sys.argv[6])
    input_matrices_filename=sys.argv[7]
    #total_bits=int(sys.argv[6])
    #bit_node_deg=int(sys.argv[7])
    #check_node_deg=int(sys.argv[8])
    #assignment_type=int(sys.argv[9])
    #print('make HGP code')
    # generate HGP code
    #H1 = generate_random_H_matrix(total_bits=total_bits,bit_node_deg=bit_node_deg,check_node_deg=check_node_deg)
    #H2 = generate_random_H_matrix(total_bits=total_bits,bit_node_deg=bit_node_deg,check_node_deg=check_node_deg)
    
    input_matrices = import_matrices(input_matrices_filename)
    code = HGP_code(input_matrices[0],input_matrices[1])
    
    dt_start = datetime.datetime.now()
    print('start simulation')
    print(dt_start)
    #dt_start = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    
    res, assignment, success_rates, success_errors, failure_rates, failure_errors, DFrates, DFerrors, nonDFLErates, nonDFLEerrors, erasure_rates = run_decoder_with_assignment_with_DFLE(
        code=code,
        num_multiplexing=num_multiplexing,
        assignment_type = assignment_type,
        num_trials=num_trials,
        max_erasure_rate=max_erasure_rate,
        min_erasure_rate=min_erasure_rate,
        num_steps=num_steps)
    
    print('simulation finished')
    dt_finished = datetime.datetime.now()
    print(dt_finished)
    #dt_finished = datetime.datetime.now().strftime("%Y%m%d%H%M%S")


    
    save = save_results_with_DFLE(
        assignment_type=assignment_type,
        assignment=assignment,
        res=res,
        success_rates=success_rates,
        success_errors=success_errors,
        failure_rates=failure_rates,
        failure_errors=failure_errors,
        DFrates=DFrates,
        DFerrors=DFerrors,
        nonDFLErates=nonDFLErates,
        nonDFLEerrors=nonDFLEerrors,
        code=code,
        num_multiplexing=num_multiplexing,
        max_erasure_rate=max_erasure_rate,
        min_erasure_rate=min_erasure_rate,
        num_steps=num_steps,
        num_trials=num_trials,
        erasure_rates=erasure_rates,
        dt_start=dt_start,
        dt_finished=dt_finished
    )
    print('saved the results to ' + save[1])
    return 0

# main()
main_with_LE()
