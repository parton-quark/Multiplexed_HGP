# -*- coding: utf-8 -*-
import numpy as np
import random
import datetime
import matplotlib.pyplot as plt
from itertools import chain, combinations
import json 
import sys

from Hypergraph_Product_Code_Construction_v3 import HGP_code, Toric3, compute_adjacency_list
from peeling_cluster_decoder import combined_peeling_and_cluster_decoder
from utilities import generate_erasure_pattern_index_set, generate_random_H_matrix, HGP_code, generate_random_error_index_set_with_erasure_support, index_to_biindex, biindex_to_index

def randomly_assign_qubits_to_photons(num_multiplexing,num_photons):
    num_qubits = num_multiplexing * num_photons
    qubits = [i for i in range(num_qubits)]
    photons = []
    
    for j in range(num_photons):
        qubits_in_photon = []
        for k in range(num_multiplexing):
            picked_qubit = random.choice(qubits)
            qubits_in_photon.append(picked_qubit)
            qubits.remove(picked_qubit)
        photons.append(qubits_in_photon)    
    return photons

def deterministically_assign_qubits_to_photons(num_multiplexing,num_photons):
    num_qubits = num_multiplexing * num_photons
    qubits = [i for i in range(num_qubits)]
    photons = []
    
    for j in range(num_photons):
        qubits_in_photon = []
        for k in range(num_multiplexing):
            picked_qubit = qubits[0]
            qubits_in_photon.append(picked_qubit)
            qubits.remove(picked_qubit)
        photons.append(qubits_in_photon)    
    return photons

# Function to assign random, non-overlapping sets up qubits in the same stabilizer supports to photons.
# Assumes the code is LDPC (in particular, that the rows of input matrix have constant weight).
# (If combining Hx and Hz, assume these have the same weight rows, though they need not in general)
# The number of qubits per photon should be a multiple of this weight.

def photon_assigment_by_stabilizer_support(H,num_stabilizers_per_photon):
    
    # Infer the weight of the rows by looking at the first row in H
    row_weight = np.count_nonzero(H[0])
    
    # Infer matrix size.
    num_rows, num_cols = np.shape(H)
    
    # Infer number of qubits.
    num_qubits = num_cols
    
    # Ideally the stabilizer support sets evenly partition the number of qubits, but this need not be true.
    # In this case, use photons of size equal to stabilizer support weight (or a multiple of this).
    # Any remaining qubits will get assigned to a "remainder photon".
    max_photons, remainder_qubits = divmod(num_qubits,row_weight)
    num_photons, remainder_stabilizers = divmod(max_photons,num_stabilizers_per_photon)

    
    # Initialize a list of photon qubit assingments.
    list_of_qubits_per_photon = []
        
    # Identify qubit index list and shuffle; each row indicates the qubit support of a stabilizer generator.
    assignment_H_list = compute_adjacency_list(H)
    np.random.shuffle(assignment_H_list)
    
    # Ideally, we choose a subset of rows which partition the qubits, but this may not exist in general.
    # We will search for such a partition by first selecting rows at random to assign to photons.
    # Then any overlapping rows are removed from the matrix and we repeat the search with the reduced matrix.
    # Any remaining qubits are grouped together in a one "bad" photon at the end.
    new_assignment_H_list = assignment_H_list
    while( (len(list_of_qubits_per_photon)<max_photons) and (len(new_assignment_H_list)>0) ):
        # Take the first row of the matrix as our first set of qubits in a photon.
        photon_row = assignment_H_list.pop()
        list_of_qubits_per_photon.append(photon_row)
        
        # Build a new assignment H list using only those rows of the original list which do not share qubits.
        new_assignment_H_list = []
        for row in assignment_H_list:
            if (photon_row.intersection(row) == set()):
                new_assignment_H_list.append(row)
        
        # Replace the list with non intersecting rows.
        assignment_H_list = new_assignment_H_list
        
    # Because of how the stabilizers generators are chosen, there may be remaining qubits.
    # (That is, there are no remaining rows of H that do not overlap with previous choices).
    # In this case, we will just do random assignment on the remaining qubits.
    # Track the number true stabilizer-support photons; this number can change and may be of interest.
    num_stab_supp_photons = len(list_of_qubits_per_photon)
    
    # To determine the remaining qubit indices, compare with the stabilizer support photons just created.
    stab_qubit_set = set()
    for temp_qubits_in_photon in list_of_qubits_per_photon:
        stab_qubit_set.update(temp_qubits_in_photon)
    # Bassically, anything which was not already thrown in a stabilizer is a remaining qubit.
    remaining_qubits = set(range(0,num_qubits)).difference(stab_qubit_set)
    # Cast as a list so it can be randomized.
    remaining_qubits_list = list(remaining_qubits)
    random.shuffle(remaining_qubits_list)
    
    # DEBUG
    #print(remaining_qubits_list,num_stab_supp_photons)
    
    for i in range(max_photons - num_stab_supp_photons):
        temp_photon = set()
        for j in range(row_weight):
            temp_photon.add(remaining_qubits_list.pop())
        list_of_qubits_per_photon.append(temp_photon)
        
    # If the requested number of photons is fewer than the maximum, combine some photons.
    # We assume here that the number of qubits in a photon is a multiple of the row weight.
    new_list_of_qubits_per_photon = []
    if (num_photons < max_photons):
        for i in range(num_photons):
            new_photon = set()
            for j in range(num_stabilizers_per_photon):
                new_photon.update(list_of_qubits_per_photon.pop())
            new_list_of_qubits_per_photon.append(new_photon)
    else:
        new_list_of_qubits_per_photon = list_of_qubits_per_photon
    
    # If there are remaining qubits or stabilizers, group these in some remainder photons.
    remainder_photon = set()
    if (remainder_stabilizers > 0):
        # The remaining stabilizer support sets can be grouped into one "remainder photon".
        for leftover_photon in list_of_qubits_per_photon:
            remainder_photon.update(leftover_photon)
    if (remainder_qubits > 0):
        #  Any remaining qubits can also be added to this photon, it will still be large enough.
        remainder_photon.update(set(remaining_qubits_list))
    
    # If the remainder photon is non-empty, add it to the list of photons; it may be smaller than the others.
    if (remainder_photon != set()):
        new_list_of_qubits_per_photon.append(remainder_photon)
        
        
    final_num_qubits_per_photon = row_weight*num_stabilizers_per_photon
    final_num_stab_photons = num_stab_supp_photons/num_stabilizers_per_photon
    final_num_rand_photons = num_photons-(num_stab_supp_photons//num_stabilizers_per_photon)
    
    # DEBUG:
#     print("Number of photons used:",num_photons)
#     print("Stabilizers per photon:",num_stabilizers_per_photon)
#     print("Qubits per photon:",final_num_qubits_per_photon)
#     print("Number of stabilizer-support photons:", final_num_stab_photons)
#     print("Number of random photons:", final_num_rand_photons)
#     print("Ratio of stabilizer photons to random photons:", final_num_stab_photons/final_num_rand_photons)
#     print("Remainder photon size:",len(remainder_photon))

    stab_assignment_list_of_qubits_per_photon = []
    for qubit_set in new_list_of_qubits_per_photon:
        stab_assignment_list_of_qubits_per_photon.append(list(qubit_set))

    return stab_assignment_list_of_qubits_per_photon

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


def agresti_coull_intetrval(pair):
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
        error = agresti_coull_intetrval([num_success, num_fail])
        errors.append(error)
    return rates, errors

def success_DF_nonDFLE_rate_and_error_bar(results, num):
    rates = []
    errors = []
    
    DFrates = []
    DFerrors = []
    
    nonDFLErates = []
    nonDFLEerrors = []
    
    for i in range(num):
        num_success = results[i][0]
        num_fail = results[i][1]
        num_DF = results[i][2]
        num_nonDFLE = results[i][3]
        
        
        rate = num_fail / (num_success + num_fail)
        rates.append(rate)
        DFrate = num_DF / (num_success + num_fail)
        DFrates.append(DFrate)
        nonDFLErate = num_nonDFLE / (num_success + num_fail)
        nonDFLErates.append(nonDFLErate)
        
        error = agresti_coull_intetrval([num_success, num_fail])
        errors.append(error)
        DFerror = agresti_coull_intetrval([num_DF, (num_success + num_fail - num_DF)])
        DFerrors.append(DFerror)
        nonDFLEerror = agresti_coull_intetrval([num_nonDFLE, (num_success + num_fail - num_nonDFLE)])
        nonDFLEerrors.append(nonDFLEerror)
    return rates, errors, DFrates, DFerrors, nonDFLErates, nonDFLEerrors


# A function to partion the the qubits in a code into different photons of equal size.
# Assumes a HGP code and will only assign qubits in different rows/columns together.

def HGP_different_row_and_col_assign_qubits_to_photons(num_multiplexing,num_photons,HGP_code):
    # This is the Sudoku assignment strategy for HGP codes.
    
    num_qubits = num_multiplexing * num_photons
    qubits = [i for i in range(num_qubits)]
    photons = []
    photons_biindex = []
    
    # Define a number to track the number of times qubits could not be selected in a different row/col.
    # In these cases, random assignment of qubit in a photon is made.
    # This is for information purposes only and does not affect the function.
    number_of_random_assignments = 0
    
    if num_qubits != HGP_code.num_qubits:
        raise Exception("The number of qubits in HGP code does not match, cannot properly assign to photons!")
        
    # Assign qubits to photons at random, checking the following conditions hold for EACH qubit in photon:
    # If qubits are from different blocks (horizontal vs vertical), no problem.
    # If qubits are from the same block, check they are in a different row and column.
    # If this condition fails, discard this qubit and try again with another selected at random.
    # If all valid assignments meeting the above conditions have been exhausted, then assign a qubit at random.
    for photon_index in range(num_photons):
        qubits_in_photon = []
        qubits_biindex_in_photon = []
        
        # Create a copy of a the original list of qubits and create a shuffled copy of it.
        # This is way, qubits can be taken at random and also exhaustively iterated through.
        shuffled_qubit_list = qubits.copy()
        random.shuffle(shuffled_qubit_list)
        
        # Set the maximum number of iterations; it's possible no remaining qubits satisfy conditions.
        max_iterations = len(shuffled_qubit_list)
        current_iteration = 0        

        for multiplex_index in range(num_multiplexing):          
            # Loop through qubits at random and try to add these to the photon, if possible meeting conditions.
            # Loop repeats a number of times equal to the multiplexing.
            while ((len(qubits_in_photon) < num_multiplexing) and (current_iteration < max_iterations)):
                current_iteration += 1
                
                # "Pop" the first element of the shuffled qubit list; shuffled ensures qubit added at random.
                # Also, popping removes the element from the list, so it is ingnored in subsequent loops.
                picked_qubit_index = shuffled_qubit_list.pop(0)

                #DEBUG
                #print("1: picked qubit index:"+str(picked_qubit_index))

                # Infer the biindex of this qubit; it depends on information from the HGP block structure.
                # There are two cases, depending on the block
                if picked_qubit_index < HGP_code.num_v_qubits:
                    # Horizontal case
                    picked_qubit_biindex = index_to_biindex(
                        index=picked_qubit_index,num_cols=HGP_code.n2,index_shift=0)
                else:
                    # Vertical case
                    picked_qubit_biindex = index_to_biindex(
                        index=picked_qubit_index,num_cols=HGP_code.r2,index_shift=HGP_code.num_h_qubits)

                # Assign if first qubit.
                if len(qubits_in_photon)==0:
                    qubits_in_photon.append(picked_qubit_index)
                    qubits_biindex_in_photon.append(picked_qubit_biindex)
                    qubits.remove(picked_qubit_index)

                else:
                    conditions_satisfied = True
                    # Iterate through each of the previously selected qubits and check they satisfy conditions.
                    for qubit_in_photon_index in range(len(qubits_in_photon)):
                        prior_qubit_index = qubits_in_photon[qubit_in_photon_index]
                        prior_qubit_biindex = qubits_biindex_in_photon[qubit_in_photon_index]

                        #DEBUG
                        #print(" 2. Comparing qubits "+str(picked_qubit_index)+" and "+str(prior_qubit_index))
                        #print("    Candidates for inclusion in photon "+str(photon_index))

                        # We only need to check the condition for qubits in the same block.
                        if (((picked_qubit_index<HGP_code.num_h_qubits) and 
                             (prior_qubit_index<HGP_code.num_h_qubits)) 
                            or ((picked_qubit_index>=HGP_code.num_h_qubits) and 
                                (prior_qubit_index>=HGP_code.num_h_qubits))):
                            # If qubits are in the same row or column, condition is not satisfied.
                            if ((picked_qubit_biindex[0]==prior_qubit_biindex[0]) or 
                                (picked_qubit_biindex[1]==prior_qubit_biindex[1])):
                                conditions_satisfied = False
                                # DEBUG
                                #print(" ---> These qubits are in same row or column of same block")
                                break

                    # Add the new qubit if it satisfies conditions on all previously selected qubits.
                    if (conditions_satisfied==True):
                        qubits_in_photon.append(picked_qubit_index)
                        qubits_biindex_in_photon.append(picked_qubit_biindex)
                        qubits.remove(picked_qubit_index)
                            
            if ((current_iteration==max_iterations) and (len(qubits_in_photon) < num_multiplexing)):
                # If all remaining qubits have been compared and none satisfy conditions, preceding loop ends.
                # In this case, default to random assignment for the remaining qubits in this photon.
                
                # DEBUG
                #print(" Photon number "+str(photon_index)+" unable to be constructed meeting conditions.")
                #print(" Instead, assigning possible remaining qubits to this photon at random.")
                
                shuffled_qubit_list = qubits.copy()
                random.shuffle(shuffled_qubit_list)
                picked_qubit_index = shuffled_qubit_list.pop()

                if picked_qubit_index < HGP_code.num_v_qubits:
                    # Horizontal case
                    picked_qubit_biindex = index_to_biindex(
                        index=picked_qubit_index,num_cols=HGP_code.n2,index_shift=0)
                else:
                    # Vertical case
                    picked_qubit_biindex = index_to_biindex(
                        index=picked_qubit_index,num_cols=HGP_code.r2,index_shift=HGP_code.num_h_qubits)
                
                qubits_in_photon.append(picked_qubit_index)
                qubits_biindex_in_photon.append(picked_qubit_biindex)
                qubits.remove(picked_qubit_index)
                
        # DEBUG
        #print("  3. Qubits in photon number "+str(photon_index)+" (index list and biindex list)")
        #print(qubits_in_photon)
        #print(qubits_biindex_in_photon)
        
        # Add the qubits in the photon to the list of photons.
        photons.append(qubits_in_photon)
        photons_biindex.append(qubits_biindex_in_photon)
        
    # DEBUG
    # print("Number of times defaulted to random assignment: "+str(number_of_random_assignments))
    
    return photons #, photons_biindex


def run_decoder_with_assignment(
    code,num_multiplexing,assignment_type,num_trials,max_erasure_rate,min_erasure_rate,num_steps):
    num_photons = code.num_qubits//num_multiplexing
    
    if assignment_type == 0:
        # deterministic, it can also be used for the case without multiplexing
        assignment = deterministically_assign_qubits_to_photons(num_multiplexing,num_photons)
    elif assignment_type == 1:
        # random
        assignment = randomly_assign_qubits_to_photons(num_multiplexing,num_photons)
    elif assignment_type == 2:
        # random with row col constraint
        assignment=HGP_different_row_and_col_assign_qubits_to_photons(
        num_multiplexing=num_multiplexing,num_photons=num_photons,HGP_code=code)
    else:
        print("invalid assignment")
    
    res = []
    step_size = (max_erasure_rate-min_erasure_rate)/num_steps
    erasure_rates = [i*step_size+min_erasure_rate for i in range(num_steps)]
    
    for i in erasure_rates:
        res_trials = get_success_prob(
            HGP=code,assignment=assignment,num_trials=num_trials,num_photons=num_photons,erasure_rate=i)
        res.append(res_trials)
        
    rate, error = rate_and_error(res,num_steps)
        
    return res,rate, error, assignment

def run_decoder_with_assignment_with_DFLE(
    code,num_multiplexing,assignment_type,num_trials,max_erasure_rate,min_erasure_rate,num_steps):
    num_photons = code.num_qubits//num_multiplexing
    
    if assignment_type == 0:
        # deterministic, it can also be used for the case without multiplexing
        assignment = deterministically_assign_qubits_to_photons(num_multiplexing,num_photons)
    elif assignment_type == 1:
        # random
        assignment = randomly_assign_qubits_to_photons(num_multiplexing,num_photons)
    elif assignment_type == 2:
        # random with row col constraint
        assignment=HGP_different_row_and_col_assign_qubits_to_photons(
        num_multiplexing=num_multiplexing,num_photons=num_photons,HGP_code=code)
    else:
        print("invalid assignment")
    
    res = []
    step_size = (max_erasure_rate-min_erasure_rate)/num_steps
    erasure_rates = [i*step_size+min_erasure_rate for i in range(num_steps)]
    
    for i in erasure_rates:
        # res_trials = get_success_prob(HGP=code,assignment=assignment,num_trials=num_trials,num_photons=num_photons,erasure_rate=i)
        res_trials = get_DF_and_LE_and_failure_prob(HGP=code,assignment=assignment,num_trials=num_trials,num_photons=num_photons,erasure_rate=i)
        res.append(res_trials)
        
    rates, errors, DFrates, DFerrors, nonDFLErates, nonDFLEerrors = success_DF_nonDFLE_rate_and_error_bar(res, num_steps) 
    return res, rates, errors, assignment, DFrates, DFerrors, nonDFLErates, nonDFLEerrors, erasure_rates

def save_results(assignment_type,assignment,res,rate,error,code,num_multiplexing,max_erasure_rate,min_erasure_rate,num_steps,num_trials):
    
    num_photons = code.num_qubits//num_multiplexing
    dt_now = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    step_size = (max_erasure_rate-min_erasure_rate)/num_steps
    # Data to be written
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
        "time": dt_now,
        "res": res,
        "rate": rate,
        "error": error,
        "timestamp": dt_now,
        "assignment": assignment,
        "HGP.Hx": code.Hx.tolist(),
        "HGP.Hz": code.Hx.tolist(),
        "H1": code.H1.tolist(),
        "H2": code.H2.tolist()
    }

 
    # Serializing json
    json_object = json.dumps(results_dictionary, indent=4)
    
    folder_name = "results/"
    file_name_base = folder_name+"results_HGP_n"+str(code.num_qubits)+"_m="+str(num_multiplexing)+"_rate"+str(rate[0])+"-"+str(rate[-1])+"_time"+str(dt_now)
    
    if assignment_type == 0:
        # deterministic, it can also be used for the case without multiplexing
        file_name = file_name_base+"_deterministic_assignment"
    elif assignment_type == 1:
        # random
        file_name = file_name_base+"_random_assignment"
    elif assignment_type == 2:
        # random with row col constraint
        file_name = file_name_base+"_row_col_assignment"
    else:
        print("invalid assignment")
    # Writing to sample.json
    
    with open(file_name+".json", "w") as outfile:
        outfile.write(json_object)
        
    return results_dictionary

def save_results_with_DFLE(assignment_type,assignment,res,rate,error,DFrates, DFerrors, nonDFLErates, nonDFLEerrors,code,num_multiplexing,max_erasure_rate,min_erasure_rate,num_steps,num_trials,erasure_rates,dt_start, dt_finished):
    
    num_photons = code.num_qubits//num_multiplexing
    dt_now = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    step_size = (max_erasure_rate-min_erasure_rate)/num_steps
    # Data to be written
    
    failure_reasonlist = []
    num_DFlist = []
    num_nonDFLElist = []
    
    for i in res:
        num_DFlist.append(i[2])
        num_nonDFLElist.append(i[3])
        failure_reasonlist.append((i[2],i[3]))
        
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
        "dt_start": dt_start, 
        "dt_finished": dt_finished,
        "dt_now": dt_now,
        "erasure_rates": erasure_rates,
        "results": res,
        "success_rate": rate,
        "error_bar_for_success_rate": error,
        "failure_reason":failure_reasonlist,
        "num_DF":num_DFlist,
        "num_nonDFLE": num_nonDFLElist,
        "DFrate": DFrates,
        "DFerrors": DFerrors,
        "nonDFLErate": nonDFLErates,
        "assignment": assignment,
        "HGP.Hx": code.Hx.tolist(),
        "HGP.Hz": code.Hx.tolist(),
        "H1": code.H1.tolist(),
        "H2": code.H2.tolist()
    }

    # Serializing json
    json_object = json.dumps(results_dictionary, indent=4)
    
    folder_name = "results/"
    file_name_base = folder_name + "results_HGP_n" + str(code.num_qubits) + "_m="+str(num_multiplexing) + "_rate=" + str(min_erasure_rate)+"-"+str(max_erasure_rate) # +"_time"+str(dt_now)
    
    if assignment_type == 0:
        # deterministic, it can also be used for the case without multiplexing
        file_name = file_name_base+"_deterministic_assignment"
    elif assignment_type == 1:
        # random
        file_name = file_name_base+"_random_assignment"
    elif assignment_type == 2:
        # random with row col constraint
        file_name = file_name_base+"_row_col_assignment"
    else:
        print("invalid assignment")
    # Writing to json file
    
    file_name = file_name + ".json"
    
    with open(file_name, "w") as outfile:
        outfile.write(json_object)
        
    return results_dictionary, file_name

def main_without_LE():

    H1_n100 = generate_random_H_matrix(total_bits=8,bit_node_deg=3,check_node_deg=4)
    H2_n100 = generate_random_H_matrix(total_bits=8,bit_node_deg=3,check_node_deg=4)
    HGP_n100 = HGP_code(H1_n100,H2_n100)

    # QM = 1
    num_multiplexing = 1
    assignment_type = 0
    num_trials=10000
    max_erasure_rate=0.55
    min_erasure_rate=0.2
    num_steps=35
    
    dt_now = datetime.datetime.now()
    print(dt_now)
    
    
    res,rate, error, assignment = run_decoder_with_assignment(
        code=HGP_n100,
        num_multiplexing=num_multiplexing,
        assignment_type = assignment_type,
        num_trials=num_trials,
        max_erasure_rate=max_erasure_rate,
        min_erasure_rate=min_erasure_rate,
        num_steps=num_steps)

    save_results(
        assignment_type=assignment_type,
        assignment=assignment,
        res=res,rate=rate,
        error=error,
        code=HGP_n100,
        num_multiplexing=num_multiplexing,
        max_erasure_rate=max_erasure_rate,
        min_erasure_rate=min_erasure_rate,
        num_steps=num_steps,
        num_trials=num_trials
        )

    print('QM=1 finished')
    dt_now = datetime.datetime.now()
    print(dt_now)
    
    

    #num_multiplexing=2
    
    for num_multiplexing in [2,5,10]:

        for i in [0,1,2]:
            assignment_type = i
    
            res,rate, error, assignment = run_decoder_with_assignment(
                code=HGP_n100,
                num_multiplexing=num_multiplexing,
                assignment_type = assignment_type,
                num_trials=num_trials,
                max_erasure_rate=max_erasure_rate,
                min_erasure_rate=min_erasure_rate,
                num_steps=num_steps)
    
            save_results(
                assignment_type=assignment_type,
                assignment=assignment,
                res=res,rate=rate,
                error=error,
                code=HGP_n100,
                num_multiplexing=num_multiplexing,
                max_erasure_rate=max_erasure_rate,
                min_erasure_rate=min_erasure_rate,
                num_steps=num_steps,
                num_trials=num_trials
                )
            print('QM=2 with assignment' + str(i) + 'finished')
            dt_now = datetime.datetime.now()
            print(dt_now)
        return 0

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
    total_bits=int(sys.argv[6])
    bit_node_deg=int(sys.argv[7])
    check_node_deg=int(sys.argv[8])
    assignment_type=int(sys.argv[9])
    print('make HGP code')
    # generate HGP code
    H1 = generate_random_H_matrix(total_bits=total_bits,bit_node_deg=bit_node_deg,check_node_deg=check_node_deg)
    H2 = generate_random_H_matrix(total_bits=total_bits,bit_node_deg=bit_node_deg,check_node_deg=check_node_deg)
    HGP = HGP_code(H1,H2)
    code = HGP
    
    dt_start = datetime.datetime.now()
    print('start simulation')
    print(dt_start)
    dt_start = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    
    res, rates, errors, assignment, DFrates, DFerrors, nonDFLErates, nonDFLEerrors, erasure_rates = run_decoder_with_assignment_with_DFLE(
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
    dt_finished = datetime.datetime.now().strftime("%Y%m%d%H%M%S")
    
    save = save_results_with_DFLE(
        assignment_type=assignment_type,
        assignment=assignment,
        res=res,
        rate=rates,
        error=errors,
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
