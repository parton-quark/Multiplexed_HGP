import json 
import sys
import numpy as np
import random
from utilities import generate_random_H_matrix
from Hypergraph_Product_Code_Construction_v3 import HGP_code

def make_and_save_matrices(total_bits, bit_node_deg, check_node_deg, file_name, save = True):
    H1 = generate_random_H_matrix(total_bits,bit_node_deg,check_node_deg)
    H2 = generate_random_H_matrix(total_bits,bit_node_deg,check_node_deg)
    Code = HGP_code(H1,H2)
    print("Input H1 size = "+str(np.shape(Code.H1)))
    print("Input H2 size = "+str(np.shape(Code.H2)))
    print("Code length = "+str(Code.num_qubits))
    print("Code dimension = "+str(Code.dim))
    print("Block 1: "+str(Code.n1)+" by "+str(Code.n2))
    print("Block 2: "+str(Code.r1)+" by "+str(Code.r2))

    if (save):
        json_object = {'H1':H1.tolist(), 'H2':H2.tolist()}
        json_object = json.dumps(json_object, indent=4)
        default_file_name = "N"+str(Code.num_qubits)+"_K"+str(Code.dim)
        folder_name = "input_matrices/"
        if (file_name==""):
            file_name = folder_name + default_file_name + ".json"
        else:
            file_name = folder_name + file_name + ".json"
        with open(file_name, "w") as outfile:
            outfile.write(json_object)
    return None

def make_and_save_symmetric_matrices(total_bits, bit_node_deg, check_node_deg, file_name, save = True):
    H1 = generate_random_H_matrix(total_bits,bit_node_deg,check_node_deg)
    H2 = H1.T

    Code = HGP_code(H1,H2)
    print("Code length = "+str(Code.num_qubits))
    print("Code dimension = "+str(Code.dim))
    print("Block 1: "+str(Code.n1)+" by "+str(Code.n2))
    print("Block 2: "+str(Code.r1)+" by "+str(Code.r2))

    default_file_name = "N"+str(Code.num_qubits)+"_K"+str(Code.dim)

    if (save):
        json_object = {'H1':H1.tolist(), 'H2':H2.tolist()}
        json_object = json.dumps(json_object, indent=4)
        folder_name = "input_matrices/"
        if file_name == "":
            file_name = folder_name + default_file_name + "_symmetric.json"
        else:
            file_name = folder_name + file_name + "_symmetric.json" 
        with open(file_name, "w") as outfile:
            outfile.write(json_object)
    return None


def chosen_symmetric_matrices(m,file_name, save = True):
    bd = 4
    cd = 4

    make_and_save_symmetric_matrices(total_bits=m, bit_node_deg=bd, check_node_deg=cd, file_name=file_name,save=save)
    return None

def chosen_equal_block_matrices(m,file_name, save = True):
    bd = 2
    cd = 4
    tb = 2*m

    make_and_save_matrices(total_bits=tb, bit_node_deg=bd, check_node_deg=cd, file_name=file_name,save=save)
    return None

# make_and_save_matrices(total_bits=16,bit_node_deg=2,check_node_deg=4,file_name="N320")
# make_and_save_symmetric_matrices(8,3,4,"test_matrix")


def save_matrices(H1,H2,file_name):

    json_object = {'H1':H1.tolist(), 'H2':H2.tolist()}
    json_object = json.dumps(json_object, indent=4)
    folder_name = "input_matrices/"
    file_name = folder_name + file_name + ".json"
    
    with open(file_name, "w") as outfile:
        outfile.write(json_object)
    return None