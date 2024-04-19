# -*- coding: utf-8 -*-
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

def randomly_assign_qubits_to_photons(num_multiplexing,num_photons):
    # Random photon assignment strategy
    
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
    # Row-Col photon assignment strategy
    
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
    # Stabilizer photon assignment strategy
    
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


# Function to take a rectangle of size height x width and construct an ordering on the cells in this rectangle.
# Height and width are assumed to be whole numbers, and the rectangle is visualized as a grid.
# Each cell in this grid enumerated with a coordinate (i,j), starting with the top left corner as (0,0).
#
# The cells in the rectangle can be spanned by one or more diagonal lines with slope -1.
# In this picture with diagonals, we visualize an infinite 2-dimensional plane tiled by rectangles of this size
# All rectangles are identified with each other, so we can restrict to looking at just a single rectangle.
# A diagonal line in the plane thus enters and exits this same rectangle multiple times.
# Eventually, a diagonal line of slope -1 starting at cell (i,j) will connect back to its starting cell.
# The number of distinct parallel diagonal lines needed to pass through every cell in the rectangle is gcd(h,w)
# All of these diagonal lines have the same length, equal to lcm(h,w).
# Starting with the cell at (0,0), we order each cell in the rectangle along these parallel diagonal lines.
# If more than one diagonal line is require to pass through each cell, we index cells sequentially by line.
# After completing the first diagonal line, we shift down to the (1,0) cell to start the next line.
# We do this as many times as necessary to cover all of the cells.
# Note that in general the same diagonal line will enter a rectangle multiple times.
# The gap between these positions is equal to the number of distinct parallel lines, gcd(h,w).
#
# After enumerating all of the cells in the rectangle, this function computes a bijective correspondence.
# This bijection matches the (i,j) coordinate of each cell with an index between 0 and (height x width).
# The function returns a bi-directional dictionary between these coordinates and indices.
# The number of items in the dictionary is exactly the number of cells in the rectangle (height x width).

def compute_rectangle_cells_bijection(height,width):
    
    # Infer the number of parallel diagonal lines and the length of each line.
    num_parallel_lines = math.gcd(height,width)
    line_length = math.lcm(height,width)
    
    # Initialize a dictionary to map indices to coordinates (i,j).
    # This dictionary will be used to create a bidirectional dictionary.
    coord_dict = {}
    
    # Intialize a starting coordinate index of -1.
    # We choose -1 and not 0 here so that the first iteration of the for-loop properly sets this to 0.
    coordinate_index = -1
    
    for line_index in range(num_parallel_lines):
        # Initialize the starting coordinate (line_index,0); this accounts for shifts between lines.
        # We intialize this as a (mutable) list, but will save as a (immutable) tuple to store the coordinate.
        # The row index matches the diagonal line index; the column index is always taken to be 0.
        current_coord = [line_index-1,-1]
        
        for line_position_index in range(line_length):
            # Shift the coordinate index by one to move to the next entry.
            coordinate_index += 1
        
            # Update the current coordinate values.
            current_coord[0] = ((current_coord[0] + 1)%height)
            current_coord[1] = ((current_coord[1] + 1)%width)
            
            coord_dict[coordinate_index] = tuple(current_coord)
            
    # Finally, convert the dictionary into a bi-directional dicitonary before returning.
    rectangle_cell_index_bidict = bidict(coord_dict)
    
    return rectangle_cell_index_bidict


# A function to partion the the qubits in a code into different photons of equal size.
# Assumes a HGP code. Assignment strategy based on using diagonals in the HGP block structure.
# In general, this will ensure that a minimal number of qubits in the same row/col are used per photon.
# Only need to specify the number of qubits per photon (num_multiplexing); the number of photons is inferred.
# Allows for a "remainder" photon with fewer qubits if this number is not a divisor of the number of qubits.

def HGP_diagonal_assign_qubits_to_photons(num_multiplexing,HGP_code):
    # Diagonal photon assignment strategy for HGP codes.
    
    # Infer the number of photons from the number of qubits and the number of multiplexing.
    # If it is not a clean divisor, include an additional remainder qubit.
    num_photons = HGP_code.num_qubits//num_multiplexing
    num_qubits_remainder = HGP_code.num_qubits%num_multiplexing
    
    # Initialize an empty list to store the photons.
    photons = []
    
    # Compute the bi-directional dictionaries relating diagonals and coordinates for each HGP block.
    h_block_diag_bidict = compute_rectangle_cells_bijection(HGP_code.n1,HGP_code.n2)
    v_block_diag_bidict = compute_rectangle_cells_bijection(HGP_code.r1,HGP_code.r2)
    
    # Initialize a list for the original HGP qubit indices ordered by the diagonal lines.
    # These are converted from the biindeces of the cells in the rectangular grid.
    # The bindex of each cell in the grid is determined by the daigonal bi-directional dictionary.
    diag_ordered_qubit_index_list = []
    
    # Construct this list for the horizontal block.
    for diag_index in h_block_diag_bidict:
        cell_biindex = h_block_diag_bidict[diag_index]
        qubit_index = biindex_to_index(
            biindex=cell_biindex,num_cols=HGP_code.n2,index_shift=0)
        diag_ordered_qubit_index_list.append(qubit_index)
        
    # Continue constructing this list using the vertical block.
    for diag_index in v_block_diag_bidict:
        cell_biindex = v_block_diag_bidict[diag_index]
        qubit_index = biindex_to_index(
            biindex=cell_biindex,num_cols=HGP_code.r2,index_shift=HGP_code.num_h_qubits)
        diag_ordered_qubit_index_list.append(qubit_index)
        
    # For each photon, slice the ordered qubit index list into equal sized pieces.
    # Add these slices to each photon in sequence; this is done by popping qubit indices from the list.
    for photon_index in range(num_photons):
        temp_photon = []
        for i in range(num_multiplexing):
            temp_photon.append(diag_ordered_qubit_index_list.pop(0))
        photons.append(temp_photon)
    
    # Any photons remaining in the final list denote the remainder qubit.
    # Add these to the list of photons if a nonzero number of qubits remain.
    if (num_qubits_remainder != 0):
        photons.append(diag_ordered_qubit_index_list.copy())
        
    # This is the final list of qubits assigned per photon that is returned.
    return photons