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
    # Strategy (i) Random photon assignment strategy
    
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


# Function to assign random, non-overlapping sets up qubits in the same stabilizer supports to photons.
# Assumes the code is LDPC (in particular, that the rows of input matrix have constant weight).
# (If combining Hx and Hz, assume these have the same weight rows, though they need not in general)
# The number of qubits per photon should be a multiple of this weight, but need not be.

# This function partions the qubits into random sets of non-overlapping stabilizer-supports.
# If no non-overlapping stabilizers remain, the remaining qubits are organized into sets at random.
# The qubits are then ordered by these randomly selected stabilizers.
# The qubits are then divided into photons of size depending on the multiplexing number.
# The goal is for each photon to cover as much of a random stabilizer as possible.
# Depending on the multiplexing number, this could be multiple stabilizers or fractions of a stabilizer.
# This works best when the stabilizer weight is a fraction or multiple of the multiplexing number.
# However, even if this is not true, it groups stabilizers into photons as much as possible.
# A "remainder" photon is also used if the multiplexing number does not divide the number of qubits.

def photon_assigment_by_stabilizer_support(H,num_multiplexing):
    # Strategy (ii) Stabilizer assignment
    
    # Infer the weight of the rows by looking at the first row in H
    row_weight = np.count_nonzero(H[0])
    
    # Infer matrix size.
    num_rows, num_cols = np.shape(H)
    
    # Infer number of qubits.
    num_qubits = num_cols
    
    # Infer the number of photons and remaining qubits.
    num_photons = num_qubits//num_multiplexing
    num_qubits_remainder = num_qubits%num_multiplexing
    
    # Infer the number of stabilizers per photon and any additional qubits that need to be added.
    num_stabilzers_per_photon = num_multiplexing//row_weight
    num_additional_qubits_per_photon = num_multiplexing%row_weight
    
    # Ideally the stabilizer support sets evenly partition the number of qubits, but this need not be true.
    # In this case, use photons of size equal to stabilizer support weight (or a multiple of this).
    # Any remaining qubits will get assigned to a "remainder photon".
    max_stabilizers, stab_remainder_qubits = divmod(num_qubits,row_weight)
    
    # Initialize a list of stabilizer qubit assingments.
    list_of_qubits_per_stabilizer = []
        
    # Identify qubit index list and shuffle; each row indicates the qubit support of a stabilizer generator.
    assignment_H_list = compute_adjacency_list(H)
    np.random.shuffle(assignment_H_list)
    
    # Ideally, we choose a subset of rows which partition the qubits, but this may not exist in general.
    # We will search for such a partition by first selecting rows at random to assign to photons.
    # Then any overlapping rows are removed from the matrix and we repeat the search with the reduced matrix.
    # Any remaining qubits are grouped together in a one "bad" remainder set at the end.
    new_assignment_H_list = assignment_H_list
    while( (len(list_of_qubits_per_stabilizer)<max_stabilizers) and (len(new_assignment_H_list)>0) ):
        # Take the first row of the matrix as our first set of qubits in a photon.
        stabilizer_row = assignment_H_list.pop()
        list_of_qubits_per_stabilizer.append(stabilizer_row)
        
        # Build a new assignment H list using only those rows of the original list which do not share qubits.
        new_assignment_H_list = []
        for row in assignment_H_list:
            if (stabilizer_row.intersection(row) == set()):
                new_assignment_H_list.append(row)
        
        # Replace the list with non intersecting rows.
        assignment_H_list = new_assignment_H_list
        
    # Because of how the stabilizers generators are chosen, there may be remaining qubits.
    # (That is, there are no remaining rows of H that do not overlap with previous choices).
    # In this case, we will just do random assignment on the remaining qubits.
    # Track the number of true non-intersecting stabilizers; this number can change and may be of interest.
    num_non_intersecting_stabilizers = len(list_of_qubits_per_stabilizer)
    
    # To determine the remaining qubit indices, compare with the stabilizer support photons just created.
    stab_qubit_set = set()
    for temp_qubits_in_stabilizer in list_of_qubits_per_stabilizer:
        stab_qubit_set.update(temp_qubits_in_stabilizer)
    # Bassically, anything which was not already thrown in a stabilizer is a remaining qubit.
    remaining_qubits = set(range(0,num_qubits)).difference(stab_qubit_set)
    # Cast as a list so it can be randomized.
    remaining_qubits_list = list(remaining_qubits)
    random.shuffle(remaining_qubits_list)
    
    # DEBUG
    #print(remaining_qubits_list,num_stab_supp_photons)
    
    # Group the remaining qubits at random into sets of size matching the stabilizer weight.
    for i in range(max_stabilizers - num_non_intersecting_stabilizers):
        temp_stab_size_set = set()
        for j in range(row_weight):
            temp_stab_size_set.add(remaining_qubits_list.pop())
        list_of_qubits_per_stabilizer.append(temp_stab_size_set)
        
    # We now have a list of sets of size matching the row weight (size of the stabilizer support).
    # As much as possible, these sets are taken to be from non-overlapping stabilizers.
    # After exhausting randomly chosen non-overlapping stabilizer generators, the remaining sets are random.
    # There may be some non-grouped remaining qubits at the end, but these will be used in a remainder photon.
    # Finally, a list of qubit indices ordered by these chosen stabilizers will be constructed.
    stabilizer_ordered_qubit_list = []
    for stabilizer_index_set in list_of_qubits_per_stabilizer:
        stabilizer_ordered_qubit_list = stabilizer_ordered_qubit_list + list(stabilizer_index_set)
    # Also add the remainder qubits to the end of this list.
    stabilizer_ordered_qubit_list = stabilizer_ordered_qubit_list + list(remaining_qubits_list)
    
    # The photon assignment will be determined by slicing up this list into subsets.
    # If photon sizes matches or is a multiple of stabilizer size, this is perfect stabilizer assignment.
    # Otherewise, the photons contain portions of one or two stabilizers.
    list_of_qubits_per_photon = []
    for photon_index in range(num_photons):
        temp_photon = []
        for i in range(num_multiplexing):
            temp_qubit_index = stabilizer_ordered_qubit_list.pop(0)
            temp_photon.append(temp_qubit_index)
        list_of_qubits_per_photon.append(temp_photon)
    
    if (len(stabilizer_ordered_qubit_list) != 0):
        list_of_qubits_per_photon.append(stabilizer_ordered_qubit_list)
          
    return list_of_qubits_per_photon


# A function to partion the the qubits in a code into different photons of equal size.
# Assumes a HGP code and will only assign qubits in different rows/columns together.

def HGP_different_row_and_col_assign_qubits_to_photons(num_multiplexing,num_photons,HGP_code):
    # Strategy (iii) Sudoku assignment strategy for HGP codes.
    
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


def deterministically_assign_qubits_to_photons(num_multiplexing,num_photons):
    # Strategy (iv) Row-Col photon assignment strategy
    
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
    # Strategy (v) Diagonal photon assignment strategy for HGP codes.
    
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

