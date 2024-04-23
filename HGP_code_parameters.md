total_bits,bit_node_deg,check_node_deg for the codes used for the paper.

# Assymmetry codes

[[320,82]]
16 2 4

[[405, 101]]
18 2 4

[[500,122]]
20 2 4

[[1000,442]]
30 2 6

total_bits/2 should be good for independent row col strategy

# Symmmetric codes (H1 = H2)
#
# LDPC parameters
# total_bits = n            = 2m   
# total_checks = r          = m    
# check_deg/bit_deg = n/r   = 2    
# bit_deg = ?               = 2           
# check_deg = ?             = 4     
#
# Classical Matrices
# H1 = [r by n] =   [ m by 2m ]
# H2 = H1
#
# HGP code
# HGP(H1,H1.T) Hz and Hx = [ nr by 2nr ]  =   [ 2m^2 by 4m^2 ]
#
# Block Size
# Upper Left = n by n = 2m by 2m
# Lower Right = r by r = m by m


# Equal Block Sizes
- total_bits = m
- total_checks = m
- bit_deg = check_deg = 4
## Classical Matrices:
- H1 = [m by m]
- H2 = [m by m]
 
## parameters
2 2 2 
4 4 4 
8 4 4 
16 4 4 

# HGP code
# HGP(H1,H2) Hz and Hx = [ m^2 by 2m^2 ]



