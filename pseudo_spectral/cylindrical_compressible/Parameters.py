#Defining the parameters using 
#in the execution of the progra

import numpy as np 
#def dadosin():

file1       =   'taxamyp' 
file2       =   'baseflowmyp' 
	###1
N	    =   200

omegai      =   0.1
#omegaf      =   3.50
omegaf      =   2.10
#domega	    = 	0.02
domega	    = 	0.2
iomega	    =	np.int64(round((omegaf-omegai)/domega))
print(iomega)
#iomega      =   128
#maximum size r
maxr    =   50
#Jets Parameters
#h      = U_secondary/U_primary
h       =   0.70
#M0=M1, sound velocity adimensionalization, Joncas is with the U1
M0      =   0.6558#0.8
#Compressibility Parameter 
Jh      =   1.0#0.8
#M0      =   0.6558#0.8
#M0     =   0.0001#0.8
#M2     =   h*M1#0.4 
#h      = U_infty/U_primary,r=0
S_T     =   0.50
#T0=T_ref, sound velocity adimensionalization, Joncas is with the U1
T_ref   =   300#0.8
T_r     =   1.0
T_0     =   S_T*T_r 
gamma 	=   1.4
#Gamma  = R2/R1
Gamma   =   2.0 
#Radii Primary jet
R1      =   1.0
###
#theta2  =   0.14 
###
#numero de modos a encontrar if  >nm breakk 
nm	=   0
# Tolerancia para procura dos autovalores
tol	=   0.001
#Max_imaginary part 
min_imag=   0.005
max_imag=   0.8
min_real=   0.005
max_real=   5.2
#
#tol	=   0.005
