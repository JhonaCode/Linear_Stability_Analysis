#Defining the parameters using 
#in the execution of the progra

import numpy as np 
#def dadosin():

#Path of the files 
path        = './'

#Number of collocation points
N	    =   150
    
#Type of Differece Matrix 
diffM       = 'Finite_d'
diffM       = 'Chebyshev'

#Type of Baseflow 
base= 'Morris_ml'
#base= 'Hu_my'
#base= 'Joncas_jet'

#Type of mapping 
mapp='not'#1)
mapp='point'#2)
#Only work with Chebyshev
#######################
#L=0.0
#mapp='tang'#3)
#L=10.0
#mapp='root'#5)
L=20.0
mapp='cyl' #4)


#Maximum domain high, not all mappins needs its.

#Type of analysis
analysis    =   'spatial' 
#analysis    =  'temporal' 


#Experiment label 
expe= 'D_%s_%s_%s'%(diffM,mapp,N) 

file1       =   'compressible_%s'%(expe) 
file2       =   'compressible_%s'%(expe) 

#file1       =   'compressible_%s'%(expe) 
#file2       =   'compressible_%s'%(expe) 

#File to save the eigenvalues 
#eig_spatial 

omegaini    =   0.001
omegaf      =   0.48
domega	    = 	0.01

#Type of analysis
#eig_temporal 

alphaini    =   0.001
alphaf      =   0.98
dalpha	    = 	0.01

#
#Azimuthal number
m	=    0
#
nn	=    0
#Jets Parameters
#h      = U_secondary/U_primary
h       =   0.70
#M0=M1, sound velocity adimensionalization, Joncas is with the U1
#repect to oxigen ref 2
M0      =   0.6558
#Parameter to control the incompressibility, Jh parameter
Jh	=   20.0
#repect to hydrogen
#M0      =   0.1592#0.8
#M_i     =   0.6#0.8
#M_o     =   0.4#0.8
#M_0     =   0.2#0.8
#M2     =   h*M1#0.4 
#h      = U_infty/U_primary,r=0
S_T     =   1.00
#T0=T_ref, sound velocity adimensionalization, Joncas is with the U1
T_ref   =   300#0.8
T_o     =   1.0 
T_0     =   1.0
T_i     =   1.0/S_T
#Gamma  = R2/R1
Gamma   =   2.0 
#Radii Primary jet
R1      =   1.0
###
###
#numero de modos a encontrar if  >nm breakk 
nm	=   0
# Tolerancia para procura dos autovalores
tol	=   0.001
#Max_imaginary part 
min_imag=    -0.001
max_imag=    0.25
min_real=    -0.001
max_real=    1.0
#
#tol	=   0.005

print 'Pseudo-spectral code to Solde the Less-Lin equation' 
print 'type of analysis:%s'%analysis 
print 'Base flow:%s'%base 
print 'Differential Matrix:%s'%diffM 
print 'Mapping:%s'%mapp 
print 'Number of points:%s'%N 
