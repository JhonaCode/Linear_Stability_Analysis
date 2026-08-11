################################3#
# Program to generate the 
# instability curves, necesaries # reproduce the Joncas and Perreu paper 
# and the Morris spectral comparison paper.
# Create by: Jhonatan Aguirre 
# Date:08/11/2018
# working:yes  
###################################
# update:16/04/20
# Another way to run the program, 
# scr file contains all the python files 
###################################


# Mathematic functions 
import numpy        as     np

#Plot library 
import matplotlib.pyplot as plt
import matplotlib        as mpl

from   subprocess       import call

########################
#Myfolders
########################

# Contains the parameters to define 
# which kind of run 
from   %_%_%Parameters  import *

#Derivative Matrix 
from   chebPy           import *

#To tranform the chebchev domain to [-\inf,\inf] 
from   mapping          import * 

# Types of Base flow to use 
from   Baseflow         import *


# To plot the base flow 
from   plotbaseflow     import * 

#Contains the rayleight equation 
# for statial and temporal analysis. 
from   type_analysis    import * 

#To define characteristics of the plots
from   plotparameters   import * 

import  spectral as sp 



D,D2,r,deltaz,L  = sp.difference_matrix_working(N,mapp,diffM)

        
#Base flow, src/Baseflow.py
Wb,Tb,Rhob,Y1,Y2  =   Baseflow_morris_mixing(r,N,L,D,D2)


#To Calculated the derivatives 
#of the base flow, necessary in the 
#Less_lin equation. 
#The derivative matrix must be 
#modified to put the boundary conditions and 
# it will not be more useful to calculate the 
#derivatives. This is a important step. 

dWbdr   =   dot(D,Wb)
d2Wbdr  =   dot(D2,Wb)

#Remove the first and the last row, 
#due to the  use of central scheme
Wb=Wb[1:N]
d2Wbdr=d2Wbdr[1:N]

#Use difference matrix
#Remove the first and the last row, 
#becuase the use of central scheme
D = D[1:N,1:N]
D2= D2[1:N,1:N]
#Boundary Conditions
#Applying the boundary condition 
#Square matrix, central difference
#Cheb.py
if diffM !='Chebyshev': 

    diff_Dirichlet(D,D2,deltaz,N) 

# to plot base flow 
#plotbase(N,x,r,Df,Df1,D,D2)   

#Spatial Stability Analysis
if (analysis=='spatial'):                    

    #Steps for stability analysis, define in parameters. 

    nomega  =   int(round((omegaf-omegaini)/domega))
    omega       =   np.arange(0,nomega,1)*domega
    omega[0]    =   omegaini

    alphar1,alphai1= spatial(N,D,D2,Wb,d2Wbdr,omega)

#Temporal Stability Analysis
else: 

    #Steps for stability analysis, define in parameters. 

    nalpha  =   int(round((alphaf-alphaini)/dalpha))
    alpha       =   np.arange(0,nalpha,1)*dalpha
    alpha[0]    =   alphaini

    alphar1,alphai1=temporal(N,D,D2,Wb,d2Wbdr,alpha)

        
plt.show()

