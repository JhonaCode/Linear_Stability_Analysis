################################3#
# Program to find the 
# Eigenvaleus of the second order 
# Derivative matrix
# See:Less and Lin file:less_lin.tex
# Create by: Jhonatan Aguirre 
    # Date:05/01/2020
    # working: no
    ###################################


import numpy      	as     np
#import matplotlib 	as     pl
import matplotlib.pyplot as plt
import matplotlib        as mpl
from   Parameters  	import *
from   chebPy           import *
from   mapping          import * 
from   Baseflow         import *
from   Boundary         import *
from   subprocess       import call
from   scipy.linalg     import eig,block_diag
from   plotbaseflow     import * 

#To define characteristics of the plots
from   plotparameters   import * 

#N=150
#omegai=0.01
#iomega=100
#domega=0.01


omega	=	np.zeros(iomega)
omega	=	np.arange(0,iomega,1)*domega+omegai

alphai  =   	np.zeros((10000,10000)) 
alphar  =   	np.zeros((10000,10000))  



#Chebichev diferential Matrix
Dc ,z	=   cheb(N) 
# Second derivative
Dc2	=   dot(Dc,Dc)

#Mapping to change the Less-Lin equation to chebichev domain[z=-1:1] 
r,dzdr,Dc,Dc2       =   mapping_squate_root(N,z,Dc,Dc2)

#Finite Diference  diferential Matrix
Df,Df2,deltax,x	=   difference(N) 
# Second derivative

#Use difference matrix
#Remove the first and the last row, 
#becuase the use of central scheme

#Which derivative is use
D =np.zeros((N-1,N-1))
D2=np.zeros((N-1,N-1))

D = Df[1:N,1:N]
D2= Df2[1:N,1:N]

#Applying the boundary condition 
#Cheb.py
diff_Dirichlet(D,D2,deltax,N) 


#Calculing the eigenvalues and the 
#Eigenvector of D2
eigvals, eigvecs = eig(D2); 

eigvals=np.sqrt(-eigvals)
eigvals=sorted(eigvals); 
eigvals=np.real(eigvals);

#Maximum distance
L=10

#Analytic solutiion: Eigenvalues
k=np.arange(1,N)
k=np.pi*k/L

#Relative erro between the analytic and the calculated with D2
erro=100.0*(k[0:10]-eigvals[0:10])/k[0:10]

#To make a file with the results 
g1       =   open('Results%d.dat'%N,'w+')

for i in range(0,10):

    print k[i],eigvals[i],erro[i]
    g1.write("%f\t%f\t%f\n"%(k[i],eigvals[i],erro[i])); 

g1.close()

