#matplotlib inline
#config InlineBackend.figure_format='svg'
from   source_lst.chebPy        import *
from   source_lst.mapping       import * 
from   source_lst.Baseflow      import *

from   scipy.linalg  import block_diag
import numpy 	         as np 
import matplotlib.pyplot as pl
#import diffusion         as df
import sys 
import os 
import decimal


"""
In: integer NN     : number of collocation points
In: character diffM: matrix of differentiation 
In: character mapp : mapping to used

Out : D     : Firt derivative Matriz 
Out : D2    : Second derivative Matriz 
Out : r     : coordinate  
Out : L     : size of the domain -L to L  

"""

def difference_matrix_working(N,mapp,diffM,L,**kwargs):

    #Mapping to change the Less-Lin equation to chebichev domain[z=-1:1] 
    #Type of difference Matrix
    #See Parameters.py

    D       = None
    D2      = None
    r       = None
    deltaz  = None

    if diffM=='Chebyshev': 
    
        #Chebichev diferential Matrix
        D ,z    =   cheb(N) 
        # Second derivative
        D2  =   dot(D,D)
    
        #To displace the mixing layer
        L  = 10.0
    
        if mapp=='cylyn': 
    
            D,D2,r = mapping_cylindrical(N,z,D,D2,L)
    
            #Warning 
            #it's not working 
            D2  =   dot(D,D)
    
        if mapp=='sqrtm': 
    
            D,D2,r,dzdr =   mapping_squate_root(N,z,D,D2)
    
        if mapp=="tanmp":
    
            D,D2,r      =  mappingtan_cartesian(N,z,D,D2) 
    
    
    if diffM=='Fnitediff': 
    
        #To work with the mapping[0:1]computational domain 
    
        #The domain begins in -1 to 1, to work as chebyshev work  
        # and use the same map
        Min= 1.0
        Max=-1.0
    
        #Finite Diference  diferential Matrix
        D,D2,deltaz,z   =   difference(N,Max,Min) 
    
        #To work with the mapping[Min:Max]computational domain 
        ##Heigh of the phisical domain 
        #To displace the mixing layer
    
        #Mapping cylindrical work to chevichev as finite difference matrix 
        if mapp=='cylyn': 
    
            ##Heigh of the phisical domain 
            D,D2,r = mapping_cylindrical(N,z,D,D2,L)
    
            #Warning 
            #it's not working 
            D2  =   dot(D,D)
    
        if mapp=='sqrtm': 
        
            #To displace the mixing layer
    
            D,D2,r,dzdr       =   mapping_squate_root(N,z,D,D2)
    
        if mapp=="tanmp":
    
            D,D2,r             =  mappingtan_cartesian(N,z,D,D2) 
    
        if  mapp=='point': 
    
            #To work with the mapping[0:1]computational domain 
            Min=0.0
            Max=1.0
    
            D,D2,deltaz,z   =   difference(N,Max,Min) 
    
            D,D2,r = mapping_point(N,z,D,D2,L)
    
    
        if  mapp=='__not': 
    
            Min=0.0
            Max=L
    
            #Define the difference matrix domain
            D,D2,deltaz,r   =   difference(N,Max,Min) 
            #Displacement mixing layer 
            #L  = (Max-Min)/2.0

    return D,D2,r,deltaz,L 


def difference_matrix(NN):

	# number of collocation nodes=N+1, because the cheb func


    #Type of difference Matrix
    #See Parameters.py
    if diffM=='Chebyshev': 

        #Chebichev,Collocation points  diferential Matrix
        
        D,z	= cheb(NN)
        D2	= np.matmul(D,D)

        
        if mapp =='root':
        
            #Square root matrix 
            D,D2,r,dzdr   = mapping_squate_root(NN,z,D,D2)
        
        elif mapp =='cyl': 
        
            #Cylindrical Mapping 
            #Concentrate the points in L/2
            D,D2,r = mapping_cylindrical(NN,z,D,D2,L)
        
            #Warning
            D2=dot(D,D) 
        
        elif mapp=='tang': 
        
            #only work with L=0
            #Concentrate the points near zero 
            #Tangent mapping
            D,D2,r  = mappingtan_cartesian(NN,z,D,D2)

        #Another way to calculated the second 
        #derivative with the mapping 
        #D2	= np.matmul(D,D)
        
        #Central Finite Diference  diferential Matrix
    else: 

        #To work with the mapping[0:1]computational domain 
        if (mapp=='cyl' or mapp=='root'or mapp=='tang') : 
        
            Min=-1.0
            Max= 1.0
            #z=computational domain
            #dr=spatial step computational domain
            D,D2,deltaz,z	=   difference(NN,Max,Min) 
        
            if mapp=='cyl': 
                D,D2,r = mapping_cylindrical(NN,z,D,D2,L)
        
            elif mapp=='root': 
                D,D2,r,dzdr = mapping_squate_root(NN,z,D,D2)
        
            elif mapp=='tang': 
                #only work with L=0
                #Concentrate the points near zero 
                #Tangent mapping
                D,D2,r  = mappingtan_cartesian(NN,z,D,D2)
        
        elif  mapp=='point': 
        
        #To work with the mapping[0:1]computational domain 
            Min=0.0
            Max=1.0
        
            D,D2,deltaz,z	=   difference(N,Max,Min) 
        
            D,D2,r = mapping_point(N,z,D,D2,L)
        
        elif  mapp=='not': 
        
        #Define the difference matrix domain
            ##Heigh of the phisical domain 
            Min=0.0
            Max= L
            D,D2,deltaz,r	=   difference(NN,Max,Min) 

    return D,D2,r


def p_spectral_matrices(NN,omega):

    D,D2,r=difference_matrix(NN)
        
    #Type of Base Flow
    if (base=='Joncas_jet'): 
        #  Base flow Joncas, Coaxial cylindrical jets 
        Wb,Tb,Rhob,Y1,Y2  =   Baseflow_joncas(r,NN)
        
    elif (base=='hu_ml'): 
        #Base flow Hu, planar mixing layer  
        Wb,Tb,Rhob,Y1,Y2  =   Baseflow_hu(r,NN)
    
    elif (base=='Morris_ml'): 
        #Base flow Morris, planar mixing layer  
        #Displacement mixing layer L/2 
        Wb,Tb,Rhob,Y1,Y2  =   Baseflow_morris_mixing(r,NN,L,D,D2)


    #Fortran program diffusion.f90 
    #   This program calculate the diffusion properties of  
    #   the differens species and the mixture of gases. 
    
    #To run the fortran subroutine.
    #f2py -m diffusion -c  global.f90 diffusion.f90 
    (mib,kb,cpb,cpb1,cpb2,D1m,Mw,Ru,Rmix,gammab,gamma_ref,gamma_j,gamma_ratio,P_r,a_ref,a,a_1,a_2)=\
    df.properties.init_diffusion(Y1,Y2,Tb,T_i,T_0,T_ref,NN)

    #Non dimensional gas ideal equation \bar{p}=1/{gamma_ref}
    #barP   =   gamma_ref*gammab*(1.0/(gamma_ref))*Jh
    #Rhob   =  Pb*gamma_ref/(Tb*Rmix) 

    #Rhob   =  1.0 

    #speed of Sound  
    #a  =\gammaRT=(\bar{p}{gamma_ref})/(\gamma\rho) 
    #To make an incompressible case 
    #a   =  Jh/(gammab*Rhob) 
    a   =  Jh/(Rhob) 

    #Data of the base flow
    #g   	=   open('baseflow_%s.dat'%base,'w+') 
	#
	#for j in range (0,NN+1):
	#
	#            g.write("%d\t%f\t%f\t%f\t%f\t%f\t%f\t%f\n"%(j,r[j],Wb[j],Tb[j],Rhob[j],a[j],Y1[j],Y2[j]));

    #g.close()

    #Derivatives 
    #calculate the derivatives before of modified the 
    #difference matrix 

    dWbdr       =   dot(D,Wb)
    d2Wbdr      =   dot(D2,Wb)

	#Density derivative 
    dRhobdr     =   dot(D,Rhob)

    #print dRhobdr

    #Boundary conditions
    #Remove the dependence of the function in the bounds
    #Use difference matrix
    #Remove the first and the last row, to impost 
    #Dirichlet boundary condition. 

    #Which derivative is use
    #D =np.zeros((N-1,N-1))
    #D2=np.zeros((N-1,N-1))
    
    D       =  D[1:NN,1:NN]
    D2      = D2[1:NN,1:NN]

    #Boundary Conditions
    #Applying the boundary condition 
    #Square matrix, central difference
    #Cheb.py

    if diffM== 'Finite_d':

        diff_Dirichlet(D,D2,deltaz,NN) 

    #If was done to D must be done to
    #the other variables 
    r       =       r[1:NN]
    a       =       a[1:NN]
    Wb      =      Wb[1:NN]
    dWbdr   =   dWbdr[1:NN]
    d2Wbdr  =  d2Wbdr[1:NN]
    Rhob    =    Rhob[1:NN]
    dRhobdr = dRhobdr[1:NN]



	## Especial Matrices
    Ov      =   np.zeros((1,(3*(NN-1))))
    Z   	=   np.zeros((NN-1,NN-1))
    I   	=   np.identity(NN-1) 
    DD      =   block_diag(D,D,D); 


	##########
	#Reference pressure 
	#nn= control de geometric cooridantes 
	

    #LESS_LIN EQUATION 
    Right1 =   np.diag(Wb*(1.0-(Wb/a)**2.0))

    L_2    =   np.diag(omega*(3.0*(Wb/a)**2.0-1.0))

    if(nn==0):

            L1_0 = Z 
    else:
            L1_0 = -np.matmul(diag((Wb*nn)/r),D)

    L_1    = - np.matmul(np.diag(Wb),D2)                    \
             + L1_0                                         \
           	 + np.matmul(np.diag(2.0*dWbdr),D)              \
             + np.matmul(np.diag(Wb/Rhob*dRhobdr),D)        \
     - np.diag(3.0*omega**2.0*Wb/(a**2.0)) 	        \
             + np.diag(Wb*(m/(r**nn))**2.0)

    if(nn==0):

            L0_1 = Z 
    else:
            L0_1 = np.multiply(omega*nn/r,D)

    L_0     =   np.multiply(omega,D2)                              \
              + L0_1                                               \
              - np.matmul(diag(omega*1/Rhob*dRhobdr),D)            \
              + np.multiply(omega**3.0/(a**2.0),I)                 \
              - np.multiply((omega*(m/(r**(nn)))**2.0),I)


    #Incompressigle Rayleight equation

    #Right1  =   np.diag(Wb)

    #L_2     =  -np.multiply(omega,I)

    #L_1     =  -np.matmul(np.diag(Wb),D2)+np.diag(d2Wbdr)

    #L_0     =  +np.multiply(omega,D2)
 

    A0   	=   np.block([[  I ,  Z ,  Z   ],       \
                               [  Z ,  I ,  Z   ],       \
                               [ -L_1, -L_2, -L_0  ]]).astype(np.float64)

    B0   	=   np.block([[Z,   Z  ,I],             \
                           [I,   Z  ,Z],             \
                           [Z,Right1,Z]]).astype(np.float64)


     #python form for generalized eigenvalue problem  
     #Ap = kBp

    return A0,B0,DD




