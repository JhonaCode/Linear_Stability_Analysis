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





