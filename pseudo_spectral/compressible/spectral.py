#matplotlib inline
#config InlineBackend.figure_format='svg'
from   chebPy        import *
from   mapping       import *
from   Baseflow      import *
from   Boundary      import *
from   scipy.linalg  import block_diag
from   Parameters    import *
import numpy 	         as np 
import matplotlib.pyplot as pl
#import diffusion         as df
import sys 
import os 

def p_spectral_matrices(NN,maxr):

	# number of collocation nodes=N+1, because the cheb func 
	D,n1	= cheb(NN)
	#D2	= dot(D,D)
	#D2[1:N,1:N]
	
	# Mappint to concentrate the point in 0. Withou using 0
	# r= Radial coordinate
	r   	= mappingtan(NN,n1,maxr,D)
	D2	= dot(D,D)
	#Base Flow

        
	Wb,Tb,Rhob =   Baseflow_geometric(r,Gamma,NN)
	#Wb,Tb,Rhob =   Baseflow_joncas(r,Gamma,NN)
	#Wb,Tb,Rhob =   Baseflow_morris(r,NN)


        #Fortran program diffusion.f90 
        #   This program calculate the diffusion properties of  
        #   the differens species and the mixture of gases. 
        
        #To run the fortran subroutine.
        #f2py -m diffusion -c  global.f90 diffusion.f90 
        
        #(mib,kb,cpb,cpb1,cpb2,D1m,Mw,Ru,Rmix,gammab,gamma_j,P_r)=\
        #df.properties.init_diffusion(Y1,Y2,Tb,T_r,T_0,T_ref,NN) 

        g	    =   open('baseflow.dat','w+')

        #####
        #Non dimensional gas ideal equation \bar{p}=1/{gamma_f}
        #Rhob  =  1.0/(Tb*Rmix) 

        #plt.figure(1)
        #plt.plot(Wb,r,'*-c')
        #plt.plot(Tb,r,'*-b')
        #plt.ylim((0,4))
        #plt.figure(2)
        #plt.plot(Rhob,r,'*-m')
        #plt.ylim((0,4))
        #plt.figure(3)
        #plt.plot(Y1,r,'o-y')
        #plt.plot(Y2,r,'o-r')
        #plt.ylim((0,4))
        #plt.show()

        #sys.exit()

        for j in range (0,NN+1):

	    g.write("%f\t%f\t%f\t%f\n"%(r[j],Wb[j],Tb[j],Rhob[j])); 
    
    	#dWbdr   =   np.gradient(Wb, dr)

    
        g.close()
	
	dWbdr      =   dot(D,Wb)
	dRhobdr    =   dot(D,Rhob)
	
	# coordinate selectors
	#rr      =   np.arange(0,N,1)#0:N+1      
	#uu      =   rr+N
	#vv      =   uu+N
	#pp      =   vv+N
	
	
	## Especial Matrices

	Ov      =   np.zeros((1,(3*(NN+1))))
	Z   	=   np.zeros((NN+1,NN+1))
	I   	=   np.identity(NN+1) 
	II      =   np.identity(3*(NN+1))
	Ee 	=   block_diag(I,I,I)
	DD      =   block_diag(D,D,D);
	

	#Cartesian
	#A   =   np.block([[diag(M),diag(Rho),Z,Z],[Z,diag(M),Z,diag(1./Rho)],[Z,Z,diag(M),Z],[Z,I,Z,diag(M)]])
	#B   =   np.block([[Z,Z,diag(Rho),Z],[Z,Z,Z,Z],[Z,Z,Z,diag(1./Rho)],[Z,Z,I,Z]])
	#C   =   np.block([[Z,Z,diag(Rhoy),Z],[Z,Z,diag(My),Z],[Z,Z,Z,Z],[Z,Z,Z,Z]])

        #Mj      =   0.5
        gamma   =   1.4

	#P base flow 
	#joncas
        #barP   =   gamma*(1.0/(gamma*M0**2.0))
        #barP   =   gamma*(1.0/(gamma*0.001**2.0))

	#jhonatan Parameter compressibility  
        #barP   =   gamma*(1.0/(gamma))
        barP   =   gamma*(1.0/(gamma))*Jh

	#Cylindrical
	Ae   =   np.block([
                          [Z         ,Z,          diag(1.0/Rhob)],
                          [Z         ,Z,          Z		],
                          [I*barP    ,Z,          Z		]
                          ])

	Be   =   np.block([
                          [diag(Wb)   ,Z          ,Z		],
                          [Z          ,diag(Wb)   ,diag(1.0/Rhob)],
                          [Z          ,I*barP     ,diag(Wb)] 
			  ])

        Cr      =   np.zeros(NN+1)

	for j in range(0,NN+1): 
	
	        Cr[j]   = barP/r[j] 

	Ce   =   np.block([
                          [Z          ,Z,Z],
                          [diag(dWbdr),Z,Z],
                          [diag(Cr)   ,Z,Z]
                          ])
	
	
	A0      =   -np.matmul(Ae,DD)-Ce;

	B0      =   np.multiply(1j,Be);
	

	return A0,Ee,B0




