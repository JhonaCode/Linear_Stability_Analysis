#################################3#
# Program to generate the 
# instability curves, necesaries 
# reproduce the Joncas and Perreu paper 
#
# Create by: Jhonatan Aguirre 
# Date:08/11/2018
# working: no
###################################


import numpy      	as     np
import matplotlib 	as     pl
from   Parameters  	import *
from   scipy.linalg  	import eig
from   scipy.special 	import airy
from   spectral   	import *
from   Boundary   	import *
from   Eigenvalues   	import *
from   path   		import pathfinder
from   subprocess       import call


# Creates a folder in the current directory called data
def createFolder(directory):
    try:
        if not os.path.exists(directory):
            os.makedirs(directory)
    except OSError:
        print ('Error: Creating directory. ' +  directory)
# Example
#createFolder('./data/')

#file1       =   'incompressible' 
#
#createFolder('./%s/'%file1)

g1       =   open('%s1.dat'%file1,'w+')
g2       =   open('%s2.dat'%file1,'w+')
g3       =   open('allspec.dat','w+')

#domega	=	(omegaf-omegai)/iomega 
omega	=	np.zeros(iomega)
omega	=	np.arange(0,iomega,1)*domega+omegai


alphai  =   	np.zeros((10000,10000)) 
alphar  =   	np.zeros((10000,10000))  


N1	=	N
		
for iomega in range(0,iomega): 

		print(omega[iomega],'omega')

		# Definition of psedo_spectral matices 

		A1,Ee1,B1	= 	p_spectral_matrices(N1,maxr)

		#temporal part or A0, change with omega
		A1      	=	np.multiply((1j*omega[iomega]),Ee1)+A1

		alphar1,alphai1	=	eigenvalues_f(A1,B1,N1,min_imag,max_imag,min_real,max_real)

	        N11		=	alphai1.shape[0]	

		i = 0

		while i < 10:

			i += 1;

		        #plt.plot(alphar1,alphai1,'*b')

			N2		=	N1+50;
		
                        print('%d.Comparando com %d pontos'%(i,N2))
			A2,Ee2,B2	= 	p_spectral_matrices(N2,maxr);

			#temporal part or A0, change with omega
			A2      	=	np.multiply((1j*omega[iomega]),Ee2)+A2;

			alphar2,alphai2	=	eigenvalues_f(A2,B2,N2,min_imag,max_imag,min_real,max_real);

		        N22		=	alphai2.shape[0];	

			#plt.plot(alphar2,alphai2,'.r')
			#plt.show()
			#os.pause

			cont 		= 	0
			
			#if (omega[iomega]>0.6 and omega[iomega]<2.0 ):

				#acoustic	= 	0.01
				#nm=1
			#else:
			acoustic	= 	0.001
				#nm=0
#
			for j in range (0,N11):

            			for k in range (0,N22):

					if (alphai2[k]<max_imag) and (alphar2[k]<max_real):

                				if (  (np.abs(alphai1[j]-alphai2[k])<tol)):
                					if ( (np.abs(alphar1[j]-alphar2[k])<tol)):
	        		#       if (alphar2[k]> 0.0 and (np.abs(alphar1[j]-alphar2[k])<0.005)):
								print 'Found Eingevalue(alphar,alphai)',alphar2[k],alphai2[k]
                    
							        alphai[iomega,cont]=alphai2[k]
							        alphar[iomega,cont]=alphar2[k]
       		                				cont=   cont+1
							       #       		             alphar[cont2,cont]=alphar2[k]
			# # de modos a encontrar
			if (N2>1100):
				cont=-2
				break	

			if (cont>nm):

				N1		=   N2-50
				break	

			alphar1,alphai1 =   alphar2,alphai2	
		        N11		=   alphai2.shape[0];	
			N1		=   N2


		if (cont==1):
	        	g1.write("%f\t%f\t%f\n"%(omega[iomega],alphar[iomega,0],alphai[iomega,0])); 
		elif (cont==2):
	        	g2.write("%f\t%f\t%f\t%f\t%f\n"%(omega[iomega],alphar[iomega,0],alphai[iomega,0],alphar[iomega,1],alphai[iomega,1])); 

		if (cont>0):

			k=cont

			for k in range (0,k):

	        		g3.write("%f\t%f\t%f\n"%(omega[iomega],alphar[iomega,k],alphai[iomega,k])); 

		elif (cont==-2):
			break

		else:

	        	alphai[iomega,cont]=-1.0
	            	alphar[iomega,cont]=-1.0

		#print 'hola2',omega[iomega],alphar[iomega,0],alphai[iomega,0]
	       ##g.write("%f\t%f\t%f\t%f\t%f\n"%(omega[j],alphai[j,0],alphar[j,0],alphai[j,1],alphar[j,1])); 
						
                plt.plot(alphar[iomega,0],alphai[iomega,0],'*b')
                plt.plot(alphar[iomega,1],alphai[iomega,1],'*b')
		
	##plt.plot(alphar,alphai,'*b')
g1.close()
g2.close()
g3.close()
plt.show()

