#################################3#
# Program to generate the 
# instability curves, necesaries # reproduce the Joncas and Perreu paper 
#
# Create by: Jhonatan Aguirre 
# Date:08/11/2018
# working: no
###################################
# Update by: Jhonatan Aguirre 
# Date:22/01/2020
# Anothes mapping functions.
# working: no
###################################


import numpy      	as     np
import matplotlib 	as     pl
from   Parameters  	import N,file1,file2
from   spectral   	import * 
from   Eigenvalues   	import *
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

g1       =   open('%s.dat'%(file1),'w+')
g2       =   open('%s.dat'%(file2),'w+')


nomega	 =   np.int(round((omegaf-omegaini)/domega))
omega    =   np.arange(0,nomega,1)*domega
omega[0] =   omegaini


alphai  =   	np.zeros((10000,10000)) 
alphar  =   	np.zeros((10000,10000))  


N1	=	N

		
for omegai in omega:

                print omegai
		# Definition of psedo_spectral matices 
		A1,B1,DD	= 	p_spectral_matrices(N1,omegai)

		alphar1,alphai1,nalpha=	eigenvalues_f(A1,B1,DD,N1,min_imag,max_imag,min_real,max_real)



		for k in range(0,nalpha): 

	        	g2.write("%f\t%f\t%f\n"%(omegai,alphar1[k],alphai1[k])); 

		if(nalpha>0):

                        #Write the maximum eigenvalue in a file
	                g1.write("%f\t%f\t%f\n"%(omegai,alphar1[0],alphai1[0])); 


			#alphai_1     =   np.max(alphai1)
    			pt           =   np.argmax(alphai1)
    			#alphai1[pt]  =   np.inf
			#alphar_1     =   alphar1[pt]

			#print omega[iomega],(alphar_1,alphai_1)

			#g1.write("%f\t%f\t%f\n"%(omega[iomega],alphar_1,alphai_1)); 

    			#alphai_2     =   np.min(alphai1)
    			#pt           =   np.argmin(alphai1)
    			#alphai1[pt]  =   np.inf
			#alphar_2     =   alphar1[pt]

			#print omega[iomega],(alphar_1,alphai_1),(alphar_2,alphai_2)

			#g2.write("%f\t%f\t%f\t%f\t%f\n"%(omega[iomega],alphar_1,alphai_1,alphar_2,alphai_2)); 

		else: 
			break

		plt.plot(alphar1,alphai1,'r*')
		#plt.show()
                #plt.plot(alphar[iomega,0],alphai[iomega,0],'*b')
                #plt.plot(alphar[iomega,1],alphai[iomega,1],'*b')
		
g1.close()
g2.close()
plt.show()

