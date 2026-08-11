#Function to calculate the 
#eigenvalues, -alpha__i, for 
#spatial stability analisis.
#and the eiga34envalues omega_i for 
#the temporal analysis.

#Created: Jhonatan A A Manco
#Date:07/01/2020

import numpy         as np

import matplotlib.pyplot as plt

import matplotlib        as mpl

from   scipy.linalg     import block_diag

from   Parameters_teste  import min_imag, \
                               max_imag, \
                               min_real, \
                               max_real, \
                               path,     \
                               file1,      \
                               file2

from   Eigenvalues      import * 


def temporal(N,D,D2,Wb,d2Wbdr,alpha): 

    fl1      =   open('%s/eig_%s.dat'%(path,file1),'w+')
    fl2      =   open('%s/allspec_%s.dat'%(path,file2),'w+')

    #Auxiliary Matrices
    DD      =   block_diag(D,D,D); 
    Z       =   np.zeros((N-1,N-1))
    I       =   np.identity(N-1) 

    #for i in range(0,nalpha): 
    for ialpha in alpha: 

        AB0     =  -(ialpha)**2.0*I
        A_10    =   ialpha*np.diag(Wb)  
        A_20    =   np.matmul(A_10,D2)  
        A_30    =   np.matmul(A_10,AB0)  
        A_40    =  -ialpha*np.diag(d2Wbdr)  

        A0      =   A_20+A_30+A_40
        B0      =   D2+AB0

        ##neig= number of eigenvalues found 
        omegar,omegai,neig= eigenvalues_f(A0,B0,DD,N,min_imag,max_imag,min_real,max_real)

        if neig>0:

            ##Print the most unstable eigenvalue
            print("alpha,omega_i",ialpha,omegai[0])
        fl1.write("%f\t%f\t%f\n"%(ialpha,omegar[0],omegai[0])); 

    for k in range(0,neig): 

        fl2.write("%f\t%f\t%f\n"%(ialpha,omegar[k],omegai[k])); 

        plt.plot(omegar,omegai,'r*')

    return omegar,omegar


def spatial(N,D,D2,Wb,d2Wbdr,omega): 

    #Spatial, Working 

    fl1      =   open('%s/eig_%s.dat'%(path,file1),'w+')
    fl2      =   open('%s/allspec_%s.dat'%(path,file2),'w+')

    #Auxiliary Matrices
    Z       =   np.zeros((N-1,N-1))
    I       =   np.identity(N-1) 
    DD      =   block_diag(D,D,D); 
    
    
    #Matrices that does depent on \omega
    
    Right1  =   np.diag(Wb)
    L_10    =  -np.dot(np.diag(Wb),D2)
    L_20    =   np.diag(d2Wbdr)
    L_1     =   L_10+L_20   
    
    #Evaluating several \omega
    
    for iomega in omega: 
    
        L_2 =  -iomega*I
    
        L_0 =  iomega*D2
    
        A0  =   np.block([[-L_0, -L_1, -L_2  ],  \
                     [  Z ,  I  ,  Z    ],       \
                     [  Z ,  Z  ,  I    ]])
    
        B0  =   np.block([[  Z ,  Z  , Right1],       \
                              [  I ,  Z  ,  Z    ],   \
                              [  Z ,  I  ,  Z    ]])
    
        # Function to calculate the eigenvalues of the 
        # generalize eigenvalue problem. 
    
        alphar,alphai,neig    = eigenvalues_f(A0,B0,DD,N,min_imag,max_imag,min_real,max_real)


        if neig>0:

            ##Print the most unstable eigenvalue
            print("omega,alpha_i,"  ,iomega,alphai[0])

            fl1.write("%f\t%f\t%f\n"%(iomega,alphar[0],alphai[0])); 

        for k in range(0,neig): 

            fl2.write("%f\t%f\t%f\n"%(iomega,alphar[k],alphai[k])); 


        plt.plot(alphar,alphai,'b*')

    return alphar,alphai



#It is the same as dot multiplication 
#L_10   =   np.zeros((N+1,N+1))
#for j in range(0,N+1): 
#
#    L_10[:,j]    =  -Wb[:]*D2[:,j]
#L_10    =   -np.matmul(np.diag(Wb),D2)

