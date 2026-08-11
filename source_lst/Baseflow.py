#! /usr/bin/python 
# Base flow to spectral methods # create by: jhonatan
# date: 11-05-2018
# Jets Parameters
# Paper: D Perrault-Joncas and S.A Maslowe,
#       "Linear Stability of compressible coaxial jet 
#        with continuos velocity and temperature profile "
#module for symbolic albegra 
import numpy as np
import math  as mt
#from math import exp, expm1
import matplotlib.pyplot as plt
#from sci to make fft
import scipy.fftpack
################################
# letura de aquivos para obter as transformadas
# number of frequency files 

def Baseflow_joncas(r,NN):

    #number of spectral points
    #N1      =   np.shape(r)

    #N       =   N1[0]
    #print(N)

    dr      =   r[1]-r[0]
    
    Wb      =   np.zeros(NN+1)

    Tb      =   np.zeros(NN+1)

    Rhob    =   np.zeros(NN+1)

    M1      =   np.zeros(NN+1)

    M2      =   np.zeros(NN+1)

    M3      =   np.zeros(NN+1)

    T1      =   np.zeros(NN+1)

    T2      =   np.zeros(NN+1)

    T3      =   np.zeros(NN+1)


    #Mass fraction profile
    Yb      =   np.zeros(NN+1)

    Y1      =   np.zeros(NN+1)

    Y2      =   np.zeros(NN+1)

    #Diameter Primary 
    D1      =   2.0*R1 
    #Radii Secondary jet
    R2      =   R1*Gamma
    #Diameter Primary 
    D2      =   2.0*R2

    #Momemtum thickness Primary 
    theta1  =   3.0/100.0*(D1+2.0/3.0*D1) 
    ##Momemtum thickness Secondary
    theta2  =   3.0/100.0*(D1+2.0/3.0*D2) 

    #print 'x',theta1,theta2
    
    #Parameters bases flow

    #Gloor
    #theta1= 0.05
    #theta2= 0.05

    b1      =   R1/(4.0*theta1)
    b2      =   R2/(4.0*theta2)
    
    
    for j in range (0,NN+1):
    
        #Primary Stream
        M1[j]   =   0.5*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1))) 
        #T1[j]   =   0.5*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1))) 

        #Secondary Stream
        M2[j]   =   0.5*(1.0+np.tanh(b2*(R2/r[j]-r[j]/R2))) 
        #T2[j]   =   0.5*(1.0+np.tanh(b2*(R2/r[j]-r[j]/R2))) 

        Wb[j]    = ((1.0-h)*M1[j]+h*M2[j])*M0
        Tb[j]    =  1.0

        #Tb[j]    =  1/(S_T)*(0.5*(1-S_T)*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1)))+S_T)

        Rhob[j]  = 1.0/(Tb[j]) 

        #Primary Stream
        Y1[j]   =   0.5*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1))) 
        #Secondary Stream
        Y2[j]   =   1-(Y1[j])


        #dWbdr   =   np.gradient(Wb, dr)


    return Wb,Tb,Rhob,Y1,Y2

def Baseflow_hu(r,NN):

    #NN = number of spectral points
    #r  = Vertical coordinate. 

    u1b     = 0.80
    u2b     = 0.20
    T1b     = 1.00
    T2b     = 0.80

    Y1b     = 1.00
    Y2b     = 0.00

    delta   = 0.40

    gamma   = 1.40

    Wb      =   np.zeros(NN+1)

    Tb      =   np.zeros(NN+1)

    Rhob    =   np.zeros(NN+1)

    Y1      =   np.zeros(NN+1)

    Y2      =   np.zeros(NN+1)


    for j in range (0,NN+1):

            Wb[j]       = 0.50*(u1b+u2b+(u1b-u2b)*np.tanh(2.0*(r[j])/delta))
                
            Tb[j]       = T1b*(Wb[j]-u2b)/(u1b-u2b)+T2b*(u1b-Wb[j])/(u1b-u2b)+(gamma-1.0)/2.0*(u1b-Wb[j])*(Wb[j]-u2b)

            Rhob[j]     = 1.0/(Tb[j])

            #Primary Stream
            Y1[j]   =   0.50*(Y1b+Y2b+(Y1b-Y2b)*np.tanh(2.0*(r[j])/delta)) 
            #Secondary Stream
            Y2[j]   =   1-(Y1[j])

    return Wb,Tb,Rhob,Y1,Y2

def Baseflow_morris_mixing(r,D,D2,L,N,**kwargs):

    #NN = number of spectral points
    #r  = Vertical coordinate. 

    #L  = Vertival displacement, put on the
    #center of the domain
    L   =L/2.0 

    #Displacement

    Wb      =   np.zeros(N+1)

    Tb      =   np.zeros(N+1)

    Rhob    =   np.zeros(N+1)

    Y1      =   np.zeros(N+1)

    Y2      =   np.zeros(N+1)

    #Analytic Derivative, Morris base flow 
    dwa         =0.5*1.0/((np.cosh(r-L))**2.0)

    d2wa        =np.tanh(r-L)/(-(np.cosh(r-L))**2.0)


    for j in range (0,N+1):

            Wb[j]   = 0.50*(1.0+np.tanh(r[j]-L))  
            Rhob[j] = 1.0                       
            Tb[j]   = 1.0                        
            Y1[j]   = 0.50*(1.0+np.tanh(r[j]-L))  
            Y2[j]   = 1-(Y1[j])                 

    #To Calculated the derivatives 
    #of the base flow, necessary in the 
    #Less_lin equation. 
    #The derivative matrix must be                                                                                                      
    #modified to put the boundary conditions and 
    # it will not be more useful to calculate the 
    #derivatives. This is a important step. 

    #To Calculated the derivatives 
    #of the base flow, necessary in the 
    #Less_lin equation. 

    
    dWbdr   =   np.dot(D,Wb)
    d2Wbdr  =   np.dot(D2,Wb)

    #plt.plot(Wb,r,'b*-')
    #plt.plot(np.dot(D,Wb),r,'*r',dwa,r,'-b')
    #plt.plot(np.dot(D2,Wb),r,'*r',d2wa,r,'-b')
    #plt.axis([-2.0, 2.0, -10.0, 10.0])
    #plt.show()

    #exit()


    #pl.figure()

    #pl.figure()
    #pl.plot(Wb,r,'*-',a,r,'*-')
    #pl.plot(Rhob,r,'*-',a,r,'*-')

    return Wb,dWbdr,d2Wbdr,Tb,Rhob,Y1,Y2
