#! /usr/bin/python 
# Base flow to spectral methods 
# create by: jhonatan
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
from  Parameters import * 
################################
# letura de aquivos para obter as transformadas
# number of frequency files 

def Baseflow_geometric(r,gamma,NN):

    #Array initicialization 	
    dr      =   r[1]-r[0]
    
    Wb      =   np.zeros(NN+1)

    Tb      =   np.zeros(NN+1)

    Rhob    =   np.zeros(NN+1)

    u1      =   np.zeros(NN+1)

    u2      =   np.zeros(NN+1)

    u3      =   np.zeros(NN+1)

    #Diameter Primary 
    D1      =   2.0*R1 
    #Momemtum thickness Primary 
    theta1  =   3.0/100.0*(D1+2.0/3.0*D1) 

    #Radii Secondary jet
    R2      =   R1*Gamma
    #Diameter Primary 
    D2      =   2.0*R2
    ##Momemtum thickness Secondary

    theta2  =   3.0/100.0*(D1+2.0/3.0*D2) 
    #theta2  =   0.14 
    
    #Parameters bases flow

    b1      =   R1/(4.0*theta1)
    b2      =   R2/(4.0*theta2)
    
    
    for j in range (0,NN+1):
    
        #Primary Stream
        u1[j]   =   0.5*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1))) 
        #Secondary Stream
        u2[j]   =   0.5*(1.0+np.tanh(b2*(R2/r[j]-r[j]/R2))) 

        #Axial Velocity Stream
        Wb[j]    = ((1.0-h)*u1[j]+h*u2[j])*M0

        Tb[j]    = 1.0
        #Tb[j]   =  (gamma-1)*Wb[j]*(M0-Wb[j])/2.0 

        Rhob[j]  = 1.0/(Tb[j]) 

    
    return Wb,Tb,Rhob

def Baseflow_joncas(r,Gamma,NN):

    #number of spectral points
    #N1      =   np.shape(r)

    #N       =   N1[0]
    #print(N)

    dr      =   r[1]-r[0]
    
    Wb      =   np.zeros(NN+1)

    Tb      =   np.zeros(NN+1)

    Rhob    =   np.zeros(NN+1)

    u1      =   np.zeros(NN+1)

    u2      =   np.zeros(NN+1)

    u3      =   np.zeros(NN+1)

    
    #Diameter Primary 
    D1      =   2.0*R1 
    #Momemtum thickness Primary 
    theta1  =   3.0/100.0*(D1+2.0/3.0*D1) 

    #Radii Secondary jet
    R2      =   R1*Gamma
    #Diameter Primary 
    D2      =   2.0*R2
    ##Momemtum thickness Secondary

    theta2  =   3.0/100.0*(D1+2.0/3.0*D2) 
    #theta2  =   0.14 
    
    #Parameters bases flow

    b1      =   R1/(4.0*theta1)
    b2      =   R2/(4.0*theta2)
    
    
    for j in range (0,NN+1):
    
        #Primary Stream
        u1[j]   =   0.5*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1))) 
        #Secondary Stream
        u2[j]   =   0.5*(1.0+np.tanh(b2*(R2/r[j]-r[j]/R2))) 

        Wb[j]    = ((1.0-h)*u1[j]+h*u2[j])

	#Crocco-Busemann
        #Tb[j]   =   1.0
        #Tb[j]   =  T_inf+(1.0-T_inf)*(gamma-1)*M0**2*Wb[j](1-Wb[j])/2 
	#hot
        Tb[j]   =   S_T*(1.0-S_T)*(1.0+np.tanh(b1*(R1/r[j]-r[j]/R1)))+S_T
	#cold
        #Tb[j]   =   -0.0592*Wb[j]**2.0-0.1032*Wb[j]+1.1624 

        #Rhob[j]  = 1.0/(Tb[j]) 

    
    return Wb,Tb,Rhob

def Baseflow_morris(r,NN):

    Wb      =   np.zeros(NN+1)

    Tb      =   np.zeros(NN+1)

    Rhob    =   np.zeros(NN+1)

    Mc   = 0.50
    #*********Thickness of the jet, base flow 
    bba  = 0.25
    #*********internal radio of the jet, base flow,solution morrismeanflow.f90 
    hba  = 0.8871800
    #*********Morris, The instability of hight speed jets 
    #*********Reservoir temperature
    TR   = 0.60
    ##*********Ambient temperature
    T0   = 0.60
    #*********Heat capacity relation 
    gamma= 1.40

    for j in range (0,NN+1):
    
        if(r[j]<hba):
            Wb[j]   = Mc  
        else: 
            Wb[j]   = Mc*np.exp(-np.log(2.0)*(r[j]-hba)**2.0/bba**2.0)  
        
    
    for j in range (0,NN+1):
        #Development Jet
        #wb[j]    = Mc*exp(-log(2.d0)*n[j]**2.d0/bba**2.d0)
        Tb[j]    = T0/TR*(1.0+(gamma-1.0)/2.0*Mc**2.0)*(1.0+(Wb[j])*(TR/T0-1.0))-(gamma-1.0)/2.0*Mc**2.0*(Wb[j]**2.0)
        Rhob[j]  = 1.0/(Tb[j])


    #plt.plot(Wb,r,Tb,r,Rhob,r)
    #plt.show()
    

    return Wb,Tb,Rhob
