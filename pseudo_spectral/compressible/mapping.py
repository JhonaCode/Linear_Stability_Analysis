from numpy import pi,tan,arange,ones,tile,dot,eye,diag,linspace,sqrt

def mappingtan(N,n,m,D1):


    #Cylindrical Coordinates
    #Paremeters of mapping
    c = 0.9
    #Max R
    #m = 50

    #radial coordinate 
    ri	= 1.e-6
    r   = c * (1-n)/(1-n**2 + 2*c/m) 
    r[0]= ri


#
    nr  = -1*(-c**2/(2*r**3) + c/r**2)/(2 *sqrt(1 + (2 *c)/m + c**2/(4 *r**2) - c/r)) - c/(2*r**2)
    nrr = (-c**2/(2 *r**3) + c/r**2)**2/(4 *(1 + (2 *c)/m + c**2/(4 *r**2) - c/r)**(3/2)) - ((3*c**2)/(2*r**4) - (2*c)/r**3)/(2 *sqrt(1 + (2 *c)/m + c**2/(4 *r**2) - c/r)) + c/r**3
#        
#
    
    for ind in range(0,N+1):

        D1[ind,:]  = D1[ind,:]*nr[ind]
        #D2[ind,:]  = D2[ind,:]*  (nr[ind]**2.0)+D1[ind,:]*nrr[ind] 


    return r 

