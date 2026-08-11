from numpy import pi,tan,arange,ones,tile,dot,eye,diag,linspace,sqrt,zeros,log,cos,sin,exp,sinh,arcsinh,diff
import matplotlib.pyplot as plt

def mapping_cylindrical(N,n,D1,D2,L):


    #Cylindrical Coordinates
    #Paremeters of mapping
    #Mapping in this point 
    #concentrated around r=rc / 2
    rc = L/2.0
    #Max R= L

    #radial coordinate 
    ri	= 1.e-6

    r   = rc * (1.0-n)/(1.0-n**2.0 + 2.0*rc/L) 
    r[0]= ri


    nr  = -1*(-rc**2/(2*r**3) + rc/r**2)/(2 *sqrt(1 + (2 *rc)/L + rc**2/(4 *r**2) - rc/r)) - rc/(2*r**2)

    nrr = (-rc**2/(2 *r**3) + rc/r**2)**2/(4 *(1 + (2 *rc)/L + rc**2/(4 *r**2) - rc/r)**(3/2)) - ((3*rc**2)/(2*r**4) - (2*rc)/r**3)/(2 *sqrt(1 + (2 *rc)/L + rc**2/(4 *r**2) - rc/r)) + rc/r**3

    for ind in range(0,N+1):

        D1[ind,:]  = D1[ind,:]*nr[ind]

        D2[ind,:]  = D2[ind,:]*(nr[ind]**2.0)+D1[ind,:]*nrr[ind] 

    #plt.plot(r[1:N],diff(r)[1:N],'b*-')
    #plt.plot(n,r,'b*')
    #plt.show()

    #exit()


    return D1,D2,r 


def mapping_squate_root(N,n,D,D2):

    #Scale factor  
    # The scaling factor r=Beta controls the distribution of grid points. Increasing r decreases
    # the number of grid points clustered around y =0, [0,5]
    # Beta=1, without mapping 
    
    #A Comparison of Numerical Methods
    #for the Rayleigh Equationin
    #Unbounded Domains
    #P. J.Morris and W.W. Liou

    D1_o     = zeros((N+1,N+1))
    D2_o     = zeros((N+1,N+1))

    #WARNIGGGGGGGGGGGGGGGGG
    #IGUALANDO MATRICES EM PYTHON, SE NAO SE CREA UM 
    #PONTEIRO, REFERENCIANDO A MESMA LOCALIZACAO. 

    #D1_o[:,:]    = D[:,:]#zeros((NN+1,NN+1))
    #D2_o[:,:]    = D2[:,:]#zeros((NN+1,NN+1))

    Beta    = 5.0
    
    #dzdr    = zeros(N-1) 
    #d2zdr   = zeros(N-1) 

    dzdr    = zeros(N+1) 
    d2zdr   = zeros(N+1) 

    
    for ind in range(0,N+1):
    #for ind in range(0,N-1):

        #First derivative metric 
        dzdr[ind]  = ((1.0-n[ind]**2.0)**(1.5))/Beta

        #Second derivative metric 
        d2zdr[ind] = -3.0*n[ind]*(1.0-n[ind]**2.0)**(1.0/2.0)/(Beta)
    
        # It is not necessary because the fuction was changed to            # the z domain [-1,1]
        D[ind,:]   = D[ind,:]*dzdr[ind]
    
        D2[ind,:]  = D2[ind,:]*(dzdr[ind]**2.0)+D[ind,:]*d2zdr[ind]

        #print ind,n[ind],r[ind]
    
    #z        =  n/((Beta**2.0+n**2)**(0.5)) 
    r       =  n*Beta/((1.0-n**2.0)**(0.5)) 
    
    
    return D,D2,r,dzdr


def mappingtan_cartesian(N,n,D1,D2):

        #Concentrate the poits near the zero 
        #Variable change 
        Beta    = 2.0
        #concentrte more points in the temperature profile between 
        #[0.6 and 1,6], if continued decresing beta you will concentrare the point more to zero.IF Beta>3 desconcnetred the point for 
        #the important region.

        D1_o    = D1[:,:]#zeros((NN+1,NN+1))
        D2_o    = D2[:,:]#zeros((NN+1,NN+1))

        #tan1= infinity

        n[0]=  0.99999
        n[N]= -0.99999

        y       = Beta*tan(pi/2.0*n)

        ny      =  2.0/(Beta*pi)*((cos(pi/2.0*n)*cos(pi/2.0*n)))
        nyy     = -2.0/(Beta)*(cos(pi/2.0*n)*sin(pi/2.0*n))

        for ind in range(0,N+1):

            D1[ind,:]  = D1[ind,:]*ny[ind]
            D2[ind,:]  = D2[ind,:]*(ny[ind]**2.0)+D1[ind,:]*nyy[ind]

        return D1,D2,y 

def mapping_point(N,ybar,D1,D2,L):
        # Computational Fluid Mechanics and Heat Transfer
        # pag:335
        # Tranformation 3
        #Variable change 
        #tau[0...Large values] mesh will refined in yc

        #Streching Parameters 
        tau   = 3.001
        h     = L
        yc    = L/2.0

        ##Computational domain x[0:1], equallly space  
        #Was defined in diff.py
        #Min=  0.0
        #Max=  1.0
        #delta  = Max/(N)  
        #ybar   = Min+arange(0,N+1)*delta


        b0  =   yc/h
        b1  =   1.0+( exp(tau)-1.0)*b0 
        b2  =   1.0+(exp(-tau)-1.0)*b0 
        B   =   1.0/(2.0*tau)*log(b1/b2)

        shtb=    sinh(tau*B)  
        #
        #Phisical domain.
        y   =   yc*(1.0+sinh(tau*(ybar-B))/shtb) 

        #Computational domain.
        #ybar=   B+1.0/tau*arcsinh(((y/yc-1.0))*sinh(tau*B))

        dybardy=zeros(N+1)

        for i in range(0,N+1):

            #First derivative metric 
            d1          =  sqrt(1.0+(y[i]/yc-1.0)**2.0*shtb*shtb)
            dybardy[i]  =  shtb/(tau*yc*d1)

            D1[i,:]     =  D1[i,:]*dybardy[i]


        #To not calcaulated the second derivative of the mapping
        D2  = dot(D1,D1)

        return D1,D2,y      


def mappingtan_boudarylayer(N,n,D1,D2,L):
        # Computational Fluid Mechanics and Heat Transfer
        # pag:335
        # Tranformation 2
        #Variable change 
        #alpha=0 mesh will refined equallly  near y=h
        #alpha=0.5 mesh will refined equallly  near y=0 and y=0

        alpha   = 0.5

        #Beta 1.001-1.5
        beta    = 4.0

        h       = 2*L

        #Computational domain x, equallly space  
        Min=  0.0
        Max=  1.0

        delta  = Max/(N)  

        r      = Min+arange(0,N+1)*delta


        i1      = (beta+1.0)
        i2      = (beta-1.0)
        i3      = (beta+2.0*alpha)
        i4      = (r-alpha)
        i5      = (1.0-alpha)
        i6      = (2.0*alpha+1.0)
        i7      = (i1/i2)**(i4/i5)
        i8      = -beta+2.0*alpha

        #Phisical domain.
        r_f     = (h)*((i3)*i7+i8)/(i6*(1.0+i7))

        #a1      = n*(2.0*alpha+1.0)/h 
        #a2      = beta+a1-2.0*alpha
        #a3      = beta-a1+2.0*alpha

        #r       = alpha+(1.0-alpha)*\
        #          log(a2/a3)/(log((beta+1.0)/(beta-1.0)))


        return r,r_f      


def simpler_mappin_chebchev(N,n,D1,D2):

#Simpler mapping 
    #LL    =   6.0
    scal  =   LL
    r     =   n*scal
    D     =   D/scal
    D2    =   D2/(scal**2.0)
    dzdr  =   1.0/scal 

    return r,D,D2
