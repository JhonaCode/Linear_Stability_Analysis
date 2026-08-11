from numpy import pi,cos,arange,ones,tile,dot,zeros,eye,diag,linspace

def cheb(N):
    #Chebushev polynomial differentiation matrix.
    # Ref.: Trefethen's 'Spectral Methods in MATLAB' book.'''
     
    x      = cos(pi*arange(0,N+1)/N)

    #for ind in range(0,N+1):

    #    print(ind,x[ind])

     
    if N%2 == 0:
    
        x[N/2] = 0.0 # only when N is even!
    
    c      = ones(N+1); c[0] = 2.0; c[N] = 2.0
    
    c      = c * (-1.0)**arange(0,N+1)
    
    c      = c.reshape(N+1,1)  
    
    X      = tile(x.reshape(N+1,1), (1,N+1))      

    dX     = X - X.T
    
    D      = dot(c, 1.0/c.T) / (dX+eye(N+1))

    D      = D - diag( D.sum(axis=1) )
    
    return D,x

def difference(N,Max,Min):

    #Min= first point of the interval  
    #L  = last  point of the interval 

    #size
    L=float(Max-Min) 

    delta  = L/N

    x      = Min+arange(0,N+1)*delta

    D      = zeros((N+1,N+1))


    for i in range(0,N+1): 
        for j in range(0,N+1): 

                if i==j-1:

                        D[i][j]=1.0            
                if i==j+1:

                        D[i][j]=-1.0 


    #Calculate the derivative at the boundaries 
    
    #FORWARD DERIVATIVE
    D[0,0]   =    -2.0
    D[0,1]   =     2.0

    #BACKWARD DERIVATIVE
    D[N,N-1] =     2.0 
    D[N,N]   =    -2.0


    D   =    1.0/(2.0*delta)*D

    #With Dirichelet boudary condition 

    D2=zeros((N+1,N+1))
    
    for i in range(0,N+1): 
        for j in range(0,N+1): 
    
    
                if i==j-1:
                        D2[i][j]=1.0            
                if i==j:
                        D2[i][j]=-2.0            
                if i==j+1:
                        D2[i][j]=1.0 

    #FORWARD DERIVATIVE
    D2[0,0]   =      1.0
    D2[0,1]   =     -2.0
    D2[0,2]   =      1.0

    #BACKWARD DERIVATIVE
    D2[N,N-2] =     1.0 
    D2[N,N-1] =    -2.0 
    D2[N,N]   =     1.0


    D2   =    1.0/(delta**2.0)*D2

    return D,D2,delta,x

def diff_Dirichlet(D,D2,delta,N):

    #Dirichlet Boundary Conditions
    #Central scheme
    #D1
    #Remember that python begins to count at 0
    D[0,0]   =     0.0
    D[0,1]   =     1.0/delta
    #
    D[N-2,N-3] =  -1.0/delta
    D[N-2,N-2] =   0.0

    #Dirichlet Boundary Conditions
    #Central scheme
    #D1
    #Remember that python begins to count at 0
    #and that the first and last row were removed
    D2[0,0]   =     -2.0/delta**2.0

    D2[0,1]   =      1.0/delta**2.0

    D2[0,2]   =      0.0

    #BACKWARD DERIVATIVE
    D2[N-2,N-4] =     0.0 
    D2[N-2,N-3] =     1.0/delta**2.0
    D2[N-2,N-2] =    -2.0/delta**2.0


    return D,D2

def diff_Dirichlet_minha(D,D2,deltax,N):

    #Dirichlet Boundary Conditions
    #Central scheme
    #D1
    D[0,0]   =     0.0
    D[0,1]   =     0.5
    #
    D[N,N]   =     0.0
    D[N,N-1] =    -0.5

    #Dirichlet Boundary Conditions
    #Central scheme
    #D1
    D2[0,0]   =     -0.0
    D2[0,1]   =     -1.0
    D2[0,2]   =      0.0

    #BACKWARD DERIVATIVE
    D2[N,N-2] =     0.0 
    D2[N,N-1] =     1.0 
    D2[N,N]   =    -0.0

    return D,D2


