import numpy        as     np

def boundary_condition_diff(A0,B0,DD,N):
    #[0,N]=N+1 PONTOS
    Ov      =   np.zeros((1,(3*(N+1))))

    # q1
    #A0[0,:]                = DD[0,:]#Ov
    A0[0,:]             = Ov
    B0[0,:]             = Ov
    A0[0,0]             = 1 

    A0[N,:]             = Ov 
    B0[N,:]             = Ov 
    A0[N,N]             = 1  
    
    # q2
    #A0[N+1,:]          = DD[0,:] 
    A0[N+1,:]           = Ov#DD[0,:] 
    B0[N+1,:]           = Ov 
    A0[(N+1),(N+1)]         = 1 
    
    A0[2*N+1,:]         = Ov
    B0[2*N+1,:]         = Ov
    A0[2*N+1,2*N+1]     = 1 

    #p 
    A0[2*(N+1),:]           = DD[0,:]
    B0[2*(N+1),:]           = Ov
    #A0[2*(N+1),:]              = 1 
    
    A0[3*N+2,:]             = Ov
    B0[3*N+2,:]             = Ov
    A0[3*N+2,3*N+2]         = 1

    return
    

def boundary_condition(A0,B0,N):
    #[0,N]=N+1 PONTOS

    Ov      =   np.zeros((1,(3*(N+1))))

    # q1
    A0[0,:]             = Ov
    A0[0,0]             = 1
    B0[0,:]             = Ov

    A0[N,:]             = Ov 
    A0[N,N]             = 1 
    B0[N,:]             = Ov 

    
    # q2
    A0[N+1,:]           = Ov 
    A0[N+1,N+1]         = 1
    B0[N+1,:]           = Ov 
    
    A0[2*N+1,:]         = Ov
    A0[2*N+1,2*N+1]         = 1 
    B0[2*N+1,:]         = Ov

    #  p
    A0[2*(N+1),:]           = Ov
    A0[2*(N+1),2*(N+1)]     = 1.0
    B0[2*(N+1),:]           = Ov
    
    A0[3*N+2,:]             = Ov
    A0[3*N+2,3*N+2]         = 1.0
    B0[3*N+2,:]             = Ov


    return A0,B0

def boundary_condition_temporal(A0,B0,N):
    #[0,N]=N+1 PONTOS

    Ov      =   np.zeros((1,(1*(N+1))))

    # u1
    A0[0,:]             = Ov
    A0[0,0]             = 1
    B0[0,:]             = Ov

    A0[N,:]             = Ov 
    A0[N,N]             = 1 
    B0[N,:]             = Ov 

    
    return A0,B0
