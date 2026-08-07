import numpy      	as     np

def boundary_condition(A0,B0,N):
	#[0,N]=N+1 PONTOS
	Ov      =   np.zeros((1,(3*(N+1))))
	# rho
	A0[0,:]        	 	= Ov
	B0[0,:]        	 	= Ov
	A0[0,0]        	 	= 1#A0[1,0]

	A0[N,:]        	 	= Ov 
	B0[N,:]        	 	= Ov 
	#A0[N,N]        	 	= 1  
	
	#u
	A0[N+1,:]    		= Ov 
	B0[N+1,:]    		= Ov 
	#A0[N+1,N+1]    		= 1#A0[N+1+1,0]  
	
	A0[2*N+1,:]     	= Ov
	A0[2*N+1,:]     	= Ov
	A0[2*N+1,2*N+1]  	= 1 

	#v
	A0[2*(N+1),:]          	= Ov
	B0[2*(N+1),:]          	= Ov
	#A0[2*(N+1),2*(N+1)]     = 1.0#A0[2*(N+1)+1,0]
	
	A0[3*N+2,:]             = Ov
	B0[3*N+2,:]             = Ov
	A0[3*N+2,3*N+2]         = 1
	
	#p
	#A0[3*(N+1),:]          	= Ov
	#B0[3*(N+1),:]          	= Ov
	##nao posso colocar qualquer coisa na pressoa
	##interpolation in the pressure
	#A0[3*(N+1),3*(N+1)]    	= A0[3*(N+1)+1,0]
	#
	#
	#A0[4*N+3,:]             = Ov
	#B0[4*N+3,:]             = Ov
	#A0[4*N+3,4*N+3]         = 1

	# rho
	#A0[0,:]               = Ov
	#B0[0,:]               = Ov
	#B0[0,0]               = 1 # A0[N,:]               = Ov 
	#B0[N,:]               = Ov 
	#B0[N,N]               = 1  
	#
	##u
	#A0[N+1,:]             = Ov 
	#B0[N+1,:]             = Ov 
	#B0[N+1,0]             = 1  
	#
	#A0[2*N,:]             = Ov
	#B0[2*N,:]             = Ov
	#B0[2*N,2*N]           = 1
	#
	##v
	#A0[2*N+1,:]           = Ov
	#B0[2*N+1,:]           = Ov
	#B0[2*N+1,N+1]         = 1
	#
	#A0[3*N,:]             = Ov
	#B0[3*N,:]             = Ov
	#B0[3*N,3*N ]          = 1
	#
	##p
	#A0[3*N+1,:]           = Ov
	#B0[3*N+1,:]           = Ov
	#B0[3*N+1,3*N+1]       = 1
	#
	#
	#A0[4*N,:]             = Ov
	#B0[4*N,:]             = Ov
	#B0[4*N,4*N]           = 1
