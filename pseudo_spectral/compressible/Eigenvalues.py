import numpy      	as     np
import matplotlib 	as     pl
from   scipy.linalg  	import eig
from   scipy.special 	import airy
from   spectral   	import *
from   Boundary   	import *

def eigenvalues_f(A0,B0,N0,min_imag,max_imag,min_real,max_real):

	#Boundary Conditions
	boundary_condition(A0,B0,N0)	

	#Eingenvalues and Eingenvectors 
	#print(max_imag)
	eigvals, eigvecs = eig(A0, B0); 

	S   =   list(filter(lambda x: (min_imag<np.imag(x)<max_imag) and (min_real< np.real(x) <max_real), eigvals));

	#Ssort= np.sorted(S,key=getKey)


	Ssort=sorted(S, key=lambda x: x.imag,reverse=True); 
	#Ssort=  np.sort(S,lambda x:x.real) 
	#print(Ssort[0])
	#print(Ssort[1])

	#print(S,'ff')
	#S= eigvals 

	Lamreal =   np.real(Ssort);
	Lamimag =   np.imag(Ssort);

	#Lamreal =   np.real(S);
	#Lamimag =   np.imag(S);

	#alphai1 =   Lamimag[0]
	#alphai2 =   Lamimag[1]

	return Lamreal,Lamimag
