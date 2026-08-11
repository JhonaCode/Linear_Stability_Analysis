import numpy         as     np
from   scipy.linalg      import eig
from   source_lst.Boundary      import *

def eigenvalues_f(A0,B0,DD,N0,min_imag,max_imag,min_real,max_real):

    #Boundary Conditions
    #boundary_condition_diff(A0,B0,DD,N0)   
    #boundary_condition(A0,B0,N0)   

    #Eingenvalues and Eingenvectors 
    #print(max_imag)
        #Ap = kBp

    eigvals, eigvecs = eig(A0,B0); 

    S   =   list(filter(lambda x: (min_imag<np.imag(x)<max_imag) and (min_real< np.real(x) <max_real), eigvals));

    #Ssort= np.sorted(S,key=getKey)

    Ssort=sorted(S, key=lambda x: x.imag,reverse=True); 
    #Ssort=  np.sort(S,lambda x:x.real) 
    #print(Ssort[0])
    #print(Ssort[1])

    #print(S,'ff')
    #S= eigvals 

    Lambdareal =   np.real(Ssort);
    Lambdaimag =   np.imag(Ssort);

    nLambda    =   Lambdaimag.shape[0]

    return Lambdareal,Lambdaimag,nLambda
