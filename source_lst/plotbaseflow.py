import matplotlib.pyplot as plt
import matplotlib        as mpl
import numpy             as np 

#To define characteristics of the plots
from   source_lst.plotparameters   import * 

def plotbase(N,x,r,Df,Df2,D,D2):   

        #To plot the base flow figures 
        Wb      =   np.zeros(N+1)
        Wbr     =   np.zeros(N+1)
        
        for j in range (0,N+1):
        
                Wbr[j]   = 0.50*(1.0+np.tanh(r[j]))
                Wb[j]    = 0.50*(1.0+np.tanh(x[j]))

        dwa         =0.5*1.0/((np.cosh(x))**2.0)
        dwar        =0.5*1.0/((np.cosh(r))**2.0)
        
        d2wa        =np.tanh(x)/(-(np.cosh(x))**2.0)
        d2war       =np.tanh(r)/(-(np.cosh(r))**2.0)
                
        mpl.rcParams.update(params)

        fig, ax = plt.subplots(1)
        
        plt.text(0.35,1.0,r'$\mathrm{\overline{w}=\dfrac{1}{2}(1+tanh(r))}$') 
        
        plt.axis([-0.01, 1.01, -6.0, 6.0])
        
        plt.plot(Wb,x,'o-m') 
        
        plt.xlabel(r' $\mathrm{\overline{W}}$')
        plt.ylabel(r' r')
        
        # Set min and max of the axis 
        
        #plt.savefig('baseflowmixingmorris.pdf', format='pdf', dpi=1000)
        plt.savefig('/home/inct1/repositories/paper_1_doc/paper_low_mach_overleaf/figures/baseflowmixingmorris.pdf', format='pdf', dpi=1000)

        plt.show()
        
        
        fig1, (ax1, ax2) = plt.subplots(1,2)
        
        ax1.set_xlabel(r' $\mathrm{D\overline{W}}$')
        ax2.set_xlabel(r' Erro ')
        
        ax1.axis([-0.01, 1.02, -4.0, 4.0])
        ax2.axis([-0.0045, 0.004, -4.0, 4.0])
        
        
        ax1.plot(np.dot(D,Wbr),r,'*b',label = 'Pseudo Spectral Matrix') 
        ax1.plot(np.dot(Df,Wb),x,'om',label = 'Finite Diference Matrix') 
        ax1.plot(dwa,x,'-k',label = 'Analytic Derivative')
        
        ax2.plot(np.dot(D,Wbr)-dwar,r,'-*b',label = 'Pseudo Spectral Matrix') 
        ax2.plot(np.dot(Df,Wb)-dwa,x,'-om',label = 'Finite Diference Matrix') 
        
        ax1.legend()
        ax1.legend(frameon=False)
        #ax2.legend()
        #ax2.legend(frameon=False)
        
        plt.savefig('baseflowd1.pdf', format='pdf', dpi=1000)
        
        
        fig1, (ax1, ax2) = plt.subplots(1,2)
        
        ax1.set_xlabel(r' $\mathrm{D^2\overline{W}}$')
        ax2.set_xlabel(r' Erro ')
        
        
        ax1.axis([-0.6, 2.5, -4.0, 4.0])
        ax2.axis([-0.02, 0.02, -4.0, 4.0])
        
        ax1.plot(np.dot(D2 ,Wbr),r,'*b',label = 'Pseudo Spectral Matrix') 
        ax1.plot(np.dot(Df2,Wb),x,'om',label = 'Finite Diference Matrix') 
        ax1.plot(d2wa,x,'-k',label = 'Analytic Derivative')
        
        ax2.plot(np.dot(D2,Wbr)-d2war,r,'-*b',label = 'Pseudo Spectral Matrix') 
        ax2.plot(np.dot(Df2,Wb)-d2wa,x,'-om',label = 'Finite Diference Matrix') 
        
        ax1.legend()
        ax1.legend(frameon=False)
        
        plt.savefig('baseflowd2.pdf', format='pdf', dpi=1000)
        
        
        #ax.xaxis.set_major_locator(plt.MultipleLocator(0.5))
        #y Axis tic 
        #ax.yaxis.set_major_locator(plt.MultipleLocator(0.1))
        
        
        plt.show()



        
