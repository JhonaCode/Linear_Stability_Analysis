#Program to read the files 
#of the eig program, compared 
#the eignevalues and make a graph 
#alphai omega

import 	numpy 			as 	np 
import 	matplotlib.pyplot 	as 	plt 
import  matplotlib as mpl
from   	plotparameters 		import 	*
from    Parameters 		import  path,     \
                                        analysis

exp= 'Df%s'%(50) 

M	=	np.loadtxt('%s/eig_%s_%s.dat'%(path,analysis,exp),unpack=True)
omega0_1=	M[0,:]
x0_1	=	M[1,:]
y0_1	=	M[2,:]
phase0_1=	(omega0_1)/(x0_1)

exp= 'Df%s'%(100) 

M	=	np.loadtxt('%s/eig_%s_%s.dat'%(path,analysis,exp),unpack=True)
omega1_1=	M[0,:]
x1_1	=	M[1,:]
y1_1	=	M[2,:]
phase1_1=	(omega1_1)/(x1_1)

exp= 'Df%s'%(150) 

M	=	np.loadtxt('%s/eig_%s_%s.dat'%(path,analysis,exp),unpack=True)
x2_1	=	M[1,:]
y2_1	=	M[2,:]
omega2_1=	M[0,:]
phase2_1=	(omega2_1)/(x2_1)

exp= 'Df%s'%(200) 

M	=	np.loadtxt('%s/eig_%s_%s.dat'%(path,analysis,exp),unpack=True)
x3_1	=	M[1,:]
y3_1	=	M[2,:]
omega3_1=	M[0,:]
phase3_1=	(omega3_1)/(x3_1)


#Load plot definitions 
mpl.rcParams.update(params)

#Autospectro 
#################################################
#New Figure
fig = plt.figure()
#New axis  
ax  = plt.axes()
ax.legend()
ax.legend(frameon=False)

#You must select the correct size of the plot in advance
#fig.set_size_inches(3.54,3.54) 

#plt.xlabel(r' $\mathrm{k_r}$')
#plt.ylabel(r'-$\mathrm{k_i}$')
#plt.xlabel(r'$\alpha_r$')
#plt.ylabel(r'$-\alpha_i$')
plt.xlabel(r'$\omega_r$')
plt.ylabel(r'$-\omega_i$')


plt.plot(x0_1,y0_1,color='k',dashes=[2, 1],label = 'D_f, N=50' )
##########x####################         ########
plt.plot(x1_1,y1_1,color='slateblue',dashes=[1, 0],label = 'D_f, N=100' )
##########x####################         ########
plt.plot(x2_1,y2_1,'g',dashes=[1, 1],label = 'D_f, N=150' )
##########x####################         ########
plt.plot(x3_1,y3_1,'m',dashes=[1, 2],label = 'D_f, N=200' )

#ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
#ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
#plt.axis([0.0, 1.0, 0.01, 0.250])

#plt.text(4.0, 0.1,r'ModeI')
#plt.text(0.8, 0.75,r'ModeII')

plt.legend()
ax.legend(frameon=False)
plt.savefig('alphari_%s.pdf'%(exp), format='pdf', dpi=1000)

#Frequencies   
###############################################33
mpl.rcParams.update(params)

fig = plt.figure()
ax  = plt.axes()

#plt.xlabel(r'$\mathrm{\omega}$')
#plt.ylabel(r'-$\mathrm{k_i}$')

plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\omega_i$')


plt.plot(omega0_1,y0_1,color='k',dashes=[2, 1],label = 'D_f, N=50' )

###############################################33
plt.plot(omega1_1,y1_1,color='slateblue',dashes=[1, 0],label = 'D_f, N=100')
#############
plt.plot(omega2_1,y2_1,color='green',dashes=[1, 1],label = 'D_f, N=150')
#############
plt.plot(omega3_1,y3_1,'m',dashes=[1, 2],label = 'D_f, N=200' )
#############
###############################################33


#ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
#ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
#plt.axis([0.0, 0.50, 0.01, 0.25])

plt.text(2.0, 0.15,r'ModeI')
#plt.text(0.5, 0.7,r'ModeII')

plt.legend()
ax.legend(frameon=False)
plt.savefig('walphai_%s.pdf'%(exp), format='pdf', dpi=1000)


#Phase Velocity 
###############################################33
mpl.rcParams.update(params)

fig = plt.figure()
ax  = plt.axes()

plt.xlabel(r'$\alpha$')
plt.ylabel(r'$\mathrm{C_p}$')

plt.ylabel(r'$\omega_i$')


plt.plot(omega0_1,phase0_1,color='k',dashes=[2, 1],label = 'D_f, N=50' )

plt.plot(omega1_1,phase1_1,color='slateblue',dashes=[1, 0],label ='D_f, N=100')
#################################                                 
plt.plot(omega2_1,phase2_1,'g',dashes=[1, 1],label = 'D_f, N=150' )
#################################                                 
plt.plot(omega3_1,phase3_1,'m',dashes=[1, 2],label = 'D_f, N=200' )
#################################


#ax.xaxis.set_major_locator(plt.MultipleLocator(0.1))
#ax.yaxis.set_major_locator(plt.MultipleLocator(0.05))
#plt.axis([0.0, 0.5, 0.45, 0.80])

#plt.text(2.0, 0.5,r'ModeI')
#plt.text(0.5, 0.3,r'ModeII')

plt.legend(loc=2)

ax.legend(frameon=False)
plt.savefig('phase_%s.pdf'%(exp), format='pdf', dpi=1000)

plt.show()


