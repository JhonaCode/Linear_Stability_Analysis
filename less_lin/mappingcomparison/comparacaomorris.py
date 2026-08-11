#Program to read the files 
#of the eig program, compared 
#the eignevalues and make a graph 
#alphai omega

import 	numpy 			as 	np 
import matplotlib
#matplotlib.use('agg')
import 	matplotlib.pyplot 	as 	plt 
import  matplotlib as mpl
#matplotlib.pyplot.switch_backend('agg')
from    plotparameters  import * 


#Path of the files 
path1        = '../morris_points_ml.txt'

path2        = 'eig_Icomp_Chebyshev_sqrtm_150_0.dat'
path3        = 'eig_Icomp_Chebyshev_tanmp_150_0.dat'
path4        = 'eig_Icomp_Chebyshev_cylyn_150_0.dat'

path5        = 'eig_Icomp_Fnitediff___not_150_0.dat'
path6        = 'eig_Icomp_Fnitediff_point_150_0.dat'
path7        = 'eig_Icomp_Fnitediff_sqrtm_150_0.dat'
path8        = 'eig_Icomp_Fnitediff_tanmp_150_0.dat'
path9        = 'eig_Icomp_Fnitediff_cylyn_150_0.dat'

M	=	np.loadtxt('%s'%(path1),unpack=True)
omega1  =	M[0,:]
alpha1  =	M[1,:]

M	=	np.loadtxt('%s'%(path2),unpack=True)
omega2  =	M[0,:]
alpha2	=	M[2,:]

M	=	np.loadtxt('%s'%(path3),unpack=True)
omega3  =	M[0,:]
alpha3	=	M[2,:]

M	=	np.loadtxt('%s'%(path4),unpack=True)
omega4  =	M[0,:]
alpha4	=	M[2,:]

M	=	np.loadtxt('%s'%(path5),unpack=True)
omega5  =	M[0,:]
alpha5	=	M[2,:]

M	=	np.loadtxt('%s'%(path6),unpack=True)
omega6  =	M[0,:]
alpha6	=	M[2,:]

M	=	np.loadtxt('%s'%(path7),unpack=True)
omega7  =	M[0,:]
alpha7	=	M[2,:]

M	=	np.loadtxt('%s'%(path8),unpack=True)
omega8  =	M[0,:]
alpha8	=	M[2,:]

M	=	np.loadtxt('%s'%(path9),unpack=True)
omega9  =	M[0,:]
alpha9	=	M[2,:]


#################################################
#Load plot definitions 
mpl.rcParams.update(params)

#New Figure
fig = plt.figure()
#New axis  
ax  = plt.axes()

plt.ylabel(r'-$\mathrm{\alpha_i}$')
plt.xlabel(r'$\mathrm{\omega_i}$')


plt.plot(omega1,alpha1,color='k',dashes=[1, 0],label = "%s"%(path1[10:20]) )

plt.plot(omega2,alpha2,'-*',color='olive',label = '%s'%(path2[20:38]) )

plt.plot(omega3,alpha3,'-+',color='r'    ,label = '%s'%(path3[20:38]) )
plt.plot(omega4,alpha4,'-.',color='g'    ,label = '%s'%(path4[20:38]) )
plt.plot(omega5,alpha5,'-+',color='b'    ,label = '%s'%(path5[10:19]) )
plt.plot(omega6,alpha6,'-*',color='gray' ,label = '%s'%(path6[20:38]) )
plt.plot(omega7,alpha7,'-*',color='magenta',label = '*%s'%(path7[20:38]) )
plt.plot(omega8,alpha8,'-*',color='cyan' ,label = '*%s'%(path8[20:38]) )
plt.plot(omega9,alpha9,'-*',color='gold' ,label = '%s'%(path9[20:38]) )
#
plt.legend(frameon=False)

plt.savefig('walphaiG2.pdf', format='pdf', dpi=1000)

plt.show()


