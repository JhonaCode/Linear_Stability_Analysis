import numpy as np
import matplotlib.pyplot as plt

def taylor(X,Y,i0,n): # Returns the Taylor polynomial for given data
# Required arguments:
# 	X,Y: array of coordinates;
#	i0: index of the point around which the series will be expanded;
#	n: number of terms in the polynomial.

	a = [Y[i0]] # First coefficient
	for i in range(n):
		Y = np.diff(Y)/np.diff(X) # Derivative of Y
		a += [Y[i0]/np.math.factorial(i+1)] # i+1th coefficient
		X = X[0:-1] # Eliminates last element of X
	return a[-1::-1] # Return backwards array, suitable to be used with numpy.polyval function

def pathfinder(x,y,i0,polydeg=3,distweight=1,stopweight=1,isextremity=True): # Finds a smooth curve from given data
# Required arguments:
#	x,y: arrays of coordinates;
#	i0: initial point of the curve to be searched. Must be far from other curves.
# Optional arguments:
#	polydeg: degree of Taylor polynomial (set to 3). Larger degrees may lead to issues related to the accuracy of the derivatives;
#	distweight: relative weight of the distance to the new point in relation to the deflection from the Taylor polynomial expansion (set to 1);
#	stopweight: relative weight of the last aberration in the stop criterium (set to 1);
#	isextremity: True if i0 is on one extremity of the curve, False otherwise (set to True).
#
# The function identifies the next point to be added to the curve by minimizing the 'aberration'. The aberration is defined by the weighted sum of the distance between the last added point to the next candidate and the deviation of the next candidate to the extrapolation of the previous added points. The function stops searching for new points if the aberration has a sudden increase, which means that the next candidate is a point belonging to another curve. If the provided initial point is not in one extremity of the curve, the search occurs in two steps: first the curve is searched in a direction chosen automatically and after the search stops the opposite direction is searched.

	# Normalization of the data:
	X = x/max(abs(x))
	Y = y/max(abs(y))

	# Finding the closest point to the initial point:
	dist = (X-X[i0])**2 + (Y-Y[i0])**2 # Distance of initial point found to every point on list
	dist[i0] = np.inf # Avoid initial point
	n = np.where(dist==min(dist))[0][0] # Closest point to initial point

	# Initialization of the output with the two initial points:
	xo = np.array([X[i0],X[n]])
	yo = np.array([Y[i0],Y[n]])

	N = np.array([i0,n]) # List of indexes of selected points
	minaber = np.array([0]) # Array of minimum aberrations

	# Main loop:
	for i in range(len(X)):
		dist = (X-xo[-1])**2 + (Y-yo[-1])**2 # Euclidean istance of last point found to every point on list
		A = taylor(xo,yo,-1,min(len(xo)-1,polydeg)) # Taylor polynomial coefficients based on selected points
		defl = abs(Y-np.polyval(A,X-xo[-1])) # Deflection of every point on list to the polynomial
		aber = distweight*dist + defl # Aberration of every point on list
		aber[N] = np.inf # Avoid previous selected points
		n = np.where(aber==min(aber))[0][0] # Finds point of minimum aberration
	
		if i>1 and aber[n]>=stopweight*4*minaber[-1]: # Stoping criterion
			break
		
		minaber = np.append(minaber,aber[n]) # Add last aberration to aberration array
		N = np.append(N,n) # Add new selected point
		xo = np.append(xo,X[n])
		yo = np.append(yo,Y[n])
	
	if not isextremity: # If the initial point is not an extremity, find the curve on the opposite direction
		dist = (X-X[i0])**2 + (Y-Y[i0])**2 # Distance of initial point found to every point on list
		dist[N] = np.inf # Avoid previously selected point
		n = np.where(dist==min(dist))[0][0] # Closest point to initial point

		# Two initial points
		xo2 = np.array([X[i0],X[n]])
		yo2 = np.array([Y[i0],Y[n]])

		N2 = np.array([i0,n]) # List of indexes of selected points
		minaber = np.array([np.inf]) # Array of minimum aberrations

		# Main loop:
		for i in range(len(X)):
			dist = (X-xo2[-1])**2 + (Y-yo2[-1])**2 # Distance of last point found to every point on list
			A = taylor(xo2,yo2,-1,min(len(xo2)-1,polydeg)) # Taylor polynomial coefficients based on selected points
			defl = abs(Y-np.polyval(A,X-xo2[-1])) # Deflection of every point on list to the polynomial
			aber = distweight*dist + defl # Aberration of every point on list
			aber[N2] = np.inf # Avoid previous selected points
			n = np.where(aber==min(aber))[0][0] # Finds point of minimum aberration

			if i>1 and aber[n]>=stopweight*4*minaber[-1]: # Stoping criterion
				break
	
			minaber = np.append(minaber,aber[n]) # Add last aberration to aberration array
			N2 = np.append(N2,n) # Add new selected point
			xo2 = np.append(xo2,X[n])
			yo2 = np.append(yo2,Y[n])
		
		# Adding the points of the opposite side:
		N = np.append(N,N2[1::])
		xo = np.append(xo,xo2[1::])
		yo = np.append(yo,yo2[1::])
		
		# Sorting the output
		ind = np.argsort(xo)
		xo = xo[ind]
		yo = yo[ind]
		N = N[ind]
		
	return (xo*max(abs(x)),yo*max(abs(y)),N) # Reescale the coordinates and return the output

# Test case:

##Curve 1:
#x1 = np.linspace(0,10,50)
#y1 = x1**2/5
#
##Curve 2:
#x2 = np.linspace(2,8,20)
#y2 = x1**2/4+3
#
## Merging both curves:
#ind = np.argsort(np.concatenate([x1,x2])) # Return the index of the sorted concatenated x arrays
#X = np.concatenate([x1,x2])[ind]
#Y = np.concatenate([y1,y2])[ind]
#
## Finding the curves:
#i1 = 40
#(X1,Y1,N) = pathfinder(X,Y,i1,isextremity=False)
#
#i2 = 41
#(X2,Y2,N) = pathfinder(X,Y,i2,isextremity=False)
#
## Ploting the curves:
#plt.plot(X,Y,'k.',label='Raw data')
#plt.plot(X1,Y1,'b--',label='Curve 1')
#plt.plot(X2,Y2,'r--',label='Curve 2')
#plt.plot(X[i1],Y[i1],'bo',label='Initial point of curve 1')
#plt.plot(X[i2],Y[i2],'ro',label='Initial point of curve 2')
#plt.legend()
#plt.axis([min(X),max(X),min(Y),max(Y)])
#plt.xlabel('x')
#plt.ylabel('y')
#plt.show()
