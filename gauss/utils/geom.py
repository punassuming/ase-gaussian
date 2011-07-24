#!/usr/bin/env python

import numpy as np
import math
from matplotlib.pyplot import scatter, show

#import plotsize as ps

def degrees(rad_angle):
	"""Converts any angle in radians to degrees.

	If the input is None, then it returns None.
	For numerical input, the output is mapped to [-180,180]
	"""
	if rad_angle is None :
		return None
	angle = rad_angle * 180 / math.pi
	while angle > 180 :
		angle = angle - 360
	while angle < -180 :
		angle = angle + 360
	return angle

def make_sphere(r, output=None, trans=[0.0,0.0,0.0]):
	"""
	make np.array based on sphere
	trans = [x,y,z]
	"""

	x = []
	y = []
	z = []

	for theta in np.arange(0,360,5):
		for phi in np.arange(0,180,5):
			x.append(r * np.sin(theta * np.pi / 180) * np.sin(phi * np.pi / 180))
			y.append(r * np.cos(theta * np.pi / 180) * np.sin(phi * np.pi / 180))
			z.append(r * np.cos(phi * np.pi / 180))

	x = [i+trans[0] for i in x]
	y = [i+trans[1] for i in y]
	z = [i+trans[2] for i in z]

	if output is None:
		scatter(x,y)
		show()
	elif output == 'std':
		for i in range(0,len(x)):
			print '%0.10f %0.10f %0.10f' % (x[i], y[i], z[i])
	elif output == 'ret':
		return x,y,z
	else:
		f = open(output, 'w')
		f.write('3\n')
		f.write('%i\n' % len(x))
		for i in range(0,len(x)):
			f.write('%0.10f %0.10f %0.10f\n' % (x[i], y[i], z[i]))
		f.close

def make_3d_space(anos, acoords, indices=None, fname=None):
	"""
	For a sequence of atoms, develop their relative radii based on tabulated
	values and then extend an array of each atom with the volumes
	This information will then be used in the QHULL application to make a 3D
	convex hull to determine atomic volumes
	"""

	# Atomic diameter computed using quantum mechanical calculations,
	# Periodic Chart of the Atoms (1979), Sargent-Welch

	atomsize = [0,1.58,0.98,4.10,2.80,2.34,1.82,1.50,1.30,1.14,1.02,4.46,3.44,3.64,2.92,2.46,2.18,1.94,1.76,5.54,4.46,4.18,4.00,3.84,3.70,3.58,3.44,3.34,3.24,3.14,3.06,3.62,3.04,2.66,2.44,2.24,2.06,5.96,4.90,4.54,4.32,4.16,4.02,3.90,3.78,3.66,3.58,3.50,3.42,4.00,3.44,3.06,2.84,2.64,2.48,6.68,5.56,5.48,4.32,4.18,4.04,3.94,3.84,3.74,3.66,3.58,3.52,4.16,3.62,3.26,3.06,2.86,2.68]

	d3x = []
	d3y = []
	d3z = []
#	color = []
	if indices is not None:
		an=[]
		ac=[]
		if isinstance(indices,int):
			anos = anos[indices]
			acoords = acoords[indices]
		else:
			for i in indices:
				an.append(anos[i])
				ac.append(acoords[i])
			acoords = ac
			anos = an

	if isinstance(anos, int):
		r = atomsize[anos]/2
		trans = acoords
		print r
		print trans
		x, y, z = make_sphere(r,output='ret',trans=trans)

		d3x.extend(x)
		d3y.extend(y)
		d3z.extend(z)

	else:

		for i in range(0,len(anos)):
			r = atomsize[anos[i]]/2
			trans = acoords[i]
			print r
			print trans
			x, y, z = make_sphere(r,output='ret',trans=trans)

			d3x.extend(x)
			d3y.extend(y)
			d3z.extend(z)
	#gain = [d*5 for d in d3z]
	#ps.set_mode('present')
	#ps.set_figsize_in(6,8.5)
	#scatter(d3x,d3y,s=gain, c=color)
	#title('Projected Atoms in 3D Space')
	#xlabel('X Coord ($\AA$)')
	#ylabel('Y Coord ($\AA$)')
	#show()


	if fname is None:
		for i in range(0,len(d3x)):
			print '%0.10f %0.10f %0.10f' % (d3x[i], d3y[i], d3z[i])

	else:
		f = open(fname, 'w')
		f.write('3\n')
		f.write('%i\n' % len(d3x))
		for i in range(0,len(d3x)):
			f.write('%0.10f %0.10f %0.10f\n' % (d3x[i], d3y[i], d3z[i]))
		f.close
