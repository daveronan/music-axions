import numpy as np
import pylab as pl

def dTheta(phase,amplitude):
	"""
	Correctly compute dTheta for a complex phase.
	z=a+1j*b
	"""
	a = amplitude*np.cos(phase)
	b = amplitude*np.sin(phase)
	return (a**2./(b**2.+a**2.))*(np.gradient(b)/a-b*np.gradient(a)/a**2.)
	


def polar(field,r,theta,phi,x0,y0,z0):
	'''
	Return a scalar field defined on (x,y,z) cartesian with origin [x0,y0,z0] in spherical polars at (r,theta,phi)
	'''
	x=r*np.sin(theta)*np.cos(phi)+x0
	y=r*np.sin(theta)*np.sin(phi)+y0
	z=r*np.cos(theta)+z0
	return field([x,y,z])
	
polar=pl.vectorize(polar)


def cartesian(r,theta,phi,x0,y0,z0):
	'''
	Return x,y,z from spherical polars r,theta,phi with origina x0,y0,z0
	'''
	x=r*np.sin(theta)*np.cos(phi)+x0
	y=r*np.sin(theta)*np.sin(phi)+y0
	z=r*np.cos(theta)+z0
	return x,y,z
cartesian = pl.vectorize(cartesian)


def Curl(v):
	"""
	curl of discrete vector field v in R3
	v is a numpy array of shape (3,:,:,:)
	"""	
	vx,vy,vz = v[0],v[1],v[2]
	dv_x,dv_y,dv_z = np.gradient(vx),np.gradient(vy),np.gradient(vz)
	curl_x = dv_z[1]-dv_y[2]
	curl_y = dv_x[2]-dv_z[0]
	curl_z = dv_y[0]-dv_x[1]
	return np.asarray([curl_x,curl_y,curl_z])
