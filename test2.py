import numpy as np
from utilities import *
import matplotlib.pyplot as plt
from scipy.interpolate import RegularGridInterpolator


# Test that the polar co-ordinates work and are centred on the maximum density

root = "/Volumes/External/JanData/"

nspeakers=8
R=10 # arbitrary radius for now
Theta = np.pi/2 # In plane
phivec=np.linspace(0.,(2.*np.pi)*(1.-1./nspeakers),nspeakers)

imax=359
#imin=36
imin=340
tsteps = imax-imin

for i in range(imin,imax):
	print i
	density = np.load(root+"AxDensity4_%03d.npy"%i)	

	
	maxdens = density.max()
	print maxdens
	centre = np.where(density==maxdens)
	x0,y0,z0 = centre[0][0],centre[1][0],centre[2][0]
	cube=np.shape(density)
	Lx,Ly,Lz=cube[0],cube[1],cube[2]
	x,y,z=np.arange(Lx),np.arange(Ly),np.arange(Lz)
	cart = cartesian(R,Theta,phivec,x0,y0,z0)
	xS,yS,zS = cart[0],cart[1],cart[2]
	
	rho = RegularGridInterpolator((x, y, z), density)
	print polar(rho,0.,Theta,0.,x0,y0,z0)
	