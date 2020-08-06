import numpy as np
import matplotlib.pyplot as plt
from utilities import *
from scipy.ndimage.filters import gaussian_filter


root = "/Volumes/External/JanData/"
# range of snapshots
imax=359
imin=36


cmD='Purples'
cmS='Reds'
cmC='magma'

for i in range(imin,imax):
	print i
	phase = np.load(root+"AxPhase_%03d.npy"%i)	
	density = np.load(root+"AxDensity4_%03d.npy"%i)	

	shape = np.shape(phase)
	# Gradient of phase -> velocity
	vel = dTheta(phase,np.sqrt(density))
	vx,vy,vz = vel[0],vel[1],vel[2]
	v_empty = np.ones(shape)
	# Soliton location and velocity
	maxdens = density.max()
	centre = np.where(density==maxdens)
	x0,y0,z0 = centre[0][0],centre[1][0],centre[2][0]
	vx_sol = vx[centre]
	vy_sol = vy[centre]
	vz_sol = vz[centre]

	# The mean, mass weighted motion
	#vx_av = np.average(vx*density)/np.average(density)
	#vy_av = np.average(vy*density)/np.average(density)
	#vz_av = np.average(vz*density)/np.average(density)
		
	# Subtract soliton velocity
	vx = vx-vx_sol*v_empty
	vy = vy-vy_sol*v_empty
	vz = vz-vz_sol*v_empty
	
	####
	# Curl of v
	####
	
	curl_v = Curl(vel)
	curlx,curly,curlz = curl_v[0],curl_v[1],curl_v[2]
	
	zval=z0
	cube=np.shape(phase)
	Lx,Ly,Lz=cube[0],cube[1],cube[2]
	x,y,z=np.arange(Lx),np.arange(Ly),np.arange(Lz)

	##########
	# Speed
	########
	
	sig = 1
	
	fig,ax=plt.subplots()
	speed=np.sqrt(vx[:,zval,:]**2.+vy[:,zval,:]**2.+vz[:,zval,:]**2.)
	val = gaussian_filter(speed, sigma=sig)
	val=np.log10(val)
	
	pmesh = ax.pcolormesh(x, z, val, 
    	cmap = cmS,vmin=-2,vmax=1)
	plt.axis([65, 135, 65, 135])
	#plt.axis([x.min(), x.max(), y.min(), y.max()])
	cbar = fig.colorbar(pmesh)
	cbar.ax.set_ylabel('Speed')
	plt.savefig('Plots/SpeedSliceMovie/BlurSpeed_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()
	
	##########
	# Speed + standard dev
	########
	fig,ax=plt.subplots()
	speed=np.sqrt(vx[:,zval,:]**2.+vy[:,zval,:]**2.+vz[:,zval,:]**2.)
	av_speed = gaussian_filter(speed, sigma=sig)
	var_kern = (speed-av_speed)**2.
	var_speed = gaussian_filter(var_kern,sigma=sig)
	var_speed = np.sqrt(var_speed)
	val=np.log10(var_speed+av_speed)
	
	pmesh = ax.pcolormesh(x, z, val, 
    	cmap = cmS,vmin=-2,vmax=1)
	plt.axis([65, 135, 65, 135])
	#plt.axis([x.min(), x.max(), y.min(), y.max()])
	cbar = fig.colorbar(pmesh)
	cbar.ax.set_ylabel('Speed')
	plt.savefig('Plots/SpeedSliceMovie/AltSpeed_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()
	
	
	
	##########
	# Curl
	########
	fig,ax=plt.subplots()
	speed=np.sqrt(curlx[:,zval,:]**2.+curly[:,zval,:]**2.+curlz[:,zval,:]**2.)
	pmesh = ax.pcolormesh(x, z, np.log10(speed), 
    	cmap = cmC,vmin=-2,vmax=1)
	plt.axis([65, 135, 65, 135])
	#plt.axis([x.min(), x.max(), y.min(), y.max()])
	cbar = fig.colorbar(pmesh)
	cbar.ax.set_ylabel('Curl')
	plt.savefig('Plots/CurlSliceMovie/Curl_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()
	
	
	#########
	# Density
	########
	fig,ax=plt.subplots()
	pmesh = ax.pcolormesh(x, y, np.log10(density[:,:,zval]), 
    	cmap = cmD,vmin=0,vmax=7)
	plt.axis([65, 135, 65, 135])
	cbar = fig.colorbar(pmesh)
	cbar.ax.set_ylabel('Density')
	plt.savefig('Plots/DensitySliceMovie/Den_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()

	########
	# Quiver velocity 
	########
	
	Qscale=9.  # Smaller number = larger arrows
	Qangle='xy'
	Qwidth=0.02 

	fig,ax=plt.subplots()
	X,Y = np.meshgrid(np.arange(Lx),np.arange(Ly))
	Q = plt.quiver(X[::3, ::3], Y[::3, ::3], vx[::3, ::3,zval], vy[::3,::3,zval],
	               pivot='mid', units='inches',angles=Qangle,scale=Qscale,width=Qwidth)
	plt.axis([65, 135, 65, 135])
	
	plt.savefig('Plots/QuiverSliceMovie/Quiver_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()

	########
	# Mix
	########

	fig,ax=plt.subplots()
	X,Z = np.meshgrid(np.arange(Lx),np.arange(Lz))
	
	pmesh = ax.pcolormesh(x, z, np.log10(density[:,zval,:]), 
			   cmap = cmD,vmin=0,vmax=7)
	
	X,Y = np.meshgrid(np.arange(Lx),np.arange(Ly))
	Q = plt.quiver(X[::3, ::3], Y[::3, ::3], vx[::3, ::3,zval], vy[::3,::3,zval],
	               pivot='mid', units='inches',angles=Qangle,scale=Qscale,width=Qwidth)
	
	plt.axis([65, 135, 65, 135])
	cbar = fig.colorbar(pmesh)
	cbar.ax.set_ylabel('Density')
	plt.savefig('Plots/MixSliceMovie/Mix_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()
