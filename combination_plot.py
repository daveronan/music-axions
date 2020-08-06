import numpy as np
import matplotlib.pyplot as plt
from utilities import *

# cmap for scatter speaker points
cm = plt.get_cmap('jet')

root = "/Volumes/External/JanData/"

# Speaker channels	
# First index: speaker
# Second index: variable (see "speaker_channels.py" for what is used)
# Third index: time step
channels = np.load("Data/channels_R10.npy")

# Speaker parameters, these should match what is given in "speaker_channels.py"
nspeakers=8
R=10 # arbitrary radius for now
Theta = np.pi/2 # In plane
phivec=np.linspace(0.,(2.*np.pi)*(1.-1./nspeakers),nspeakers)
# range of snapshots
imax=359
imin=36
tsteps = imax-imin

# Make plots

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
	# Subtract soliton velocity
	vx = vx-vx_sol*v_empty
	vy = vy-vy_sol*v_empty
	vz = vz-vz_sol*v_empty
	
	cube=np.shape(phase)
	Lx,Ly,Lz=cube[0],cube[1],cube[2]
	x,y,z=np.arange(Lx),np.arange(Ly),np.arange(Lz)
	
	zval = z0
	zind = np.where(z==z0)
	zind = zind[0][0]
# Speaker locations
	cart = cartesian(R,Theta,phivec,x0,y0,z0)
	xS,yS,zS = cart[0],cart[1],cart[2]

	

#################
# Plot
###############

	fig, ((ax1, var1), (var2, var3)) = pl.subplots(2, 2,figsize=(9,7))
	
	# in ax1 show the "raw" data of density and velocity arrows, overlay speaker locations
	# in var1, var2, var3, show three of the variables and trace them over time
	
	# Density mesh
	pmesh = ax1.pcolormesh(np.log10(density[:,:,zind]), 
    	cmap = 'Purples',vmin=0,vmax=7)
		
	# Velocity quiver
	Qscale=7.  # Smaller number = larger arrows
	Qangle='xy'
	Qwidth=0.02
	
	X,Y = np.meshgrid(np.arange(Lx),np.arange(Ly))	
	Q = ax1.quiver(X[::3, ::3], Y[::3, ::3], vx[::3, ::3,zind], vy[::3, ::3,zind], 
		pivot='mid', units='inches',angles=Qangle,scale=Qscale,width=Qwidth)
	
	

	for j in range(nspeakers):
		cmap_num=j*(1./nspeakers)
		# Scatter speaker locations on density field
		# The density field has x and y labels reversed in location for this plot
		# I have verified this is just a plotting problem and not a problem for the speaker channels
		ax1.scatter(yS[j],xS[j],color=cm(cmap_num))
		ax1.scatter(y0,x0,marker="*",facecolor='r')
		x_scat = 1.*(i-imin)
		##########
		# Plot the three variables sent to each speaker
		# Line of whole time series, and moving coloured dot
		# Offest y-axis by fixed amount
		##########
		# Log10 Density	
		offset = 3.	
		var1.plot(channels[j,0,:]+offset*j,'-k',zorder=1)	
		y1_scat = channels[j,0,i-imin]+offset*j
		var1.scatter(x_scat,y1_scat,color=cm(cmap_num),zorder=2)
		# Log10 Speed	
		offset = 3.		
		var2.plot(channels[j,1,:]+offset*j,'-k',zorder=1)	
		y2_scat = channels[j,1,i-imin]+offset*j
		var2.scatter(x_scat,y2_scat,color=cm(cmap_num),zorder=2)
		# Log10 Curl
		offset = 3.			
		var3.plot(channels[j,2,:]+offset*j,'-k',zorder=1)	
		y3_scat = channels[j,2,i-imin]+offset*j
		var3.scatter(x_scat,y3_scat,color=cm(cmap_num),zorder=2)
	
	# Set axes
	ax1.axis([65, 135, 65, 135])
	ax1.set_yticks([])
	ax1.set_xticks([])
	
	#var1.axis([0,tsteps,0,30])
	var1.set_yticks([])
	var1.set_xticks([])
	
	#var2.axis([0,tsteps,-10,90])
	var2.set_yticks([])
	var2.set_xticks([])
	
	#var3.axis([0,tsteps,0,30])
	var3.set_yticks([])
	var3.set_xticks([])
	
	plt.savefig('Plots/CombinationMovie/Combination_%03d.png'%i,bbox_inches='tight')
	plt.clf()
	plt.close()