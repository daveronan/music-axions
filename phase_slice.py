import numpy as np
import matplotlib.pyplot as plt
from phase_colormap import *

root = "/Volumes/External/JanData/"

const=np.random.uniform(-np.pi,np.pi) # phase defined up to arbitrary constant, don't be visually fooled!

for i in range(36,365):
	print i
	phase = np.load(root+"AxPhase_%03d.npy"%i)

	fig,ax=plt.subplots()

	cube=np.shape(phase)
	Lx,Ly,Lz=cube[0],cube[1],cube[2]
	x,y,z=np.arange(Lx),np.arange(Ly),np.arange(Lz)

	zval=100
	
	cm=hpluv_anglemap
	# there is something weird with plotting if I use x,y which cuts for some x??
	pmesh = ax.pcolormesh(x, z, (phase[:,zval,:]+const)/np.pi, 
	    cmap = cm, vmin=-1, vmax=1)
	plt.axis([x.min(), x.max(), y.min(), y.max()])
	cbar = fig.colorbar(pmesh)
	cbar.ax.set_ylabel('Phase [pi]')


	#plt.contourf(phase[:,:,100])
	plt.savefig('Plots/PhaseSliceMovie/Phase_slice%03d.png'%i,bbox_inches='tight')





