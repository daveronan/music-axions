from utilities import *
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.filters import gaussian_filter

nameOUTR1 = "extra_channels_S10_R1"
nameOUTR2 = "extra_channels_S10_R10"
nameOUTR3 = "extra_channels_S10_R30"
nspeakers=10
ntheta=2
Rvec = np.asarray([1.,10.,30.])



########################
n_var=3
dphi = (2.*np.pi)*(1.-1./nspeakers)
# Generate two sets of channels at displaced theta values
Theta1=np.pi/2.
Theta2=np.pi/2.+dphi/10.
phivec=np.linspace(0.,dphi,nspeakers)

folder = "/Volumes/External/JanData/"
DRoot = "AxDensity4_"
TRoot = "AxPhase_"

# range of snapshots
imax=359
#imax=37
imin=36
tsteps = imax-imin

# generate empty output arrays
# First index: speaker
# Second index: velocity
# Third index: curl velocity
# Doubled for extra theta values
outputR1 = np.zeros((nspeakers+1,2*n_var,tsteps))
outputR2 = np.zeros((nspeakers+1,2*n_var,tsteps))
outputR3 = np.zeros((nspeakers+1,2*n_var,tsteps))

# Loop over snapshots
for i in range(imin,imax):
	print "working on snapshot  ", i-imin+1, "  of  ", imax-imin+1
	
	# Density, maximum defines co-ordinate origin
	density = np.load(folder+DRoot+"%03d.npy"%i)
	# Define co-ordinates
	maxdens = density.max()
	centre = np.where(density==maxdens)
	x0,y0,z0=centre[0][0],centre[1][0],centre[2][0]	
	cube=np.shape(density)
	Lx,Ly,Lz=cube[0],cube[1],cube[2]
	x,y,z=np.arange(Lx),np.arange(Ly),np.arange(Lz)

	# Gradient of phase -> velocity
	phase = np.load(folder+TRoot+"%03d.npy"%i)
	shape = np.shape(phase)
	vel = dTheta(phase,np.sqrt(density))
	vx,vy,vz = vel[0],vel[1],vel[2]
	v_empty = np.ones(shape)
	vx_sol = vx[centre]
	vy_sol = vy[centre]
	vz_sol = vz[centre]
	# Subtract soliton velocity
	vx = vx-vx_sol*v_empty
	vy = vy-vy_sol*v_empty
	vz = vz-vz_sol*v_empty
	
	# Scalar speed
	speed = np.sqrt(vx**2.+vy**2.+vz**2.)
	
	####
	# Curl of v
	####	
	curl_v = Curl(vel)
	curlx,curly,curlz = curl_v[0],curl_v[1],curl_v[2]
	# Scalar curl magnitude
	curl = np.sqrt(curlx**2.+curly**2.+curlz**2.)

	
	
	
	# Loop over speakers in the snapshot
	# Output is log of scalar fields interpolated to speaker locations
	
	
	rho = RegularGridInterpolator((x, y, z), density)
	speed = RegularGridInterpolator((x, y, z), speed)
	curl = RegularGridInterpolator((x, y, z), curl)
	av_speed = speed
 
	

	
	
	####################################################
	
	for j in range(nspeakers+1):
		if j==0:
			# 0th entry is centre of halo
			outputR1[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputR1[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputR1[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputR1[j,3,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputR1[j,4,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputR1[j,5,i-imin]=polar(curl,0,0,0,x0,y0,z0)



			outputR2[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputR2[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputR2[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputR2[j,3,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputR2[j,4,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputR2[j,5,i-imin]=polar(curl,0,0,0,x0,y0,z0)


			outputR3[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputR3[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputR3[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputR3[j,3,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputR3[j,4,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputR3[j,5,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			
		else:
			outputR1[j,0,i-imin]=polar(rho,Rvec[0],Theta1,phivec[j-1],x0,y0,z0)
			outputR1[j,1,i-imin]=polar(av_speed,Rvec[0],Theta1,phivec[j-1],x0,y0,z0)
			outputR1[j,2,i-imin]=polar(curl,Rvec[0],Theta1,phivec[j-1],x0,y0,z0)
			outputR1[j,3,i-imin]=polar(rho,Rvec[0],Theta2,phivec[j-1],x0,y0,z0)
			outputR1[j,4,i-imin]=polar(av_speed,Rvec[0],Theta2,phivec[j-1],x0,y0,z0)
			outputR1[j,5,i-imin]=polar(curl,Rvec[0],Theta2,phivec[j-1],x0,y0,z0)
		
			outputR2[j,0,i-imin]=polar(rho,Rvec[1],Theta1,phivec[j-1],x0,y0,z0)
			outputR2[j,1,i-imin]=polar(av_speed,Rvec[1],Theta1,phivec[j-1],x0,y0,z0)
			outputR2[j,2,i-imin]=polar(curl,Rvec[1],Theta1,phivec[j-1],x0,y0,z0)
			outputR2[j,3,i-imin]=polar(rho,Rvec[1],Theta2,phivec[j-1],x0,y0,z0)
			outputR2[j,4,i-imin]=polar(av_speed,Rvec[1],Theta2,phivec[j-1],x0,y0,z0)
			outputR2[j,5,i-imin]=polar(curl,Rvec[1],Theta2,phivec[j-1],x0,y0,z0)
		
			outputR3[j,0,i-imin]=polar(rho,Rvec[2],Theta1,phivec[j-1],x0,y0,z0)
			outputR3[j,1,i-imin]=polar(av_speed,Rvec[2],Theta1,phivec[j-1],x0,y0,z0)
			outputR3[j,2,i-imin]=polar(curl,Rvec[2],Theta1,phivec[j-1],x0,y0,z0)
			outputR3[j,3,i-imin]=polar(rho,Rvec[2],Theta2,phivec[j-1],x0,y0,z0)
			outputR3[j,4,i-imin]=polar(av_speed,Rvec[2],Theta2,phivec[j-1],x0,y0,z0)
			outputR3[j,5,i-imin]=polar(curl,Rvec[2],Theta2,phivec[j-1],x0,y0,z0)
			
			
		
				
				
				
##################################
# Save				

# Save CSV arrays
for j in range(nspeakers+1):
	nameCSV = nameOUTR1+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputR1[j,:,:].T, delimiter=",")
	nameCSV = nameOUTR2+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputR2[j,:,:].T, delimiter=",")
	nameCSV = nameOUTR3+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputR3[j,:,:].T, delimiter=",")
	nameCSV = nameOUTR1+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputR1[j,:,:]).T, delimiter=",")
	nameCSV = nameOUTR2+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputR2[j,:,:]).T, delimiter=",")
	nameCSV = nameOUTR3+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputR3[j,:,:]).T, delimiter=",")




# Save full array	
np.save("Data/"+nameOUTR1+"_LN"+".npy",outputR1)
np.save("Data/"+nameOUTR2+"_LN"+".npy",outputR2)
np.save("Data/"+nameOUTR3+"_LN"+".npy",outputR3)

np.save("Data/"+nameOUTR1+"_LY"+".npy",np.log10(outputR1))
np.save("Data/"+nameOUTR2+"_LY"+".npy",np.log10(outputR2))
np.save("Data/"+nameOUTR3+"_LY"+".npy",np.log10(outputR3))
