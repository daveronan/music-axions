from utilities import *
import numpy as np
from scipy.interpolate import RegularGridInterpolator
from scipy.ndimage.filters import gaussian_filter

nameOUTN1R1 = "blur_channels_S8_R1"
nameOUTN1R2 = "blur_channels_S8_R10"
nameOUTN1R3 = "blur_channels_S8_R30"
nameOUTN2R1 = "blur_channels_S10_R1"
nameOUTN2R2 = "blur_channels_S10_R10"
nameOUTN2R3 = "blur_channels_S10_R30"
nspeakers1=8
nspeakers2=10
Rvec = np.asarray([1.,10.,30.])



########################
n_var=4
Theta=np.pi/2.
phivec1=np.linspace(0.,(2.*np.pi)*(1.-1./nspeakers1),nspeakers1)
phivec2=np.linspace(0.,(2.*np.pi)*(1.-1./nspeakers2),nspeakers2)


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
outputN1R1 = np.zeros((nspeakers1+1,n_var,tsteps))
outputN1R2 = np.zeros((nspeakers1+1,n_var,tsteps))
outputN1R3 = np.zeros((nspeakers1+1,n_var,tsteps))
outputN2R1 = np.zeros((nspeakers2+1,n_var,tsteps))
outputN2R2 = np.zeros((nspeakers2+1,n_var,tsteps))
outputN2R3 = np.zeros((nspeakers2+1,n_var,tsteps))

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

	
	av_speed = gaussian_filter(speed, sigma=1)
	var_kern = (speed-av_speed)**2.
	var_speed = gaussian_filter(var_kern,sigma=1)
	var_speed = np.sqrt(var_speed)
	var_speed = av_speed + var_speed
	
	# Loop over speakers in the snapshot
	# Output is log of scalar fields interpolated to speaker locations
	
	
	rho = RegularGridInterpolator((x, y, z), density)
	speed = RegularGridInterpolator((x, y, z), speed)
	curl = RegularGridInterpolator((x, y, z), curl)
	
 
	av_speed = RegularGridInterpolator((x, y, z), av_speed)
	var_speed = RegularGridInterpolator((x, y, z), var_speed)

	
	
	####################################################
	
	for j in range(nspeakers1+1):
		# Main output
		if j==0:
			# 0th entry is centre of halo
			outputN1R1[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputN1R1[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputN1R1[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputN1R1[j,3,i-imin]=polar(var_speed,0,0,0,x0,y0,z0)
			
			
			outputN1R2[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputN1R2[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputN1R2[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputN1R2[j,3,i-imin]=polar(var_speed,0,0,0,x0,y0,z0)
			
			
			outputN1R3[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputN1R3[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputN1R3[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputN1R3[j,3,i-imin]=polar(var_speed,0,0,0,x0,y0,z0)
			
			
			
				
				
		else:
			outputN1R1[j,0,i-imin]=polar(rho,Rvec[0],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R1[j,1,i-imin]=polar(av_speed,Rvec[0],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R1[j,2,i-imin]=polar(curl,Rvec[0],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R1[j,3,i-imin]=polar(var_speed,Rvec[0],Theta,phivec1[j-1],x0,y0,z0)
			
			
			outputN1R2[j,0,i-imin]=polar(rho,Rvec[1],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R2[j,1,i-imin]=polar(av_speed,Rvec[1],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R2[j,2,i-imin]=polar(curl,Rvec[1],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R2[j,3,i-imin]=polar(var_speed,Rvec[1],Theta,phivec1[j-1],x0,y0,z0)
			
			
			outputN1R3[j,0,i-imin]=polar(rho,Rvec[2],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R3[j,1,i-imin]=polar(av_speed,Rvec[2],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R3[j,2,i-imin]=polar(curl,Rvec[2],Theta,phivec1[j-1],x0,y0,z0)
			outputN1R3[j,2,i-imin]=polar(var_speed,Rvec[2],Theta,phivec1[j-1],x0,y0,z0)
			
			

# Repeat for nspeakers2			
			
	for j in range(nspeakers2+1):
		# Main output
		if j==0:
			# 0th entry is centre of halo
			outputN2R1[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputN2R1[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputN2R1[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputN2R1[j,3,i-imin]=polar(var_speed,0,0,0,x0,y0,z0)
			

			outputN2R2[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputN2R2[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputN2R2[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputN2R2[j,3,i-imin]=polar(var_speed,0,0,0,x0,y0,z0)
			

			outputN2R3[j,0,i-imin]=polar(rho,0,0,0,x0,y0,z0)
			outputN2R3[j,1,i-imin]=polar(av_speed,0,0,0,x0,y0,z0)
			outputN2R3[j,2,i-imin]=polar(curl,0,0,0,x0,y0,z0)
			outputN2R3[j,3,i-imin]=polar(var_speed,0,0,0,x0,y0,z0)
			




		else:
			outputN2R1[j,0,i-imin]=polar(rho,Rvec[0],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R1[j,1,i-imin]=polar(av_speed,Rvec[0],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R1[j,2,i-imin]=polar(curl,Rvec[0],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R1[j,3,i-imin]=polar(var_speed,Rvec[0],Theta,phivec2[j-1],x0,y0,z0)
			

			outputN2R2[j,0,i-imin]=polar(rho,Rvec[1],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R2[j,1,i-imin]=polar(av_speed,Rvec[1],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R2[j,2,i-imin]=polar(curl,Rvec[1],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R2[j,3,i-imin]=polar(var_speed,Rvec[1],Theta,phivec2[j-1],x0,y0,z0)
			

			outputN2R3[j,0,i-imin]=polar(rho,Rvec[2],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R3[j,1,i-imin]=polar(av_speed,Rvec[2],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R3[j,2,i-imin]=polar(curl,Rvec[2],Theta,phivec2[j-1],x0,y0,z0)
			outputN2R3[j,3,i-imin]=polar(var_speed,Rvec[2],Theta,phivec2[j-1],x0,y0,z0)
			
		
				
				
				
##################################
# Save				

# Save CSV arrays
for j in range(nspeakers1+1):
	nameCSV = nameOUTN1R1+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputN1R1[j,:,:].T, delimiter=",")
	nameCSV = nameOUTN1R2+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputN1R2[j,:,:].T, delimiter=",")
	nameCSV = nameOUTN1R3+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputN1R3[j,:,:].T, delimiter=",")
	nameCSV = nameOUTN1R1+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputN1R1[j,:,:]).T, delimiter=",")
	nameCSV = nameOUTN1R2+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputN1R2[j,:,:]).T, delimiter=",")
	nameCSV = nameOUTN1R3+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputN1R3[j,:,:]).T, delimiter=",")


for j in range(nspeakers2+1):
	nameCSV = nameOUTN2R1+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputN2R1[j,:,:].T, delimiter=",")
	nameCSV = nameOUTN2R2+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputN2R2[j,:,:].T, delimiter=",")
	nameCSV = nameOUTN1R1+"_LN"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, outputN2R3[j,:,:].T, delimiter=",")
	nameCSV = nameOUTN2R1+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputN2R1[j,:,:]).T, delimiter=",")
	nameCSV = nameOUTN2R2+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputN2R2[j,:,:]).T, delimiter=",")
	nameCSV = nameOUTN1R1+"_LY"+"_"+str(j)+".csv"
	np.savetxt("Data/"+nameCSV, np.log10(outputN2R3[j,:,:]).T, delimiter=",")


# Save full array	
np.save("Data/"+nameOUTN1R1+"_LN"+".npy",outputN1R1)
np.save("Data/"+nameOUTN1R2+"_LN"+".npy",outputN1R2)
np.save("Data/"+nameOUTN1R3+"_LN"+".npy",outputN1R3)
np.save("Data/"+nameOUTN2R1+"_LN"+".npy",outputN2R1)
np.save("Data/"+nameOUTN2R2+"_LN"+".npy",outputN2R2)
np.save("Data/"+nameOUTN2R3+"_LN"+".npy",outputN2R3)	
np.save("Data/"+nameOUTN1R1+"_LY"+".npy",np.log10(outputN1R1))
np.save("Data/"+nameOUTN1R2+"_LY"+".npy",np.log10(outputN1R2))
np.save("Data/"+nameOUTN1R3+"_LY"+".npy",np.log10(outputN1R3))
np.save("Data/"+nameOUTN2R1+"_LY"+".npy",np.log10(outputN2R1))
np.save("Data/"+nameOUTN2R2+"_LY"+".npy",np.log10(outputN2R2))
np.save("Data/"+nameOUTN2R3+"_LY"+".npy",np.log10(outputN2R3))