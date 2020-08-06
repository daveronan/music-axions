import numpy as np

name_root = ["channels_S8_R1_","channels_S8_R10_","channels_S8_R30_","channels_S10_R1_","channels_S10_R10_","channels_S10_R30_"]
nspeakers1=8
nspeakers2=10


for i in range(6):

		
		name = name_root[i]+'LY'
		dat=np.load("Data/"+name+'.npy')
		dens = dat[1:,0,:]
		vel = dat[1:,1,:]
		cur = dat[1:,2,:]
		dmin,dmax = dens.min(),dens.max()
		vmin,vmax = vel.min(),vel.max()
		cmin,cmax = cur.min(),cur.max()
		name = 'min_max_'+name+'.csv'
		vals = np.asarray([dmin,dmax,vmin,vmax,cmin,cmax])
		print name
		print vals
		np.savetxt("Data/"+name, vals, delimiter=",")
		
		name = name_root[i]+'LN'
		dat=np.load("Data/"+name+'.npy')
		dens = dat[1:,0,:]
		vel = dat[1:,1,:]
		cur = dat[1:,2,:]
		dmin,dmax = dens.min(),dens.max()
		vmin,vmax = vel.min(),vel.max()
		cmin,cmax = cur.min(),cur.max()
		name = 'min_max_'+name+'.csv'
		vals = np.asarray([dmin,dmax,vmin,vmax,cmin,cmax])
		print name
		print vals
		np.savetxt("Data/"+name, vals, delimiter=",")
			
			
