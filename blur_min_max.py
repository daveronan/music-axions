import numpy as np

name_root = ["channels_S8_R1_","channels_S8_R10_","channels_S8_R30_","channels_S10_R1_","channels_S10_R10_","channels_S10_R30_"]


for i in range(6):

		
		name = name_root[i]+'LY'
		dat=np.load("Data/"+"blur_"+name+'.npy')
		dens = dat[1:,0,:]
		vel = dat[1:,1,:]
		cur = dat[1:,2,:]
		avel = dat[1:,3,:]
		dmin,dmax = dens.min(),dens.max()
		vmin,vmax = vel.min(),vel.max()
		cmin,cmax = cur.min(),cur.max()
		avmin,avmax = avel.min(),avel.max()
		name = 'blur_min_max_'+name+'.csv'
		vals = np.asarray([dmin,dmax,vmin,vmax,cmin,cmax,avmin,avmax])
		print name
		print vals
		np.savetxt("Data/"+name, vals, delimiter=",")
		
		name = name_root[i]+'LN'
		dat=np.load("Data/"+"blur_"+name+'.npy')
		dens = dat[1:,0,:]
		vel = dat[1:,1,:]
		cur = dat[1:,2,:]
		avel = dat[1:,3,:]
		dmin,dmax = dens.min(),dens.max()
		vmin,vmax = vel.min(),vel.max()
		cmin,cmax = cur.min(),cur.max()
		avmin,avmax = avel.min(),avel.max()
		name = 'blur_min_max_'+name+'.csv'
		vals = np.asarray([dmin,dmax,vmin,vmax,cmin,cmax,avmin,avmax])
		print name
		print vals
		np.savetxt("Data/"+name, vals, delimiter=",")
			
			
