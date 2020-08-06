import numpy as np

name_root = ["extra_channels_S10_R1_","extra_channels_S10_R10_","extra_channels_S10_R30_"]
nspeakers=10


for i in range(3):

		
		name = name_root[i]+'LY'
		dat=np.load("Data/"+name+'.npy')
		
		dens = dat[1:,0,:]
		vel = dat[1:,1,:]
		cur = dat[1:,2,:]
		
		
		dens2 = dat[1:,3,:]
		vel2 = dat[1:,4,:]
		cur2 = dat[1:,5,:]
		
		dens = np.append(dens,dens2)
		vel = np.append(vel,vel2)
		cur = np.append(cur,cur2)
		
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
		
		
		dens2 = dat[1:,3,:]
		vel2 = dat[1:,4,:]
		cur2 = dat[1:,5,:]
		
		dens = np.append(dens,dens2)
		vel = np.append(vel,vel2)
		cur = np.append(cur,cur2)
		
		dmin,dmax = dens.min(),dens.max()
		vmin,vmax = vel.min(),vel.max()
		cmin,cmax = cur.min(),cur.max()
		
		
		name = 'min_max_'+name+'.csv'
		
		vals = np.asarray([dmin,dmax,vmin,vmax,cmin,cmax])
		print name
		print vals
		np.savetxt("Data/"+name, vals, delimiter=",")
			
			
