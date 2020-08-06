import numpy as np
import matplotlib.pyplot as plt

root = 'blur_channels_S10_R10_LY'

alldat = np.load('Data/'+root+'.npy')

fig,((ax1),(ax2),(ax3),(ax4)) = plt.subplots(4,1)

for i in range(1,11):
	
	d = alldat[i,0,:]
	v = alldat[i,1,:]
	c = alldat[i,2,:]
	av = alldat[i,3,:]
	ax1.plot(d,label=str(i))
	
	ax2.plot(v)
	ax3.plot(c)
	ax4.plot(av)
	ax1.legend()
	
plt.savefig('Data/Plots/'+root+'_plot.pdf',bbox_inches='tight')


plt.close()


v1 = alldat[1,1,:]
av1 = alldat[1,3,:]
v2 = alldat[2,1,:]
av2 = alldat[2,3,:]
plt.plot(v1,'-k')
plt.plot(av1,'--k')
#plt.plot(v2,'-r')
#plt.plot(av2,'--r')
plt.show()