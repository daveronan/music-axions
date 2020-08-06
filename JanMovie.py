# Not saure about indentation of the for loop

import numpy as np
from scipy import stats
import numexpr as ne
import yt
import math
root = "/Volumes/External/JanData/"
array = np.load(root+"AxDensity4_036.npy")
maxdens = array.max()
for i in range(36,365):
	print i
	array = np.load(root+"AxDensity4_%03d.npy"%i)
	data = dict(density=array)
	bbox = np.array([[0., 1.0], [0, 1.0], [0., 1.0]])
	ds = yt.load_uniform_grid(data, array.shape,length_unit="cm",bbox=bbox)
	ad = ds.all_data()
	sc = yt.create_scene(ds, field="density")
	sc.camera.width = (0.3,'cm')
	alpha = i/(365.0-36.0)*2*math.pi
	sc.camera.position = np.array([0.5,0.5,0.5])+np.array([np.cos(alpha),0,np.sin(alpha)])
	sc.camera.focus = np.array([0.5,0.5,0.5])
	sc.camera.north_vector = np.array([0,1.0,0])
# sc.camera.focus = np.array(np.unravel_index(array.argmax(),array.shape))/406.0

	source =sc[0]
	bounds = (maxdens/1e3, maxdens)

	tf = yt.ColorTransferFunction(np.log10(bounds))
	tf.sample_colormap(np.log10(maxdens*0.75), w=.005,colormap='RAINBOW')
	tf.sample_colormap(np.log10(maxdens/2.0e1), w=.005, alpha =
	0.2,colormap='RAINBOW')
	tf.sample_colormap(np.log10(maxdens/1.0e2), w=.005, alpha =
	0.2,colormap='RAINBOW')

	source.tfh.tf = tf
	source.tfh.bounds = bounds


	sc.save(root+'DensityMovie/rotrendering%03d.png'%i,sigma_clip=4.0)