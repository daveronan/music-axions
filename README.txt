Remove phase gradients:

For the gradient you can problably use numpy.gradient()
then you have to account for transitions from 2*pi to 0 which you can do with numpy.mod()
then multiply the gradients with the density array
then take the mean() of that array divided by mean() of the density in order to get the mass weighted mean gradient
the phases you want to substract you can construct with numpy.fromfunction (with gradient*index ) using the mean gradients
substract the constructed array from the Phase array
take numpy.mod() again to have the range 0 to 2*pi

------------------------

ssh marsh@login.gwdg.de

---------------------------
Jan Data: in numpy files

# Jan densities and phases

/home/uni09/cosmo/jveltma/Doddy

AxDensity4_xxx.npy ,  xxx = [036,365],[36,99]
AxPhase_xxx.npy ,  xxx=[035,359]

Copied into Volumes/External/JanData

----------------------
Bodo Data: best accessed with yt

/home/uni09/cosmo/bschwab/multimerge_mass9_halos10_64kpc/ # appears corrupted. Hope we didn't fuck it up

/home/uni09/cosmo/bschwab/multimerge_mass10_halos10_8kpc # another run


* get the whole repository

scp -r marsh@login.gwdg.de:/home/uni09/cosmo/bschwab/multimerge_mass9_halos10_64kpc path_to_save

* get a single snapshot

scp -r marsh@login.gwdg.de:/home/uni09/cosmo/bschwab/multimerge_mass9_halos10_64kpc/plot##### path_to_save

plot11801

plot##### is a directory. yt read the whole directory as a HDF5 file


* Find fields

ds = yt.load("IsolatedGalaxy/galaxy0030/galaxy0030")
ds.field_list

For example, AxDens axion density field

