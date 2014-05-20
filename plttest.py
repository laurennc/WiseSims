import matplotlib
matplotlib.use('Agg')

from plotting_routines import *

#plot_factors(0,'fillFactors','fillFactorsAgain.png','fillFactorsAgainDiff.png')
#plot_factors(0,'massFactors','massFactorsAgain.png','massFactorsAgainDiff.png')

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
field = 'SolarMetals'
dim = 'x'
data = cPickle.load(open('clump_dict.cpkl','rb'))

#for i in range(len(data['halonum'])):
#for i in range(4):
#	tegradii = [data['avgRadius'][i],data['zMinRadius'][i],data['zMaxRadius'][i]]
#	fileout = 'halo'+str(data['halonum'][i])+'allradii'
#	width = tegradii[2]*2.0+4.0
#	print 'i is ', str(i)
#	plot_tegmark_radii(pf,data,i,field,dim,width,tegradii,fileout)


timesteps = read_in_timesteps('/u/10/l/lnc2115/vega/data/Wise/time_steps.dat')

halonum = 3

data_source = pf.h.sphere(data['centers'][halonum],data['rvirs'][halonum]/pf['cm'])
sm, ct, metal = sfr_quants(data_source)


plot_stellar_t_vs_Z(pf,ct,metal,'testing.png',mark_timesteps=True,timesteps=timesteps)




