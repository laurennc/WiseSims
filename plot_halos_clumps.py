import matplotlib
matplotlib.use('Agg')

from yt.mods import *
from myanyl import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt
from lauren import *

#halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]
#halos_to_run = [8,11,13]
halos_to_run = [0,3]

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
field = 'Metallicity'

data = cPickle.load(open('clump_dict_Z6.cpkl','rb'))
counter = 0

while counter < len(halos_to_run):
	halonum = halos_to_run[counter]
	print 'Analyzing halo: ',halonum
	
	filein = "/u/10/l/lnc2115/vega/data/Wise/pickles/clumps/halo"+str(halonum)+"_clumps_Z6.cpkl"
	data2 = cPickle.load(open(filein,'rb'))
	master_clump = data2[1]
	all_clumps = get_lowest_clumps(master_clump)

	idx = np.where(data['halonum'] == halonum)[0][0]
	
	center = data['centers'][idx]
        radius = data['rvirs'][idx]
        virial_sphere = pf.h.sphere(center,radius/pf['cm'])
        index = int(data['clump_index'][idx])

	fileout = 'clump_plots/halo'+str(halonum)+'_allhalos_June16'

	clump_contours(pf,center,radius,'Metallicity','x',all_clumps[index],25.,fileout,'slice')

	#pc = PlotCollection(pf,center=center)
	#pc.add_slice('Metallicity','x')
	#pc.set_width(25.,'kpc')
	#pc.plots[0].modify['point'](center,'o',text_args={'color':'blue'})
	#pc.plots[0].modify['sphere'](center,radius,text_args={'color':'blue','linestyle':'--','linewidth':2.5})
	
	#for i in range(len(data['halonum'])):
   #		pc.plots[0].modify['point'](data['centers'][i],str(data['halonum'][i]))
   #		pc.plots[0].modify['sphere'](data['centers'][i],data['rvirs'][i])

#	pc.set_cmap('binary')
#	pc.save('clump_plots/halo'+str(halonum)+'_allhalos')
	
	counter = counter + 1




