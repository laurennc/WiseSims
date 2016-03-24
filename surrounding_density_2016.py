import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import cPickle
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *

YEAR = 3.155693e+7 #s/yr

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')

count = 0
clump_data = cPickle.load(open('clump_dict_Z6.cpkl','rb'))

dens_maxR = []
dens_minR = []

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')

count = 0
halos_to_run = [0,1,2,3,4,5,6,7,8,9,11,13,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]

while (count < len(halos_to_run)):
        halonum = halos_to_run[count]
        print 'Analyzing halo number: ',halonum

        radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
        filein = '/u/10/l/lnc2115/vega/data/Wise/pickles/clumps/halo'+str(halonum)+'_clumps_Z6.cpkl'
        data = cPickle.load(open(filein,'rb'))
        master_clump = data[1]
	
	#need to loop through the lowest clumps to find which has the smallest offset from the center of the halo
	num = 0
        comparing = 1e9
        index = 1000
	all_clumps = get_lowest_clumps(master_clump)

	index = int(clump_data['clump_index'][count])

	com = all_clumps[index].quantities['CenterOfMass']()
	radii = np.sqrt((com[0]-all_clumps[index]['x'])**2.0+(com[1]-all_clumps[index]['y'])**2.0+(com[2]-all_clumps[index]['z'])**2.0)
	
	imin = np.where(radii == radii.min())[0]
	imax = np.where(radii == radii.max())[0]

	dens_maxR = np.append(dens_maxR,np.log10(all_clumps[index]['Density'][imin]))
	dens_minR = np.append(dens_minR,np.log10(all_clumps[index]['Density'][imax]))

	count = count + 1


##the problem with this is that the minimum is just going to be the cell closest to the COM and not the minimum edge of the surface....
##so still need to be able to find the surface....


