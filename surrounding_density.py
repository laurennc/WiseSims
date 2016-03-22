from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import matplotlib.pyplot as plt
from lauren import *
import numpy as np

YEAR = 3.155693e+7 #s/yr

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')

count = 0
clump_data = cPickle.load(open('clump_dict_Z6.cpkl','rb'))

avgDens_maxR = []
avgDens_avgR = []

while (count < len(clump_data['halonum'])):
	halonum = clump_data['halonum'][count]
	print 'Analyzing halo number: ',halonum

	#stellar mass of the halo within the virial radius
	data_source = pf.h.sphere(clump_data['centers'][count],clump_data['maxRadius'][count]/pf['kpc'])
	
	#use a tolerance of 5 pc within the radii
	radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],clump_data['centers'][count])*pf['kpc']
	
	idx = np.where((radii > (clump_data['avgRadius'][count]-0.05)) & (radii < (clump_data['avgRadius'][count]+0.05)))
	avgDens_avgR = np.append(avgDens_avgR,np.average(data_source['Density'][idx]))

	idx = np.where((radii > (clump_data['maxRadius'][count]-0.05)) & (radii < (clump_data['maxRadius'][count]+0.05)))
        avgDens_maxR = np.append(avgDens_maxR,np.average(data_source['Density'][idx]))     
	
	
	count = count + 1

clump_data['avgDens_maxR'] = avgDens_maxR
clump_data['avgDens_avgR'] = avgDens_avgR

cPickle.dump(clump_data,open('clump_dict_Z6.cpkl','wb'),protocol=-1)

fig,ax = plt.subplots(1,2,sharey=True)
ax[0].plot(clump_data['star_masses'],np.log10(clump_data['avgDens_avgR']),'bo')
ax[1].plot(clump_data['star_masses'],np.log10(clump_data['avgDens_maxR']),'ro')

plt.savefig('density_outside_bubble_estimate.png')
plt.close()


