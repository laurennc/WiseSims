#import matplotlib as mpl
#mpl.use('agg')
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from myanyl import *
import cPickle

#halos_to_run = [59, 188, 362, 900]
#halos_to_run = [0,3]

halos_to_run = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]
count = 0

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')

while (count < len(halos_to_run)):
	halonum = halos_to_run[count]

	radius, mass, center = r200(pf, halos["pos"][halonum,:],halos["mass"][halonum],verbose=False)
	field = "Metallicity"
	step = 2.0
	#changed from 25 to see if this is where that limit was coming from
	data_source = pf.h.sphere(center,30.0/pf['kpc'])
	#c_min = 1.0e-17
	#c_max = 1.0e-15
	c_min = 5.623e-7
	c_max = 1.778e-6
	function = 'self.data["Metallicity"].size > 5e+3'
	#function = 'self.data["Metallicity"].size > 3e+2'

	master_clump = Clump(data_source,None, field, function=function)
	find_clumps(master_clump,c_min,c_max,step)
	#all_clumps = get_lowest_clumps(master_clump)

	fileout = 'halo'+str(halonum)+'_clumps_Z6.cpkl'
	cPickle.dump(master_clump,open(fileout,"wb"),protocol=-1)

	count = count + 1
	#halonum = halonum+1



