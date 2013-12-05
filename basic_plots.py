import matplotlib
matplotlib.use('Agg')

from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import matplotlib.pyplot as plt
from lauren import *

halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]
#halos_to_run = [22]
#counter = 13
counter = 0

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
field = 'Metallicity'

from numpy import *

data = cPickle.load(open('clump_dict.cpkl','rb'))

#while ( counter < 14) :
while (counter < len(halos_to_run)):
	halonum = halos_to_run[counter]
	#halonum = 22
	print 'Analyzing halo number: ',halonum

	#radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	filein = '../WiseSimsData/pickles/halo'+str(halonum)+'_clumps.cpkl'
	data2 = cPickle.load(open(filein,'rb'))
	master_clump = data2[1]
	all_clumps = get_lowest_clumps(master_clump)	

	#properties of virial halo
	center = data['centers'][counter]
	radius = data['rvirs'][counter]
	virial_sphere = pf.h.sphere(center,radius/pf['cm'])
	halo_metals = virial_sphere['Metallicity']

	index = int(data['clump_index'][counter])	
	#clump metallicities
	metallicities = all_clumps[index]['Metallicity']


	#Projection Plots
	fileout = 'analysis_plots/halonum'+str(halonum)
	for dim in "xyz":
		clump_contours(pf,center,radius,field,dim,all_clumps[index],25.0,fileout,'slice')
		clump_contours(pf,center,radius,'Density',dim,all_clumps[index],25.0,fileout,'slice')

	#SFH plots
	sfr = sfr_halo(pf,virial_sphere)
	
	fileout2 = fileout+'_sfh.png'
	sfh_plot(sfr,fileout2)

	print 'Clump Metallicities:'
	print metallicities
	print 'Halo Metallicities:'
	print halo_metals
	
	#Clump Metal Distribution
	fileout3 = fileout+'_zclump.png'
	if absolute(log(metallicities.min()) - log(metallicities.max())) > 20:	
		n, bins, patches = metallicity_hist(metallicities,25,fileout3)
	else:
		print 'Clump metallicities length = 0'
	#Halo Metal Distribution
	fileout4 = fileout+'_zhalo.png'
	#print halo_metals
	print absolute(log(halo_metals.min()) - log(halo_metals.max()))
	if absolute(log(halo_metals.min()) - log(halo_metals.max())) > 20:
		n2, bins2, patches = metallicity_hist(halo_metals,25,fileout4)
	else:
		print 'Halo metallicities length = 0'

	counter = counter + 1


