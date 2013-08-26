import matplotlib as mpl
mpl.use('agg')
from lauren import *
import triangle
from matplotlib.artist import *

halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,
                36,42,48,59,188,362,900]
count = 0

pf = load('/u/10/l/lnc2115/home/WiseSimsData/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/home/WiseSimsData/groups_02797.dat')
data = cPickle.load(open('clump_dict.cpkl','rb'))

while (count < len(halos_to_run)):
	halonum = halos_to_run[count]
	print 'Analayzing halo ',halonum
	if halonum != 22:
		master_clump, radius, mass, center = readin_masterclump(pf,halos,halonum)
		data_source = pf.h.sphere(center,radius/pf['cm'])
		metals = make_solar_metallicity(data_source['Metallicity'])
		dists = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)
		
		datain = np.zeros((len(dists),2))
		datain[:,0] = dists*pf['kpc']
		datain[:,1] = metals
	
	#from matplotlib.artist import *
	#from matplotlib.pyplot import *

	#figure = triangle.corner(datain,labels=['Distance from Center (kpc)','Metallicity'])
		fileout = 'analysis_plots/metal_tri/metal_dist_halo'+str(halonum)+'_all.png'
	#figure.savefig(fileout)
		triangle.corner(datain,labels=['Distance from Center (kpc)','Metallicity']).savefig(fileout)

		idx = (metals > -5)
		metals = metals[idx]
		dists = dists[idx]
		datain = np.zeros((len(metals),2))
		datain[:,0] = dists*pf['kpc']
		datain[:,1] = metals
	
		figure = triangle.corner(datain,labels=['Distance from Center (kpc)','Metallicity'])
		fileout = 'analysis_plots/metal_tri/metal_dist_halo'+str(halonum)+'_gt5.png'
		figure.savefig(fileout)

	#sm,ct,metal = sfr_quants(data_source)	
	
	count = count + 1
