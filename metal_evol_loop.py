import matplotlib as mpl
mpl.use('agg')

from lauren import *


YEAR = 3.155693e+7 #s/yr

pf = load('/u/10/l/lnc2115/home/WiseSimsData/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/home/WiseSimsData/groups_02797.dat')

count = 0
halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]

while (count < len(halos_to_run)):
	halonum = halos_to_run[count]
	print 'Halonum is ',halonum
	radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	data_source = pf.h.sphere(center,radius/pf['cm'])
	
	fileout = '/u/10/l/lnc2115/home/WiseSims/analysis_plots/MDF/halo'+str(halonum)+'_stellarMDF.png'

	sm, ct, metals = sfr_quants(data_source)
	plot_stellar_t_vs_Z(pf,ct,metals,fileout)
	
	fileout2 = '/u/10/l/lnc2115/home/WiseSims/analysis_plots/MDF/halo'+str(halonum)+'_gasMDF.png'
	if halonum != 22:
		halo_metals = data_source['Metallicity']
		metallicity_hist(metals,25,fileout2)

	count = count + 1


