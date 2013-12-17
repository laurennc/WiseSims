from lauren import *

def create_structures(halonum):

	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
	
	radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	data_source = pf.h.sphere(center,radius/pf['cm'])

	return pf, data_source, radius, center

def slices_with_clumps_and_stars(halonum,clump_index,pf,center,radius)
	filein = '/u/10/l/lnc2115/vega/data/Wise/pickles/halo'+str(halonum)+'_clumps.cpkl'
        data = cPickle.load(open(filein,'rb'))

        master_clump = data[1]
        all_clumps = get_lowest_clumps(master_clump)


	clump_contours(pf,center,radius,'SolarMetals','z',all_clumps['clump_index'],30.,'halo'+str(halonum),'slice')
	#clump index for halo 0 is 3!
	#clump index for halo 3 is 5!
	#they're both the last clump in the list all_clumps so len-1

	sp = SlicePlot(pf,'z','SolarMetals',center,width=(30.,'kpc'),axes_unit=["kpc","kpc"])
	sp.annotate_particles(0.05,p_size=2.5,stars_only=True)
	sp.set_zlim('all',-18,0)
	sp.set_cmap('all','spectral')
	sp.save('halo'+str(halonum)+'_stars')
	return


def 


