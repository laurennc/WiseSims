from lauren import *

def create_structures(halonum):

	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
	
	radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	data_source = pf.h.sphere(center,radius/pf['cm'])

	return pf, data_source, radius, center

def slices_with_clumps_and_stars(halonum,clump_index,pf,center,radius):
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


def sf_time_properties(pf,data_sphere,halonum,fileout):
	sfr = sfr_halo(pf,data_sphere)
	sm, ct, metals = sfr_quants(data_sphere)
	metals = np.log10(metals)-np.log10(0.02)
        #t = np.log10(pf['Time']*ct/YEAR)	
	t = pf['Time']*ct/YEAR

	fig,ax = plt.subplots(3,1,sharex=True)
	fig.set_size_inches(5,11)
	fig.subplots_adjust(hspace=0.1,wspace=0.45)

	ax[0].plot(sfr.time,np.log10(sfr.Msol_yr),'g-',linewidth=1.3)
	ax[0].set_ylabel('log(SFR) (Msol/yr)')
	ax[0].set_ylim(-5.0)
	#ax[0].set_xlim(8,9)
	ax[0].set_xscale('log')

	ax[1].plot(t,metals,'bo')
	ax[1].set_ylabel('Metallicity')
#        ax[1].set_xlim(8,9)
	ax[1].set_xscale('log')
	ax[2].plot(t,metals,'bo')
        ax[2].set_ylabel('Metallicity')
	ax[2].set_xlabel('log(time) (yrs)')
	ax[2].set_ylim(-4,1)
	ax[2].set_xlim(1e8,1e9)	
	ax[2].set_xscale('log')
	plt.savefig(fileout)
	plt.close()
	return 'plotted'

def plot_filling_factors(halonum,fileout):
	data = cPickle.load(open('halo'+str(halonum)+'.cpkl','rb'))
	want = 6
	count = 0
	while count < len(data['redshift']):
		plt.plot(np.log10(np.zeros(want)+data['time'][count]),data['fillFactors'][count][:6],'s')
		count = count + 1
	plt.ylabel('Filling Factors')
	plt.xlabel('log(time) (yrs)')
	#plt.xscale('log')
	#plt.xlim(7e8,8e8)
	plt.xlim(8.85,8.875)
	plt.savefig(fileout)
	plt.close()
	return 'plotted'

def basic_slice_plot(pf,center,radius,dim,fileout):
	sp = SlicePlot(pf,dim,'SolarMetals',center,width=(radius,'cm'),axes_unit=["kpc","kpc"])
        sp.annotate_particles(0.05,p_size=3.0,stars_only=True)
        sp.set_zlim('all',-18,0)
        sp.set_cmap('all','spectral')
        sp.save(fileout)
	return



