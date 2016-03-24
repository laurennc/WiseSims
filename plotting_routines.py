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


	clump_contours(pf,center,radius,'SolarMetals','x',all_clumps[clump_index],30.,'halo'+str(halonum),'slice')
	#clump index for halo 0 is 3!
	#clump index for halo 3 is 5!
	#they're both the last clump in the list all_clumps so len-1

	sp = SlicePlot(pf,'x','SolarMetals',center,width=(30.,'kpc'),axes_unit=["kpc","kpc"])
	sp.annotate_particles(0.05,p_size=2.5,stars_only=True)
	sp.set_zlim('all',-18,0)
	#sp.set_cmap('all','spectral')
	sp.set_cmap('all','cool')
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

def plot_factors(halonum,factorname,fileout,fileout2):
	#data = cPickle.load(open('halo'+str(halonum)+'.cpkl','rb'))
	data0 = cPickle.load(open('halo0.cpkl','rb'))
	data3 = cPickle.load(open('halo3.cpkl','rb'))
	want = 6
	count = 0
	while count < len(data0['redshift']):
		plt.plot(np.log10(np.zeros(want)+data0['time'][count]),data0[factorname][count][:6],'s')
		plt.plot(np.log10(np.zeros(want)+data3['time'][count]),data3[factorname][count][:6],'^')
		count = count + 1
	plt.ylabel(factorname)
	plt.xlabel('log(time) (yrs)')
	#plt.xscale('log')
	#plt.xlim(7e8,8e8)
	plt.xlim(8.85,8.875)
	plt.savefig(fileout)
	plt.close()

	count = 0
	while count < len(data0['redshift']):
                plt.plot(np.log10(np.zeros(want)+data0['time'][count]),(data0[factorname][count][:6]-data3[factorname][count][:6]),'o')
                count = count + 1
        plt.ylabel(factorname+' Difference')
        plt.xlabel('log(time) (yrs)')
        #plt.xscale('log')
        #plt.xlim(7e8,8e8)
        plt.xlim(8.85,8.875)
        plt.savefig(fileout2)
        plt.close()
        

	return 'plotted'

def basic_slice_plot(pf,center,radius,dim,fileout):
	sp = SlicePlot(pf,dim,'SolarMetals',center,width=(radius,'cm'),axes_unit=["kpc","kpc"])
        sp.annotate_particles(0.05,p_size=3.0,stars_only=True)
        #sp.set_zlim('all',-18,0)
        #sp.set_cmap('all','spectral')
        sp.set_cmap('all','cool')
	sp.save(fileout)
	return

def plot_tegmark_radii(pf,data,idx,field,dim,width,tegradii,fileout):
	radius = data['rvirs'][idx]/pf['cm']
	center = data['centers'][idx]
	pc = PlotCollection(pf,center=center)
	pc.add_slice(field,dim)
	pc.set_width(width,'kpc')

	pc.plots[0].modify['point'](center,'.')
	pc.plots[0].modify['point'](data['CofMs'][idx],'x')
	
	pc.plots[0].modify['sphere'](center,radius,circle_args={'lw':2.0,'color':'white'})
	pc.plots[0].modify['sphere'](center,tegradii[0]/pf['kpc'],circle_args={'ls':'dashed','lw':2.0,'color':'white'})
	pc.plots[0].modify['sphere'](center,tegradii[1]/pf['kpc'],circle_args={'ls':'dotted','lw':2.0,'color':'white'})
	pc.plots[0].modify['sphere'](center,tegradii[2]/pf['kpc'],circle_args={'ls':'dotted','lw':2.0,'color':'white'})
	
	dist = distance_from_center(data['centers'][:,0],data['centers'][:,1],data['centers'][:,2],center)
	idc = np.where(dist*pf['kpc'] <= tegradii[2])[0]
	for val in idc:
		#pc.plots[0].modify['sphere'](data['centers'][val],data['rvirs'][val]/pf['kpc'],circle_args={'ls':'dashdot','color':'green','lw':2.0})	
		pc.plots[0].modify['point'](data['centers'][val],str(data['halonum'][val]))

	print 'all marked'
	pc.set_zlim(-6,0)
	pc.set_cmap('spectral')
	pc.save(fileout)
	return

def mark_nearby_halos(data,plotter,halonum,maxrad):
	#want to mark nearby halos and their virial radii on the slice plots
	#If I do it this way, I'm going to be plotting the host halo twice which I don't think should be a problem but noted here in case
	dist = distance_from_center(data['centers'][:][0],data['centers'][:][1],data['centers'][:][2],data['centers'][halonum])
	idx = np.where(dist <= maxrad)[0]
	for i in idx:
		plotter.modify['sphere'](data['centers'][i],data['rvirs'][i],circle_args={'ls':'dashdot','color':'green','lw':2.0})
		plotter.modify['point'](data['centers'][i],'c')
	return


def plotting_surrounding_densities():
	clump_data = cPickle.load(open('clump_dict_Z6.cpkl','rb'))
	fig,ax = plt.subplots(1,2,sharey=True)
	ax[0].plot(np.log10(clump_data['star_masses']),np.log10(clump_data['avgDens_avgR']),'bo')
	ax[1].plot(np.log10(clump_data['star_masses']),np.log10(clump_data['avgDens_maxR']),'ro')
	
	plt.savefig('density_outside_bubble_estimate.png')
	plt.close()
	return
