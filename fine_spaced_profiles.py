from lauren import *

#ok what exactly do I want to do? 
#I want some way of showing how the metallicity changes at different radii...
#can look at some sort of radial profile for sure!
#instead of storing just one value want to build a table of them!
#instead of doing percents in each radius bin, can use a universally normalized quantity --- either just summing the mass directly or normalizing it to the total mass instead of the mass within the bin


def build_radial_profile_matrix(pf,data_source,center,rvirKPC,bins):
	#metals = np.arange(101)*-6./100.
	metals = np.linspace(0,-6,num=bins)
	radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)*pf['kpc']
	metallicities = make_solar_metallicity(data_source['Metallicity'])

	#WANT TO CALCULATE METALLICITY THAT WOULD BE USED IN SAMS
	sam_Z = make_solar_metallicity(np.sum(data_source['Metal_MassMsun'])/np.sum(data_source['CellMassMsun']))

	#metal_masses = data_source['Metal_MassMsun']
	#DOING VOLUMES BUT CALLED METAL_MASSES FOR EASE!
	#metal_masses = data_source['CellVolume']/(pf['cm']**3.0)*(pf['pc']**3.0)
	#DOING DENSITY BUT CALLED METAL_MASSES FOR EASE!
#	metal_masses = data_source['Density']
	#DOING TOTAL GAS MASS BUT CALLED METAL_MASSES FOR EASE!
	metal_masses = data_source['CellMassMsun']
	dr = rvirKPC/bins
	rp_r = np.arange(bins)*dr + dr/2.0
	rp_vals = np.zeros((len(metals),bins))
	total_mass = np.sum(metal_masses)
	for irad in range(int(bins)):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (radii>=minrad) * (radii<maxrad)
		idx = np.where(metallicities[thisindex] > 0.)[0]
		#LETS NORMALIZE BASED OFF THE TOTAL QUANTITY IN THAT ANNULUS
       		#total_mass = np.sum(metal_masses[thisindex])
		#LETS NOT NORMALIZE
		#total_mass = 1.0
		rp_vals[0,irad] = metal_masses[thisindex][idx].sum()
		rp_vals[0,irad] = rp_vals[0,irad]/total_mass
		i = 1
		while i < (len(metals)-1):
			idx = np.where( (metallicities[thisindex] > metals[i+1]) & (metallicities[thisindex] <= metals[i]) )[0]
			rp_vals[i,irad] = metal_masses[thisindex][idx].sum()
			rp_vals[i,irad] = rp_vals[i,irad]/total_mass
			i = i + 1 
		i = len(metals)-1
		idx = np.where(metallicities[thisindex] < -6.)[0]
		rp_vals[i,irad] = metal_masses[thisindex][idx].sum()
		rp_vals[i,irad] = rp_vals[i,irad]/total_mass

	return rp_vals, sam_Z

def plot_rp_matrix(rp_vals,rvirKPC,sam_Z,fileout):
	#plt.imshow(np.log10(rp_vals),cmap='spectral',interpolation='none',extent=(0,rvirKPC,-6,0),aspect='auto')
	plt.imshow(rp_vals,cmap='spectral',interpolation='none',extent=(0,rvirKPC,-6,0),aspect='auto')#,vmax=0.3)
	plt.colorbar()
	plt.axhline(y=sam_Z,color='w',linewidth=2.5)
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Metallicity') 
	plt.savefig(fileout)
	plt.close()

def run_many_halos(len):
	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	data = cPickle.load(open('clump_dict.cpkl','rb'))
	bins = 25.
	i = 0
	while i < len:
		print 'i is ',i
		data_source = pf.h.sphere(data['centers'][i],data['rvirs'][i]/pf['cm'])
		rvirKPC = data['rvirs'][i]/pf['cm']*pf['kpc']

		rpvals, sam_Z = build_radial_profile_matrix(pf,data_source,data['centers'][i],rvirKPC,bins)
		
		#fileout = 'MassMetalProfiles/halo'+str(data['halonum'][i])+'_massmetalprofile_normal.png'
		#fileout = '/u/10/l/lnc2115/vega/data/Wise/Plots/VolMetalProfiles/bins25/halo'+str(data['halonum'][i])+'_volmetalprofile_25totalmass.png'	
		#fileout = 'DensMetalProfiles/halo'+str(data['halonum'][i])+'_densmetalprofile_normal.png'
		fileout = '/u/10/l/lnc2115/vega/data/Wise/Plots/GasMetalProfiles/bins25/halo'+str(data['halonum'][i])+'_gasmetalprofile_25totalmass.png'
	
		plot_rp_matrix(rpvals,rvirKPC,sam_Z,fileout)
		i = i + 1

def gauss_function(x,a,mu,sigma):
	return a*np.exp(-(x-mu)**2.0/(2.0*sigma**2.0))




