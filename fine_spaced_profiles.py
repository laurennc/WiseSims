from lauren import *

#ok what exactly do I want to do? 
#I want some way of showing how the metallicity changes at different radii...
#can look at some sort of radial profile for sure!
#instead of storing just one value want to build a table of them!
#instead of doing percents in each radius bin, can use a universally normalized quantity --- either just summing the mass directly or normalizing it to the total mass instead of the mass within the bin


def build_radial_profile_matrix(pf,data_source,center,rvirKPC,bins):
	metals = np.arange(101)*-6./100.
	radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)*pf['kpc']
	metallicities = make_solar_metallicity(data_source['Metallicity'])
	metal_masses = data_source['Metal_MassMsun']
	#DOING VOLUMES BUT CALLED METAL_MASSES FOR EASE!
	#metal_masses = data_source['CellVolume']/(pf['cm']**3.0)*(pf['pc']**3.0)
	dr = rvirKPC/bins
	rp_r = np.arange(bins)*dr + dr/2.0
	rp_vals = np.zeros((len(metals),bins))
	for irad in range(int(bins)):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (radii>=minrad) * (radii<maxrad)
		idx = np.where(metallicities[thisindex] > 0.)[0]
		#LETS ALSO NORMALIZE BASED OFF THE TOTAL QUANTITY IN THAT ANNULUS
       		total_mass = np.sum(metal_masses[thisindex])
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

	return rp_vals

def plot_rp_matrix(rp_vals,rvirKPC,fileout):
	#plt.imshow(np.log10(rp_vals),cmap='spectral',interpolation='none',extent=(0,rvirKPC,-6,0),aspect='auto')
	plt.imshow(rp_vals,cmap='spectral',interpolation='none',extent=(0,rvirKPC,-6,0),aspect='auto')
	plt.colorbar()
	plt.xlabel('Radius (kpc)')
	plt.ylabel('Metallicity') 
	plt.savefig(fileout)
	plt.close()

def run_many_halos(len):
	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	data = cPickle.load(open('clump_dict.cpkl','rb'))
	i = 0
	while i < len:
		print 'i is ',i
		data_source = pf.h.sphere(data['centers'][i],data['rvirs'][i]/pf['cm'])
		rvirKPC = data['rvirs'][i]/pf['cm']*pf['kpc']

		rpvals = build_radial_profile_matrix(pf,data_source,data['centers'][i],rvirKPC,100.)
		
		fileout = 'MassMetalProfiles/halo'+str(data['halonum'][i])+'_massmetalprofile_normal.png'
		#fileout = 'VolMetalProfiles/halo'+str(data['halonum'][i])+'_volmetalprofile_normal.png'		

		plot_rp_matrix(rpvals,rvirKPC,fileout)
		i = i + 1



