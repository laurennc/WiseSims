from lauren import *

#ok what exactly do I want to do? 
#I want some way of showing how the metallicity changes at different radii...
#can look at some sort of radial profile for sure!
#instead of storing just one value want to build a table of them!
#instead of doing percents in each radius bin, can use a universally normalized quantity --- either just summing the mass directly or normalizing it to the total mass instead of the mass within the bin


def build_radial_profile_matrix(pf,data_source,center,rvirKPC,metal_bins,radii_bins):
	#metals = np.arange(101)*-6./100.
	metals = np.linspace(0,-6,num=metal_bins)
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
	dr = rvirKPC/radii_bins
	rp_r = np.arange(radii_bins)*dr + dr/2.0
	rp_vals = np.zeros((metal_bins,radii_bins))
#	total_mass = np.sum(metal_masses)
	for irad in range(int(radii_bins)):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (radii>=minrad) * (radii<maxrad)
		idx = np.where(metallicities[thisindex] > 0.)[0]
		#LETS NORMALIZE BASED OFF THE TOTAL QUANTITY IN THAT ANNULUS
       		total_mass = np.sum(metal_masses[thisindex])
		#LETS NOT NORMALIZE
		#total_mass = 1.0
		rp_vals[0,irad] = metal_masses[thisindex][idx].sum()
		rp_vals[0,irad] = rp_vals[0,irad]/total_mass
		i = 1
		while i < (metal_bins-1):
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
	radii_bins = 25.
	metal_bins = 500.
	i = 0
	while i < len:
		print 'i is ',i
		data_source = pf.h.sphere(data['centers'][i],data['rvirs'][i]/pf['cm'])
		rvirKPC = data['rvirs'][i]/pf['cm']*pf['kpc']

		rpvals, sam_Z = build_radial_profile_matrix(pf,data_source,data['centers'][i],rvirKPC,mtal_bins,radii_bins)
		
	#fileout = 'MassMetalProfiles/halo'+str(data['halonum'][i])+'_massmetalprofile_normal.png'
		#fileout = '/u/10/l/lnc2115/vega/data/Wise/Plots/VolMetalProfiles/bins25/halo'+str(data['halonum'][i])+'_volmetalprofile_25totalmass.png'	
		#fileout = 'DensMetalProfiles/halo'+str(data['halonum'][i])+'_densmetalprofile_normal.png'
		fileout = '/u/10/l/lnc2115/vega/data/Wise/Plots/GasMetalProfiles/bins25/halo'+str(data['halonum'][i])+'_gasmetalprofile_25totalmass.png'
	
		plot_rp_matrix(rpvals,rvirKPC,sam_Z,fileout)
		i = i + 1


def find_percentile_values(rp,lower,upper,metal_bins,radii_bins,tolerance):
	#What I Want This Code to Do!
	#Need the normalized rp values because then the percentages make sense!
	# For each radial bin, I want to 
	metals = np.linspace(0,-6,num=metal_bins)
	per_vals = np.zeros((3,radii_bins))
	for idr in range(int(radii_bins)):
		sum = 0.0
		for i in range(int(metal_bins)-1):
			i = i + 1
			sum = sum + rp[i,idr]
			if (sum <= lower+tolerance) and (sum >= lower-tolerance):
				#print 'lower',i, sum
				per_vals[0,idr] = metals[i]
			elif (sum <= 0.5+tolerance) and (sum >= 0.5-tolerance):
				#print 'median',i, sum
				per_vals[1,idr] = metals[i]
			elif (sum <= upper+tolerance) and (sum >= upper-tolerance):
				#print 'upper',i, sum
				per_vals[2,idr] = metals[i]
			else:
				sum = sum
	return per_vals
	

def plot_percentile_values(per_vals,rvirKPC,radii_bins,fileout):
	#NEED TO FILL IN HOW WHAT PARAMETERS I NEED TO PAS BUT FOR NOW.....
	#LET'S ASSUME THAT I HAVE PER_VALS
	#x = rp_r
	#y_low = per_vals[0,*], y_med=per_vals[1,*],y_high=per_vals[2,*]
	dr = rvirKPC/radii_bins
	rp_r = np.arange(radii_bins)*dr + dr/2.0
	plt.plot(rp_r,per_vals[0,:],'k')
	plt.plot(rp_r,per_vals[2,:],'k')
	plt.plot(rp_r,per_vals[1,:],'k',lw=2.5)
	plt.fill_between(rp_r,per_vals[0,:],per_vals[2,:],color='grey',alpha='0.5')
	plt.save(fileout)
	return


	


