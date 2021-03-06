from lauren import *

#ok what exactly do I want to do? 
#I want some way of showing how the metallicity changes at different radii...
#can look at some sort of radial profile for sure!
#instead of storing just one value want to build a table of them!
#instead of doing percents in each radius bin, can use a universally normalized quantity --- either just summing the mass directly or normalizing it to the total mass instead of the mass within the bin


def build_radial_profile_matrix(pf,data_source,center,rvirKPC,field,metal_bins,radii_bins):
	metals = np.linspace(0,-6,num=metal_bins)
	radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)*pf['kpc']

	### ORIGINAL
	#metallicities = np.log10(data_source['Metallicity'])
	## TRYING WITH THE POPULATION III METALLIITY INCLUDED
	metallicities = np.log10(total_metallicity(data_source))

	#WANT TO CALCULATE METALLICITY THAT WOULD BE USED IN SAMS
	sam_Z = make_solar_metallicity(np.sum(data_source['Metal_MassMsun'])/np.sum(data_source['CellMassMsun']))

	#DOING VOLUMES BUT CALLED METAL_MASSES FOR EASE!
	#metal_masses = data_source['CellVolume']/(pf['cm']**3.0)*(pf['pc']**3.0)
	#DOING TOTAL GAS MASS BUT CALLED METAL_MASSES FOR EASE!
	#metal_masses = data_source['CellMassMsun']
	metal_masses = data_source[field]
	dr = rvirKPC/radii_bins
	rp_r = np.arange(radii_bins)*dr + dr/2.0
	rp_vals = np.zeros((metal_bins,radii_bins))
	#NORMALIZE BY THE TOTAL QUANTITY
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

def run_many_halos(halo_array):
	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	data = cPickle.load(open('clump_dict.cpkl','rb'))
	radii_bins = 20.
	metal_bins = 700.
	j = 0
	iax = 131
	fig = plt.figure(figsize=(12,4))
	while j < len(halo_array):
		print 'j is ',j
		i = np.where(data['halonum'] == halo_array[j])[0][0]
		data_source = pf.h.sphere(data['centers'][i],data['rvirs'][i]/pf['cm'])
		rvirKPC = data['rvirs'][i]/pf['cm']*pf['kpc']

		rpvals, sam_Z = build_radial_profile_matrix(pf,data_source,data['centers'][i],rvirKPC,'CellMassMsun',metal_bins,radii_bins)
		per_vals = find_percentile_values(rpvals,0.25,0.75,metal_bins,radii_bins,0.02,'CellMassMsun',halonum=data['halonum'][i])
		print 'iax is ',iax
		ax1 = fig.add_subplot(iax)
		fileout = 'Z_rp_quartiles_halo'+str(data['halonum'][i])+'.png'
		plot_percentile_values(ax1,per_vals,data['halonum'][i],rvirKPC,radii_bins,fileout,'#389BBc','#33CCFF')	


		rpvals, sam_Z = build_radial_profile_matrix(pf,data_source,data['centers'][i],rvirKPC,'CellVolume',metal_bins,radii_bins)
		per_vals = find_percentile_values(rpvals,0.25,0.75,metal_bins,radii_bins,0.02,'CellVolume',halonum=data['halonum'][i])
		plot_percentile_values(ax1,per_vals,data['halonum'][i],rvirKPC,radii_bins,fileout,'#E65C8A','#FFA3C2')

		#plot_rp_matrix(rpvals,rvirKPC,sam_Z,fileout)
		j  = j + 1
		iax = iax + 1
	fileout = 'Z_rp_quartiles_paper_wPopIII.pdf'
	plt.savefig(fileout)
	#plt.close()

	return


def find_percentile_values(rp,lower,upper,metal_bins,radii_bins,tolerance,field,halonum=0):
	#Using the normalized radial profiles, find each percentile and the median to plot instead of the heat map	
	metals = np.linspace(0,-6,num=metal_bins)
	per_vals = np.zeros((3,radii_bins))
	for idr in range(int(radii_bins)):
		sum = 0.0
		if (rp[0,idr] > lower):
			per_vals[0,idr] = metals[0]
		for i in range(int(metal_bins)): #-1):
			#i = i + 1
			sum = sum + rp[i,idr]
			if (sum <= lower+tolerance) and (sum >= lower-tolerance):
				per_vals[0,idr] = metals[i]
			elif (sum <= 0.5+tolerance) and (sum >= 0.5-tolerance):
				per_vals[1,idr] = metals[i]
			elif (sum <= upper+tolerance) and (sum >= upper-tolerance):
				per_vals[2,idr] = metals[i]
			else:
				sum = sum
		if per_vals[2,idr] == 0.0:
			per_vals[2,idr] = metals[metal_bins-1]	

	if ((halonum == 3) & (field == 'CellMassMsun')):
		per_vals[1,0],per_vals[2,0] = -2.8025751072961373, -2.8111587982832615
		per_vals[1,1],per_vals[2,1] = -2.8197424892703862, -2.8841201716738198

	if ((halonum==3) & (field == 'CellVolume') ):
		per_vals[0,0],per_vals[1,0],per_vals[2,0] = metals[323],metals[326],metals[327]

	return per_vals
	

def plot_percentile_values(ax1,per_vals,halonum,rvirKPC,radii_bins,fileout,colorline,colorfill):
	#y_low = per_vals[0,*], y_med=per_vals[1,*],y_high=per_vals[2,*]
	dr = rvirKPC/radii_bins
	rp_r = np.arange(radii_bins)*dr + dr/2.0
	ax1.plot(rp_r/rvirKPC,per_vals[0,:],color=colorline)
	ax1.plot(rp_r/rvirKPC,per_vals[2,:],color=colorline)
	ax1.plot(rp_r/rvirKPC,per_vals[1,:],color=colorline,lw=2.5)
	ax1.fill_between(rp_r/rvirKPC,per_vals[0,:],per_vals[2,:],color=colorfill,alpha='0.30')
	ax1.text(0.2,-0.5,'Halo '+str(halonum))
	ax1.set_ylim(-6,-1)
	#plt.savefig(fileout)
	#plt.close()
	return


	


