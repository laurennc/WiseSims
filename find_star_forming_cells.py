from lauren import *

proton_mass = 1.67e-24 #g

def 

	dens_cut = 400.

	radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)*pf['kpc']
	number_density = data_source['Density']/proton_mass
	metallicities = make_solar_metallicity(data_source['Metallicity'])

	thisindex = (number_density > dens_cut) * (data_source['Temperature'] < 1000.)
	idx = np.where(thisindex == True)[0]

	plt.plot(radii[idx],metallicities[idx],'rs')#,'bo')

	plt.xlabel('Radius (kpc)')
	plt.ylabel('Metallicity')
	plt.xlim(0.,rvirKPC)
	plt.savefig('halo0_starformcells2.png')

