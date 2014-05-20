from lauren import *

#The idea of this file is to accumulate functions that I would like to be able to run on all of my halos and be able to save and output the results!

def radius_containing_stars(data_source,center):
	ct = data_source['creation_time']
	sm = data_source['ParticleMassMsun']
	stars = ((ct > 0) & (sm > 1))
	dists = distance_from_center(data_source['particle_position_x'][stars],data_source['particle_position_y'][stars],data_source['particle_position_z'][stars],center)
	return dists.max(),np.average(dists)

def center_by_densest_points(data_source,npoints,center):
	idSort = np.argsort(data_source['Density'])[::-1]
	idSort = idSort[0:npoints]
	center_dens = [np.average(data_source['x'][idSort]),np.average(data_source['y'][idSort]),np.average(data_source['z'][idSort])]
	dist = distance_from_center(center_dens[0],center_dens[1],center_dens[2],center)
	return center_dens, dist

def run_through_halos():
	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
	data = cPickle.load(open('0062.cpkl','rb'))
	
	#data I'm collecting this time...
	centers = np.array([0.,0.,0.])
	dist_apart = []
	star_dists1 = []
	star_dists2 = []
	for i in range(len(data['halonum'])):
		center = data['centers'][i]
		rvir = data['rvirs'][i]
		data_source = pf.h.sphere(center,rvir/pf['cm'])
		center_dens, dist = center_by_densest_points(data_source,100,center)
		centers = np.vstack((centers,center_dens))
		dist_apart = np.append(dist_apart,dist)
		starDistMax, starDistAvg = radius_containing_stars(data_source,center)
		star_dists1 = np.append(star_dists1,starDistMax/(rvir/pf['cm']))
		star_dists2 = np.append(star_dists2,starDistAvg/(rvir/pf['cm']))

		#starDistMax, starDistAvg = radius_containing_stars(data_source,center_dens)
		#star_dists2 = np.append(star_dists2,starDist)

	centers = np.delete(centers,0,0)
	return centers, dist_apart, star_dists1,star_dists2
 
