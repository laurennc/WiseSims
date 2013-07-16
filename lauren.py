
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt



def clump_contours(pf,center,radius,field,dim,clump,width,fileout,plottype):
	#all_clumps = get_lowest_clumps(master_clump)
	#num_tot = len(master_clump.children)-1
	radius = radius/pf['cm']	

	pc = PlotCollection(pf,center=center)
	if (plottype == 'proj'):
		pc.add_projection(field,dim)
	elif (plottype == 'slice'):
		pc.add_slice(field,dim)
	else:
		return false
	pc.set_width(width,'kpc')
	clump_plot = [clump]
	pc.plots[0].modify['clumps'](clump_plot)
	pc.plots[0].modify['point'](center,'o')
	#pc.plots[0].modify['point'](all_clumps[1].quantities['CenterOfMass'](),'x')
	pc.plots[0].modify['point'](clump.quantities['CenterOfMass'](),'x')
	pc.plots[0].modify['sphere'](center,radius)
	#pc.set_cmap('cool')
	pc.save(fileout)
	return True

def sfr_halo(pf,virial_sphere):
	sm = virial_sphere['ParticleMassMsun']
	ct = virial_sphere['creation_time']
	stars = (ct > 0)
	ct = ct[stars]
	sm = sm[stars]
	total_volume = virial_sphere.quantities['TotalQuantity']('CellVolume')
	sfr = StarFormationRate(pf,star_mass=sm,star_creation_time=ct,volume=total_volume)
	return sfr

def sfr_clumps(pf,master_clump):
	num_tot = len(master_clump.children) - 1
	sm = master_clump.children[num_tot]['ParticleMassMsun']
	ct = master_clump.children[num_tot]['creation_time']
	stars = (ct > 0)
	ct = ct[stars]
	sm = sm[stars]
	total_volume = master_clump.children[num_tot].quantities['TotalQuantity']('CellVolume')

	sfr = StarFormationRate(pf, star_mass=sm, star_creation_time=ct, volume=total_volume)
	return sfr

def sfh_plot(sfr,fileout):
	#plt.figure(fig_unit)
	plt.plot(sfr.lookback_time,sfr.Msol_yr)
	plt.xscale('log')
	plt.yscale('log')
	plt.xlabel('Lookback Time (yrs)')
	plt.ylabel('SFR (Msol/yr)')
	plt.savefig(fileout)
	plt.close()
	return True


def metallicity_hist(metallicities,nbins,fileout):
	#plt.figure(fig_unit)
	indicies = np.log(metallicities) > -10
	log_metals = np.log(metallicities[indicies])-np.log(0.02)
	n, bins, patches = plt.hist(log_metals,nbins,facecolor='g')
	plt.xlabel('Metallicity (solar)')
	plt.ylabel('Number of Cells')
	plt.savefig(fileout)
	plt.close()	
	return n, bins, patches
