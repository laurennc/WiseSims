
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt

YEAR = 3.155693e7 #sec / yr


def clump_contours(pf,center,radius,field,dim,clump,width,fileout,plottype):
	#all_clumps = get_lowest_clumps(master_clump)

	#num_tot = len(master_clump.children)-1
	radius = radius/pf['cm']	

	pc = PlotCollection(pf,center=center)
	
	#zlim = {"Density": (1e-27, 1e-22), "Metallicity":()}	
	if field == 'Density':
		pc.set_zlim(1e-27,1e-22)

	if (plottype == 'proj'):
		pc.add_projection(field,dim)
	elif (plottype == 'slice'):
		pc.add_slice(field,dim)
	else:
		return false

	pc.set_width(width,'kpc')
	clump_plot = [clump]
	pc.plots[0].modify['clumps'](clump_plot,plot_args={'colors':'DarkOrange'})#,'linewidth':2.5})
	pc.plots[0].modify['point'](center,'o',text_args={'color':'Purple'})
	pc.plots[0].modify['point'](clump.quantities['CenterOfMass'](),'x',text_args={'color':'DarkOrange'})
	pc.plots[0].modify['sphere'](center,radius,circle_args={'color':'Purple'})#,'linestyle':'--','linewidth':2.5})
	pc.set_cmap('Blues_r')
	pc.set_zlim(10**-8,1)
#	pc.set_cmap('spectral')
#	pc.set_cmap('binary')
	pc.save(fileout)
	return True

def sfr_halo(pf,virial_sphere):
	sm = virial_sphere['ParticleMassMsun']
	ct = virial_sphere['creation_time']
	stars = (ct > 0)
	ct = ct[stars]
	sm = sm[stars]
	total_volume = virial_sphere.quantities['TotalQuantity']('CellVolume')
	sfr = StarFormationRate(pf,data_source=virial_sphere,star_mass=sm,star_creation_time=ct,volume=total_volume)
	return sfr

def sfr_quants(region):
	sm = region['ParticleMassMsun']
	ct = region['creation_time']
	metals = region['metallicity_fraction']
	metal = region['metallicity_fraction']
	stars = (ct > 0)
	sm = sm[stars]
	ct = ct[stars]
	metal = metal[stars]
	return sm, ct, metal

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
	#plt.plot(sfr.lookback_time,sfr.Msol_yr,'g*')
	plt.plot(np.log10(sfr.time),np.log10(sfr.Msol_yr),'g-',linewidth=1.3)
	#plt.xscale('log')
	#plt.yscale('log')
	#plt.xlabel('Lookback Time (yrs)')
	plt.xlabel('log(Time) (yrs)')
	plt.ylabel('log(SFR) (Msol/yr)')
	#plt.ylim(1.0e-5,1.0)
	plt.ylim(-5,0)
	plt.savefig(fileout)
	plt.close()
	return True

def plot_stellar_t_vs_Z(pf,t,metals,fileout,mark_timesteps=False,timesteps=0.):
	#t and metals should be the values from the SIMULATIONS 
	#otherwise, these corrections don't make sense
	#t = t * pf['Time'] / YEAR
	metals = np.log10(metals)-np.log10(0.02)	
	#t = np.log10(pf['Time']*t/YEAR)

	fig,axes = plt.subplots(2,1,sharex=True)
	fig.subplots_adjust(hspace=0.1,wspace=0.1)
	axes[0].plot(t,metals,'bo')
	#axes[0].set_xlabel('Particle Creation Time (yrs)')
	axes[0].set_ylabel('Metallicity')
	axes[1].plot(t,metals,'bo')
	axes[1].set_xlabel('Particle Creation Time (yrs)')
	axes[1].set_ylabel('Metallicity')
	axes[1].set_ylim(-4,1)

	if mark_timesteps:
		idx = np.where((timesteps['sim_time'] <= t.max()) & (timesteps['sim_time'] >= t.min()))[0]
		for val in idx:
			axes[0].axvline(x=timesteps['sim_time'][val],linestyle='dashed',color='red')
			axes[1].axvline(x=timesteps['sim_time'][val],linestyle='dashed',color='red')
			axes[0].text(timesteps['sim_time'][val],-2,str(time_steps['stamp'][val]))
			axes[1].text(timesteps['sim_time'][val],-2,str(time_steps['stamp'][val]))

	fig.savefig(fileout)
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


def dilution_model(tsn,tnow,t_dil,M_dil):
	return M_dil*(1.0 - np.exp((tsn-tnow)/t_dil))

def return_snwinds_scaled_radii(t_start,mass_in):
	L_sun = 6.3706e-21 #solar mass kpc^2 yr^-3
	G = 4.4986e-24 #kpc^3 yr^-2 solar mass^-1
	H_0 = 7.48609e-11 #yr^-1
	little_h = 0.732
	omega_b = 0.0416
	f_sn = 0.01

	L_scale = 1.2*L_sun*mass_in
	return (L_scale**0.2)*(G**0.2)*(t_start)

def make_solar_metallicity(metals):
	return np.log10(metals) - np.log10(0.02)

def distance_from_center(x,y,z,center):
	return ((x-center[0])**2.0+(y-center[1])**2.0+(z-center[2])**2.0)**0.5

def readin_masterclump(pf,halos,halonum):
	radius,mass,center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	filein = '/hmt/hpcstorage1/hpc/astro/users/lnc2115/WiseSimsData/pickles/halo'+str(halonum)+'_clumps.cpkl'
	data = cPickle.load(open(filein,'rb'))
	master_clump = data[1]
	return master_clump,radius,mass,center

def output_file_readin(fileinput):
	return np.genfromtxt(fileinput,dtype=[("x",float),("y",float),("z",float),("N_part",float),("m_halo",float),("m_FOF",float),("m_star",float),("m_gas",float),("rvir",float),("v_bulk",float),("v_rms",float),("lambda_DM",float),("lambda_gas",float)])

def write_virial_values(pf,halos,fileout):
	f = open(fileout,'wb')
	for i in range(len(halos['mass'])):
		radius, mass, center = r200(pf, halos['pos'][i,:],halos['mass'][i],verbose=False)
		s = str(i)+'\t'+str(center)+'\t'+str(radius)+'\t'+str(mass)+'\n'
		f.write(s)
	f.close()

def make_radial_profile(pf,data_source,center,rvirKPC,bins):
	radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)*pf['kpc']
	dr = rvirKPC/bins
	rp_r = np.arange(bins)*dr + dr/2.0
	rp_vals = []
	for irad in range(int(bins)):
		minrad = irad*dr
		maxrad = minrad + dr
		thisindex = (radii>=minrad) * (radii<maxrad)
		rp_vals = np.append(rp_vals,np.average(make_solar_metallicity(data_source['Metallicity'][thisindex])))
	return rp_r,rp_vals

def gauss_function(x,a,mu,sigma):
	return a*np.exp(-(x-mu)**2.0/(2.0*sigma**2.0))

def load_tree(hid, fn):
    if not os.path.exists(fn):
        print "Cannot find %s" % fn
        sys.exit()
    lines = open(fn, "r").readlines()

    # Search for tree
    start_line = None
    rank = -hid
    num = 0
    if hid <= 0:
        for l in lines:
            if l.startswith("#tree"):
                if num == -hid:
                    hid = int(l.split()[1])
                    break
                num += 1
        print "Halo %d (by mass) has ID %d" % (rank, hid)

    header = "#tree %d\n" % hid
    for i,l in enumerate(lines):
        if l.startswith("#Omega_M"):
            h = float(l.split("=")[3])
        if l == header:
            start_line = i+1
        elif l.startswith("#tree") and start_line != None:
            end_line = i
            break

    # Load data from strings
    nleaves = end_line - start_line
    nfields = len(lines[start_line].split())
    data = na.empty((nleaves, nfields))
    for i in range(start_line, end_line):
        data[i-start_line,:] = map(float, lines[i].split())
    del lines
    return hid, data

def read_in_timesteps(filein):
	data = np.genfromtxt(filein,delimiter=' ')
	timesteps = {}
	timesteps['stamp'] = data[:,0]
	timesteps['sim_time'] = data[:,1]
	return timesteps
