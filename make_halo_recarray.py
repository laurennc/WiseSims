from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import matplotlib.pyplot as plt

halonum = 0
nhalos = 12
YEAR = 3.155693e+7 #s/yr

pf = load('/u/10/l/lnc2115/home/WiseSimsData/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/home/WiseSimsData/groups_02797.dat')
field = 'Metallicity'
step = 2.0
c_min = 1.0e-17
c_max = 1.0e-15
function = 'self.data["Metallicity"].size > 5e+2'
#function = 'self.data["Metallicity"].size > 5e+3'

from numpy import *

#Lists that I'll need
#make a list of the centers
centers = [0.,0.,0.]
#Virial Radii in cm
rvirs = [0.]
#Stellar Mass within halo's virial radius
star_masses=[0.]
#need a list of the CofM of the main clumps
CofMs = [0.,0.,0.]
#list of the offset
offsets = [0.]
#list of the masses
masses = [0.]
#list of the filling factors
metallicityList = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.0]
fillFactors = array([0.,0.,0.,0.,0.,0.,0.])
#list of the radii
avgRadius = [0.]
maxRadius = [0.]
minRadius = [0.]
zcutRadius = [0.]
#SFR quantities
sfrAvg = [0.]
sfrMax = [0.]
#List of the index of the clump with smallest offset from halo 
clump_index = [0.]
t_start = [0.0]

while (halonum < nhalos):
	print 'Analyzing halo number: ',halonum

	radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	filein = '../WiseSimsData/halo'+str(halonum)+'_clumps30.cpkl'
	data = cPickle.load(open(filein,'rb'))
	master_clump = data[1]
	
	#stellar mass of the halo within the virial radius
	data_source = pf.h.sphere(center,radius/pf['cm'])
	star_masses = append(star_masses,data_source.quantities['TotalQuantity']('StarMassMsun'))

	#should calculate SFR inside virial radius only!
	sm = data_source['ParticleMassMsun']
	ct = data_source['creation_time']
	stars = (ct > 0)
        ct = ct[stars]
        sm = sm[stars]
	total_volume = data_source.quantities['TotalQuantity']('CellVolume')

        if (len(stars[where(stars==True)]) == 0):
               sfrAvg = append(sfrAvg,0.0)
               sfrMax = append(sfrMax,0.0)
               t_start = append(t_start,0.0)
        else:
               sfr = StarFormationRate(pf,star_mass=sm,star_creation_time=ct,volume=total_volume)
               sfrAvg = append(sfrAvg,average(sfr.Msol_yr))
               sfrMax = append(sfrMax,sfr.Msol_yr.max())
	       t_start_here = min(ct)*pf['Time']/YEAR
               t_start = append(t_start,t_start_here)
	

	#num_tot = len(master_clump.children)-1
	centers = vstack((centers,center))
	rvirs = append(rvirs,radius)
	#need to loop through the lowest clumps to find which has the smallest offset from the center of the halo
	num = 0
	comparing = 1e9
	index = 1000
	all_clumps = get_lowest_clumps(master_clump)
	while (num < len(all_clumps)):
		holder = all_clumps[num].quantities['CenterOfMass']()
		offsets_temp = (((holder[0]-center[0])**2.0+(holder[1]-center[1])**2.0+(holder[2]-center[2])**2.0)**0.5)*pf['kpc']
		if (offsets_temp < comparing):
			comparing = copy(offsets_temp)
			index = copy(num)
		num = num + 1

	clump_index = append(clump_index,index)
	masses = append(masses,all_clumps[index].quantities['TotalQuantity']('CellMassMsun'))
	holder = all_clumps[index].quantities['CenterOfMass']()
	CofMs = vstack((CofMs,holder))
	offsets_temp = (((holder[0]-center[0])**2.0+(holder[1]-center[1])**2.0+(holder[2]-center[2])**2.0)**0.5)*pf['kpc']
	offsets = append(offsets,offsets_temp)

	avgRadius = append(avgRadius,average(all_clumps[index]['Radius'])/pf['cm']*pf['kpc'])	
	maxRadius = append(maxRadius,all_clumps[index]['Radius'].max()/pf['cm']*pf['kpc'])
	minRadius = append(minRadius,all_clumps[index]['Radius'].min()/pf['cm']*pf['kpc'])
	#Example of old code: minRadius = append(maxRadius,master_clump.children[num_tot]['Radius'].max()/pf['cm']*pf['kpc'])

	metallicities = all_clumps[index]['Metallicity']

	volumes = all_clumps[index]['CellVolume']
	total_volume = all_clumps[index].quantities['TotalQuantity']('CellVolume')

	loopthrough = 0
	fillFactorsTemp = array([0.,0.,0.,0.,0.,0.,0.])
	while (loopthrough < 7):
		goodvalue = metallicityList[loopthrough]
		indices = where(metallicities >= goodvalue) 
		fillFactorsTemp[loopthrough] = volumes[indices].sum()/total_volume
		loopthrough = loopthrough + 1
	fillFactors = vstack((fillFactors,fillFactorsTemp))

	goodvalue = 1e-6
	indices = where(all_clumps[index]['Metallicity'] >= goodvalue)
	zcutRadius = append(zcutRadius,all_clumps[index]['Radius'][indices].max()*pf['kpc']/pf['cm'])


	halonum = halonum + 1
	loopthrough = 0


#Need to delete the zeroes from the beginning of the arrays!
centers = delete(centers,0,0)
rvirs = delete(rvirs,0)
star_masses = delete(star_masses, 0)
CofMs = delete(CofMs,0,0)
masses = delete(masses,0)
offsets = delete(offsets,0)
avgRadius = delete(avgRadius,0)
maxRadius = delete(maxRadius,0)
minRadius = delete(minRadius,0)
fillFactors = delete(fillFactors,0,0)
sfrAvg = delete(sfrAvg,0)
sfrMax = delete(sfrMax,0)
clump_index = delete(clump_index,0)
t_start = delete(t_start,0)
zcutRadius = delete(zcutRadius,0)

#create and pickle the recarray for easy loading
data = {}
data['centers'] = centers
data['rvirs'] = rvirs
data['star_masses'] = star_masses
data['CofMs'] = CofMs
data['masses'] = masses
data['offsets'] = offsets
data['avgRadius'] = avgRadius
data['maxRadius'] = maxRadius
data['minRadius'] = minRadius
data['fillFactors'] = fillFactors
data['sfrAvg'] = sfrAvg
data['sfrMax'] = sfrMax
data['clump_index'] = clump_index
data['t_start'] = t_start
data['zcutRadius'] = zcutRadius

fileout = 'clump_dict30.cpkl'
cPickle.dump(data,open(fileout,'wb'),protocol=-1)



