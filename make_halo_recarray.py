from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import matplotlib.pyplot as plt
from lauren import *

YEAR = 3.155693e+7 #s/yr

pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat')
field = 'Metallicity'
step = 2.0
c_min = 1.0e-17
c_max = 1.0e-15
function = 'self.data["Metallicity"].size > 5e+2'
#function = 'self.data["Metallicity"].size > 5e+3'

from numpy import *

#Lists that I'll need
#halonum to identify the halos!
halonumber = [0]
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
mvirs = [0.]
#list of the filling factors
#metallicityList = [1e-6,1e-5,1e-4,1e-3,1e-2,1e-1,1.0]
metallicityList = [-6,-5,-4,-3,-2,-1,0]
fillFactors = array([0.,0.,0.,0.,0.,0.,0.])
#list of the radii
avgRadius = [0.]
maxRadius = [0.]
minRadius = [0.]
zMaxRadius = [0.]
zMinRadius = [0.]
massRadius = [0.]
#SFR quantities
sfrAvg = [0.]
sfrMax = [0.]
#List of the index of the clump with smallest offset from halo 
clump_index = [0.]
t_start = [0.0]
#Quantities needed for Dilution Mass Models
dilmass6 = [0.]
dilmass17 = [0.]
dilmassvir = [0.]

count = 0
halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]

while (count < len(halos_to_run)):
	halonum = halos_to_run[count]
	halonumber = append(halonumber,halonum)
	print 'Analyzing halo number: ',halonum

	radius, mass, center = r200(pf, halos['pos'][halonum,:],halos['mass'][halonum],verbose=False)
	filein = '/u/10/l/lnc2115/vega/data/Wise/pickles/halo'+str(halonum)+'_clumps.cpkl'
	data = cPickle.load(open(filein,'rb'))
	master_clump = data[1]
	mvirs = append(mvirs,mass)	

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
	

	num_tot = len(master_clump.children)-1
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

	if (halonum == 15):
		print 'changing for halo 15'
		index = 2
	clump_index = append(clump_index,index)
	masses = append(masses,all_clumps[index].quantities['TotalQuantity']('CellMassMsun'))
	holder = all_clumps[index].quantities['CenterOfMass']()
	CofMs = vstack((CofMs,holder))
	offsets_temp = (((holder[0]-center[0])**2.0+(holder[1]-center[1])**2.0+(holder[2]-center[2])**2.0)**0.5)*pf['kpc']
	offsets = append(offsets,offsets_temp)

	avgRadius = append(avgRadius,average(all_clumps[index]['Radius'])/pf['cm']*pf['kpc'])	
	maxRadius = append(maxRadius,all_clumps[index]['Radius'].max()/pf['cm']*pf['kpc'])
	minRadius = append(minRadius,all_clumps[index]['Radius'].min()/pf['cm']*pf['kpc'])
	massRadius = append(massRadius,all_clumps[index].quantities['WeightedAverageQuantity']('Radius','CellMassMsun')/pf['cm']*pf['kpc'])

	metallicities = make_solar_metallicity(data_source['Metallicity'])
	volumes = data_source['CellVolume']
	total_volume = volumes.sum()

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
	if (len(indices) > 0):
		print all_clumps[index]['Radius'][indices].max()*pf['kpc']/pf['cm']
		zMaxRadius = append(zMaxRadius,all_clumps[index]['Radius'][indices].max()*pf['kpc']/pf['cm'])
		zMinRadius = append(zMinRadius,all_clumps[index]['Radius'][indices].min()*pf['kpc']/pf['cm'])
	else:
		zMaxRadius = append(zMaxRadius,0.0)
		zMinRadius = append(zMinRadius,0.0)

	dist = ((all_clumps[index]['x']-center[0])**2.0+(all_clumps[index]['y']-center[1])**2.0+(all_clumps[index]['z']-center[2])**2.0)**0.5
	idx = where(all_clumps[index]['Metallicity'] >= 1e-6)
	temp = where(dist[idx] > radius/pf['cm'])
	dilmass6 = append(dilmass6, all_clumps[index]['CellMassMsun'][temp].sum())
	
	idx = where(all_clumps[index]['Metallicity'] >= 1e-17)
	temp = where(dist[idx] > radius/pf['cm'])
	dilmass17 = append(dilmass17, all_clumps[index]['CellMassMsun'][temp].sum())

	temp = where(dist < radius/pf['cm'])
	dilmassvir = append(dilmassvir,all_clumps[index]['CellMassMsun'][temp].sum())

	count = count + 1
	loopthrough = 0


#Need to delete the zeroes from the beginning of the arrays!
halonumber = delete(halonumber,0)
centers = delete(centers,0,0)
rvirs = delete(rvirs,0)
star_masses = delete(star_masses, 0)
CofMs = delete(CofMs,0,0)
masses = delete(masses,0)
offsets = delete(offsets,0)
avgRadius = delete(avgRadius,0)
maxRadius = delete(maxRadius,0)
minRadius = delete(minRadius,0)
massRadius = delete(massRadius,0)
fillFactors = delete(fillFactors,0,0)
sfrAvg = delete(sfrAvg,0)
sfrMax = delete(sfrMax,0)
clump_index = delete(clump_index,0)
t_start = delete(t_start,0)
zMaxRadius = delete(zMaxRadius,0)
zMinRadius = delete(zMinRadius,0)
dilmass6 = delete(dilmass6,0)
dilmass17 = delete(dilmass17,0)
dilmassvir = delete(dilmassvir,0)
mvirs = delete(mvirs,0)

#create and pickle the recarray for easy loading
data = {}
data['halonum'] = halonumber
data['centers'] = centers
data['rvirs'] = rvirs
data['mvirs'] = mvirs
data['star_masses'] = star_masses
data['CofMs'] = CofMs
data['masses'] = masses
data['offsets'] = offsets
data['avgRadius'] = avgRadius
data['maxRadius'] = maxRadius
data['minRadius'] = minRadius
data['massRadius'] = massRadius
data['fillFactors'] = fillFactors
data['sfrAvg'] = sfrAvg
data['sfrMax'] = sfrMax
data['clump_index'] = clump_index
data['t_start'] = t_start
data['zMaxRadius'] = zMaxRadius
data['zMinRadius'] = zMinRadius
data['dilmass6'] = dilmass6
data['dilmass17'] = dilmass17
data['dilmassvir'] = dilmassvir

fileout = 'clump_dict.cpkl'
cPickle.dump(data,open(fileout,'wb'),protocol=-1)



