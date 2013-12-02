#from /hmt/hpcstorage1/hpc/astro/users/lnc2115/WiseSims/lauren import *
from lauren import *

YEAR = 3.155693e+7 #s/yr

#yt datasets for each of the halos
#pfL = load('/u/10/l/lnc2115/home/WiseSimsData/DD0062/output_0062')
#halosL = readhalos(foffile='/u/10/l/lnc2115/home/WiseSimsData/groups_02797.dat')
#pfL = load('/hpc/astrostats/astro/users/lnc2115/DD0135/output_0135')
#halosL = readhalos(foffile='/hpc/astrostats/astro/users/lnc2115/groups_02789.dat')
#pfL = load('/hpc/astrostats/astro/users/lnc2115/DD0134/output_0134')
#halosL = readhalos(foffile='/hpc/astrostats/astro/users/lnc2115/groups_02784.dat')
#pfL = load('/hpc/astrostats/astro/users/lnc2115/DD0133/output_0133')
#halosL = readhalos(foffile='/hpc/astrostats/astro/users/lnc2115/groups_02779.dat')
#pfL = load('/hpc/astrostats/astro/users/lnc2115/DD0059/output_0059')
#halosL = readhalos(foffile='/hpc/astrostats/astro/users/lnc2115/groups_02685.dat')
pfL = load('/hpc/astrostats/astro/users/lnc2115/DD0037/output_0037')
halosL = readhalos(foffile='/hpc/astrostats/astro/users/lnc2115/groups_01260.dat')
pfE = load('/hpc/astrostats/astro/users/lnc2115/DD0036/output_0036')
halosE = readhalos(foffile='/hpc/astrostats/astro/users/lnc2115/groups_01172.dat')

count = 0
#halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]
#halos_to_run = [0,3,4,5,6,7,10,16,18,27,28,29,33,34,35,43,51,57]
#halos_to_run = [0,3,4,5,6,7,10,16,19,28,29,30,34,35,36,41,57]
#halos_to_run = [0,3,4,5,6,7,10,16,18,20,29,30,35,37,38,44,59]
#halos_to_run = [0,2,5,6,7,15,17,20,27,29,33,36,39,42,55]
halos_to_run = [0,16,2,6,5,19,51,9,13,10,15,38,30,27,60,247,42,228]
#halos_to_run = [0]

#output_file = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0062-halos.dat')
#output_file = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0135-halos.dat')
#output_file = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0134-halos.dat')
#output_file = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0133-halos.dat')
#output_file = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0059-halos.dat')
output_file = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0037-halos.dat')
output_fileE = output_file_readin('/hpc/astrostats/astro/users/lnc2115/output_0036-halos.dat')

#idstars = np.where(output_fileE['m_star'] > 0.0)[0]
#idstars = np.where(np.log(output_fileE['m_halo']) >= 6.0 )[0]
#all_ids = np.where(output_fileE['m_halo'] > 0.0)[0]

tnow = pfL.current_time*pfL['Time']/YEAR
tthen = pfE.current_time*pfE['Time']/YEAR

f = open('dd0037_traceback_all.out','wb')
#f = open('testing_all.out','wb')

while (count < len(halos_to_run)):
	halonum = halos_to_run[count]
	rvirL,massL,centerL = r200(pfL,halosL['pos'][halonum,:],halosL['mass'][halonum],verbose=False)
	haloL = pfL.h.sphere(centerL,rvirL/pfL['cm'])
	particlesL = haloL['particle_index']

	#r = search_radius(pfL,5.0,output_file['v_bulk'][count],tnow,tthen)
	#let's assume that code unit positions remain the same at different time steps
	#check for particles in halos that are within this distance, r
	r = 25./pfE['kpc']
	#all of these values should be in code units
	dists = distance_from_center(halosE['pos'][:,0],halosE['pos'][:,1],halosE['pos'][:,2],centerL)
	idxE = np.where(dists <  r)[0]

	percents = []
	indices = []
	for i in idxE:
	#for i in all_ids:
	#for i in idstars:
		rvirE,massE,centerE = r200(pfE,halosE['pos'][i,:],halosE['mass'][i],verbose=False)
		if rvirE > 0:
			print i
			haloE = pfE.h.sphere(centerE,rvirE/pfE['cm'])
			particlesE = haloE['particle_index']
			c = list(set(particlesE).intersection(particlesL))
			percents = np.append(percents,float(len(c))/float(len(particlesL)))
			indices = np.append(indices,i)

	#ok now that I have all the halos that I want to search...
	#perhaps write them to some sort of file to be looked at....
	#really -- just want the halo with the highest percentage for each one but should examine to make sure it's a sizeable amount and not something like a 10 percent match... shouldn't be too hard to do for all of the individually... but still printing to a file will make the run time easier!
	#maxper = max(percents)
	#maxind = indices[np.argmax(percents)]
	#print percents
	#print indices
	maxind = np.where(percents > 0.1)[0]
	maxpers = percents[maxind]
	maxind = idxE[maxind]
	s = str(halonum)+'\t'+str(maxpers)+'\t'+str(maxind)+'\n'
	f.write(s)
	f.write('\n')

	print percents.sum()
	count = count + 1

f.close()

