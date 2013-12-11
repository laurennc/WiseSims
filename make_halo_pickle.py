from lauren import *

halo0 = [0,0,0,0,0]
halo3 = [3,5,6,6,6]


timepickles = ['0062.cpkl','0135.cpkl','0134.cpkl','0133.cpkl','0059.cpkl']

count = 0

centers = np.array([0.,0.,0.])
rvirs = []
mvirs = []
fillFactors = np.array([0.,0.,0.,0.,0.,0.,0.])
ravg = np.array([0.,0.,0.,0.,0.,0.,0.])
rmax = np.array([0.,0.,0.,0.,0.,0.,0.])
redshift = [] 
time = []

listappends  = [rvirs,mvirs,redshift,time]
appendstrings = ['rvirs','mvirs','redshift','time']
listvstacks = [centers,fillFactors,ravg,rmax]
vstackstrings = ['centers','fillFactors','ravg','rmax']

while count < 5:
	data = cPickle.load(open(timepickles[count],'rb'))
	halonum = halo0[count]
	idx = np.where(data['halonum'] == halonum)[0]
	
	i = 0
	while i < len(listvstacks):
		listvstacks[i] = np.append(listvstacks[i],data[vstackstrings[i]][idx])
		i = i + 1
	i = 0
	while i < len(listappends):
		listappends[i] = np.append(listappends[i],data[appendstrings[i]][idx])
		i = i + 1

	count = count + 1

fillFactors = np.delete(fillFactors,0,0)
ravg = np.delete(ravg,0,0)
rmax = np.delete(rmax,0,0)
centers = np.delete(centers,0,0)

data = {}
data['halonum'] = halo0
data['centers'] = centers
data['rvirs'] = rvirs
data['mvirs'] = mvirs
data['fillFactors'] = fillFactors
data['ravg'] = ravg
data['rmax'] = rmax
data['redshift'] = redshift
data['time'] = time

fileout = 'halo0.cpkl'
cPickle.dump(data,open(fileout,'wb'),protocol=-1)


