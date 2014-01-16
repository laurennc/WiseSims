from lauren import *

halo0 = [0,0,0,0,0]
halo3 = [3,5,6,6,6]


timepickles = ['0062.cpkl','0135.cpkl','0134.cpkl','0133.cpkl','0059.cpkl']

count = 0

centers = np.array([0.,0.,0.])
rvirs = []
mvirs = []
fillFactors = np.array([0.,0.,0.,0.,0.,0.,0.])
massFactors = np.array([0.,0.,0.,0.,0.,0.,0.])
ravg = np.array([0.,0.,0.,0.,0.,0.,0.])
rmax = np.array([0.,0.,0.,0.,0.,0.,0.])
redshift = [] 
time = []

appendstrings = ['rvirs','mvirs','redshift','time']
vstackstrings = ['centers','fillFactors','massFactors','ravg','rmax']

while count < 5:
	data = cPickle.load(open(timepickles[count],'rb'))
	halonum = halo3[count]
	idx = np.where(np.array(data['halonum']) == halonum)[0]

	rvirs = np.append(rvirs,data['rvirs'][idx])
	mvirs = np.append(mvirs,data['mvirs'][idx])
	redshift = np.append(redshift,data['redshift'])
	time = np.append(time,data['time'])
	
	centers = np.vstack((centers,data['centers'][idx]))
	fillFactors = np.vstack((fillFactors,data['fillFactors'][idx]))
	massFactors = np.vstack((massFactors,data['massFactors'][idx]))
	ravg = np.vstack((ravg,data['ravg'][idx]))
	rmax = np.vstack((rmax,data['rmax'][idx]))

	count = count + 1

fillFactors = np.delete(fillFactors,0,0)
massFactors = np.delete(massFactors,0,0)
ravg = np.delete(ravg,0,0)
rmax = np.delete(rmax,0,0)
centers = np.delete(centers,0,0)

data = {}
data['halonum'] = halo3
data['centers'] = centers
data['rvirs'] = rvirs
data['mvirs'] = mvirs
data['fillFactors'] = fillFactors
data['massFactors'] = massFactors
data['ravg'] = ravg
data['rmax'] = rmax
data['redshift'] = redshift
data['time'] = time

fileout = 'halo3.cpkl'
cPickle.dump(data,open(fileout,'wb'),protocol=-1)


