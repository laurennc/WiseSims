from lauren import *

list62 = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59]
list135 = [0,3,4,5,6,7,10,16,18,27,28,29,33,34,35,43,51,57]
list134 = [0,3,4,5,6,7,10,16,19,28,29,30,34,35,36,41,57]
list133 = [0,3,4,5,6,7,10,16,18,20,29,30,35,37,38,44,59]
list59 = [0,2,5,6,7,15,17,20,27,29,33,36,39,42,55]
list37 = [0,2,5,6,9,10,13,15,16,19,27,30,38,42,51,60,228,247]
list36 = [] 
list31 = []

#list62 = [0,26,3,9]
#list135 = [0,29,5,10]
#list134 = [0,30,6,10]
#list133 = [0,30,6,10]
#list59 = [0,29,6]
#list37 = [0,10,6]

master_index_list = [list62,list135,list134,list133,list59,list37,list36,list31]

datasets = ['/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062','/u/10/l/lnc2115/vega/data/Wise/DD0135/output_0135','/u/10/l/lnc2115/vega/data/Wise/DD0134/output_0134','/u/10/l/lnc2115/vega/data/Wise/DD0133/output_0133','/u/10/l/lnc2115/vega/data/Wise/DD0059/output_0059','/u/10/l/lnc2115/vega/data/Wise/DD0037/output_0037','/u/10/l/lnc2115/vega/data/Wise/DD0036/output_0036','/u/10/l/lnc2115/vega/data/Wise/DD0031/output_0031']

halosets = ['/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02789.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02784.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02779.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02685.dat','/u/10/l/lnc2115/vega/data/Wise/groups_01260.dat','/u/10/l/lnc2115/vega/data/Wise/groups_01172.dat','/u/10/l/lnc2115/vega/data/Wise/groups_00751.dat','/u/10/l/lnc2115/vega/data/Wise/groups_00399.dat']

timesteps=['0062','0135','0134','0133','0059','0037','0036','0031']

count = 0

metallicityList=[-6,-5,-4,-3,-2,-1]

#need to assemble some sort of dictionary that links timesteps to a dictionary of values!


#while count < len(datasets):
while count < 5:
	timehere = timesteps[count]
	pf = load(datasets[count])
	halos = readhalos(foffile=halosets[count])


	mvirs=[]
	centers = np.array([0.,0.,0.])
	fillFactors=np.array([0.,0.,0.,0.,0.,0.,0.])	
	massFactors=np.array([0.,0.,0.,0.,0.,0.,0.])
	rmax = np.array([0.,0.,0.,0.,0.,0.,0.])
	ravg = np.array([0.,0.,0.,0.,0.,0.,0.])
	rvirs = []

	for i in master_index_list[count]:
		print timehere, i
		rvir,mvir,center = r200(pf,halos['pos'][i,:],halos['mass'][i],verbose=False)
		mvirs = np.append(mvirs,mvir)
		rvirs = np.append(rvirs,rvir)
		centers = np.vstack((centers,center))		

		fields = ["Density","SolarMetals"]
		zlim = {"Density":(), "SolarMetals":(-6,1)}

		#for dim in "xyz":
	                #pc = make_halo_plot(pf,center,rvir,15.,'slice',dim,'SolarMetals')
		#	sp = SlicePlot(pf,dim,fields,center,width=(10.,'kpc'),axes_unit=["kpc","kpc"])
		#	sp.annotate_particles(0.05,p_size=2.5,stars_only=True)
		#	fileout='starparticles/'+str(timehere)+'_halo'+str(i)
		#	sp.save(fileout)

		#for dim in "xyz":
			#pc = make_halo_plot(pf,center,rvir,15.,'slice',dim,'Density')
		#	fileout='fillfactors/'+str(timehere)+'_halo'+str(i)
                 #       pc.save(fileout)

	
		data_source = pf.h.sphere(center,rvir/pf['cm'])
		fillFactorsTemp = np.array([0.,0.,0.,0.,0.,0.,0.])
		massFactorsTemp = np.array([0.,0.,0.,0.,0.,0.,0.])
		rmaxTemp = np.array([0.,0.,0.,0.,0.,0.,0.])
		ravgTemp = np.array([0.,0.,0.,0.,0.,0.,0.])
		
		metallicities = data_source['SolarMetals']
		volumes = data_source['CellVolume']
		cellmasses = data_source['CellMassMsun']
		total_mass = cellmasses.sum()
		total_volume = volumes.sum()
		radii = distance_from_center(data_source['x'],data_source['y'],data_source['z'],center)*pf['kpc']

		for loopthrough in range(6):
			#print loopthrough
			goodval = metallicityList[loopthrough]
			indices = np.where(metallicities >= goodval)
			fillFactorsTemp[loopthrough] = volumes[indices].sum()/total_volume
			massFactorsTemp[loopthrough] = cellmasses[indices].sum()/total_mass
			if len(indices[0]) == 0:
				print 'ERROR: INDICES IN FILL FACTORS LOOP ARE ZERO'
			elif i== 22:
				print '22 is a pain' 
			else:	
				rmaxTemp[loopthrough] = np.max(radii[indices])
				ravgTemp[loopthrough] = np.average(radii[indices])

		fillFactors = np.vstack((fillFactors,fillFactorsTemp))
		massFactors = np.vstack((massFactors,massFactorsTemp))
		rmax = np.vstack((rmax,rmaxTemp))
		ravg = np.vstack((ravg,ravgTemp))

		
	#here make a pickle for each timestep for now
	fillFactors = np.delete(fillFactors,0,0)
	massFactors = np.delete(massFactors,0,0)
	ravg = np.delete(ravg,0,0)
	rmax = np.delete(rmax,0,0)
	centers = np.delete(centers,0,0)

	data = {}
	data['halonum'] = master_index_list[count]
	data['centers'] = centers
	data['rvirs'] = rvirs
	data['mvirs'] = mvirs
	data['fillFactors'] = fillFactors
	data['massFactors'] = massFactors
	data['ravg'] = ravg
	data['rmax'] = rmax
	data['redshift'] = pf.current_redshift
	data['time'] = pf.current_time*pf['Time']/(3.15569*10**7.0)

	fileout = timehere+'.cpkl'
	cPickle.dump(data,open(fileout,'wb'),protocol=-1)	
	
	count = count + 1



