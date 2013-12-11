import matplotlib as mpl
mpl.use('agg')
from lauren import *

#need to be able to to loop through each halo and plot! woot
datasets = ['/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062','/u/10/l/lnc2115/vega/data/Wise/DD0135/output_0135','/u/10/l/lnc2115/vega/data/Wise/DD0134/output_0134','/u/10/l/lnc2115/vega/data/Wise/DD0133/output_0133','/u/10/l/lnc2115/vega/data/Wise/DD0059/output_0059','/u/10/l/lnc2115/vega/data/Wise/DD0037/output_0037','/u/10/l/lnc2115/vega/data/Wise/DD0036/output_0036','/u/10/l/lnc2115/vega/data/Wise/DD0031/output_0031']

halosets = ['/u/10/l/lnc2115/vega/data/Wise/groups_02797.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02789.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02784.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02779.dat','/u/10/l/lnc2115/vega/data/Wise/groups_02685.dat','/u/10/l/lnc2115/vega/data/Wise/groups_01260.dat','/u/10/l/lnc2115/vega/data/Wise/groups_01172.dat','/u/10/l/lnc2115/vega/data/Wise/groups_00751.dat','/u/10/l/lnc2115/vega/data/Wise/groups_00399.dat']

timesteps=['0062','0135','0134','0133','0059','0037','0036','0031']

halo0 = [0,0,0,0,0]
halo3 = [3,5,6,6,6]
halolist = [halo3]#, halo3]

for haloids in halolist:
	count = 0
	while count < 5:
		print timesteps[count]
		halonum = haloids[count]
		picklehere = timesteps[count]+'.cpkl'
		data = cPickle.load(open(picklehere,'rb'))
		pf = load(datasets[count])
		idx = np.where(np.array(data['halonum']) == halonum)[0][0]
		print idx
		data_source = pf.h.sphere(data['centers'][idx],data['rvirs'][idx]/pf['cm'])
		rp_r,rp_vals = make_radial_profile(pf,data_source,data['centers'][idx],data['rvirs'][idx]/pf['cm']*pf['kpc'],25.)
		#plt.plot(rp_r,rp_vals,linewidth=2.,label=timesteps[count])
		
		#######################################
		######### found in yt cookbook ########
		#######################################

		profile = BinnedProfile1D(data_source,25,"Radiuskpc",0.1,data['rvirs'][idx]/pf['cm']*pf['kpc'])
		profile.add_fields('SolarMetals')
		plt.loglog(profile['Radiuskpc'],profile['SolarMetals'],label='Mean')
		plt.loglog(profile['Radiuskpc'],profile['SolarMetals_std'],label='std')


		count = count + 1
	plt.xlabel('Radius (kpc)')
	splt.ylabel('Metallicity')
	plt.legend()
	plt.xlim(0,3.5)
	plt.ylim(-3.5,0)
	#outfile = 'halo'+str(haloids[0])+'rp_same.png'
	outfile = 'halo'+str(haloids[0])+'ytrp.png'
	plt.savefig(outfile)
	plt.close()




