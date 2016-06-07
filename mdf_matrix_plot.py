from yt.mods import *
import numpy as np
import cPickle
import matplotlib.pyplot as plt
from lauren import *

##The goal of this plot is to show a series of specifc MDFs to make some points in the paper. So three mass bins! and two classifications -- burst/flat and isolated/group!

fileout = 'mdf_matrix_paper.pdf'
xlen,ylen = 3,2
fig,ax = plt.subplots(ylen,xlen,sharex=True,sharey=True)
fig.set_size_inches(12,6)
ax = ax.flat
i = 0

pf = load('/media/caldisk/lauren/data/Wise/DD0062/output_0062')
lastStepData = cPickle.load(open('/media/caldisk/lauren/WiseAnalysis/halo_pckls/0062.cpkl','rb'))

green_hist_nums = [362,16,2,11,3,9]
purple_hist_nums = [59,42,4,25,8,1]

while i < len(ax):
	idx = np.where(lastStepData['halonum'] == green_hist_nums[i])[0][0]
	data_region = pf.h.sphere(lastStepData['centers'][idx],lastStepData['rvirs'][idx]/pf['cm'])
	sm,ct,metal = sfr_quants(data_region)
	metal = make_solar_metallicity(metal)
	wanted = (metal > -4)
	metal = metal[wanted]
	ax[i].hist(metal,histtype='step',color='MediumSpringGreen',label=str(green_hist_nums[i]),normed=True,range=(-4,1),bins=11)


	idx = np.where(lastStepData['halonum'] == purple_hist_nums[i])[0][0]
        data_region = pf.h.sphere(lastStepData['centers'][idx],lastStepData['rvirs'][idx]/pf['cm'])
        sm,ct,metal = sfr_quants(data_region)
        metal = make_solar_metallicity(metal)
        wanted = (metal > -4)
        metal = metal[wanted]
        ax[i].hist(metal,histtype='step',color='MediumVioletRed',label=str(purple_hist_nums[i]),normed=True,range=(-4,1),bins=11)
	
	i = i + 1


plt.savefig(fileout)
plt.close()


