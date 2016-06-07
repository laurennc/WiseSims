##want to find the time difference between the last Pop III and the first Pop II
#might be an indicator of the Population III triggering a starburst 

import matplotlib
matplotlib.use('Agg')

import numpy as np
import matplotlib.pyplot as plt
import cPickle
from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *

YEAR = 3.155693e7 #sec / yr

def find_time_values():
	pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
	lastStepData = cPickle.load(open('/u/10/l/lnc2115/vega/data/Wise/0062.cpkl','rb'))
	
	first_popII = []
	last_popIII = []
	time_diff = []
	
	
	for i in range(len(lastStepData['rvirs'])):
		data_source = pf.h.sphere(lastStepData['centers'][i],lastStepData['rvirs'][i]/pf['cm'])
	
		sm = data_source['ParticleMassMsun']
		ct = data_source['creation_time']
		metals = data_source['metallicity_fraction']
		stars = (ct > 0)
		ct, sm, metals = ct[stars], sm[stars], metals[stars]
		
		stars  = (sm > 1.)
		sm,ct,metals = sm[stars],ct[stars],metals[stars]
		
		popIII = np.where(np.log10(metals) < -4.)[0]
		popII = np.where(np.log10(metals) > -4.)[0]
		idx_last_popIII = np.where(ct[popIII] == ct[popIII].max())[0]
		idx_first_popII = np.where(ct[popII] == ct[popII].min())[0]
		
		time_difference_here = (ct[popII][idx_first_popII] - ct[popIII][idx_last_popIII])*pf['Time']/YEAR

		t = np.log10(pf['Time']*t/YEAR)

		last_popIII = np.append(last_popIII,ct[popIII][idx_last_popIII]*pf['Time']/YEAR)
		first_popII = np.append(first_popII,ct[popII][idx_first_popII]*pf['Time']/YEAR)
		time_diff = np.append(time_diff,time_difference_here)

	return last_popIII,first_popII,time_diff

def plot_time_diff():
        pf = load('/u/10/l/lnc2115/vega/data/Wise/DD0062/output_0062')
        lastStepData = cPickle.load(open('/u/10/l/lnc2115/vega/data/Wise/0062.cpkl','rb'))

	last_popIII,first_popII,time_diff = find_time_values()

	plt.plot(np.log10(lastStepData['mvirs']),np.log10(time_diff),'bo')
	plt.xlabel('Halo Mass [log(Msun)]')
	plt.ylabel('Time between Pop III/II [log(yr)]')
	plt.savefig('time_diff_popIII_popII.png')
	plt.close()

	return

