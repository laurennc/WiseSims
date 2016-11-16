from lauren import *
import itertools
from matplotlib import rc
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import norm
import matplotlib.mlab as mlab
import math
import cPickle

def build_values(pf,lastStepData):
	avg_logZ_stars_mass_recent = []
	avg_logZ_gas_weightvol_50r = []
	halonums = []
	predicted_Z_total = []
	predicted_Z_50r = []	

	for i in range(len(lastStepData['halonum'])):
		print lastStepData['halonum'][i]
		data_region = pf.h.sphere(lastStepData['centers'][i],lastStepData['rvirs'][i]/pf['cm'])
		gas_mass_total,stellar_mass_total = np.sum(data_region['CellMassMsun']),np.sum()
		idx = np.where(make_solar_metallicity(metal) > -4.)[0]
		sm,ct,metal = sm[idx],ct[idx],metal[idx]

		data_regsmall = pf.h.sphere(lastStepData['centers'][i],lastStepData['rvirs'][i]*0.5/pf['cm'])
		sm_small,ct_small,metal_small = sfr_quants(data_regsmall)
		gas_mass_small,stellar_mass_small = np.sum(data_regsmall['CellMassMsun']),np.sum(sm_small)

	
		recent_cut = pf.current_time - 50e6*YEAR/pf['Time']
		recent_ids = np.where(ct >= recent_cut)[0]
		metal = make_solar_metallicity(metal)

		if len(recent_ids > 0):
			avg_logZ_stars = np.append(avg_logZ_stars,np.average(metal[recent_ids]))
			avg_logZ_gas_weightvol_50r = np.append(avg_logZ_gas_weightvol_50r, np.average(np.log10(data_regsmall['Total_Metallicity']),weights=data_regsmall['CellVolume'] ))
			predicted_Z_totalTEMP = np.log10(-0.01*np.log(gas_mass_total/(gas_mass_total+stellar_mass_total)))-np.log10(Z_solar)
			predicted_Z_50rTEMP = np.log10(-0.01*np.log(gas_mass_small/(gas_mass_small+stellar_mass_small)))-np.log10(Z_solar)
			predicted_Z_total = np.append(predicted_Z_total,predicted_Z_totalTEMP)
			predicted_Z_50r = np.append(predicted_Z_50r,predicted_Z_50rTEMP)
			halonums = np.append(halonums,lastStepData['halonum'][i])


return avg_logZ_stars_mass_recent,avg_logZ_gas_weightvol_50r,predicted_Z_total,predicted_Z_50r,halonums


pf = load('/media/caldisk/lauren/data/Wise/DD0062/output_0062')
lastStepData = cPickle.load(open('/media/caldisk/lauren/WiseAnalysis/halo_pckls/0062.cpkl','rb'))

avg_logZ_stars_mass_recent,avg_logZ_gas_weightvol_50r,predicted_Z_total,predicted_Z_50r,halonums = build_values(pf,lastStepData)

data = {}

data['avg_logZ_stars_mass_recent'] = avg_logZ_stars_mass_recent
data['avg_logZ_gas_weightvol_50r'] = avg_logZ_gas_weightvol_50r
data['predicted_Z_total'] = predicted_Z_total
data['predicted_Z_50r'] = predicted_Z_50r
data['halonums'] = halonums


fileout = 'recentstars_halo_metal_total.cpkl'
cPickle.dump(data,open(fileout,'wb'),protocol=-1)



