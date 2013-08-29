import matplotlib as mpl
mpl.use('agg')

from lauren import *

data = cPickle.load(open('clump_dict.cpkl','rb'))

grp1 = [0,1,14,15,19,20,23]
grp2 = [4,6,7,11,21]
grp3 = [5,13]
#independent
indp = [2,3,8,9,10,16,17,22]
#these ones are weird
skip = [12,18]
#for the plots -- since all the grps have the same clump, should just be able to use one of them to get the clump quantities that I need
good_idx = [0,4,5,2,3,8,9,10,16,17,22]

#tnow = pf.current_time*pf['Time']/YEAR
tnow = 745464315.387

r_scale_grp = np.zeros(3)
r_scale_grp[0] = return_snwinds_scaled_radii(data['t_start'][grp1].min(),data['star_masses'][grp1].sum())
r_scale_grp[1] = return_snwinds_scaled_radii(data['t_start'][grp2].min(),data['star_masses'][grp2].sum())
r_scale_grp[2] = return_snwinds_scaled_radii(data['t_start'][grp3].min(),data['star_masses'][grp3].sum())

r_scale_indp = return_snwinds_scaled_radii(data['t_start'][indp],data['star_masses'][indp])

r_scale_avg = np.zeros(len(good_idx))
r_scale_avg[0] = data['avgRadius'][0]/r_scale_grp[0]
r_scale_avg[1] = data['avgRadius'][4]/r_scale_grp[1]
r_scale_avg[2] = data['avgRadius'][5]/r_scale_grp[2]
r_scale_avg[3:] = data['avgRadius'][indp]/r_scale_indp

r_scale_max = np.zeros(len(good_idx))
r_scale_max[0] = data['maxRadius'][0]/r_scale_grp[0]
r_scale_max[1] = data['maxRadius'][4]/r_scale_grp[1]
r_scale_max[2] = data['maxRadius'][5]/r_scale_grp[2]
r_scale_max[3:] = data['maxRadius'][indp]/r_scale_indp


pltmasses = np.zeros(len(good_idx))
pltmasses[0:3] = [data['star_masses'][grp1].sum(),data['star_masses'][grp2].sum(),data['star_masses'][grp3].sum()]
pltmasses[3:] = data['star_masses'][indp]

plt.plot(pltmasses[0:3],r_scale_avg[0:3],'go',label='Avg Groups')
plt.plot(pltmasses[3:],r_scale_avg[3:],'g^',label='Avg Indp')
plt.plot(pltmasses[0:3],r_scale_max[0:3],'ro',label='Max Groups')
plt.plot(pltmasses[3:],r_scale_max[3:],'r^',label='Max Indp')

r_tegmark = [1.49896,1.32790,1.10524,1.43802,0.840769,0.155246,0.150740,0.675077,0.443687,0.386675,0.488499]

plt.plot(pltmasses[0:3],r_tegmark[0:3],'bo',label='Tegmark Groups')
plt.plot(pltmasses[3:],r_tegmark[3:],'b^',label='Tegmark Indp')

plt.xscale('log')
plt.xlabel('log(Stellar Mass)')
plt.ylabel('Scaled Radius')
plt.legend(loc='upper left')
plt.savefig('snwinds_model_addedmass.png')



#need to figure out which plots I want to make
#now have the scaled data -- need the model predictions
#these should be the values that I need for the wechs.pro code on macha
#print data['t_start'][grp1].min(),data['star_masses'][grp1].sum()
#print data['t_start'][grp2].min(),data['star_masses'][grp2].sum()
#print data['t_start'][grp3].min(),data['star_masses'][grp3].sum()
#print data['t_start'][indp],data['star_masses'][indp]

