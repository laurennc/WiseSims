import matplotlib as mpl
mpl.use('agg')
from clump_manager import *
from lauren import *

YEAR = 3.155693e+7 #s/yr
all_halos = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,36,42,48,59,188,362,900]
good_halos = [2,3,4,9,16,36,362]

data = cPickle.load(open('clump_dict.cpkl','rb'))
good_idx = [2,3,4,8,10,17,22]

#tnow = pf.current_time*pf['Time']/YEAR
tnow = 745464315.387

#First I'll worry about the dilution model since that's so much easier
yDilMod = dilution_model(data['t_start'],tnow,5.0e9,5.0e8)
yDilModVir = dilution_model(data['t_start'],tnow,1.0e7,1.0e6)

plt.plot(data['t_start'][good_idx],yDilMod[good_idx],'ro',label='Jasons Model')
plt.plot(data['t_start'][good_idx],data['dilmass6'][good_idx],'bs',label='1e-6 Metal Cut')#,markersize=4)
plt.plot(data['t_start'][good_idx],data['dilmass17'][good_idx],'g^',label='1e-17 Metal Cut')#,markersize=8)
plt.xlabel('Start Time of SF')
plt.ylabel('Dilution Mass')
plt.legend(loc='lower left')
plt.xscale('log')
plt.yscale('log')
plt.savefig('dilution_model_good.png')

plt.plot(data['t_start'],yDilMod,'ro',label='Jasons Model')
plt.plot(data['t_start'],data['dilmass6'],'bs',label='1e-6 Metal Cut')#,markersize=4)
plt.plot(data['t_start'],data['dilmass17'],'g^',label='1e-17 Metal Cut')#,markersize=8)
plt.xlabel('Start Time of SF')
plt.ylabel('Dilution Mass')
plt.legend(loc='lower left')
plt.xscale('log')
plt.yscale('log')
plt.savefig('dilution_model_all.png')

fig,axes = plt.subplots(2,1)
axes[0].plot(data['t_start'],yDilModVir,'ro',label='Jasons Model')
axes[0].plot(data['t_start'],data['dilmassvir'],'bs',label='Inside Rvir')
axes[0].set_xlabel('Start Time of SF')
axes[0].set_ylabel('Dilution Mass')
axes[0].legend(loc='lower left')
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].plot(data['t_start'][good_idx],yDilModVir[good_idx],'ro',label='Jasons Model')
axes[1].plot(data['t_start'][good_idx],data['dilmassvir'][good_idx],'bs',label='Inside Rvir')
axes[1].set_xlabel('Start Time of SF')
axes[1].set_ylabel('Dilution Mass')
axes[1].legend(loc='lower left')
axes[1].set_xscale('log')
axes[1].set_yscale('log')
fig.savefig('dilution_model_vir.png')


#Remake the Tegmark Plot
lauren_r = [1.47619,1.49896,1.43802,0.840769,1.17984,1.10524,0.788054,0.788054,0.155246,0.150740,0.675077,1.32790    0.0565202,0.899243,1.40328,0.710530,0.443687,0.386675,0.163982,0.483630,0.308490,0.788054,0.488499,0.292885]

tegmark_r = [2.46459, 2.49453, 2.41438, 1.63200, 2.07518, 1.97735, 1.56361, 1.56361,0.805612,0.801571, 1.41791, 2.26963,0.732799, 1.70808, 2.36871, 1.46348, 1.12637, 1.05719,0.813593, 1.17566,0.965485, 1.56361,1.18171,0.947766]

lauren_r = np.array(lauren_r)
tegmark_r = np.array(tegmark_r)

r_scaled = return_snwinds_scaled_radii(data['t_start'],data['star_masses'])
r_scaled = data['avgRadius']/r_scaled
plt.plot(data['star_masses'][good_idx],r_scaled[good_idx],'g^',label='Average Radius')
plt.plot(data['star_masses'][good_idx],lauren_r[good_idx],'bs',label='Lauren SFH')
plt.plot(data['star_masses'][good_idx],tegmark_r[good_idx],'ro',label='Tegmark SFH')
plt.xscale('log')
plt.xlabel('log(Stellar Mass)')
plt.ylabel('Scaled Radius')
plt.legend(loc='upper left')
plt.savefig('snwinds_model_good.png')


