import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import numpy as np
import cPickle 

data = cPickle.load(open('global_halo_metal_total.cpkl','rb'))

xlen,ylen = 2,1
fig,ax = plt.subplots(ylen,xlen)#,sharey=True)
fig.set_size_inches(6,3)
ax = ax.flat

fileout = 'closedbox_paper_total_weightmvol_small.pdf'

x = np.arange(-3.5,0.5,0.001)

#ax[0].plot(data['avg_logZ_gas_weightvol'],data['predicted_Z'],'ko')
ax[0].plot(data['predicted_Z'],data['avg_logZ_gas_weightvol_small'],'ko')
ax[0].plot(x,x,'--',color='Grey',linewidth=2)
ax[0].set_ylim(-6,0)
ax[0].set_xlim(-3.5,0)
ax[0].set_xlabel('Predicted Z')
ax[0].set_ylabel('<Z gas>')


#ax[1].plot(data['avg_logZ_stars_mass'],data['predicted_Z'],'ko')
ax[1].plot(data['predicted_Z'],data['avg_logZ_stars_mass'],'ko')
ax[1].plot(x,x,'--',color='Grey',linewidth=2)
ax[1].set_ylim(-6,0)
ax[1].set_xlim(-3.5,0)
ax[1].set_ylabel('<Z stars>')
ax[1].set_xlabel('Predicted Z')

#ax[2].plot(data['avg_logZ_gas_weightvol'],data['avg_logZ_stars_mass'],'ko')
#ax[2].plot(x,x,'--',color='Grey',linewidth=2)
#ax[2].set_xlim(-6,0)
#ax[2].set_ylim(-3.5,0)
#ax[2].set_xlabel('<Z gas>')
#ax[2].set_ylabel('<Z stars>')

## Actual Bursty
#bursty = [3,9,11,15,16,17,26,33,42,59,362]
## Really flat but named bursty for ease
#bursty = [2,8,25,35,36,48]
#idx_burst = []
#for i in range(len(bursty)):
#	wanted = np.where(data['halonums'] == bursty[i])[0][0]
#	idx_burst = np.append(idx_burst,int(wanted))
#print idx_burst

#idx_burst = [int(i) for i in idx_burst]
#print idx_burst
#print len(data['avg_logZ_gas'])

#ax[0].plot(data['avg_logZ_gas'][idx_burst],data['predicted_Z'][idx_burst],'bo')

#ax[1].plot(data['avg_logZ_stars'][idx_burst],data['predicted_Z'][idx_burst],'bo')

#ax[2].plot(data['avg_logZ_gas'][idx_burst],data['avg_logZ_stars'][idx_burst],'bo')


plt.tight_layout()
plt.savefig(fileout)
plt.close()




