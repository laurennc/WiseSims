import matplotlib as mpl
mpl.use('agg')

from lauren import *

metallicityList = [-6,-5,-4,-3,-2,-1,0]
indp = [2,3,8,9,10,16,17,22]
grps = [0,1,14,15,19,20,23,4,6,7,11,21,5,13]

data = cPickle.load(open('clump_dict.cpkl','rb'))

i  = 0
colors = ['m','b','c','k','g','y','r']
fillingfactors = np.array(data['fillFactors'])

fig,axes = plt.subplots(2,1)

while i < 6:#7:
	label_string = str(metallicityList[i])
	axes[0].plot(data['star_masses'][indp],fillingfactors[indp,i],
		colors[i],marker='o',label=label_string,linestyle='')
	axes[1].plot(data['star_masses'][grps],fillingfactors[grps,i],
		colors[i],marker='s',label=label_string,linestyle='')
	i = i + 1

for j in [0,1]:
	axes[j].set_xscale('log')  
	axes[j].set_xlabel('log(Stellar Mass)')
	axes[j].set_ylabel('Filling Factor')
	axes[j].legend()
fig.savefig('filling_factors.png')
 
