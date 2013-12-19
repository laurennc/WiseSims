from lauren import *

data = cPickle.load(open('/u/10/l/lnc2115/vega/repos/WiseSims/clump_dict.cpkl','rb'))

grp1 = [0,1,14,15,19,20,23]
grp2 = [4,6,7,11,21]
grp3 = [5,13]
maxs = [0,4,5]
#independent
indp = [2,3,8,9,10,16,17,22]
#these ones are weird
skip = [12,18]
#for the plots -- since all the grps have the same clump, should just be able to use one of them to get the clump quantities that I need
good_idx = [0,4,5,2,3,8,9,10,16,17,22]

#tnow = pf.current_time*pf['Time']/YEAR
tnow = 745464315.387

r_scale = return_snwinds_scaled_radii(data['t_start'],data['star_masses'])

lauren_r = [1.47619,1.49896,1.43802,0.840769,1.17984,1.10524,0.788054,0.788054,0.155246,0.150740,0.675077,1.32790,0.0565202,0.899243,1.40328,0.710530,0.443687,0.386675,0.163982,0.483630,0.308490,0.788054,0.488499,0.292885]

ravg = data['avgRadius']/r_scale
rmax = data['maxRadius']/r_scale

ravgcomp = ravg/lauren_r
rmaxcomp = rmax/lauren_r

plt.plot(np.log10(data['star_masses'])[indp],ravgcomp[indp],'rs',label='Indp')
plt.plot(np.log10(data['star_masses'])[grp1],ravgcomp[grp1],'bo',label='Groups')
plt.plot(np.log10(data['star_masses'])[grp2],ravgcomp[grp2],'bo')
plt.plot(np.log10(data['star_masses'])[grp3],ravgcomp[grp3],'bo')
plt.plot(np.log10(data['star_masses'])[maxs],ravgcomp[maxs],'g*',markersize=13.5,label='Max Grp Mem')
#plt.plot(np.arange(6)+2.5,np.zeros(6)+1,'k--',linewidth=1.5)
plt.axhline(y=1,linewidth=1.5,linestyle='--',color='k')
plt.xlabel('log(Stellar Mass)')
plt.ylabel('R Avg / Tegmark')
plt.legend()
plt.savefig('best_tegmark_compare_avg.eps')
plt.close()

plt.plot(np.log10(data['star_masses'])[indp],rmaxcomp[indp],'rs',label='Indp')
plt.plot(np.log10(data['star_masses'])[grp1],rmaxcomp[grp1],'bo',label='Groups')
plt.plot(np.log10(data['star_masses'])[grp2],rmaxcomp[grp2],'bo')
plt.plot(np.log10(data['star_masses'])[grp3],rmaxcomp[grp3],'bo')
plt.plot(np.log10(data['star_masses'])[maxs],rmaxcomp[maxs],'g*',markersize=13.5,label='Max Grp Mem')
#plt.plot(np.arange(6)+2.5,np.zeros(6)+1,'k--',linewidth=1.5)
plt.axhline(y=1,linewidth=1.5,linestyle='--',color='k')
plt.xlabel('log(Stellar Mass)')
plt.ylabel('R Max / Tegmark')
plt.savefig('best_tegmark_compare_max.eps')
plt.close()

