from lauren import *

data = cPickle.load(open('/u/10/l/lnc2115/vega/repos/WiseSims/clump_dict_Z6.cpkl','rb'))

#grp1 = [0,1,14,15,19,20,23]
#grp2 = [4,6,7,11,21]
#grp3 = [5,13]
grps = [1,14,15,19,20,23,4,6,7,11,21,13]
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

#THESE OLD VALUES CAME FROM USING THE WRONG MASS IN THE SNRATES FILE
#lauren_r = [1.47619,1.49896,1.43802,0.840769,1.17984,1.10524,0.788054,0.788054,0.155246,0.150740,0.675077,1.32790,0.0565202,0.899243,1.40328,0.710530,0.443687,0.386675,0.163982,0.483630,0.308490,0.788054,0.488499,0.292885]

#THESE ARE THE VALUES I GET WHEN I BELIEVE I'VE CORRECTED IT
lauren_r = [2.78808,2.83108,2.71597,1.58790,2.22833,2.08742,1.48833,1.48833,0.292994,0.284480,1.27494,2.50799,0.106331,1.69834,2.65037,1.34190,0.837878,0.730189,0.309501,0.913325,0.582502,1.48833,0.922523,0.553023]

ravg = data['avgRadius']/r_scale
#rmin = data['zMinRadius']/r_scale
#rmax = data['zMaxRadius']/r_scale
rmin = data['minRadius']/r_scale
rmax = data['maxRadius']/r_scale

ravgcomp = ravg/lauren_r
rmaxcomp = rmax/lauren_r
rmincomp = rmin/lauren_r

fig,ax = plt.subplots(nrows=1,ncols=2)
fig.set_size_inches(11,5)

ax[0].plot(np.log10(data['star_masses'])[indp],ravgcomp[indp],'rs',label='Indp')
ax[0].errorbar(np.log10(data['star_masses'])[indp],ravgcomp[indp],yerr=[rmincomp[indp],rmaxcomp[indp]],color='red',fmt='.')
ax[0].plot(np.log10(data['star_masses'])[grps],ravgcomp[grps],'bo',label='Groups')
ax[0].errorbar(np.log10(data['star_masses'])[grps],ravgcomp[grps],yerr=[rmincomp[grps],rmaxcomp[grps]],color='blue',fmt='.')
ax[0].plot(np.log10(data['star_masses'])[maxs],ravgcomp[maxs],'g*',markersize=13.5,label='Max Grp Mem')
ax[0].errorbar(np.log10(data['star_masses'])[maxs],ravgcomp[maxs],yerr=[rmincomp[maxs],rmaxcomp[maxs]],color='green',fmt='.')
ax[0].axhline(y=1,linewidth=1.5,linestyle='--',color='k')
ax[0].set_xlabel('log(Stellar Mass)')
ax[0].set_ylabel('R Avg / Tegmark')
ax[0].legend()

ax[1].plot(np.log10(data['star_masses'])[indp],ravgcomp[indp],'rs',label='Indp')
ax[1].errorbar(np.log10(data['star_masses'])[indp],ravgcomp[indp],yerr=[rmincomp[indp],rmaxcomp[indp]],color='red',fmt='.')
ax[1].plot(np.log10(data['star_masses'])[grps],ravgcomp[grps],'bo',label='Groups')
ax[1].errorbar(np.log10(data['star_masses'])[grps],ravgcomp[grps],yerr=[rmincomp[grps],rmaxcomp[grps]],color='blue',fmt='.')
ax[1].plot(np.log10(data['star_masses'])[maxs],ravgcomp[maxs],'g*',markersize=13.5,label='Max Grp Mem')
ax[1].errorbar(np.log10(data['star_masses'])[maxs],ravgcomp[maxs],yerr=[rmincomp[maxs],rmaxcomp[maxs]],color='green',fmt='.')
ax[1].axhline(y=1,linewidth=1.5,linestyle='--',color='k')
ax[1].set_xlabel('log(Stellar Mass)')
ax[1].set_ylabel('R Avg / Tegmark')
ax[1].set_yscale('log')
#ax[1].legend()

#plt.savefig('best_tegmark_compare_avg_fixed_werrs.eps')
plt.savefig('best_tegmark_compar_Z6.eps')
plt.close()


