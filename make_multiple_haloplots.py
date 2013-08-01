import matplotlib as mpl
mpl.use('agg')
from clump_manager import *

#13 was deleted?
halos_to_run = [0,1,2,3,4,5,6,7,9,15,16,17,22,25,26,33,35,
                36,42,48,59,188,362,900]
#halos_to_run = [0]
count = 0

pf = load('/u/10/l/lnc2115/home/WiseSimsData/DD0062/output_0062')
halos = readhalos(foffile='/u/10/l/lnc2115/home/WiseSimsData/groups_02797.dat')
data = cPickle.load(open('clump_dict.cpkl','rb'))

while (count < len(halos_to_run)):
	halonum = halos_to_run[count]
	print 'Analyzing halo: ',str(halonum)
	cm = ClumpManager(pf,halonum,data['centers'][count],data['rvirs'][count],25.,data)
	fileout = '/u/10/l/lnc2115/home/WiseSims/analysis_plots/all_marked/halo'+str(halonum)+'_allmarked'
	num = 0
	for dim in "xyz":
		cm.make_main_plot('slice',dim,'Density',num)
		cm.mark_all_halos(data,num)
		cm.save_plot(fileout)
		#num = num + 1

	for dim in "xyz":
		cm.make_main_plot('slice',dim,'Metallicity',num)
		cm.mark_all_halos(data,num)
		cm.save_plot(fileout)
		#num = num + 1
	cm.save_plot(fileout)
	count = count + 1



