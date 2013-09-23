from lauren import *

data1 = cPickle.load(open('0062.cpkl','rb'))
data2 = cPickle.load(open('0135.cpkl','rb'))
data3 = cPickle.load(open('0134.cpkl','rb'))
data4 = cPickle.load(open('0133.cpkl','rb'))
data5 = cPickle.load(open('0059.cpkl','rb'))
data6 = cPickle.load(open('0037.cpkl','rb'))

datas = [data1,data2,data3,data4,data5]
colors = ['r','g','b']

for j in range(5):
	datanow = datas[j]
	for i in range(3):
		plt.plot(datanow['redshift'],datanow['fillFactors'][i][4],colors[i],marker='o',linestyle='')

plt.xlabel('Redshift')
plt.ylabel('Filling Factor of -2')
plt.savefig('fillfactor_vs_t.png')

plt.close()

for j in range(5):
        datanow = datas[j]
        for i in range(3):
                #rnow = datanow['rmax'[i][4]
		#rs = [rnow,rnow,rnow,rnow,rnow]
		plt.plot(datanow['ravg'][i,:],datanow['fillFactors'][i,:],colors[i],marker='o',linestyle='')
	plt.xlabel('Avg Radius (kpc)')
	plt.ylabel('Filling Factor')
	fileout = 'fillfactor_vs_rmax'+str(j)+'.png'
	plt.savefig(fileout)
	plt.close()
