from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt


class Clump_Manager:
	def __init__(self,pf, halonumber, center, radius, width, halos):
		self.pf = pf
		self.halonum = halonumber
		self.center = center
		self.radius = radius
		self.width = width

		self.index = np.where(halos['halonum'] == halonumber)
		self.index = self.index[0]
		
		#add a plot collection so all the methods can simply edit the plots
		self.pc = PlotCollection(self.pf,center=self.center)
	
		#do the clump finding in here!
		filein = '../WiseSimsData/pickles/halo'+str(halonum)+'_clumps.cpkl'
		data2 = cPickle.load(open(filein,'rb'))
		master_clump = data2[1]
		self.all_clumps = get_lowest_clumps(master_clump)
		self.clump = self.all_clumps[halos['clump_index'][self.index]]
	
	def make_main_plot(self, plottype, dim, field):
		if field=='Density':
			self.pc.set_zlim(1e-27,1e-22)
		
		if (plottype == 'proj'):
			self.pc.add_projection(self.field,dim)
		elif (plottype == 'slice'):
			self.pc.add_projection(self.field,dim)
		else:
			print 'Plot type should be either proj or slice'
			return False
		
		self.pc.set_width(self.width,'kpc')
		clump_plot=[self.clump]
		rad = self.radius/self.pf['cm']
		
		self.pc.plots[0].modify['clumps'](clump_plot)
		self.pc.plots[0].modify['point'](self.center,'o')
		self.pc.plots[0].modify['point'](self.clump.quantities['CenterOfMass'](),'x')
		self.pc.plots[0].modify['sphere'](self.center,rad)
		self.pc.set_cmap('spectral')
		return True
		
	def save_plot(self,fileout):
		self.pc.save(fileout)
		return True
	
	def mark_all_halos(self,halos):
		i = 0
		#print self.center
		while i < len(halos['masses']):
			offset = (((self.center[0]-halos['centers'][i][0])**2.0+(self.center[1]-halos['centers'][i][1])**2.0+
						(self.center[2]-halos['centers'][i][2])**2.0)**0.5)*self.pf['kpc']
			#if (self.halonum == halos['halonum'][i]):
			#	pass
			if (offset < self.width):
				#print 'adding halonumber ',halos['halonum'][i]
				#print halos['centers'][i],halos['halonum'][i]
				self.pc.plots[0].modify['point'](halos['centers'][i],str(halos['halonum'][i]))
				self.pc.plots[0].modify['sphere'](halos['centers'][i],halos['rvirs'][i]/self.pf['cm'])
			#else:
			#	pass
			i = i + 1
		return True
	
	
