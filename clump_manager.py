from yt.mods import *
from yt.analysis_modules.level_sets.api import *
from yt.analysis_modules.star_analysis.api import *
from myanyl import *
import cPickle
import numpy as np
import matplotlib.pyplot as plt


class ClumpManager:
	def __init__(self,pf, halonum, center, radius, width, halos):
		self.pf = pf
		self.halonum = halonum
		self.center = center
		self.radius = radius
		self.width = width

		self.cindex = np.where(halos['halonum'] == self.halonum)
		self.cindex = self.cindex[0][0]
		#add a plot collection so all the methods can simply edit the plots
		self.pc = 'woot'#PlotCollection(self.pf,center=self.center)
	
		#do the clump finding in here!
		filein = '../WiseSimsData/pickles/halo'+str(self.halonum)+'_clumps.cpkl'
		data2 = cPickle.load(open(filein,'rb'))
		master_clump = data2[1]
		self.all_clumps = get_lowest_clumps(master_clump)
		self.clump = self.all_clumps[int(halos['clump_index'][self.cindex])]
	
	def make_main_plot(self, plottype, dim, field, unit):
		self.pc = PlotCollection(self.pf,center=self.center)
		if field=='Density':
			self.pc.set_zlim(1e-27,1e-22)
		
		if (plottype == 'proj'):
			self.pc.add_projection(field,dim)
		elif (plottype == 'slice'):
			self.pc.add_slice(field,dim)
		else:
			print 'Plot type should be either proj or slice'
			return 
		
		self.pc.set_width(self.width,'kpc')
		clump_plot=[self.clump]
		rad = self.radius/self.pf['cm']
		
		self.pc.plots[unit].modify['clumps'](clump_plot)
		self.pc.plots[unit].modify['point'](self.center,'o')
		self.pc.plots[unit].modify['point'](self.clump.quantities['CenterOfMass'](),'x')
		self.pc.plots[unit].modify['sphere'](self.center,rad)
		self.pc.set_cmap('spectral')
		return 
		
	def save_plot(self,fileout):
		self.pc.save(fileout)
		return 
	
	def mark_all_halos(self,halos,unit):
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
				self.pc.plots[unit].modify['point'](halos['centers'][i],str(halos['halonum'][i]))
				self.pc.plots[unit].modify['sphere'](halos['centers'][i],halos['rvirs'][i]/self.pf['cm'])
			#else:
			#	pass
			i = i + 1
		return 
	
	def where_outside_rvir(metal_limit):
		rad = self.radius/self.pf['cm']
		dist = ((self.clump['x']-self.center[0])**2.0+(self.clump['y']-self.center[1])**2.0+(self.clump['z']-self.center[2])**2.0)**0.5
		idx = np.where(self.clump['Metallicity'] >= metal_limit)
		return np.where(dist[idx] > rad)

	def total_clump_quantity(indices,keyin):
		quantity_wnted = self.clump[keyin]
		return quantity_wanted[indices].sum()
	
