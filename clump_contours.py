class Clump_Contours:
	def __init__(self,pf, center, radius, field, clump, width)
		self.pf = pf
		self.center = center
		self.radius = radius
		self.field = field
		self.clump = clump
		self.width = width
		
		self.pc = PlotCollection(self.pf,center=self.center)
	
		#everything will have to do with plots so make the plotcollection object here too!
		
	
	def make_main_plot(self, plottype, dim):
		if self.field='Density':
			self.pc.set_zlim(1e-27,1e-22)
		
		if (plottpe = 'proj'):
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
		self.pc.plots[0].modify['point'](self.center,'c')
		self.pc.plots[0].modify['point'](self.clump.quantities['CenterOfMass'](),'x')
		self.pc.plots[0].modify['sphere'](self.center,rad)
		self.pc.set_cmap('spectral')
		return True
		
	def save_plot(self,fileout):
		self.pc.save(fileout)
		return True
	
	def mark_all_halos(self,halos):
		i = 0
		while i lt len(halos['masses']):
			offset = (((self.center[0]-halos['centers'][i][0])**2.0+(self.center[1]-halos['centers'][i][1])**2.0+
						(self.center[2]-halos['centers'][i][2])**2.0)**0.5)*self.pf['kpc']
			if (offset < self.width):
				self.pc.plots[0].modify['point'](halos['centers'][i],halos['halonum'][i])
				self.pc.plots[0].modify['sphere'](halos['centers'][i],halos['rvirs'][i]/self.pf['cm'])
			i = i + 1
		return True
	
	