from yt.mods import *
from myanyl import *

halonum = 5

pf = load('DD0062/output_0062')
halos = readhalos(foffile='groups_02797.dat')

# Calculate the virial mass, radius, and halo center
radius, mass, center = r200(pf, halos["pos"][halonum,:], halos["mass"][halonum],
                            verbose=True)

# Plot a slice of density, temperature, electron fraction, and
# metallicity (Pop II and III separately) through the halo center with
# a field of view of 10*rvir
#fields = ["Metallicity", "Metallicity3"]
fields = ["Density", "Temperature", "Electron_Fraction", "Metallicity", "Metallicity3"]
#fields = ["Density", "Temperature", "Electron_Fraction"]
zlim = {"Density": (), 
        "Temperature": (), 
        "Electron_Fraction": (), 
        "Metallicity": (1e-6,1), 
        "Metallicity3": (1e-6,1)}
for dim in "xyz":
    SlicePlot(pf, dim, fields, center, width = (10*radius, "cm"),axes_unit=["kpc","kpc"]).save()
    ProjectionPlot(pf, dim, fields, center, weight_field="Density",width = (10*radius, "cm"), axes_unit=["kpc","kpc"]).save()

# Plot a radial profile for the same quantities
pc = PlotCollection(pf, center=center)
for f in fields:
    pc.add_profile_sphere(radius, "cm", ["Radiuspc", f], weight="CellMass")
pc.save()

