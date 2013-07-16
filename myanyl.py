import os, sys
import numpy as na
from yt.mods import *
from yt.utilities.cosmology import *

def readhalos(name=None, cycle=None, foffile=None, maxhalos=-1):
    if name != None:
        lines = open(name).readlines()
        for l in lines:
            if l.startswith("InitialCycle"):
                cycle_num = int(l.split("=")[1])
                break
        foffile = "FOF/groups_%5.5d.dat" % cycle_num
    elif cycle != None:
        foffile = "FOF/groups_%5.5d.dat" % cycle
    elif foffile is None:
        print "Need either to specify parameter file or FOF filename!"
        return None
    print foffile
    if not os.path.exists(foffile):
        print "FOF file %s not found" % foffile
        sys.exit()
    lines = open(foffile).readlines()

    halos_in_file = 0
    for l in lines:
        if not (l.startswith("#") or l.startswith("datavar") or len(l) < 2):
            halos_in_file += 1

    if maxhalos > 0:
        nhalos = min(maxhalos, halos_in_file)
    else:
        nhalos = halos_in_file
    halos = {}
    halos['pos'] = na.zeros((nhalos,3))
    halos['vel'] = na.zeros((nhalos,3))
    halos['npart'] = na.zeros(nhalos, dtype='int')
    halos['mvir'] = na.zeros(nhalos)
    halos['mass'] = na.zeros(nhalos)
    halos['rvir'] = na.zeros(nhalos)
    halos['vc'] = na.zeros(nhalos)
    halos['spin'] = na.zeros(nhalos)
    #added by me
    halos['mstar'] = na.zeros(nhalos)
    halos['halonum'] = na.zeros(nhalos)

    halonum = 0
    for l in lines:
        if not (l.startswith("#") or l.startswith("datavar") or len(l) < 2):
            values = l.split()
            halos['pos'][halonum,:] = map(float, values[0:3])
            halos['vel'][halonum,:] = map(float, values[9:12])
            halos['npart'][halonum] = int(values[4])
            halos['mass'][halonum] = float(values[5])
            halos['mvir'][halonum] = float(values[6])
            #added by me
	    #halos['mstar'][halonum] = float(values[7])
	    #halos['halonum'][halonum] = float(values[3])
	    #
	    halos['rvir'][halonum] = float(values[8])
            halos['vc'][halonum] = float(values[12])
            halos['spin'][halonum] = float(values[16])
	    halonum += 1
            if halonum >= nhalos:
                break
    del lines
    halos['NumberOfHalos'] = len(halos['mass'])
    return halos

def r200(pf, center, mass, verbose=False, delta=200.0):
    # Obtain cosmology parameters
    mh = 1.673e-24
    omegam = float(pf.get_parameter("CosmologyOmegaMatterNow"))
    omegal = float(pf.get_parameter("CosmologyOmegaLambdaNow"))
    zz = pf["CosmologyCurrentRedshift"]
    hh = pf["CosmologyHubbleConstantNow"]
    cosmo = Cosmology()
    ez = na.sqrt(omegal + omegam * (1+zz)**3)
    omegamz = omegam * (1+zz)**3 / ez**2
    hz = (hh * 100.0 / 3.086e19) * ez
    rhoc = cosmo.CriticalDensity(pf["CosmologyCurrentRedshift"])

    # Calculate virial radius from given mass, then find r200
    rvir = (6.673e-8 * mass * 1.989e33 /
            (100.0 * omegamz * hz**2))**(1.0/3)
    sp = pf.h.sphere(center, 3*rvir/pf["cm"])
    rmin = max(pf.h.get_smallest_dx() * pf["cm"], 0.1*rvir)
    prof = BinnedProfile1D(sp, 100, "Radius", rmin, 5*rvir)
    prof.add_fields(["TotalMassMsun", "CellVolume"], accumulation=True,
                    weight=None)
    avg_overden = 1.989e33*prof["TotalMassMsun"] / prof["CellVolume"] / rhoc
    avg_overden = avg_overden[-1:0:-1]
    rr = prof['Radius'][-1:0:-1]
    r200 = na.interp(delta, avg_overden, rr)
    del prof
    del sp
    if na.isnan(r200):
        return None, None, None
    
    # Find center of mass, using this r200.  Then find r200 again.
    sp = pf.h.sphere(center, r200/pf["cm"])
    com = sp.quantities["CenterOfMass"](use_particles=True)
    del sp
    sp = pf.h.sphere(com, 1.5*r200/pf["cm"])
    rmin = max(pf.h.get_smallest_dx() * pf["cm"], 0.1*r200)
    prof = BinnedProfile1D(sp, 100, "Radius", rmin, 1.5*r200)
    prof.add_fields(["TotalMassMsun", "CellVolume"], accumulation=True,
                    weight=None)
    avg_overden = 1.989e33*prof["TotalMassMsun"] / prof["CellVolume"] / rhoc
    avg_overden = avg_overden[-1:0:-1]
    rr = prof['Radius'][-1:0:-1]
    r200_2 = na.interp(delta, avg_overden, rr)
    M200 = delta * rhoc * (4.0/3) * na.pi * r200_2**3 / 1.989e33
    del prof
    del sp
    if na.isnan(r200_2):
        return None, None, None

    if verbose:
        print "given mass = %g Msun" % mass
        print "rvir(mass) = %g kpc" % (rvir / 3.086e21)
        print "r%03d     = %g kpc" % (int(delta), r200_2 / 3.086e21)
        print "M%03d     = %g Msun" % (int(delta), M200)
        print "original center = %s" % center
        print "center of mass  = %s" % com

    return r200_2, M200, com

def com_particle_iter(position, mass, factor=0.05):
    N = mass.size
    inside = na.ones(N, dtype="bool")
    com = na.zeros(3)
    rmax = None
    it = 0
    while N > 2:
        for dim in range(3):
            com[dim] = (position[inside, dim] * mass[inside]).sum() / mass[inside].sum()
        dr = position - com
        r = na.sqrt((dr**2).sum(1))
        if rmax == None: rmax = r.max()
        rmax *= 1.0 - factor
        inside = r < rmax
        N = na.count_nonzero(inside)
        it += 1
    return com
