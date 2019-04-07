"""
To set up a simplified system with randomly diffusion insulin secretory granules.

"""
from __future__ import print_function, division
import IMP.atom
import IMP.algebra
import IMP.rmf
import IMP.core
import RMF
import IMP.container
import IMP.display
import sys
import math
import random
import numpy
from random import choice
from random import uniform
import GranuleFactory

def convert_time_ns_to_frames(time_ns, step_size_fs):
    '''
    Given time in nanoseconds time_ns and step size in femtosecond
    step_size_fs, return an integer number of frames greater or equal
    to 1, such that time_ns*step_size_fs is as close as possible to
    time_ns.
    '''
    FS_PER_NS= 1E6
    time_fs= time_ns * FS_PER_NS
    n_frames_float= (time_fs+0.0) / step_size_fs
    n_frames= int(round(n_frames_float))
    return max(n_frames, 1)

def create_nucleus(m, R):
    '''
    Generate a coarse-grained spherical nuclear envelope
    of radius R in model m
    '''
    p= IMP.Particle(m, "nucleus")
    xyzr = IMP.core.XYZR.setup_particle(p)
    xyzr.set_coordinates_are_optimized(True)
    xyzr.set_coordinates([0,0,0])
    xyzr.set_radius(R)
    IMP.display.Colored.setup_particle(p,
                                       IMP.display.get_display_color(2))
    IMP.atom.Hierarchy.setup_particle(p)
    return p

#---------- Simulation parameters ----------

# I. Parts parameters
L = 50000 # Length of our bounding box, A
R = 20000 # PBC radius, A
R_NUCLEUS = 10000 # NE radius, A
N_GRANULES = 50 # Number of granules
R_GRANULES = 1500 # Radius of granules, A

# II. Interaction parameters
K_BB = 0.1  # Strength of the harmonic boundary box in kcal/mol/A^2
K_EXCLUDED=0.1 # Strength of lower-harmonic excluded volume score in kcal/mol/A^2

# III. Time parameters
BD_STEP_SIZE_SEC= 10E-8 # a time step of 10 ns
SIM_TIME_SEC= 0.050 # a total simulation time of 5 ms
bd_step_size_fs= BD_STEP_SIZE_SEC * 1E+15
sim_time_ns= SIM_TIME_SEC * 1E+9
RMF_DUMP_INTERVAL_NS= sim_time_ns / 1000.0
sim_time_frames= convert_time_ns_to_frames(sim_time_ns, bd_step_size_fs)
rmf_dump_interval_frames= convert_time_ns_to_frames(RMF_DUMP_INTERVAL_NS, bd_step_size_fs)
print("Simulation time {:.1e} ns / {} frames; "
      "RMF dump interval {:.1e} ns / {} frames".format(sim_time_ns,
                                                      sim_time_frames,
                                                       RMF_DUMP_INTERVAL_NS,
                                                      rmf_dump_interval_frames))

# -------- I. System parts --------
# Model:
m = IMP.Model()
# Root of parts hierarchy:
p_root= IMP.Particle(m, "root")
h_root = IMP.atom.Hierarchy.setup_particle(p_root)
# Outer bounding box for simulation:
bb = IMP.algebra.BoundingBox3D(IMP.algebra.Vector3D(-L/2, -L/2, -L/2), IMP.algebra.Vector3D(L/2, L/2, L/2))
# PBC cytoplasm bounding sphere:
pbc_sphere= IMP.algebra.Sphere3D([0,0,0], R)
# Nucleus:
p_nucleus= create_nucleus(m, R_NUCLEUS)
IMP.atom.Mass.setup_particle(p_nucleus, 1.0) # fake mass
h_nucleus= IMP.atom.Hierarchy(p_nucleus)
h_root.add_child(h_nucleus)
# Granules hierarchy root:
p_granules_root= IMP.Particle(m, "Granules")
IMP.atom.Mass.setup_particle(p_granules_root, 1.0) # fake mass
h_granules_root= IMP.atom.Hierarchy.setup_particle(p_granules_root)
h_root.add_child(h_granules_root)
# Actual granules:
nucleus_sphere= IMP.core.XYZR(p_nucleus).get_sphere()
gf=GranuleFactory.GranuleFactory\
    (model= m,
     default_R= R_GRANULES,
     cell_sphere= pbc_sphere,
     nucleus_sphere= nucleus_sphere)
for i in range(N_GRANULES):
    granule= gf.create_simple_granule("Granule_{}".format(i))
    h_granule= IMP.atom.Hierarchy(granule)
    h_granules_root.add_child(h_granule)
    IMP.atom.Diffusion(granule).set_diffusion_coefficient(1.5E-10) # A^2/fs
    #print(IMP.atom.Diffusion(granule).get_diffusion_coefficient())
    #print(IMP.core.XYZ(granule))
#print(h_root.get_children())

# ----- II. System interactions: -----
# Add enclosing spheres for pbc and outer simulation box
bb_harmonic= IMP.core.HarmonicUpperBound(0, K_BB)
pbc_bsss = IMP.core.BoundingSphere3DSingletonScore(bb_harmonic,
                                                   pbc_sphere)
outer_bbss = IMP.core.BoundingBox3DSingletonScore(bb_harmonic,
                                                  bb)
# Restraints - match score with particles:
rs = []
rs.append(IMP.container.SingletonsRestraint(pbc_bsss,
                                            h_granules_root.get_children()))
rs.append(IMP.container.SingletonsRestraint(outer_bbss,
                                            h_granules_root.get_children()))
# Add excluded volume restraints among all (close pairs of) particles:
ev = IMP.core.ExcludedVolumeRestraint(IMP.atom.get_leaves(h_root),
                                      K_EXCLUDED,
                                      10, # slack affects speed only
                                          # (slack of close pairs finder)
                                      "EV")
rs.append(ev)
# Scoring Function from restraints
sf = IMP.core.RestraintsScoringFunction(rs, "SF")
#print(h_root.get_children())

# -------- III. System dynamic: --------
bd = IMP.atom.BrownianDynamics(m)
bd.set_log_level(IMP.SILENT)
bd.set_scoring_function(sf)
bd.set_maximum_time_step(bd_step_size_fs) # in femtoseconds
bd.set_temperature(310.15)

# -------- Add RMF visualization --------
rmf = RMF.create_rmf_file("random.rmf")
rmf.set_description("Brownian dynamics trajectory with {}fs timestep.\n"\
                    .format(bd_step_size_fs))
IMP.rmf.add_hierarchy(rmf, h_root)
IMP.rmf.add_restraints(rmf, rs)
IMP.rmf.add_geometry(rmf, IMP.display.BoundingBoxGeometry(bb))
IMP.rmf.add_geometry(rmf, IMP.display.SphereGeometry(pbc_sphere))
    
# Pair RMF with model using an OptimizerState ("listener")
sos = IMP.rmf.SaveOptimizerState(m, rmf)
sos.set_log_level(IMP.SILENT)
sos.set_simulator(bd)
sos.set_period(rmf_dump_interval_frames)
bd.add_optimizer_state(sos)
# Dump initial frame to RMF
sos.update_always("initial conformation")

# -------- Run simulation ---------
print("Running simulation")
#m.update()
print("Score before: {:f}".format(sf.evaluate(True)))
bd.optimize(sim_time_frames)
print("Run finished succesfully")
print("Score ater: {:f}".format(sf.evaluate(True)))
