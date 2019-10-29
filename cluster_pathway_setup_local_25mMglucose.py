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
import GlucoseFactory
#import IMP.npctransport

########TODO###########3
#1, dffusion coefficien of ISG

# record the simulation and ISG distance
f1=open('pathway_ISG_record_25mMglucose.txt', 'w')
f2=open('pathway_ISG-NE_distance_25mMglucose.xvg', 'w')

#def time_step_estimation(D_GLUCOSE):
    
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

def do_stats_after_optimize(granules):
    '''write out ISG coordinates long time'''
    distance = []
    for g in granules:
        distance.append(round(IMP.core.get_distance(IMP.core.XYZ(h_nucleus),IMP.core.XYZ(g)),4))
    return distance

def get_random_vector_in_cytoplasm(outer_sphere, inner_sphere):
    '''
    Return a random vector inside cell (outer) sphere and
    outside the nuclear envelope (=inner) sphere
    '''
    R_inner= inner_sphere.get_radius()
    while True:
        vector = IMP.algebra.get_random_vector_in(outer_sphere)
        if vector.get_magnitude() > R_inner:
            return vector

#---------- Simulation parameters ----------

# I. Parts parameters
L = 140000 # Length of our bounding box, A
R = 65351 # PBC radius, A
R_NUCLEUS = 42550 # NE radius, A
N_GRANULES = 306 # Number of granules
R_GRANULES = 1289 # Radius of granules, A
N_PATCHES = 6 # Number of patches on the surface of each ISG
D_GRANULES = 3.2E-5 # A^2/fs, 100,000 times of the real diffusion coeff.
N_GLUCOSE = 203704   #50926944 #
R_GLUCOSE = 100 # A, 20 times of the real radius
D_GLUCOSE = 9.6E-5 # A^2/fs
R_CLUSTER = 4000 # Radius of the cluster, A
R_CLUSTER_CUTOFF = 2000 # Radius of the cluster cutoff, A
N_GRANULES_IN_CLUSTER = 40 # number of granules in the cluster

# II. Interaction parameters
K_BB = 0.1  # Strength of the harmonic boundary box in kcal/mol/A^2
K_EXCLUDED=0.1 # Strength of lower-harmonic excluded volume score in kcal/mol/A^2
K_BB_CLUSTER = 0.5

# III. Time parameters
#BD_STEP_SIZE_SEC= 10E-8 # a time step of 10 ns
BD_STEP_SIZE_SEC= 1E-10 # a time step of 0.001 ns
#SIM_TIME_SEC= 0.050 # a total simulation time of 5 ms
SIM_TIME_SEC= 1E-9 # a total simulation time of 100 ns 1ns = 1E-9s
bd_step_size_fs= BD_STEP_SIZE_SEC * 1E+15
sim_time_ns= SIM_TIME_SEC * 1E+9
#RMF_DUMP_INTERVAL_NS= sim_time_ns / 1000.0
RMF_DUMP_INTERVAL_NS= sim_time_ns / 100.0
sim_time_frames= convert_time_ns_to_frames(sim_time_ns, bd_step_size_fs)
rmf_dump_interval_frames= convert_time_ns_to_frames(RMF_DUMP_INTERVAL_NS, bd_step_size_fs)
print("Simulation time {:.1e} ns / {} frames; "
      "RMF dump interval {:.1e} ns / {} frames".format(sim_time_ns,
                                                      sim_time_frames,
                                                       RMF_DUMP_INTERVAL_NS,
                                                      rmf_dump_interval_frames), file = f1)

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
grnf=GranuleFactory.GranuleFactory\
    (model= m,
     default_R= R_GRANULES,
     cell_sphere= pbc_sphere,
     nucleus_sphere= nucleus_sphere)
for i in range(N_GRANULES):
    granule= grnf.create_granule_with_interaction_patches("Granule_{}".format(i), N_PATCHES)
    h_granule= IMP.atom.Hierarchy(granule)
    #print(h_granule.get_children())
    h_granules_root.add_child(h_granule)
    IMP.atom.Diffusion(granule).set_diffusion_coefficient(D_GRANULES)
    #print(IMP.atom.Diffusion(granule).get_diffusion_coefficient())
    #print(IMP.core.XYZ(granule))
#print(h_root.get_children())

# Glucose hierarchy root:
p_glucose_root= IMP.Particle(m, "Glucose")
IMP.atom.Mass.setup_particle(p_glucose_root, 1.0) # fake mass
h_glucose_root= IMP.atom.Hierarchy.setup_particle(p_glucose_root)
h_root.add_child(h_glucose_root)
# Actual glucose:
gluf=GlucoseFactory.GlucoseFactory\
    (model= m,
     default_R= R_GLUCOSE,
     cell_sphere= pbc_sphere,
     nucleus_sphere= nucleus_sphere)

for i in range(N_GLUCOSE):
    glucose= gluf.create_simple_glucose("Glucose_{}".format(i))
    h_glucose= IMP.atom.Hierarchy(glucose)
    h_glucose_root.add_child(h_glucose)
    IMP.atom.Diffusion(glucose).set_diffusion_coefficient(D_GLUCOSE)
    #print(IMP.atom.Diffusion(glucose).get_diffusion_coefficient())
    #print(IMP.core.XYZ(glucose))

# Cluster hierarchy root:
p_cluster_root= IMP.Particle(m, "Granule_in_cluster")
IMP.atom.Mass.setup_particle(p_cluster_root, 1.0) # fake mass
h_cluster_root= IMP.atom.Hierarchy.setup_particle(p_cluster_root)
h_root.add_child(h_cluster_root)
# Actual granules in cluster:
cluster_sphere = IMP.algebra.Sphere3D(get_random_vector_in_cytoplasm(
    IMP.algebra.Sphere3D([0,0,0], R-2*(R_CLUSTER+R_CLUSTER_CUTOFF)),
    IMP.algebra.Sphere3D([0,0,0], R_NUCLEUS+R_CLUSTER_CUTOFF)), R_CLUSTER)
cluster_inner_sphere = IMP.algebra.Sphere3D(cluster_sphere.get_center(),1E-9)
gfc = GranuleFactory.GranuleFactory(
    model = m,
    default_R = R_GRANULES,
    cell_sphere = cluster_sphere,
    nucleus_sphere = cluster_inner_sphere)
for i in range(N_GRANULES_IN_CLUSTER):
    granule = gfc.create_granule_with_interaction_patches("Granule_in_cluster_{}".format(i), N_PATCHES)
    h_granule= IMP.atom.Hierarchy(granule)
    h_cluster_root.add_child(h_granule)
    IMP.atom.Diffusion(granule).set_diffusion_coefficient(D_GRANULES) # A^2/fs

# ----- II. System interactions: -----
#print(IMP.core.RigidBody(h_granules_root.get_children()[1]).get_rigid_members()[1:])

# Add granule activition
#gaos= IMP.npctransport.GranuleActivationOptimizerState(h_granules_root.get_children(), h_glucose_root.get_children(), 500, 100, 1) #granule and patches, glucose, contact_range, slack, frame interval
#gaos.set_log_level(IMP.TERSE)
#print("Model", gaos.get_model())

# Add enclosing spheres for pbc and outer simulation box
bb_harmonic= IMP.core.HarmonicUpperBound(0, K_BB)
bb_harmonic_cluster = IMP.core.HarmonicUpperBound(0,K_BB_CLUSTER)
pbc_bsss = IMP.core.BoundingSphere3DSingletonScore(bb_harmonic,
                                                   pbc_sphere)
outer_bbss = IMP.core.BoundingBox3DSingletonScore(bb_harmonic,
                                                  bb)
cluster_bbss = IMP.core.BoundingSphere3DSingletonScore(bb_harmonic_cluster,
                                                    cluster_sphere)
# Restraints - match score with particles:
rs = []
rs.append(IMP.container.SingletonsRestraint(pbc_bsss,
                                            h_granules_root.get_children()))
rs.append(IMP.container.SingletonsRestraint(outer_bbss,
                                            h_granules_root.get_children()))
rs.append(IMP.container.SingletonsRestraint(pbc_bsss,
                                            h_glucose_root.get_children()))
rs.append(IMP.container.SingletonsRestraint(outer_bbss,
                                            h_glucose_root.get_children()))
rs.append(IMP.container.SingletonsRestraint(cluster_bbss,
                                            h_cluster_root.get_children()))                                          
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
rmf = RMF.create_rmf_file("pathway_25mMglucose.rmf")
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
#bd.add_optimizer_state(gaos)

# Dump initial frame to RMF
sos.update_always("initial conformation")

# -------- Run simulation ---------
print("Running simulation", file = f1)
print("Score before: {:f}".format(sf.evaluate(True)), file = f1)
n_frames_left=sim_time_frames

#frames_per_cycle=1000 # can be 10, or 100, dependeing on the frequency of the optimizer state and how much ISG moves.
frames_per_cycle=1000 # can be 10, or 100, dependeing on the frequency of the optimizer state and how much ISG moves.
while n_frames_left>0:
    cur_n_frames=min(frames_per_cycle, n_frames_left)
    bd.optimize(cur_n_frames)
    print(*do_stats_after_optimize(h_granules_root.get_children()), sep=" ", file = f2)
    print(sim_time_frames - n_frames_left, sep=" ", file = f1)
    #print(n_frames_left)
    #print(cur_n_frames)
    n_frames_left = n_frames_left - cur_n_frames

print("Run finished succesfully", file = f1)
print("Score ater: {:f}".format(sf.evaluate(True)), file = f1)
