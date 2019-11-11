"""
Factory for generating granules in the cluster
"""

from __future__ import print_function, division
import IMP.atom
import IMP.algebra
import IMP.core
import IMP.display
import numpy as np
import IMP.insulinsecretion

def generate_given_vectors(default_R,center):
    '''
    Return a list of all the given coords for the ISGs in the cluster
    '''
    xs,ys,zs = [],[],[]
    f = open('/Users/tracy/Lab/Salilab/Code/IMP_code/granule_cluster/outcome/cluster_ISG_coords_070510','r')
    coords = f.readlines()
    for item in coords:
        xs.append(float(item.split(',')[0].split('(')[1]))
        ys.append(float(item.split(',')[1]))
        zs.append(float(item.split(',')[2].split(')')[0]))
    vectors = []
    for i in range(len(xs)):
        coord = (xs[i],ys[i],zs[i])
        vector = np.array(coord) - [0,0,0]
        distance = np.linalg.norm(vector)/14267 * default_R
        vector_new = center + vector/np.linalg.norm(vector) * distance
        vectors.append(vector_new)
    return vectors

class ClusterGranuleFactory:
    '''
    A class for generating insulin granules in the cluster
    '''
    def __init__(self, model, default_R, cluster_sphere, number):
        self.model= model
        self.default_R= default_R
        self.cluster_sphere= cluster_sphere
        self.number= number

    def _create_granule_cores(self, name):
        '''
        Create granules in model m with name "Granule_in_cluster_i" (todo)
        and radius R at given positions in the cluster
        '''
        ps = []
        center = self.cluster_sphere.get_center()
        vectors = generate_given_vectors(self.default_R,center)
        for coord in vectors:
            p= IMP.Particle(self.model, name)
            xyzr= IMP.core.XYZR.setup_particle(p)
            xyzr.set_coordinates_are_optimized(True)
            xyzr.set_coordinates(coord)
            xyzr.set_radius(self.default_R)
            # Setup (fake) mass and hierearchy
            IMP.atom.Mass.setup_particle(p, 1.0)   
            IMP.atom.Hierarchy.setup_particle(p) 
            IMP.display.Colored.setup_particle(p,
                                                IMP.display.get_display_color(0))
            ps.append(p)
        return ps

    def _create_granule_interaction_patch(self,
                                          coords,
                                          name,
                                          R= 100.0):
        '''
        Creates a granule interaction_patch at a random
        surface location on p_granule_core
        '''
        p=IMP.Particle(self.model, name)
        xyzr= IMP.core.XYZR.setup_particle(p)
        xyzr.set_coordinates(coords)
        xyzr.set_coordinates_are_optimized(True)
        xyzr.set_radius(R)
        # Setup (fake) mass and hierearchy
        IMP.atom.Mass.setup_particle(p, 1)   
        IMP.atom.Hierarchy.setup_particle(p)
        IMP.display.Colored.setup_particle(p,
                                           IMP.display.get_display_color(2))
        return p

    def create_simple_granules(self, name):
        '''
        Create simple granules with specified name - just a diffusive
        sphere decorated with mass, color and hierarchy
        '''
        ps= self._create_granule_cores(name)
        for p in ps:
            IMP.atom.Diffusion.setup_particle(p)
        return ps

    def create_granule_with_interaction_patches(self,
                                                name,
                                                n_patches):
        '''
        Create a granule with n_patches interaction patches
        and specified name
        '''
        ps = []
        # Create core particle and add as rigid body core
        p_granule_cores= self._create_granule_cores(name + "_core")
        for p_granule_core in p_granule_cores:
            p= IMP.Particle(self.model,
                        name)
            # Setup (fake) mass and hierearchy
            IMP.atom.Mass.setup_particle(p, 1.0)
            IMP.insulinsecretion.SecretionCounter.setup_particle(p, 0)
            IMP.insulinsecretion.Maturation.setup_particle(p, 0)
            h= IMP.atom.Hierarchy.setup_particle(p)
            rb= IMP.core.RigidBody.setup_particle(p,
                                              [p_granule_core])
            h.add_child(IMP.atom.Hierarchy(p_granule_core))
            rb.set_coordinates_are_optimized(True)
            # compute patch coordinates on core surface
            sphere_granule_core= IMP.core.XYZR(p_granule_core).get_sphere()
            patch_coords_list= \
                IMP.algebra.get_uniform_surface_cover(sphere_granule_core,
                                                        n_patches)
            for i, patch_coords in enumerate(patch_coords_list):
                patch_name= "{}_patch{}".format(name, i)
                p_patch= self._create_granule_interaction_patch(patch_coords,
                                                                patch_name)
                rb.add_member(p_patch)
                h.add_child(IMP.atom.Hierarchy(p_patch))

            # Set default diffusion coefficient according to the
            # Stokes radius of the core particle:
            rbd= IMP.atom.RigidBodyDiffusion.setup_particle(p)
            core_R= IMP.core.XYZR(p_granule_core).get_radius()
            D= IMP.atom.get_einstein_diffusion_coefficient(core_R)
            rbd.set_diffusion_coefficient(D)
            Drot= IMP.atom.get_einstein_rotational_diffusion_coefficient(core_R)
            rbd.set_rotational_diffusion_coefficient(Drot)
            ps.append(p)
        return ps
