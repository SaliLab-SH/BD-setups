"""
Factory for generating granules
"""

from __future__ import print_function, division
import IMP.atom
import IMP.algebra
import IMP.core
import IMP.display
import IMP.insulinsecretion

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

def get_random_vector_in_cytoplasm_except_cluster(outer_sphere, inner_sphere, cluster_sphere):
    '''
    Return a random vector inside cell (outer) sphere and
    outside the nuclear envelope (=inner) sphere
    '''
    R_inner= inner_sphere.get_radius()
    while True:
        vector = IMP.algebra.get_random_vector_in(outer_sphere)
        if vector.get_magnitude() > R_inner and all(map(lambda x:(vector - x.get_center()).get_magnitude() > x.get_radius(),cluster_sphere)):
            return vector

class GranuleFactory:
    '''
    A class for generating insulin granules
    '''
    def __init__(self, model, default_R, cell_sphere, nucleus_sphere, cluster_sphere ):
        self.model= model
        self.default_R= default_R
        self.cell_sphere= cell_sphere
        self.nucleus_sphere= nucleus_sphere
        self.cluster_sphere= cluster_sphere

    def _create_granule_core(self, name):
        '''
        Create a granule in model m with name "Granule_i"
        and radius R at a random position in the cytoplasm
        (within cell_sphere and outside ne_sphere)
        '''
        p= IMP.Particle(self.model, name)
        xyzr= IMP.core.XYZR.setup_particle(p)
        xyzr.set_coordinates_are_optimized(True)
        v= get_random_vector_in_cytoplasm_except_cluster(self.cell_sphere,
                                                         self.nucleus_sphere,
                                                         self.cluster_sphere)
        xyzr.set_coordinates(v)
        xyzr.set_radius(self.default_R)
        # Setup (fake) mass and hierearchy
        IMP.atom.Mass.setup_particle(p, 1)
        # Setup the decorator for counting the number of secretion events
        IMP.insulinsecretion.SecretionCounter.setup_particle(p, 0)
        IMP.insulinsecretion.Maturation.setup_particle(p, 0)
        # set the coordinate values
        IMP.atom.Hierarchy.setup_particle(p) 
        IMP.display.Colored.setup_particle(p,
                                           IMP.display.get_display_color(0))
        return p

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

    def create_simple_granule(self, name):
        '''
        Create a simple granule with specified name - just a diffusive
        sphere decorated with mass, color and hierarchy
        '''
        p= self._create_granule_core(name)
        IMP.atom.Diffusion.setup_particle(p)
        return p

    def create_granule_with_interaction_patches(self,
                                                name,
                                                n_patches):
        '''
        Create a granule with n_patches interaction patches
        and specified name
        '''
        p= IMP.Particle(self.model,
                        name)
        # Setup (fake) mass and hierearchy
        IMP.atom.Mass.setup_particle(p, 1.0)
        IMP.insulinsecretion.SecretionCounter.setup_particle(p, 0)
        IMP.insulinsecretion.Maturation.setup_particle(p, 0)
        h= IMP.atom.Hierarchy.setup_particle(p)
        # Create core particle and add as rigid body core
        p_granule_core= self._create_granule_core(name + "_core")
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
        return p
