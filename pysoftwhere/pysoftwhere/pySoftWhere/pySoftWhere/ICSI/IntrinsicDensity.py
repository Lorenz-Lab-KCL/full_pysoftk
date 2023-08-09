import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances
import scipy.stats as stats
import numpy as np
from itertools import product as product
import .icsi as isci
import ../clustering.cluster as psw_cluster

class IntrinsicDensity:


    def __init__(self,
                    u,
                    cluster_resids,
                    cluster_atoms_positions,
                    core_sel_atoms_positions,
                    shell_sel_atoms_positions,
                    frame=-1,
                    no_bins=31,
                    no_random_points=10000,
                    normalisation_run=False):

        self.u = u
        self.frame = frame
        self.no_bins = no_bins
        self.cluster_resids = cluster_resids
        self.cluster_atoms_positions=cluster_atoms_positions
        self.core_sel_atoms_positions=core_sel_atoms_positions
        self.shell_sel_atoms_positions=shell_sel_atoms_positions
        self.no_random_points=no_random_points
        self.normalisation_run=normalisation_run
        
        self._calc_id()
        
    def _calc_id(self):

        self.intrinsic_r,self.spherical_r, self.interface_vals = intrinsic_r, spherical_r, interface_vals
       


 

    def __call__(self):
        """
        .
        """     

        return self.intrinsic_r,self.spherical_r, self.interface_vals


