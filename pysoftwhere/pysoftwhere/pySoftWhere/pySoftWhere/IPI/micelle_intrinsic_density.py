
import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances
import cluster


class spherical_intrinsic_density:
    """Returns the intrinsic density of given species in a quasispherical nanostructure.  """ 

    
    def __init__(self, u, cluster_atoms, cluster_cutoff, micelle_selection, core_selection, measurement_selection, frame):

        self.u = u
        self.cluster_atoms = cluster_atoms
        self.cluster_cutoff = cluster_cutoff
        self.micelle_selection = micelle_selection
        self.core_selection = core_selection
        self_measurement_selection = measurement_selection 
        self.frame = frame
       
        self._intrinsic_density()

        
    def _intrinsic_density(self):
        
        
        """
        calculate the two profiles
        """
       
        ##will first call cluster.largest_cluster to find the ids for the micelle...      


  
        ts=self.u.trajectory[self.frame]
        
        cluster_resids = cluster.largest_cluster(self.u,self.cluster_atoms,self.cluster_cutoff,self.frame)()

        self.cluster_resids = cluster_resids


        ###now it needs to to do the intrinsic density calculation 



        return self.cluster_resids













    def __call__(self):
        
        """Return per frame: histogram of profile, normalisation histogram, core size, cluster size
        """
        return self.cluster_resids
