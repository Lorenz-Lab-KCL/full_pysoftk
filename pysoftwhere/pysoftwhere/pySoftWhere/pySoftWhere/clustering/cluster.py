import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances
import scipy.stats as stats
import numpy as np

class find_largest_cluster:
    
    def __init__(self,u,
                      frame,
                      selection,
                      cutoff_distance=10, 
                      define_clustering_atoms=False): #name of single atom
    
        self.u = u
        self.frame = frame
        self.selection = selection
        self.cutoff_distance = cutoff_distance
        self.cluster_atoms = define_clustering_atoms
        
       
        self._cluster_gt()
    
    def _cluster_gt(self):


        """
        do the clustering
        """

        ts =  self.u.trajectory[self.frame]
    
        polymers=self.u.select_atoms(self.selection)
    
        if  self.cluster_atoms==False:
            self.cluster_atoms='name '+(polymers.atoms.names[0]) 
        else:
            pass    
    
        cluster_atoms_sel = polymers.select_atoms(self.cluster_atoms)

        dist_array=MDAnalysis.analysis.distances.contact_matrix(cluster_atoms_sel.atoms.positions,
                                                            cutoff=self.cutoff_distance,
                                                            box=self.u.dimensions)

        G = nx.from_numpy_matrix(dist_array)
        clusters=[h for h in nx.connected_components(G)]
    
        for c in clusters:
            if len(c)==max([len(c_) for c_ in clusters]):
                largest_cluster_resids = [i+1 for i in  c]
            
        self.largest_cluster_resids = largest_cluster_resids
    
    def __call__(self):
        """
        .
        """     
        return self.largest_cluster_resids

class find_all_clusters:
    
    def __init__(self,u,
                      frame,
                      selection,
                      cutoff_distance=10, 
                      define_clustering_atoms=False): #name of single atom
    
        self.u = u
        self.frame = frame
        self.selection = selection
        self.cutoff_distance = cutoff_distance
        self.cluster_atoms = define_clustering_atoms
        
       
        self._cluster_gt()
    
    def _cluster_gt(self):


        """
        do the clustering
        """

        ts =  self.u.trajectory[self.frame]
    
        polymers=self.u.select_atoms(self.selection)
    
        if  self.cluster_atoms==False:
            self.cluster_atoms='name '+(polymers.atoms.names[0]) 
        else:
            pass    
    
        cluster_atoms_sel = polymers.select_atoms(self.cluster_atoms)

        dist_array=MDAnalysis.analysis.distances.contact_matrix(cluster_atoms_sel.atoms.positions,
                                                            cutoff=self.cutoff_distance,
                                                            box=self.u.dimensions)

        G = nx.from_numpy_matrix(dist_array)
        clusters=[list(h) for h in nx.connected_components(G)]
    
        self.clusters = clusters
    
    def __call__(self):
        """
        .
        """     
        return self.clusters


class make_cluster_whole:
    
    def __init__(self,u,
                       frame,
                       cluster_resids,
                       core_selection,  #need to add  -> easier to do it now
                       shell_selection, #need to add- > easier to do it now 
                       solvent=False): 
    
        self.u = u
        self.frame = frame
        self.cluster_resids = cluster_resids
        self.core_selection = core_selection
        self.shell_selection = shell_selection
        self.solvent = solvent 
       
        self._make_whole()
        
    def _make_whole(self):
        
        self.u.trajectory[self.frame]
        
        cluster_sel = self.u.select_atoms('resid '+" ".join([str(i) for i in self.cluster_resids]))
        core_sel = self.u.select_atoms(self.core_selection).select_atoms('resid '+" ".join([str(i) for i in self.cluster_resids]))
        
        if self.solvent==False:
            shell_sel = self.u.select_atoms(self.shell_selection).select_atoms('resid '+" ".join([str(i) for i in self.cluster_resids]))
        elif self.solvent==True:
            shell_sel = self.u.select_atoms(self.shell_selection)

        cluster_atoms_positions=cluster_sel.positions.copy()
        core_sel_atoms_positions=core_sel.positions.copy()
        shell_sel_atoms_positions=shell_sel.positions.copy()
        
    
    
            
        box=self.u.dimensions
    
        for dimension in range(0,3):
    
            y = stats.binned_statistic(cluster_sel.positions[:,dimension],
                                         cluster_sel.positions[:,dimension],
                                         bins=np.arange(-50,box[0]+50,10),
                                         statistic='count').statistic
    
            a = stats.binned_statistic(cluster_sel.positions[:,dimension],
                                       cluster_sel.positions[:,dimension],
                                       bins=np.arange(-50,box[0]+50,10),
                                       statistic='count').bin_edges
    
            move_above=False
            move_water=False
    
            filled=np.where(y!=0)[0]
    
            for i in range(len(filled)-1):
                if filled[i]!=filled[i+1]-1:
                    move_above = a[filled[i+1]]
                    move_lower_bound=a[filled[i]]
                    move_water=(move_lower_bound+move_above)/2
            ####
    
            if move_above!=False:
                for atom in range(len(cluster_atoms_positions)):
                    if cluster_atoms_positions[atom][dimension]>move_above:
                        cluster_atoms_positions[atom][dimension]=cluster_atoms_positions[atom][dimension]-box[dimension]
    
                for atom in range(len(core_sel_atoms_positions)):
                    if core_sel_atoms_positions[atom][dimension]>move_above:
                        core_sel_atoms_positions[atom][dimension]=core_sel_atoms_positions[atom][dimension]-box[dimension]
                
                if self.solvent==False:
    
                    for atom in range(len(shell_sel_atoms_positions)):
                        if shell_sel_atoms_positions[atom][dimension]>move_above:
                            shell_sel_atoms_positions[atom][dimension]=shell_sel_atoms_positions[atom][dimension]-box[dimension]
    
    
                if self.solvent==True:
                    for atom in range(len(shell_sel_atoms_positions)):
                        if shell_sel_atoms_positions[atom][dimension]>move_water:
                            shell_sel_atoms_positions[atom][dimension]=shell_sel_atoms_positions[atom][dimension]-box[dimension]
                
        self.cluster_atoms_positions,self.core_sel_atoms_positions,self.shell_sel_atoms_positions = cluster_atoms_positions,core_sel_atoms_positions,shell_sel_atoms_positions
     



    def __call__(self):
        """
        .
        """     

        return self.cluster_atoms_positions,self.core_sel_atoms_positions,self.shell_sel_atoms_positions
