import MDAnalysis as mda
import networkx as nx
import MDAnalysis.analysis.distances
import scipy.stats as stats
import numpy as np
from itertools import product as product


class ipi:


    def __init__(self,u,frame,
                             core_sel,
                             density_sel,
                             interface,
                             no_bins,
                             cluster=True,
                             recombine=True,
                            interpolate_interface=False):
        self.u = u
        self.frame = frame
        self.core_sel = core_sel
        self.density_sel = density_sel       
        self.interface = interface
        self.no_bins = no_bins
        self.cluster = cluster
        self.recombine = recombine
        self.interpolate_interface = interpolate_interface
 
        self._intrinsic_density()
        
    def _intrinsic_density(self):

        '''
        u:   MDAnalysis Universe
        
        frame: MDAnalysis trajectory frame (int)
           
        core_sel: the atom selection that describes the reference object (membrane, slab, etc.)
        
        density_sel: the atom selection for which the intrinsic density is to be calculated
        
        interface: "Lower" or "Upper" interface of core_sel to calculate intrinsic density from 
        
        no_bins: number of bins in each lateral dimension (larger number gives finer grid) see 
        Ziolek et al Langmuir 2019 SI for discussion of grid size effects
        
        cluster: check whether any core_sel molecules have diffused out of the main structure 
        
        recombine: slow step that recombines the core_sel structure if it's split over the z PBC
    
        interpolate_interface: simple linear interpolation of the intrinsic interface, useful 
                              for patchy monolayers, etc. (like Ziolek et al Langmuir 2019 )
        '''
    
        self.u.trajectory[self.frame] 
    
        ### clustering check  
        
        core = self.u.select_atoms(self.core_sel)
        
        if self.cluster==True:
                
            cutoff_cluster = 25
    
            cluster_atoms = core.select_atoms('name '+str(core.atoms.names[0]))
    
            dist_array=MDAnalysis.analysis.distances.contact_matrix(cluster_atoms.atoms.positions,cutoff=cutoff_cluster,box=self.u.dimensions)
    
            G = nx.from_numpy_matrix(dist_array)
            clusters=[h for h in nx.connected_components(G)]
    
            for c in clusters:
                if len(c)==max([len(c_) for c_ in clusters]):
                    slab_resids = [i+1 for i in  c] #from np array to resids
    
            ###make sure the core selection doesn't break up and save positions as slab_calc_positions
            core = core.select_atoms('resid '+str(' '.join([str(i) for i in slab_resids])))
    
        ### make whole in z direction if slab/bilayer is split
    
        if self.recombine==True:
    
            # bonds_original=core.bonds
            new_bond_atoms = core.select_atoms('name '+str(core.atoms.names[0]))
            new_bonds=[((int(new_bond_atoms[i].index)),(int(new_bond_atoms[i+1].index))) for i in range(len(new_bond_atoms)-1)]
            self.u.add_bonds(new_bonds)
    
    
        #just take x,y as is and update z using this protocol
        #then just duplicate the positions of the solvent etc by taking away the box size in  z and go for certain distance only...
        #only valid until half way ...
        
    
    
            slab_calc_postions=np.zeros_like(core.positions)
            
            slab_calc_postions[:,0:2]=core.positions[:,0:2] #x and y positions can be from initial wrapped traj
    
            core.unwrap()
    
            slab_calc_postions[:,2]=core.positions[:,2]  #z positions can be from the newly unwrapped traj
        
        elif self.recombine==False:
            
            slab_calc_postions=np.zeros_like(core.positions)
            slab_calc_postions=core.positions
            
        if self.interface=='Lower':
            
            interface_grid = stats.binned_statistic_2d(slab_calc_postions[:,0]%self.u.dimensions[0],
                                                         slab_calc_postions[:,1]%self.u.dimensions[1],
                                                         slab_calc_postions[:,2],
                                                          bins=[np.linspace(0,self.u.dimensions[0],self.no_bins),np.linspace(0,self.u.dimensions[1],self.no_bins)],
                                                         statistic=np.min)
            interface_grid=interface_grid.statistic
        
        else:
            
            interface_grid = stats.binned_statistic_2d(slab_calc_postions[:,0]%self.u.dimensions[0],
                                                         slab_calc_postions[:,1]%self.u.dimensions[1],
                                                         slab_calc_postions[:,2],
                                                          bins=[np.linspace(0,self.u.dimensions[0],self.no_bins),np.linspace(0,self.u.dimensions[1],self.no_bins)],
                                                         statistic=np.max)
            interface_grid=interface_grid.statistic
            
        selection=self.u.select_atoms(self.density_sel)
        sel_grid_positions = stats.binned_statistic_2d(selection.positions[:,0]%self.u.dimensions[0],
                                                     selection.positions[:,1]%self.u.dimensions[1],
                                                     selection.positions[:,2],
                                                     bins=[np.linspace(0,self.u.dimensions[0],self.no_bins),np.linspace(0,self.u.dimensions[1],self.no_bins)],
                                                     statistic='count',
                                                     expand_binnumbers=True)
        
        pc_coverage=1
        if self.interpolate_interface==True:
            box_ = self.u.dimensions
            
            pc_coverage = np.count_nonzero(~np.isnan(interface_grid)) / (len(np.linspace(0,box_[0],num=self.no_bins))-1)*(len(np.linspace(0,box_[0],num=self.no_bins))-1)

            while np.count_nonzero(~np.isnan(interface_grid))!=(len(np.linspace(0,box_[0],num=self.no_bins))-1)*(len(np.linspace(0,box_[0],num=self.no_bins))-1):
                
                for i in range(len(interface_grid)):
                    for j in range(len(interface_grid)):
    
                        if np.isnan(interface_grid[i][j]):
                            n_i=[i-1,i,i+1]
                            n_j=[j-1,j,j+1]
                            interface_grid[i][j]=np.nanmean([interface_grid[divmod(ip,len(interface_grid))[1]][divmod(jp,len(interface_grid))[1]] for ip,jp in product(n_i,n_j)])
    
        
        intrinsic_z_=[]
        
        for atom in range(len(selection)):
            interface_position = interface_grid[sel_grid_positions.binnumber[0][atom]-1,sel_grid_positions.binnumber[1][atom]-1]
            intrinsic_z_.append(selection.positions[atom,2]-interface_position)
    
        half_z = self.u.trajectory[self.frame].dimensions[2]/2
        
        new_intrinsic_z_=[]
        
        for atom in intrinsic_z_:
            
            if atom > half_z:
                new_intrinsic_z_.append(atom-half_z*2)
            elif atom < -half_z:
    
                new_intrinsic_z_.append(atom+half_z*2)
    
            else:
                new_intrinsic_z_.append(atom)
 
        self.intrinsic_z_,self.new_intrinsic_z_ ,self.interface_grid,self.pc_coverage =   intrinsic_z_,new_intrinsic_z_,interface_grid,pc_coverage

 

    def __call__(self):
        """
        .
        """     
        return self.new_intrinsic_z_,self.interface_grid,self.pc_coverage








