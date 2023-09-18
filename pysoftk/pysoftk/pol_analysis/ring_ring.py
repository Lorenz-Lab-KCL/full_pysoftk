import MDAnalysis as mda
from MDAnalysis.lib import distances 

import numpy as np
import matplotlib.pyplot as plt
from tools.utils_mda import MDA_input

# suppress some MDAnalysis warnings when writing PDB files
import warnings
warnings.filterwarnings('ignore')

from rdkit import Chem
from rdkit.Chem import AllChem

from itertools import combinations
import itertools

from scipy.spatial.distance import cdist
from pysoftk.pol_analysis.micelle_whole import micelle_whole

from tqdm import tqdm
import time
from functools import wraps


def timeit(func):
    @wraps(func)
    def timeit_wrapper(*args, **kwargs):
        import time
        start_time = time.perf_counter()
        result = func(*args, **kwargs)
        end_time = time.perf_counter()
        total_time = end_time - start_time
        print(f'Function {func.__name__} Took {total_time:.4f} seconds')
        return result
    return timeit_wrapper


class RSA(MDA_input):
    """A class used to compute a Ring Stacking Analysis (RSA).
    """
    
    def __init__(self, tpr_file, xtc_file):
      """Instantiating the MDA_input class as
         super function.
      """

      super().__init__(tpr_file, xtc_file)

      
    def calc_angle(self, u_norm, u_norm2):
        """ Function to compute the dot product between two normal 
            vectors. 

        Parameters
        -----------

        u_norm: numpy.ndarray
           An user-provided normal vector. 

        u_norm2: numpy.ndarray
           An user-provided normal vector. 

        Return
        -------

        None:
           Result of a dot product between two vectors.

        """

        return abs(np.rad2deg(np.arccos(np.dot(u_norm,u_norm2)/np.linalg.norm(u_norm)/np.linalg.norm(u_norm2))))

    def GetRingSystems(self, mol, min_1, max_1, min_2, max_2, includeSpiro=False):
        """Function to detect the atomic indexes for a provided 
           molecular complex (single, dimer, timer and so on molecules).
 

           Parameters
           -----------

           mol: MDAnalysis.Universe
              User provided mda universe

           includeSpiro: class.bool
              To check for spiro sites


           Returns
           --------
 
           ring_contact_indices: class.list
              A list of arrays with the mda atom indices of the 
              rings of the polymers.

        """

        universe = mol.select_atoms('index '+str(int(min_1))+':'+str(int(max_1))+' or index '+str(int(min_2))+':'+str(int(max_2)))
        rdkit_mol = universe.convert_to("RDKIT")

        ri = rdkit_mol.GetRingInfo()
        systems = []
        for ring in ri.AtomRings():
            ringAts = set(ring)
            nSystems = []
            for system in systems:
                nInCommon = len(ringAts.intersection(system))
                if nInCommon and (includeSpiro or nInCommon>1):
                    ringAts = ringAts.union(system)
                else:
                    nSystems.append(system)
                    
            nSystems.append(ringAts)
            systems = nSystems 

        #getting the atom indices in mdanalysis
        ring_contact_indices = [universe.atoms[list(item)].indices for item in systems]
   
        return ring_contact_indices 


    def separate_rings(self, ring_contact_indices, pol2_min):
        """Function to get the possible ring combinations between two polymers to 
           calculate the distances between them.
 

           Parameters
           -----------

           mol: MDAnalysis.Universe
              User provided mda universe

           includeSpiro: class.bool
              To check for spiro sites


           Returns
           --------
 
           ring_contact_indices: class.list
              A list of arrays with the mda atom indices of the 
              rings of the polymers.

        """

        pol1_indices = [item for item in ring_contact_indices if np.all(item < pol2_min)]
        pol2_indices = [item for item in ring_contact_indices if np.any(item >= pol2_min)]

        ring_comb = list(itertools.product(pol1_indices, pol2_indices))

        return ring_comb

    def get_polymer_resids(self, pol_indices):
        """Function to get the possible ring combinations between two polymers to 
           calculate the distances between them.
 

           Parameters
           -----------

           pol_indices: np.array
              Numpy array with the atom indices of the dimer.

           Returns
           --------
 
           ring_contact_indices: class.list
             A list of arrays with the mda atom indices of the 
             rings of the polymers.

        """

        ring_comb = list(itertools.product(pol_indices[0], pol_indices[1]))
 
        return ring_comb 


    def svd(self, coordinates):
        """Function to obtain the Single Value Decomposition (SVD) 
           of coordinates.

        Parameters
        -----------
    
        array: numpy.ndarray
           Array containing the coordinates upon which svd is going 
           to be calculated.
    
        Return
        -------
    
        None:
            A numpy.ndarray with the right unitary singular vector, orthogonal to 
            the given set of coordinates.

        """

        u_, s, vh = np.linalg.svd(coordinates - coordinates.sum(axis=0) / coordinates.shape[0])
   
        return vh[2, :] 

    
    def rings_stacking_dist(self, u, ring_comb, polymer_pos_1,
                            polymer_pos_2, pol_indices, cut_off_stack=4, par_scheme='OpenMP'):
        """Function to find if two rings are stacked. This is done by checking if 
           the distance of geometry between two rings is smaller than an 
           user-defined cut off distance.
 

            Parameters
            -----------

            ring_comb: numpy.ndarray
               Numpy array with all possible ring combinations

            polymer_pos_1: numpy.ndarray
               Numpy array with the correct positions of all atoms 
               of the polymer.

            polymer_pos_2: numpy.ndarray
                Numpy array with the correct positions of all atoms 
                of the polymer.

            pol_indices: numpy.ndarray
                Numpy array with the first atom index of each polymer. The 
                first one needs to be the one corresponding to polymer_pos_1
                and the second one to polymer_pos_2

            cut_off_stack : class.float
                Cut off distance for ring stackings


           Return
           --------
 
           vector_angle: class.float or class.boolean
                Angle between the two planes of atoms


        """

        mol1_pos = polymer_pos_1[ring_comb[0]-pol_indices[0]]
        mol2_pos = polymer_pos_2[ring_comb[1]-pol_indices[1]]
        tot_dist = mda.lib.distances.distance_array(mol1_pos,
                                                    mol2_pos,
                                                    backend=str(par_scheme),
                                                    box=u.dimensions)

        if np.min(tot_dist) < cut_off_stack:      
            u_norm = self.svd(mol1_pos)
            u_norm2 = self.svd(mol2_pos)

            #finding the angle between the two vectors
            vector_angle = abs(np.rad2deg(np.arccos(np.dot(u_norm,u_norm2)/np.linalg.norm(u_norm)/np.linalg.norm(u_norm2))))

        else:
            vector_angle = None

            
        return vector_angle

    def run_rings_stacking_dist(self, u, ring_comb, cut_off, ang_cutoff):
        """Function to find if two rings are stacked. This is done by checking 
           if the distance of geometry between two rings is smaller than a 
           user-defined cut off distance.
 

           Parameters
           -----------

           ring_comb: numpy.ndarray
              Numpy array with all possible ring combinations

           cutoff : class.float
              Cut off distance for ring stackings

           ang_cut_off : class.float
              Cut off angle for ring stackings


        Returns
        --------
 
           ring_contact_atoms: class.list
              A list of arrays with the mda atom indices of the rings of the polymers

        """
    
        resid_1=np.unique(u.atoms[ring_comb[0][0]].resids)
        resid_2=np.unique(u.atoms[ring_comb[0][1]].resids)
        len_mol_1 = len(u.select_atoms('resid '+str(resid_1[0])))

        pol_resids = [list(resid_1), list(resid_2)]
        pol_resids_f = [item for sublist in pol_resids for item in sublist]

        #number to substract to the atom indices to only do the atom positions once
        pol1_first_index = u.select_atoms('resid '+str(pol_resids_f[0]))[0].index
        pol2_first_index = u.select_atoms('resid '+str(pol_resids_f[1]))[0].index
        pol_indices_f = [pol1_first_index, pol2_first_index]

        atom_pos_list = u.select_atoms('resname UNK and resid '+str(pol_resids_f[0])+'  '+str(pol_resids_f[1])).positions
        atom_positions_tot=atom_pos_list

        polymer_1_pos = atom_positions_tot[:len_mol_1]
        polymer_2_pos = atom_positions_tot[len_mol_1:]

        angles_tot_list= list(map(lambda i: self.rings_stacking_dist(u, ring_comb[i], polymer_1_pos,
                                                                polymer_2_pos, pol_indices_f), range(len(ring_comb))))

        ring_index_cutoff = [i for i, x in enumerate(angles_tot_list)
                             if x is not None and (x < ang_cutoff or x > 180-ang_cutoff)]
        
        ring_contact_atoms = [ring_comb[i] for i in ring_index_cutoff]

        #Printing the found dimers.
        if len(ring_contact_atoms)>0:
            snapshot = u.select_atoms('resid '+str(pol_resids_f[0])+' '+str(pol_resids_f[1]))
            snapshot.positions = atom_positions_tot
            
            with mda.Writer('snapshot_'+str(pol_resids_f[0])+'_'+str(pol_resids_f[1])+'.pdb', snapshot.n_atoms) as W:              
                    W.write(snapshot)

            pol_idx = (pol_resids_f[0], pol_resids_f[1])

            return (ring_contact_atoms, pol_idx)


    def find_several_rings_stacked(self, rings_df_name):
        """Function to find several rings stacked using a provided
           pandas DataFrame from previous calculations.

        Parameters
        -----------

        rings_df_name: pandas.DataFrame
           Pandas dataframe containing the results from the ring_ring 
           analysis.


        Results
        --------

        connected_components: class.list
           List of connected polymers along the trajectory.

        """
        import networkx as nx
        import pandas as pd

        rings_df = pd.read_parquet(rings_df_name)
        connections = rings_df.iloc[:, 1].tolist()

        graph = nx.Graph()

        # Add edges to the graph based on the connections
        for connection in connections:
            elem1, elem2 = connection
            graph.add_edge(elem1, elem2)

        # Find the connected components in the graph
        connected_components = list(nx.connected_components(graph))

        return connected_components


    def nearest_neighbours(self, u, lab1, lab2):
        """Function to filter different molecules
           by closest interatomic distance.

        Parameters
        -----------

        u: MDAnalysis.object
           Universe containing the trajectory

        lab1: str
           User provided label of one molecule

        lab2: str
           User provided label for second molecule

        frame : int
            Frame for calculation

        Returns
        --------

        float: np.float
          Minimum distance between the atoms belonging to the polymer pair 
          and all the atompositions.

        """

        mol=u.select_atoms(str(lab1)) 
        mol1=u.select_atoms(str(lab2))

        tot_dist = distances.distance_array(mol.positions, mol1.positions, box=u.dimensions)

        return np.min(tot_dist) 

    
    def number_pol(self, u):
        """Number of polymers inside an user-provided box.

        Parameters
        -----------

        u: MDAnalysis.Universe
           An user-provided MDAnalysis universe.

    
        Return
        -------

        None:
          Number of polymers within a box.

        """
    
        return len(u.atoms.residues) 

    
    def pol_cutoff(self, u, cutoff):
        """Function to count the number of polymers within an 
           user-provided cutoff.

        Parameters
        -----------

        u: MDAnalysis.universe
           An user-provided universe of MDAnalysis
 
        cutoff: float
           An user provided cutoff defining closest distances.

        frames : list
           User provided frames where calculation will be carried out

        Return
        -------

        index_cutoff: list
           List of intergers with the indexes of the selected
           molecules.
     
        """

        a=[i for i in range(np.min(np.unique(u.atoms.resids)),len(np.unique(u.atoms.resids))+1)]

        #list of combinations of residues/renum of the dimers of the different polymers
        b=list(combinations(a, 2))  

        #resid same as resnum
        all_comb=[['{} {}'.format('resnum',values[0]), '{} {}'.format('resnum', values[1])]
                  for idx, values in enumerate(b)] 

        dist_tot=np.array([self.nearest_neighbours(u, str(i), str(j)) for i, j in tqdm(all_comb)])

        #returns list with the pair indices of the polymers that have at least two atoms that are in contact
        index_cutoff=np.where(dist_tot < float(cutoff))[0] 
    
        return index_cutoff,b
    
    @timeit
    def stacking_analysis(self, cut_off, ang_cut, output_name="output.parquet"):
        """Function to perform a stacking analysis
        """
        #How to check the different number of residues being sure that you
        #are taking a polymer and not a water molecule or solvent.
        import pandas as pd
        
        print('Ring Stacking analysis has started')


        #MDAnalysis Universe
        #u=mda_to_rdkit_sel("f8bt_slab_quench.tpr", "less_frames_trajectory.xtc") 
        u=super().get_mda_universe()
        
        #Get the indices of b of the polymers that are in contact
        #Dimers with smallest distances according to radius cutoff
        index_cutoff,b=self.pol_cutoff(u, float(cut_off)) 
     
        dimer_mol=[b[i] for i in index_cutoff] #list of resnum of polymers that are in contact

        #Array with atom indices of the polymers that are in contact given in only one array
        dimer_ind=[u.select_atoms('resnum {} {}'.format(values[0], values[1])).indices
                   for idx, values in enumerate(dimer_mol)]  

        #Counting of number of atoms in a single polymer
        num_atoms=len(u.select_atoms('resnum 1').residues.atoms)

        #Printing in this order: dimer_1, dimer_2: min dimer_1, max dimer_1, min dimer_2, max dimer_2
        #first and last atom indices of the polymers that are in contact.
        
        dimers_indexes=[(min(values[:num_atoms]),max(values[:num_atoms]), min(values[num_atoms:]),
                         max(values[num_atoms:])) for idx, values in enumerate(dimer_ind)] 
  
        us = [u]*len(dimers_indexes)  
        dim2 = list(tqdm(map(self.GetRingSystems, us, [x[0] for x in dimers_indexes],
                             [x[1] for x in dimers_indexes], [x[2] for x in dimers_indexes],
                             [x[3] for x in dimers_indexes]), total=len(dimers_indexes),
                             desc="Detecting atoms in Rings"))


        u4 = [u]*len(dim2)

        dim3 = tqdm(map(self.separate_rings, dim2, [x[2] for x in dimers_indexes]),
                    total=len(dim2), desc="Separating Rings")

        c = list(dim3)

        u4 = [u]*len(c)
        cuts = [float(cut_off)]*len(c)
        ang_cuts = [float(ang_cut)]*len(c)

        dim4 = list(tqdm(map(self.run_rings_stacking_dist, u4, c, cuts, ang_cuts),
                         total=len(dim3), desc="Computing stacking distances"))
    
        results = [result for result in dim4 if result is not None]
        results1, results2 = zip(*results)

        #Create a pandas DataFrame from the results
        print ("\n")
        print ("Preparing DataFrame to store the results")
    
    
        df = pd.DataFrame({'atom_index': results1,
                           'pol_resid': results2})

        # Save the DataFrame to a file
        df.to_parquet(str(output_name), index=False)
    
        #rings_df = pd.read_parquet(str(output_name))
        #results_f = self.find_several_rings_stacked(rings_df)
        #print (results_f)
        print ("Information succesfully stored in {}".format(str(output_name)))
        print ("Stacking analysis has succesfully finished!")
    


