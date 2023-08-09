Intrinsic Core-Shell Interface (ICSI)
=========


.. _icsi_density:

Calculating an intrinsic density with ICSI
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')

``pysoftwhere.ICSI.icsi()`` has been designed for maximum flexibility. While calculating the intrinsic density profile of a single micelle, centered in the middle of the simulation box over time can be acomplished by calling ``pysoftwhere.ICSI.icsi()`` immediately, there are more complicated problems that may require the use of pySoftWhere's auxiliary functions. For example, there might be multiple nanoparticles in the simulation (like small surfactant clusters), which span periodic boundary conditions (with nanoparticles greater in size than half of the simulation box size), and so on.  


.. code-block:: python3
    

    ''' create lists to save data '''

    intrinsic_r_tx114,intrinsic_r_rand_tx114 , box_size = [] , [] , []
    
    
    for ts in u.trajectory[start:stop:skip]:
        
        #find the largest cluster (this is the micelle)
        
        largest_cluster_resids=pysw_cluster.find_largest_cluster(u,
                                       frame=ts.frame,
                                       selection='name C*',
                                       cutoff_distance=10,
                                       define_clustering_atoms=False)()
    
        #unwrap the micelle coordinates
        
        cluster_atoms_positions,core_sel_atoms_positions,shell_sel_atoms_positions=pysw_cluster.make_cluster_whole(u,
                                        frame = ts.frame,
                                        cluster_resids  = largest_cluster_resids,
                                        core_selection  = micelle_core,
                                        shell_selection = micelle_shell)()
        
        #calculate the intrinsic positions of the micelle shell using ICSI
        
        intrinsic_r, spherical_r, icsi_vals = pysw_icsi.icsi(u, 
                                                            
                                        cluster_resids=largest_cluster_resids,
                                        cluster_atoms_positions=cluster_atoms_positions,
                                        core_sel_atoms_positions=core_sel_atoms_positions,
                                        shell_sel_atoms_positions=shell_sel_atoms_positions,
                                        frame=ts.frame,
                                        no_bins=11,
                                        no_random_points=n_rand_points,
                                        normalisation_run=False)()
    
        #bin and save the intrinsic positions 
        
        intrinsic_r_tx114.append(stats.binned_statistic(intrinsic_r,intrinsic_r,bins=np.arange(-40.5,150,0.5),statistic='count').statistic)
            
        #calculate the normalisation factor 
        
        intrinsic_r_rand, spherical_r_rand, icsi_vals = pysw_icsi.icsi(u, 
                                                            
                                        cluster_resids=largest_cluster_resids,
                                        cluster_atoms_positions=cluster_atoms_positions,
                                        core_sel_atoms_positions=core_sel_atoms_positions,
                                        shell_sel_atoms_positions=shell_sel_atoms_positions,
                                        frame=ts.frame,
                                        no_bins=11,
                                        no_random_points=n_rand_points,
                                        normalisation_run=True)()
        
        #bin and save the normalisation positions 
    
        intrinsic_r_rand_tx114.append(stats.binned_statistic(intrinsic_r_rand,intrinsic_r_rand,bins=np.arange(-40.5,150,0.5),statistic='count').statistic)
        
        #save the box size for each timestep
    
        box_size.append(u.dimensions[0]*u.dimensions[1]*u.dimensions[2])
        
        print(ts)
        
        
    ###calculate average count for the intrinsic distance vector - CHECK!
    intrinsic_r_tx114_profile=np.sum(np.array(intrinsic_r_tx114),axis=0) / len(np.arange(start,stop,skip))
    
    ###calculate the normalisation vector
    S_bar=np.sum(np.array(intrinsic_r_rand_tx114),axis=0) *np.mean(box_size) / (len(np.arange(start,stop,skip))*n_rand_points)


.. _ipi_tracking:


Tracking solute molecules with ICSI
----------------------------------


Now call ``pysoftwhere.ICSI.icsi()`` to calculate the intrinsic positions of solute atoms for multiple frames

.. code-block:: python3
    
    
    ''' 
    u:           an MDAnalysis Universe that contains bond information
    core_sel:    an MDAnalysis atom selection of the atoms used to construct the intrinsic interface
    density_sel: an MDAnalysis atom selection of the atoms whose intrinsic positions will be calculated
    interface:   'Lower' or 'Upper' - which interface of the core_sel to calculate with respect to
    no_bins:     number of grid edges to use in each lateral dimension
    recombine:   'True' or 'False' - set 'True' to recombine core_sel over the PBC in the vertical direction
    cluster:     'True' or 'False' - set 'True' to check if any core_sel molecules have diffused away from the main selection and remove them from the analysis
    interpolate_interface: 'True' or 'False' - set "True' to linearly interpolate the interface if it is patchy
    '''

    solute_positions = []
    
    for frame_sel in range(0,100):
    
        intrinsic_z,interface_vals  =  pysoftwhere.IPI.ipi(u,
                                                           frame=frame_sel,
                                                           core_sel='resname SLAB and prop mass >2',
                                                           density_sel='resname SOLUTE and name C1',
                                                           interface='Lower',                        
                                                           no_bins=51,
                                                           recombine=True,                          
                                                           cluster=True,                            
                                                           interpolate_interface=False)()           
        solute_positions.append(intrinsic_z)

Now we can easily access the intrinsic solute positions as a function of time. These can be used to find out when solute molecules penetrate into the ``core_sel`` atoms, which might represent a polymer slab, surfactant monolayer, lipid bilayer, and so on.
 


