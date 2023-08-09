Intrinsic Planar Interface (IPI)
=========


.. _ipi_density:

Calculating an intrinsic density with IPI
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.IPI.ipi()`` to calculate the intrinsic density for multiple frames

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

    intrinsic_z_ = []
    interface_vals_ = []
    
    for frame_sel in range(0,100):
    
        intrinsic_z,interface_vals  =  pysoftwhere.IPI.ipi(u,
                                                           frame=frame_sel,
                                                           core_sel='resname UNK and prop mass >2',
                                                           density_sel='name OW',
                                                           interface='Lower',                        
                                                           no_bins=51,
                                                           recombine=True,                          
                                                           cluster=True,                            
                                                           interpolate_interface=False)()           
        intrinsic_z_.append(intrinsic_z)
        interface_vals_.append(interface_vals)



.. _ipi_tracking:


Tracking solute molecules with IPI
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.IPI.ipi()`` to calculate the intrinsic positions of solute atoms for multiple frames

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
 


.. _ipi_profiles:

Custom Intrinsic Profiles with IPI
----------------------------------

We can use ``pysoftwhere.IPI.ipi()`` to easily create intrinsic profiles of any quantity that we may wish to calculate. For example, intrinsic profiles can be constructed for the number of hydrogen bonds in the solute, hydration number, solute mobility, orientation, and so. 


