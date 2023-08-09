Basic Use
=========

.. _importing:

Importing pySoftWhere
----------------------

.. code-block:: python3

    import MDAnalysis as mda
    import matplotlib.pyplot as plt
    import numpy as np
    import sys


Set the directory where pySoftWhere is located (if a full install has not been completed)

.. code-block:: python3

    sys.path.append('<path to>/pysoftwhere') 

Import the pySoftWhere functions

.. code-block:: python3

    import pysoftwhere.clustering.cluster
    import pysoftwhere.ICSI.icsi
    import pysoftwhere.IPI.ipi



.. _ipi_simple:

Intrinsic Planar Interface (IPI)
----------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.IPI.ipi()`` to calculate the intrinsic density for a single frame

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

    intrinsic_z,interface_vals  =  pysoftwhere.IPI.ipi(u,
                                                       frame=-1,
                                                       core_sel='resname UNK and prop mass >2',
                                                       density_sel='name OW',
                                                       interface='Lower',                        
                                                       no_bins=51,
                                                       recombine=True,                          
                                                       cluster=True,                            
                                                       interpolate_interface=False)()           


The ``pysoftwhere.IPI.ipi()`` function returns two results when called. These are ``intrinsic_z``,  a list of the intrinsic positions of atoms selected in ``density_sel``, and ``interface_vals``, an array containing the z-positions of the intrinsic interface.

 
.. _icsi_simple:

Intrinsic Core-Shell Interface (ICSI)
--------------------------------------

Firstly, define an MDAnalysis Universe

.. code-block:: python3

   u = mda.Universe('system.tpr','system.xtc')


Now call ``pysoftwhere.ICSI.icsi()`` to calculate the intrinsic density for a single frame. Note that ``pysoftwhere.ICSI.icsi()`` uses atomic positions, rather than an atom selection as used for IPI. This difference is made to allow flexibility that becomes useful in more complex use cases.

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

    
    intrinsic_r, spherical_r, icsi_vals = pysoftwhere.ICSI.icsi(u, 
                                                        
                                    cluster_resids=largest_cluster_resids,
                                    cluster_atoms_positions=cluster_atoms_positions,
                                    core_sel_atoms_positions=core_sel_atoms_positions,
                                    shell_sel_atoms_positions=shell_sel_atoms_positions,
                                    frame=ts.frame,
                                    no_bins=31,
                                    no_random_points=n_rand_points,
                                    normalisation_run=False)()

The ``pysoftwhere.IPI.ipi()`` function returns three results when called. These are ``intrinsic_r``,  a list of the intrinsic positions of atoms selected in ``density_sel``, ``spherical_r``, the radial distance of each ``density_sel`` atom to the nanoparticle center of mass,  and ``icsi_vals``, an array containing the r-positions of the intrinsic core-shell interface. 
