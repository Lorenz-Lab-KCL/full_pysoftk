import os
from functools import wraps
import numpy as np

from pysoftk.pol_analysis.tools.utils_mda import MDA_input
from pysoftk.pol_analysis.tools.utils_tools import *
from pysoftk.pol_analysis.clustering import SCP
from pysoftk.pol_analysis.make_micelle_whole import micelle_whole
from pysoftk.pol_analysis.intrinsic_density_micelle import intrinsic_density

topology=os.path.join(os.getcwd(), 'intrinsic_density_movie.tpr')
trajectory=os.path.join(os.getcwd(), 'intrinsic_density_movie.xtc')
resids = os.path.join(os.getcwd(), 'cyclic_id.parquet')

MA_names_c = ['C028', 'C024', 'C022', 'C02E', 'C02I', 'C02M', 'C02Q', 'C02U', 'C02Y', 'C00S', 'C00Q', 'C00M', 
             'C00L', 'C027',  'C023',  'C021',  'C02H',  'C02L',  'C02P',  'C02T', 'C02X', 'C00U',  'C00R',  'C00P', 
             'O00A', 'O00D', 'O00C', 'O01G', 'O01F', 'O01D', 'O01A', 'O019', 'O016', 'O015', 'O010', 'O012', 
             'O00B', 'O00E', 'O01I', 'O01H', 'O01C', 'O01B', 'O01E', 'O018', 'O017', 'O014', 'O00Z', 'O013', 
             'C00O', 'C00W', 'C031', 'C030', 'C02W', 'C02S', 'C02O', 'C02K', 'C02G', 'C02D', 'C026', 'C02B']

resids = micelle_whole(topology, trajectory).obtain_largest_micelle_resids(resids)
atom_pos = micelle_whole(topology, trajectory).running_make_cluster_whole(['LIG'], resids, 75, 10000, 1)
density_whole = intrinsic_density(topology, trajectory).run_intrinsic_density(resids, atom_pos, MA_names_c, MA_names_c, ['OW'], 75, 10000, 1)

print  (list(density_whole))
