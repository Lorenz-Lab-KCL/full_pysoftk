import pytest
import os
from functools import wraps
from  utils_mda import MDA_input
from utils_tools import *
from clustering import SCP
from make_micelle_whole import micelle_whole
from ring_ring import RSA
#from pysoftk.pol_analysis.ring_ring import RSA
import numpy as np
import filecmp
import pandas as pd




@pytest.fixture
def rootdir():
    return os.path.dirname(os.path.abspath(__file__))


def test_RSA(rootdir):

    from pandas.testing import assert_frame_equal

    topology=os.path.join(rootdir, 'data/f8bt_slab_quench.tpr')
    trajectory=os.path.join(rootdir, 'data/1_frame_traj.xtc')

    rsa_og = os.path.join(rootdir, 'data/rsa.parquet')

    rsa_og_df= pd.read_parquet(rsa_og)

    
    results='data/rsa_test.parquet'

     
    ang_c=30

    
    dist_c=5

    rsa=RSA(topology, trajectory).stacking_analysis(dist_c, ang_c, results)

    rsa_df=pd.read_parquet(results)
    
    assert_frame_equal(rsa_df, rsa_og_df)
