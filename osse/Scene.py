'''
Class that defines properties of exospheric scene and returns radiance profiles 

TODO:
- Implement call function for any length of passed LOS variables
- Implement call for spatially varying los (tail)
- Implement call for global temporal enhancements
'''
import numpy as np
from multiprocessing import Pool

from scip.data_processing.helpers import *

import radiative_transfer.rt_common as rt_common
import radiative_transfer.rt_msis as rt_msis
import radiative_transfer.rt_thermosphere as rt_thermosphere
import radiative_transfer.rt_h_dens as rt_h_dens

def point_rad(hexo, vg, hdm, sf):
    '''
    get the radiance of a single passed point
    Needs to be globally defined for multiprocessing
    '''
    return rt_common.run_forward_los(hdm, vg, hexo, sf, ncore=0)[0, 5]


class Scene:
    def __init__(self, hexo, solar_flux, solar_activity, twoD=False):
        '''
        FIXME: Explain variables
        NOTE: Tail factor should be >1
        Contains exosphere and solar activity parameters
        FIXME: Shouldn't contain instrument pointing 
        '''
        # Base exosphere parameters
        self.hexo = hexo # Mean exobase density
        self.solar_flux = solar_flux
        self.solar_activity = solar_activity

        self.twoD = twoD
        if twoD: # Define exobase varying w/solar angle
            self.A = 1 # 2D Toggle variable
            self.t = 0
        else:
            self.A = 0 # Verified: Setting A=0 causes RT output to be equivalent to 1D

    def __call__(self, los):
        '''
        FIXME Docstring
        FIXME SLOW AS HELL
        FIXME Can use 2D density distribution for geotail ==> Refactor!

        Args: 
            los (ndarray) - Pregenerated viewing geometry for some schedule
        '''
        # Thermosphere generation
        msis = rt_msis.MSIS(self.solar_activity, 0, 90) # FIXME Defaults to [0,90] tpsge coord
        exo = msis.exobase()
        mlt_saber = rt_thermosphere.MLT_SABER()
        self.exo_obj = exo
        # FIXME: need array of TPGSE coords to be accurate
        
        thm_obj = rt_thermosphere.Thermosphere(msis, mlt_saber, exo)
        hdm = rt_h_dens.HDensityBishop(thm_obj, A=self.A, hsat=0)
        # A = 0 should have same output as 1D model

        # Run RT
        out = rt_common.run_forward_los(hdm, np.array(los), self.hexo, self.solar_flux)

        return np.array(out[:,5])