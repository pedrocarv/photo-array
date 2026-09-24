# !/usr/bin/env python3

import numpy as np
import warnings
# warnings.filterwarnings("ignore")

import glide.science.radiative_transfer.rt_common as rt_common
import glide.science.radiative_transfer.rt_msis as rt_msis
import glide.science.radiative_transfer.rt_thermosphere as rt_thermosphere
import glide.science.radiative_transfer.rt_h_dens as rt_h_dens


def run_forward_albedo(h_model, h_exo, n_sza=64):
    """ Calculates albedo by running the Forward model. Albedo is the Source function (ts).

        Args:
            h_model (obj): H Model object.
            h_exo (float): H density at the exobase.
            n_source (int, optional): Number of SZA points in source function.

        Return:
            ndarray: Total source function (tsrc).

    """
    #First step: set h_exo and z
    h_model.set_h(z=h_model.Z_ALBEDO, h_exo=h_exo)
    
    exopt, cdens, thermo, brad, bchi, chin, rhot, ts = rt_common.run_forward(h_model, n_sza=n_sza)
    
    tsrc = ts[:,:,1]

    return tsrc


def create_albedo(tpgse_coord, h_exo, solar_flux, solar_activity, model = 'Chamberlain'):
    """ Invert (optimize, match) procedure. Minimize Root Mean Squared Error (RMSE). Return calculated
    H exo and indices of a profile and wedge.

        Args:
            tpgse_coord (ndarray): Tangent point (lat, lon) in GSE coordinates to generate albedo at.
            h_exo (float): H density at the exobase.
            solar_flux (float): Solar flux at line-center.
            solar_activity (SolarActivity): Solar activity object.
            model (str, optional): H density model to use. Default is Chamberlain.

        Return:
            ndarray: A source function.

    """
    # Default grids (zgrid=None and szas=None)  are defined in InputRTModel class. Grids can be changed
    # by supplying np.array zgrid and szas to the InputRTModel here.
    # NOTE: Length of zgrid must be the same as PARAMETER IKNT in global_parameters.f and lyao_rt.f
    # (lyao_los.f not used in albedo)
    # and
    # length of szas must be the same as PARAMETER JKNT in global_parameters.f and lyao_rt.f
    # (lyao_los.f not used in albedo)
    # Additionally PARAMETER MDIM in lyao_rt.f must be IKNT * JKNT
    
    #Set up H density model
    msis = rt_msis.MSIS(solar_activity, tpgse_coord[0], tpgse_coord[1])
    exo = msis.exobase()
    exo.h_exo = h_exo
    mlt_saber = rt_thermosphere.MLT_SABER()
    thm_obj = rt_thermosphere.Thermosphere(msis, mlt_saber, exo)
    
    if model == 'Bishop':
        hdm = rt_h_dens.HDensityBishop(thm_obj)
    elif model == 'Chamberlain':
        hdm = rt_h_dens.HDensityChamberlain(thm_obj)
    else:
        raise ValueError('Select Bishop or Chamberlain')
        
    tsrc = run_forward_albedo(hdm, h_exo)

    return tsrc


