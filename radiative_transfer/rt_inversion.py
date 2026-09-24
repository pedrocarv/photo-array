# !/usr/bin/env python3

import os
import numpy as np
import scipy as sp
import warnings
# warnings.filterwarnings("ignore")
import multiprocessing

import radiative_transfer.rt_common as rt_common
import radiative_transfer.rt_msis as rt_msis
import radiative_transfer.rt_thermosphere as rt_thermosphere
import radiative_transfer.rt_h_dens as rt_h_dens
#TODO: order inputs consistently

def _normalize_profiles(calculated, observed, obs_per_los=None):
    """Normalize matching profiles, optionally in equal-length segments."""
    calculated = np.asarray(calculated, dtype=float).reshape(-1)
    observed = np.asarray(observed, dtype=float).reshape(-1)

    if calculated.size != observed.size:
        raise ValueError(
            "Calculated and observed radiance profiles must have the same "
            f"length ({calculated.size} != {observed.size})."
        )
    if calculated.size == 0:
        raise ValueError("Radiance profiles must contain at least one observation.")

    if obs_per_los in (None, -1):
        obs_per_los = calculated.size
    elif not isinstance(obs_per_los, (int, np.integer)) or obs_per_los <= 0:
        raise ValueError("obs_per_LOS must be a positive integer, None, or -1.")

    if calculated.size % obs_per_los != 0:
        raise ValueError(
            "The radiance profile length must be divisible by obs_per_LOS "
            f"({calculated.size} % {obs_per_los} != 0)."
        )

    calculated_normalized = np.empty_like(calculated)
    observed_normalized = np.empty_like(observed)
    for start in range(0, calculated.size, obs_per_los):
        stop = start + obs_per_los
        calculated_scale = calculated[start]
        observed_scale = observed[start]
        if calculated_scale == 0 or observed_scale == 0:
            raise ValueError("Cannot normalize a radiance segment whose first value is zero.")
        calculated_normalized[start:stop] = calculated[start:stop] / calculated_scale
        observed_normalized[start:stop] = observed[start:stop] / observed_scale

    return calculated_normalized, observed_normalized


def cost_function(hexo, hdm, obs, orbit_info, solar_flux, obs_per_LOS=None):
    """
    Cost function optimizing hexo based on the output norm.

    Args:
        hexo (float): H density at exobase.
        hdm (HDensityModel): H density model object.
        channel (str): Type of data, instrument channel; 1x, NFI, WFI, inv6
        obs (ndarray): Observed radiance profile.
        orbit_info (ndarray): Orbit info (viewing geometry).
        solar_flux (float): Solar flux at line center.

    Returns:
        norm (float): Error metric.

    """
    lyao_los = rt_common.run_forward_los(hdm, orbit_info, hexo, solar_flux)

    # A single continuous profile is the default.  If multiple instruments or
    # profiles are concatenated, obs_per_LOS normalizes each segment separately.
    rad_rt, rad_obs = _normalize_profiles(lyao_los[:, 5], obs, obs_per_LOS)

    norm = (((rad_rt - rad_obs) ** 2).sum())
    return norm


def optimize_profile(tpgse_coord, radiance, orbit_info, solar_flux, solar_activity, model='Chamberlain', obs_per_LOS=None):
    """ Invert (optimize, match) procedure. Minimize Root Mean Squared Error (RMSE). Return calculated
    H exo and indices of a profile and wedge.

        Args:            
            tpgse_coord (list): A wedge coordinates.
            radiance (ndarray): (1x6) Radiance profile.
            orbit_info (ndarray): Orbit info.
            solar_flux (float): A solar flux.
            solar_activity (SolarActivity): Solar activity object containing date, aps, f10.7 and f10.7a.
            channel (str, optional): Type of data, instrument channel; 1x, NFI, WFI, for inversion it is always inv6
            model (str, optional): Model type to use - Chamberlain or Bishop. Default is Chamberlain.
            obs_per_LOS (int, optional): Number of observations in each
                independently normalized profile. By default, treat all
                observations as one continuous profile.

        Return:
            float: hexo
            float: sf_real, calculated solar flux.
            int: number of iterations to converge

    """
    
    #Set up H density model
    msis = rt_msis.MSIS(solar_activity, tpgse_coord[0], tpgse_coord[1])
    exo = msis.exobase()
    mlt_saber = rt_thermosphere.MLT_SABER()
    thm_obj = rt_thermosphere.Thermosphere(msis, mlt_saber, exo)
    
    if model == 'Bishop':
        hdm = rt_h_dens.HDensityBishop(thm_obj, hsat=0)
    elif model == 'Chamberlain':
        hdm = rt_h_dens.HDensityChamberlain(thm_obj)
    else:
        raise ValueError('Select Bishop or Chamberlain')

    # convert 1 x 6 two-dimensional array to a 6 element one-dimensional array.
    rad_obs = np.squeeze(radiance)
    
    #Error checks
    

    # 9.9e3,3e5 are the physical boundaries of hexo, we don't want to go bellow and above
    results = sp.optimize.fminbound(cost_function, 9.9e3, 3e5,
                                    args=(hdm, rad_obs, orbit_info, solar_flux, obs_per_LOS),
                                    full_output=True, xtol=10, maxfun=30, disp=3)
    try:
        minimizer, fval, ierr, numfunc = results

        hexo = minimizer

        los_results = rt_common.run_forward_los(hdm, orbit_info, float(hexo), solar_flux)
        rad_rt = los_results[:, 5]

        fac = np.mean(rad_obs / rad_rt)
        sf_real = fac * solar_flux

    except TypeError as err:
        sf_real = solar_flux
        print(err)

    return [hexo, sf_real, numfunc]


def inversion_init(image, orbit_info, tpgse, solar_activity, solar_flux, model='Chamberlain', cores=24):
    """ Initialize H_exo calculations (inversions). Assumes image has been binned into nr x nth curvilinear
    pixels. Each wedge will be inverted independently to return nth h_exo, solar_flux, z_exo, t_exo.

        Args:
            image (ndarray): Binned profiles to be inverted. [nr x nth] array
            orbit_info (ndarray): Orbit info consisting of r_obs, sza, los_zenith, los_azimuth in each column.
                Array is [nr x nth x 4].
            tpgse (ndarray): Tangent points for each profile (GSE lat, lon). [nth x 2] array
            solar_activity (SolarActivity): Solar activity object containing date, aps, f10.7 and f10.7a.
            solar_flux (float): Initial solar flux.
            model (str, optional): H density model to use - either Chamberlain or Bishop. Default is Chamberlain.
            cores (int, optional): cores - 2, 4, 8, 16, 24 (1 profile, default), 32, 48, 64. Set total processors
                         (cores) used. This option splits the 24 wedges into equal batches where each batch
                         of LOS is assigned to a new processor and runs in parallel.

        Return:
            hexo (ndarray): Retrieved H density at the exobase for each wedge.
            solar_flux (ndarray): Retrieved solar flux at Ly-a center for each wedge.

    """
    
    #Error check input dimensions
    nth = image.shape[1]
    assert orbit_info.shape[0] == image.shape[0], 'Orbit info must be nr x nth x 4'
    assert orbit_info.shape[1] == nth, 'Orbit info must be nr x nth x 4'
    assert orbit_info.shape[2] == 4, 'Orbit info must be nr x nth x 4'
    assert tpgse.shape[0] == nth, 'tpgse must be nth x 2'
    assert tpgse.shape[1] == 2, 'tpgse must be nth x 2'
    
    
    #Construct input list
    items = [(tpgse[ii], image[:,ii], orbit_info[:,ii,:], solar_flux, solar_activity, model) for ii in
            range(nth)]
    if cores > 0:
        # The Fortran RT code is OpenMP-threaded and its LU solve uses a threaded
        # LAPACK (Accelerate on macOS, usually OpenBLAS on Linux); give each worker
        # process an equal share of the CPUs so cores x threads doesn't
        # oversubscribe the machine. Workers inherit the environment, and thread
        # counts the user set explicitly are kept.
        per_worker = str(max(1, (os.cpu_count() or 1) // cores))
        set_vars = [v for v in ('OMP_NUM_THREADS', 'OPENBLAS_NUM_THREADS', 'VECLIB_MAXIMUM_THREADS')
                    if v not in os.environ]
        for v in set_vars:
            os.environ[v] = per_worker
        try:
            # 'spawn' (the macOS default) instead of Linux's 'fork': GNU OpenMP
            # (libgomp) is not fork-safe, and a worker forked after this process
            # has run the Fortran code hangs in its first parallel region.
            with multiprocessing.get_context('spawn').Pool(processes=cores) as pool:
                res = pool.starmap(optimize_profile, items)
        finally:
            for v in set_vars:
                del os.environ[v]
    else:
        res_all = []
        for ii in range(nth):
            res = optimize_profile(tpgse[ii], image[:,ii], orbit_info[:,ii,:], solar_flux, model)
            res_all.append(res)
    res_all = np.vstack(res)
    
    #Parse results
    hexo = res_all[:,0]
    solar_flux = res_all[:,1]
    num_iter = res_all[:,2]

    return hexo, solar_flux, num_iter
