# !/usr/bin/env python3


import numpy as np
import warnings
import multiprocessing
import importlib



class Exobase:
    """ A class for holding h_exo, z_exo, t_exo  and lat, lon location(s)""" 
    h_exo: np.ndarray
    z_exo: np.ndarray
    t_exo: np.ndarray
    lat: np.ndarray
    lon: np.ndarray
    
    def __init__(self, h_exo, z_exo, t_exo, lat, lon, szas=None):
        """
        Initialize Exobase. Inputs may be floats, lists, or arrays.

        Args:
            h_exo (ndarray): H density at exobase.
            z_exo (ndarray): Altitude (km) at exobase.
            t_exo (ndarray): Temperature (K) at exobase.
            lat (ndarray): Latitude grid points in GSE coordinates.
            lon (ndarray): Longitude grid points in GSE coordinates.
            TODO: sza array starting from here
            
        Returns:
            None.

        """
        
        
        self._h_exo = np.atleast_1d(h_exo)
        self.z_exo = np.atleast_1d(z_exo)
        self.t_exo = np.atleast_1d(t_exo)
        self.lat = np.atleast_1d(lat)
        self.lon = np.atleast_1d(lon)
        
        self._error_check(self._h_exo)
        
    
    @property
    def h_exo(self):
        return self._h_exo
    
    @h_exo.setter
    def h_exo(self, h_val):
        self._error_check(np.atleast_1d(h_val))
        self._h_exo = np.atleast_1d(h_val)
        
    
    @property
    def h_flat(self):
        return np.ravel(self.h_exo)
    
    @property
    def z_flat(self):
        return np.ravel(self.z_exo)
    
    @property
    def t_flat(self):
        return np.ravel(self.t_exo)
    
    @property
    def npt(self):
        return self.lat.size
    
    def __repr__(self):
        if self.npt == 1:
            return f'Exobase({self.h_exo[0]}, {self.z_exo[0]}, {self.t_exo[0]}, [{self.lat[0]}, {self.lon[0]}])'
        elif self.h_exo.ndim > 1:
            shape = self.h_exo.shape
            return f'Exobase({shape[0]} x {shape[1]} points)'
        else:
            return f'Exobase({self.npt} points)'
        
    def _error_check(self, h_exo):
        
        #Error check
        assert h_exo.shape == self.z_exo.shape, "H exo and Z exo must be same shape."
        assert h_exo.shape == self.t_exo.shape, "H exo and T exo must be same shape."
        
        if self.lat.shape != h_exo.shape:
            assert self.lat.size == h_exo.shape[0], "Latitudes must be same shape as H exo or same size as first dimension of H exo."
        if self.lon.shape != h_exo.shape:
            assert self.lon.size == h_exo.shape[1], "Longitudes must be same shape as H exo or same size as second dimension of H exo."


def los_parallel(losx, losy, solar_flux, orbit, exopt, cdens, thermo, brad, bchi, chin, rhot, ts,
                 rt_params_int, rt_params_real, ap_hist, thbkg, zgrid):
    """ Pool parallelization for calculating of los. Used only in Forward (ground truth) model.

        Args:
            channel: (str): channel (str): Type of data, instrument channel; 1x, NFI, WFI
            losx (obj): losx.
            losy (obj): losy.
            solar_flux (float): solar_flux.
            orbit (ndarray): orbit.
            exopt (ndarray): exopt.
            cdens (ndarray): cdens.
            thermo (ndarray): thermo.
            brad (ndarray): brad. Radial grid of source function
            bchi (ndarray): bchi. Solar angle grid of source function
            chin (ndarray): chin.
            rhot (ndarray): rhot.
            ts (ndarray): ts.
            rt_params_int (ndarray): rt_params_int.
            rt_params_real (ndarray): rt_params_real.
            ap_hist (ndarray): ap_hist.
            rtbkg (ndarray): rtbkg.
            thbkg (ndarray): thbkg.
            zgrid (ndarray): zgrid.

        Return:
            ndarray: LOS results.

    """
    los_results = None

    try:
        los = importlib.import_module(f'.los_1x', package='radiative_transfer')

        # No rtbkg and solar angle grid!
        los_results = los.ft_main(losx, losy, solar_flux, orbit, exopt, cdens, thermo, brad, bchi, chin, rhot, ts,
                                  rt_params_int, rt_params_real, ap_hist, thbkg, zgrid)
    except ModuleNotFoundError as err:
        raise SystemExit(err)
    return los_results


def run_forward(h_model, n_sza=32):
    """ Calculates initial radiance and background files by running the Forward and LOS models.

        Args:
            h_model (HDensityModel): H Density model object.
            n_sza (int, optional): Number of SZAs in source function. Default is 32.

        Return:
            ndarray: exopt. Array contining optical and exospheric parameters (e.g. exobase temps and dens, satellite info, branching ratio, etc.)
            ndarray: cdens. 'Effective Zenith Column Densities'
            ndarray: thermo. Thermosphere info
            ndarray: brad. Source Function radial grid points
            ndarray: bchi. Source Function solar angle grid points (solar angle ==> angle from sun-earth line)
            ndarray: chin. Solar angle centroids (points in the middle of defined grid boundaries)
            ndarray: rhot. Radial Centroid Points (0), H Density (1), O2 Density (2), and temperature [?] (3)
            ndarray: ts. Solar transmission ([:,:,0]) and source function ([:,:,1])

    """

    
    rt_bkg, thermo_bkg, rt_params_int, rt_params_real = h_model.to_fortran_input(ncopy=1)
    zgrid = h_model.z
    ap_hist = h_model.msis.aps
    sagrid = h_model.solar_ang
    sagridx = len(sagrid)

    # Fortran array dimensions 
    iknt = zgrid.size
    jknt = n_sza

    try:
        if iknt == 32 and jknt == 32:
            import radiative_transfer.forward as forward
        elif iknt == 76 and jknt == 64:
            import radiative_transfer.forward_hres as forward
        else:
            raise ValueError('Fortran iknt, jknt sizes must be either 32 x 32 or 76 x 64.')

        rtbkgx, rtbkgy, rtbkgz = rt_bkg.shape # (3, 1, iknt)
        thbkgx, thbkgy, thbkgz = thermo_bkg.shape # (7, 1, 61), (itherm+1)
        
        #Error check thermo_bkg size
        assert thermo_bkg.shape[0]==7 and thermo_bkg.shape[2]==61, 'thermo_bkg must be 7 x 1 x 61'

        

        # No rtbkg
        exopt, cdens, thermo, brad, bchi, chin, rhot, ts = forward.ft_main(14, 4, 5, thbkgz, iknt + 1,
                                                                           jknt + 1, jknt, iknt, 4, iknt, jknt, 2,
                                                                           rt_params_int, rt_params_real, ap_hist,
                                                                           thermo_bkg, zgrid, h_model.A, h_model.w, h_model.t,
                                                                           [5, 16, len(ap_hist),
                                                                            thbkgx, thbkgy, thbkgz, iknt])
    except ModuleNotFoundError as err:
        raise SystemExit(err)

    #EXPLAIN: What are all these outouts?
    return exopt, cdens, thermo, brad, bchi, chin, rhot, ts


def run_los(h_model, orbit_info, solar_flux, exopt, cdens, thermo, brad, bchi, chin, rhot,
            ts, ncore=1):
    """ Runs LOS code using outputs from run_forward.

        Args:
            h_model (HDensityModel): H Density Model object.
            orbit_info (ndarray): Orbit info data.
            solar_flux (float): Solar flux.
            exopt (ndarray): exopt.
            cdens (ndarray): cdens.
            thermo (ndarray): thermo.
            brad (ndarray): brad.
            bchi (ndarray): bchi.
            chin (ndarray): chin.
            rhot (ndarray): rhot.
            ts (ndarray): ts.
            parallel (bool): LOS code run in parallel with Pool module. Make sure max_os and import los match.
                             Default False

        Return:
            ndarray: Lyao source file with radiance and ratio.

    """

    rt_bkg, thermo_bkg, rt_params_int, rt_params_real = h_model.to_fortran_input(ncopy=1)
    zgrid = h_model.z
    ap_hist = h_model.msis.aps

    # output array
    #   NFI 1024*1024, WFI 512*512
    #   sequential    1x 1024     NFI 1024*1024   WFI 512*512
    #   parallel 8    1x  128     NFI  128*1024   WFI 128*512
    #   parallel 16   1x   64     NFI   64*1024   WFI  64*512
    #   parallel 32   1x   32     NFI   32*1024   WFI  32*512
    #   inversion     6 (INV6   MAXLOS = 6)
    losx = orbit_info.shape[0]
    losy = 7

    if ncore == 0:
        #Error check
        assert losx % ncore == 0, 'LOS must be divisible by number of cores requested.'
        
        max_los = losx // ncore  # 8:128, 16:64
        orbit_par = []
        for i in range(0, ncore):
            low = i * max_los
            high = (i + 1) * max_los

            orbit_par.append(orbit_info[int(low):int(high), :])

        # No rt_bkg
        items = [(max_los, 7, solar_flux, orbitp, exopt, cdens, thermo, brad, bchi, chin, rhot,
                  ts, rt_params_int, rt_params_real, ap_hist, thermo_bkg, zgrid) for orbitp in orbit_par]
        # 'spawn': GNU OpenMP (libgomp) is not fork-safe, see rt_inversion.inversion_init
        with multiprocessing.get_context('spawn').Pool(processes=ncore) as pool:
            res = pool.starmap(los_parallel, items)

        los_results = np.vstack(res)
    else:
        try:
            # FIXME: Forcing to use the 1x driver
            los = importlib.import_module(f'.los_1x', package='radiative_transfer')

            # No rtbkg
            los_results = los.ft_main(losx, losy, solar_flux, h_model.A, h_model.w, h_model.t, orbit_info, exopt, cdens, thermo, brad, bchi, chin, rhot, ts,
                                      rt_params_int, rt_params_real, ap_hist, thermo_bkg, zgrid)
        except ModuleNotFoundError as err:
            raise SystemExit(err)
    
    return los_results


def run_forward_los(h_model, orbit_info, hexo, solar_flux, ncore=1):
    """
    Calculates radiance by running H density model, forward source function, and los codes.

    Args:
        h_model (HDensityModel): H Density Model object.
        orbit_info (ndarray): npt x 4 array formatted for radiative transfer code. 
                              Columns are r_obs, SZA_sc, los_zenith, and los_azimuth.
        hexo (float): H density at exobase.
        solar_flux (float): Solar flux at 1216 A (line-center).
        ncore (int, optional): Number of cores to use. If > 1, run_los will run in parallel. Default is 1.

    Returns:
        los_results (ndarray): DESCRIPTION.

    """
    # 1st step - H density model - generate H density
    h_model.set_h(h_exo=hexo)  

    # 2nd step - Forward code
    exopt, cdens, thermo, brad, bchi, chin, rhot, ts = run_forward(h_model)

    # 3rd step - Line of Site code
    los_results = run_los(h_model, orbit_info, solar_flux, exopt, cdens, thermo, brad, bchi, chin, rhot,
                          ts, ncore=ncore)

    return los_results