#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 28 20:15:52 2025

@author: hrf
"""

import numpy as np
import scipy as sp
from radiative_transfer.rt_msis import MSIS
from astropy import constants as const
from astropy import units as unit
R_EARTH = const.R_earth.to_value(unit.km)


class HDensityModel:
    """ 
    A H density model class which constitutes the part-1 of the 3-part RT code.
    Author - Pratik P. Joshi (ppjoshi2@illinois.edu), ECE Illinois, H. Filippini (hrf@illinois.edu)
    """
    GG = 0.8862269  # Boltzmann (OR \"Barometric\") factor [unitless]
    GM = 6.6743e-8  # Gravitational constant [cm^3 / g s^2]
    HMASS = 1.6738e-24  # mass of atomic hydrogen [g] 1.67×10−27 kg
    MP = 5.972e27  # Mass of Earth [g]
    BOLTZ = 1.3806e-16  # Boltzmann constant [cm^2.g.s^-2.K^-1]
    kB = 1.38064852e-23  # Boltzmann constant [m^2.kg.s^-2.K^-1]
    PLANETR = R_EARTH * 1e5 # Mean planet radius [cm]
    TSOL = 1.10e6  # Solar ionization lifetime [s]

    # Gaussian quadrature points and weights
    GQ16 = np.asarray([[2.71524601e-02, -9.89400923e-01], [6.22535236e-02, -9.44575012e-01],
                       [9.51585099e-02, -8.65631223e-01], [1.24628969e-01, -7.55404413e-01],
                       [1.49595991e-01, -6.17876232e-01], [1.69156522e-01, -4.58016783e-01],
                       [1.82603419e-01, -2.81603545e-01], [1.89450607e-01, -9.50125083e-02],
                       [1.89450607e-01, 9.50125083e-02], [1.82603419e-01, 2.81603545e-01],
                       [1.69156522e-01, 4.58016783e-01], [1.49595991e-01, 6.17876232e-01],
                       [1.24628969e-01, 7.55404413e-01], [9.51585099e-02, 8.65631223e-01],
                       [6.22535236e-02, 9.44575012e-01], [2.71524601e-02, 9.89400923e-01]])
    
    # Default z_grid with 32 altitudes in km above the Earth surface.
    Z_GRID = np.array([
        74.5, 76.0, 78.5, 82.0, 86.5, 92.0, 98.5, 106.0, 116.5, 133.5,
        161.0, 205.5, 277.5, 394.0, 582.5, 887.5, 1381.0, 2179.5, 3471.5, 5562.0,
        8945.0, 14417.0, 23273.0, 37601.0, 60785.0, 98299.0, 158999.0, 257199.0, 416109.0, 673219.0,
        1089229.0, 1762329.0
    ])
    
    #Default 76-point albedo z-grid
    Z_ALBEDO = np.array([74.5, 75.5, 77.5, 80.5, 84.5, 89.5, 95.5, 102.5, 111.5, 125.5,
                          147.5, 183.5, 242.0, 336.0, 488.5, 735.5, 1134.5, 1780.5, 2826.0, 4517.0,
                          7253.0, 10843.0, 13073.0, 13753.0, 14469.0, 15222.0, 16013.0, 16847.0, 17723.0, 18645.0,
                          19615.0, 20635.0, 21709.0, 22838.0, 24026.0, 25275.0, 26591.0, 27974.0, 29429.0, 30960.0,
                          32571.0, 34265.0, 36047.0, 37923.0, 39895.0, 41970.0, 44153.0, 46451.0, 48867.0, 51409.0,
                          54083.0, 56897.0, 59857.0, 62971.0, 66247.0, 69693.0, 73317.0, 77131.0, 81144.0, 85365.0,
                          89806.0, 94479.0, 99389.0, 104559.0, 109999.0, 115729.0, 121749.0, 128079.0, 134739.0, 141749.0,
                          149119.0, 205049.0, 336649.0, 544669.0, 881229.0, 1425829.0])
    
    SOLAR_ANGLE = np.array([0.0000000000000000,     
        3.6439273268664917e-3,
        0.36722528005866900,     
        0.52232635511499836,     
        0.64350996086882573,     
        0.74759151082735187,     
        0.84107460873999973,     
        0.92730019731419144,     
        1.0082642681179304,     
        1.0852817075808143,     
        1.1592823782645494,     
        1.2309617646146440,     
        1.3008653679044733,     
        1.3694397612041413,     
        1.4370656305727307,     
        1.5040806219779146,     
        1.5707963267948970,     
        1.6375120316118785,     
        1.7045270230170624,     
        1.7721528923856518,     
        1.8407272856853198,     
        1.9106308889751491,     
        1.9823102753252437,     
        2.0563109460089786,     
        2.1333283854718630,     
        2.2142924562756017,     
        2.3005180448497935,     
        2.3940011427624412,     
        2.4980826927209674,     
        2.6192662984747948,     
        2.7743673735311241,     
        3.1379487262629264]) # Solar angle centroid grid default for fortran source function
    
    
    def __init__(self, thermosphere, solar_ang=SOLAR_ANGLE, A=1, w=1, t=0):
        """
        Initialize H density model with Thermosphere instance.

        Args:
            thermosphere (Thermosphere): Thermosphere object.

        Returns:
            None.

        """
        
        self.thermosphere = thermosphere
        self.solar_ang = solar_ang
        self.A = A # Acts as a toggle for turning on the varying exobase
        self.w = w # Frequency of sinusoid... should never be changed from one. 
        self.t = t # NOTE: Offset from -pi/2
        # Retrieval would be for the phase shift and hexo
    
    @property
    def z_exo(self):
        return self.thermosphere.exo.z_exo
    
    @property
    def t_exo(self):
        return self.thermosphere.exo.t_exo
    
    @property
    def h_exo(self):
        return self.thermosphere.exo.h_exo
    
    @property
    def msis(self):
        return self.thermosphere.msis
    
    @property
    def npt(self):
        return self.thermosphere.npt
    
    @property
    def r(self):
        """ Radial distance from Earth in cm """
        return (self.z + R_EARTH) * 1e5
    
    def set_h(self, z=Z_GRID, h_exo=None):
        """
        Sets H density for input z grid.

        Args:
            z (ndarray, optional): Altitude array (km). Defaults to Z_GRID.
            h_exo (ndarray, optional): Updates H density at exobase if not None. Default is None.

        Returns:
            None.

        """
        
        #Update h_exo if necessary
        if h_exo is not None:
            self.thermosphere.set_h(h_exo)
        
        self.z = z
        
        #Find exobase crossing point in z_grid
        zexo = self.thermosphere.exo.z_flat
        exo_ind = np.zeros(zexo.size, dtype=int)
        
        for ii in range(exo_ind.size):
            exo_ind[ii] = np.argmin(np.abs(z - zexo[ii]), axis=0)
        
        #Calculate h density in exosphere
        h_dens = self.get_h(self.r)
        h_shape = (self.r.size, ) + self.z_exo.shape
        h_dens = np.reshape(h_dens, (self.r.size, self.npt))

        #Fill H density with thermospheric density up to exobase
        #TODO: replace interp1d with modern interpolator
        for ii in range(self.npt):
            f = sp.interpolate.interp1d(self.thermosphere.z_flat[:,ii], self.thermosphere.h_flat[:,ii], kind='cubic', fill_value='extrapolate')
            h_dens[0:exo_ind[ii],ii] = f(z[0:exo_ind[ii]])
                

        #Reshape     
        h_dens = np.reshape(h_dens, h_shape)

        
        self.h = h_dens
    
    def to_fortran_input(self, ncopy=1):
        """
        Packages H density model attributes into Fortran radiative transfer code inputs. Array sizes much match Fortan code expectations.

        Args:
            ncopy (int, optional): Copies 1D arrays ncopy times. Defaults to 1. Must correspond to NPT variable in global_parameters.f.

        Returns:
            rt_bkg (ndarray): H, O2, and T arrays over z grid stacked into array with dimensions [3 x ncopy x nz].
            thermo_bkg (ndarray): Thermosphere background arrays (zthermo, rthermo, T, hthermo, O2, N2, O) stacked into array with dimensions [7 x ncopy x nzthermo].
            rt_params_int (ndarray): Integer parameters: MSIS day-of-year, MSIS year, IAPH, MOD_MSIS, IGEO.
            rt_params_real (ndarray): Float parameters: MSIS UT, geodetic latitude, geodetic longitude, f10.7, f10.7a, h_exo, adjt, scaln, base, top, z_exo, flux, hmeso, tsat, hsat, zmeso.

        """
        
        #FIXME: Remove unused variables from interface
        
        #Error check
        # assert self.npt==1, 'Fortran code only currently supports 1D input.'
        
        #Build rt_bkg: 3 x ncopy x nz
        #h, nO2, T
        #FIXME: rt_bkg is parsed in Fortran, but never used
        zmsis = MSIS(self.msis.solar_activity, self.msis.lat, self.msis.lon)
        zmsis.set_background(self.z)
        rt_bkg = np.dstack((self.h, zmsis.O2, zmsis.T))
        rt_bkg = np.swapaxes(rt_bkg, 0, 2)
        rt_bkg = rt_bkg * np.ones((rt_bkg.shape[0], ncopy, rt_bkg.shape[2]))
        
        #Build thermo_bkg: 7 x ncopy x nzthermo
        #zthermo (km), rthermo (cm), T, hthermo, nO2, nN2, nO
        zthermo = self.thermosphere.z
        rthermo = (self.thermosphere.z + 6371.0) * 1e5 #FIXME: Fortran code uses 6371 as Earth radius. Code breaks if astropy constant is used.
        Tn = self.msis.T
        hthermo = self.thermosphere.h
        nO2 = self.msis.O2
        nN2 = self.msis.N2
        nO = self.msis.O
        thermo_bkg = np.dstack((zthermo, rthermo, Tn, hthermo, nO2, nN2, nO))
        thermo_bkg = np.swapaxes(thermo_bkg, 0, 2)
        thermo_bkg = thermo_bkg * np.ones((thermo_bkg.shape[0], ncopy, thermo_bkg.shape[2]))
        
        #Integer RT parameters
        iaph = 0 #EXPLAIN: What is this? Parsed in Fortran but not used
        rt_params_int = [self.msis.uday, self.msis.year, iaph, self.msis.mod_msis, self.igeo]
        rt_params_int = np.array(rt_params_int, dtype=int, order="F")
        

        #Real RT parameters
        d_max = self.thermosphere.hmeso
        alt_max = self.thermosphere.zmeso
        adjt = 1 #EXPLAIN: What is this? Not used in Fortran
        scaln = 1 #EXPLAIN: What is this? Not used in Fortran
        base = self.thermosphere.z[0,0]
        top = self.z_exo[0] #EXPLAIN: This is redundant to z_exo and not used on the Fortran side
        lat_geo, lon_geo = self.msis.geo_grid
        
        if hasattr(self, 'tsat'):
            satt = self.tsat
        else:
            satt = 750.0
        if hasattr(self, 'hsat'):
            satd = self.hsat
        else:    
            satd = 7e5
        
        #msis_p2 - real parameters
        rt_params_real = [self.msis.ut, lat_geo[0], lon_geo[0], self.msis.f107, self.msis.f107a, self.h_exo[0], adjt,
             scaln, base, top, self.z_exo[0], self.thermosphere.flux, d_max, satt, satd, alt_max]
        rt_params_real = np.array(rt_params_real, dtype=float, order="F")
        
        
        return rt_bkg, thermo_bkg, rt_params_int, rt_params_real
        

class HDensityBishop(HDensityModel):
    """ A class to define the Bishop H density model"""
    
    def __init__(self, thermosphere, tsat=750.0, hsat=7e5, **kwargs):
        
        super().__init__(thermosphere, **kwargs)
        
        self.tsat = tsat
        self.hsat = hsat
        self.igeo = 1 #Fortran code switch for H density model type
    
    def get_h(self, rr):
        """
        Gets Bishop H density model defined over input radial grid.

        Args:
            rr (ndarray): Radial grid.

        Returns:
            h_dens (ndarray): H density array.

        """
        
        top = self.thermosphere.exo.z_flat
        texo = self.thermosphere.exo.t_flat
        dexo = self.thermosphere.exo.h_flat
        tsat = self.tsat
        dsat = self.hsat
        f107 = self.msis.f107

        rtpi = np.sqrt(np.pi)
        # Generic offset value [unitless]
        offset = 1e-6
        
        h_dens = np.zeros((rr.size, self.npt))
        
        #Calculate density over all lat, lon points
        for ii in range(self.npt):
        
            # Calculated parameters
            # Radial distance of the exobase [cm]
            rc = self.PLANETR + top[ii] * 1e5
            # 3.579, Ly-a line-center flux used in calculating RP [/cm^2 s]
            radpf = 2.91 * (1.0 + 0.002 * (f107 - 65.0))
            # 800*25066846121.282, Radial distance of the exopause [cm]
            rp = np.sqrt(self.GM * self.MP / 0.1774e0 / radpf)
    
            #Set satellite H and T to exobase values if dsat < 0
            if (dsat < 0):
                tsat = texo[ii]
                dsat = dexo[ii]
    
            # 0.906, Scaling factor for satellite Exobase temperature [unitless]
            ftsat = tsat / texo[ii]
            # 4.211, Scaling factor for satellite Exobase density [unitless]
            fdsat = dsat / dexo[ii]
            # intermediate parameter [cm K]
            conl = self.GM * self.HMASS * self.MP / self.BOLTZ
            ww, xx = [self.GQ16[:, 0], self.GQ16[:, 1]]
    
            cac = conl / (rc * texo[ii])
            car = conl / (rp * texo[ii])
            # cas = conl / (ftsat * texo)
            caa = conl / (rr * texo[ii])
            
            
            #Condition: caa >= (cac * (1 - offset))
            ind_1 = caa >= (cac * (1 - offset))
            h_dens[ind_1,ii] = dexo[ii] * np.exp(caa[ind_1] - cac)
            
            #Condition: caa > car and caa < (cac * (1 - offset))
            ind_2 = np.logical_and(caa > car, caa < (cac * (1 - offset)))
            y1 = (caa[ind_2] ** 2) / (cac + caa[ind_2])
            ya = caa[ind_2] - car
            yb = caa[ind_2] - car * cac / (cac + car)
            ye = (caa[ind_2] ** 2) / (car + caa[ind_2])
            ga = gamma16(ya, xx, ww)
            gb = gamma16((yb - y1), xx, ww)
            gs = gamma16((yb - y1) / ftsat, xx, ww)
            za = zamma16((ye - ya), xx, ww)
            zb = zamma16((ye - yb), xx, ww)
            zs = zamma16((ye - yb) / ftsat, xx, ww)
            acoef = 2 * np.sqrt(caa[ind_2] ** 2 - car ** 2) * np.exp(-ye) / car / rtpi
            ccoef = 2 * np.sqrt(cac ** 2 - caa[ind_2] ** 2) * np.exp(-y1) / cac / rtpi
            balst = 2 * ga / rtpi + acoef * (za - zb) - ccoef * gb
            escap = (self.GG - ga) / rtpi - acoef * (za - zb) / 2 - ccoef * (self.GG - gb) / 2
            part1 = (balst + escap) * np.exp(caa[ind_2] - cac)
            satel = acoef * zs * np.exp((ftsat - 1) * ye / ftsat) + ccoef * gs * np.exp((ftsat - 1) * y1 / ftsat)
            dk = 2.97e6 * np.arcsin(1 - rc / rr[ind_2]) / np.sqrt(rr[ind_2]) / radpf
            satel = satel * np.exp(-dk / self.TSOL)
            part2 = fdsat * satel * np.exp((caa[ind_2] - cac) / ftsat)
            h_dens[ind_2,ii] = dexo[ii] * (part1 + part2)
        
        #Reshape
        h_shape = (rr.size, ) + self.z_exo.shape
        h_dens = np.reshape(h_dens, h_shape)
        
        return h_dens

class HDensityChamberlain(HDensityModel):
    
    """ A class to define the Chamberlain H density model"""
    
    def __init__(self, thermosphere, tsat=750.0, **kwargs):
        
        super().__init__(thermosphere, **kwargs)
        
        self.tsat = tsat
        self.igeo = 0 #Fortran code switch for H density model type
    
    def get_h(self, rr):
        """
        Gets Chamberlain H density model defined over input radial grid.

        Args:
            rr (ndarray): Radial grid.

        Returns:
            h_dens (ndarray): H density array.

        """
        
        top = self.thermosphere.exo.z_flat
        texo = self.thermosphere.exo.t_flat
        dexo = self.thermosphere.exo.h_flat
        tsat = self.tsat

        rtpi = np.sqrt(np.pi)
        # Generic offset value [unitless]
        offset = 1e-6
        
        h_dens = np.zeros((rr.size, self.npt))

        ww, xx = [self.GQ16[:, 0], self.GQ16[:, 1]]
        
        #Calculate density over all lat, lon points
        for ii in range(self.npt):
            
            # Calculated parameters
            # Radial distance of the exobase [cm]
            rc = self.PLANETR + top[ii] * 1e5
            # 2*6371e5, Satellite critical radial distance [cm]
            ftsat = tsat * self.PLANETR
            # 837119542228.0166, intermediate parameter [cm K]
            conl = self.GM * self.HMASS * self.MP / self.BOLTZ

            lamb_c = conl / (rc * texo[ii]) # Absolute potential energy at exobase level [k Texo]
            lamb_sat = conl / (ftsat * texo[ii]) # Absolute potential energy at satellite critical level [k Texo]
            lamb = conl / (rr * texo[ii]) # Absolute potential energy at radial bin ctr [k Texo]
            
            # Condition: lamb >= (lamb_c * (1 - offset))
            # Condition: For radial bins below the exobase
            ind_1 = lamb >= (lamb_c * (1 - offset))
            h_dens[ind_1,ii] = dexo[ii] * np.exp(lamb[ind_1] - lamb_c)
            
    
            # Condition: lamb > lamb_s and lamb < (lamb_c * (1 - offset))
            # Condition: For radial bins above exobase and below satellite critical level
            ind_2 = np.logical_and(lamb > lamb_sat, lamb < (lamb_c * (1 - offset)))
            ya = (lamb[ind_2] ** 2) / (lamb_c + lamb[ind_2])
            ga = gamma16(lamb[ind_2], xx, ww)
            gc = gamma16((lamb[ind_2] - ya), xx, ww)
            ccoef = np.sqrt(lamb_c ** 2 - lamb[ind_2] ** 2) * np.exp(-ya) / lamb_c
            balst = (ga - ccoef * gc) * 2 / np.sqrt(np.pi)
            escap = (1 - ccoef) * self.GG / rtpi - balst / 2
            bound = 2 * ga / rtpi
            part = bound + escap
            h_dens[ind_2,ii] = part * dexo[ii] * np.exp(lamb[ind_2] - lamb_c)
            
            # Condition: lamb > 0 and lamb <= lamb_s
            # Condition: Radial bins above satellite critical level and less than infinity
            ind_3 = np.logical_and(lamb > 0, lamb <= lamb_sat)
            ya = (lamb[ind_3] ** 2) / (lamb_c + lamb[ind_3])
            ga = gamma16(lamb[ind_3], xx, ww)
            gc = gamma16((lamb[ind_3] - ya), xx, ww)
            ccoef = np.sqrt(lamb_c ** 2 - lamb[ind_3] ** 2) * np.exp(-ya) / lamb_c
            balst = (ga - ccoef * gc) * 2 / np.sqrt(np.pi)
            escap = (1 - ccoef) * self.GG / rtpi - balst / 2
            bound = 2 * ga / rtpi
            ys = (lamb[ind_3] ** 2) / (lamb_sat + lamb[ind_3])
            gs = gamma16(lamb[ind_3] - ys, xx, ww)
            scoef = np.sqrt(lamb_sat ** 2 - lamb[ind_3] ** 2) * np.exp(-ys) / lamb_sat
            satel = (ga - scoef * gs) * 2 / np.sqrt(np.pi) - balst
            part = balst + satel + escap
            h_dens[ind_3,ii] = part * dexo[ii] * np.exp(lamb[ind_3] - lamb_c)
        
        #Reshape
        h_shape = (rr.size, ) + self.z_exo.shape
        h_dens = np.reshape(h_dens, h_shape)

        return h_dens

class HDensity2DLCM(HDensityChamberlain):
    ''' 
    Describes 2D H-density model described by Vidal-Madjar and Bertaux '72, approximation using 'local' chamberlain models
    
    Density varies radially and azimuthally.

    NOTE: For now, jsut for visualization... not yet the version for forward modeling
    NOTE: Passed thermosphere must be defined at 0, 90, and 180 longitude GSE
    '''
    # NOTE: alpha/azimuthal grid 
    ALPHA_GRID = np.linspace(0, 180, 20)

    def get_h_2d(self, rr):
        '''
        Returns 2D H-density grid
        NOTE: Keeping get_h available to compare the spherically symmetric dist (w/ varying exosphere)
        to the 2D distribution
        '''
        top = self.thermosphere.exo.z_flat
        texo = self.thermosphere.exo.t_flat
        dexo = self.thermosphere.exo.h_flat
        tsat = self.tsat

        rtpi = np.sqrt(np.pi)
        # Generic offset value [unitless]
        offset = 1e-6
        
        h_dens = np.zeros((rr.size, self.npt))

        ww, xx = [self.GQ16[:, 0], self.GQ16[:, 1]]
        
        #Calculate density over all lat, lon points
        for ii in range(self.npt):
            # Calculated parameters
            # Radial distance of the exobase [cm]
            rc = self.PLANETR + top[ii] * 1e5
            # 2*6371e5, Satellite critical radial distance [cm]
            ftsat = tsat * self.PLANETR
            # 837119542228.0166, intermediate parameter [cm K]
            conl = self.GM * self.HMASS * self.MP / self.BOLTZ

            lamb_c = conl / (rc * texo[ii]) # Absolute potential energy at exobase level [k Texo]
            lamb_sat = conl / (ftsat * texo[ii]) # Absolute potential energy at satellite critical level [k Texo]
            lamb = conl / (rr * texo[ii]) # Absolute potential energy at radial bin ctr [k Texo]
            
            # Condition: lamb >= (lamb_c * (1 - offset))
            # Condition: For radial bins below the exobase
            ind_1 = lamb >= (lamb_c * (1 - offset))
            h_dens[ind_1,ii] = dexo[ii] * np.exp(lamb[ind_1] - lamb_c)
            
    
            # Condition: lamb > lamb_s and lamb < (lamb_c * (1 - offset))
            # Condition: For radial bins above exobase and below satellite critical level
            ind_2 = np.logical_and(lamb > lamb_sat, lamb < (lamb_c * (1 - offset)))
            ya = (lamb[ind_2] ** 2) / (lamb_c + lamb[ind_2])
            ga = gamma16(lamb[ind_2], xx, ww)
            gc = gamma16((lamb[ind_2] - ya), xx, ww)
            ccoef = np.sqrt(lamb_c ** 2 - lamb[ind_2] ** 2) * np.exp(-ya) / lamb_c
            balst = (ga - ccoef * gc) * 2 / np.sqrt(np.pi)
            escap = (1 - ccoef) * self.GG / rtpi - balst / 2
            bound = 2 * ga / rtpi
            part = bound + escap
            h_dens[ind_2,ii] = part * dexo[ii] * np.exp(lamb[ind_2] - lamb_c)
            
            # Condition: lamb > 0 and lamb <= lamb_s
            # Condition: Radial bins above satellite critical level and less than infinity
            ind_3 = np.logical_and(lamb > 0, lamb <= lamb_sat)
            ya = (lamb[ind_3] ** 2) / (lamb_c + lamb[ind_3])
            ga = gamma16(lamb[ind_3], xx, ww)
            gc = gamma16((lamb[ind_3] - ya), xx, ww)
            ccoef = np.sqrt(lamb_c ** 2 - lamb[ind_3] ** 2) * np.exp(-ya) / lamb_c
            balst = (ga - ccoef * gc) * 2 / np.sqrt(np.pi)
            escap = (1 - ccoef) * self.GG / rtpi - balst / 2
            bound = 2 * ga / rtpi
            ys = (lamb[ind_3] ** 2) / (lamb_sat + lamb[ind_3])
            gs = gamma16(lamb[ind_3] - ys, xx, ww)
            scoef = np.sqrt(lamb_sat ** 2 - lamb[ind_3] ** 2) * np.exp(-ys) / lamb_sat
            satel = (ga - scoef * gs) * 2 / np.sqrt(np.pi) - balst
            part = balst + satel + escap
            h_dens[ind_3,ii] = part * dexo[ii] * np.exp(lamb[ind_3] - lamb_c)

        # After chamberlain Distribution is computed, apply 2Dimensionality
        # < 1.5Re: same as regular chamberlain distribution
        # 1.5Re < z < 2Re: All densities equal to mean chamberlain dist
        # > 2Re: N = 0.95(N_mean + (N_mean - N)/3.5)
        hden_mean = np.mean(h_dens, axis=1)
        eq_dens = np.logical_and(rr > 1.5*self.PLANETR, rr <= 2*self.PLANETR)
        up_dens = rr > 2*self.PLANETR
        h_dens[eq_dens, :] = np.repeat(np.expand_dims(hden_mean[eq_dens], axis=1), self.npt, axis=1)
        up_mean = np.repeat(np.expand_dims(hden_mean[up_dens], axis=1), self.npt, axis=1)
        h_dens[up_dens, :] = 0.95 * (up_mean + (up_mean - h_dens[up_dens, :])/3.5)

        #Reshape
        h_shape = (rr.size, ) + self.z_exo.shape
        h_dens = np.reshape(h_dens, h_shape)

        return h_dens
    

# Incomplete Gamma function
def gamma16(aa, xx, ww):
    if np.any(aa) < 0:
        raise ValueError("Input must be non-negative.")

    root = np.sqrt(aa)
    erf = 0.0
    for j in range(0, 16):
        xpnt = ((1 + xx[j]) * root / 2.0) ** 2
        erf += ww[j] * np.exp(-xpnt) * root / 2.0
    gamma = erf - root * np.exp(-aa)
    return gamma


def zamma16(aa, xx, ww):
    zamma = 0.0
    for j in range(0, 16):
        xpnt = (1 + xx[j]) * aa / 2.0
        zamma += ww[j] * np.exp(xpnt) * np.sqrt(xpnt) * aa / 2.0
    return zamma
