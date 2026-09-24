#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 21 07:31:20 2025

@author: hrf
"""

from pathlib import Path
import numpy as np
import pandas as pd
from scipy import interpolate
import xarray as xr
from astropy import constants as const
from astropy import units as unit
R_EARTH = const.R_earth.to_value(unit.km)

class Thermosphere:
    """ 
    A model of the thermosphere
    Author - Pratik P. Joshi (ppjoshi2@illinois.edu), ECE Illinois; Heather Filippini (hrf@illinois.edu), CSL Illinois
    """
    
    def __init__(self, msis, mlt_saber, exo, nzthermo=61):
        """
        Initializes Thermosphere object. Currently only supports 1D_Bishop and 1D_Chamberlain models. 
        TODO: Implement 2D H density based on Vidal-Madjar [1971]

        Args:
            msis (MSIS): MSIS object.
            mlt_saber (MLT_SABER): MLT SABER object.
            exo (Exobase): Exobase object. Exobase array dimensions should match MSIS lat, lon dimensions.
            nzthermo (int, optional): Number of thermosphere altitude points. Default is 61.
            
        """
        
        # Thermospheric atomic hydrogen parameters, Fortran nomenclature
        self.flux = 2.53602354e8  # H flux from diffusion [atoms/cm2s]
        self.hmeso = 3.5e8  # 3.18043322e8
        self.zmeso = 80.0
        self.zjoint = 112.5
        
        #Set properties
        self.msis = msis
        self.exo = exo
        
        #Set thermospheric altitude grid
        self.set_z(nzthermo)
        
        #Set msis background
        self.msis.set_background(self.z)
        
        #Get parameters from mlt_saber
        self.phiH = mlt_saber.get_phiH(self.z)
        self.kzz = mlt_saber.get_k(self.z)
        
        #Set thermosphere H density
        self.set_h()
    
    @property
    def npt(self):
        return self.msis.npt
    
    @property
    def z_flat(self):
        
        return np.reshape(self.z, (self.z.shape[0], self.npt))
    
    @property
    def h_flat(self):
        
        return np.reshape(self.h, (self.z.shape[0], self.npt))
    
    @property
    def r(self):
        return (self.z + R_EARTH) * 1e5

    def set_z(self, nzthermo):
        """
        Sets z grid for thermosphere with nzthermo points using thermospheric_alt_grid.

        Args:
            nzthermo (int): Number of points in z grid.

        Returns:
            None.

        """
        
        self.z = thermospheric_alt_grid(self.exo.z_exo, nzthermo=nzthermo)
    
    def set_h(self, h_exo=None):
        """
        Sets H density profile for thermosphere.
        
        Args:
            h_exo (ndarray, optional): Updates H density at exobase if not None. Default is None.

        Returns:
            None.

        """
        
        #Update h_exo if necessary
        if h_exo is not None:
            self.exo.h_exo = h_exo
        
        # Array lengths
        nz = self.z.shape[0]
        npt = self.npt
        
        # Flatten msis arrays to be nz x npt
        flat_shape = self.z_flat.shape
        t_flat = np.reshape(self.msis.T, flat_shape)
        o_flat = np.reshape(self.msis.O, flat_shape)
        o2_flat = np.reshape(self.msis.O2, flat_shape)
        n2_flat = np.reshape(self.msis.N2, flat_shape)
        
        # dz and dT
        dz = np.ones(flat_shape)
        dz[1:] = np.diff(self.z_flat * 1e5, axis=0)
        dT = np.ones(t_flat.shape)
        dT[1:] = np.diff(t_flat, axis=0)

        # thermal diffusion constant based on Hodges [1994]
        alpha = -0.25
        
        # H scale height
        grav = 3.9898e20 / ((self.z_flat * 1e5) + 6371e5) ** 2  # g = cm/s2
        scal_p = (1.3806e-16 * t_flat) / (1.6735e-24 * grav)
        Ha = (1.3806e-16 * t_flat) / (28.964 * 1.6735e-24 * grav)

        # Molecular diffusion coefficient based on Hodges [1994]
        d_O = 4.43e17 * (t_flat ** 0.750) / o_flat
        d_O2 = 4.75e17 * (t_flat ** 0.711) / o2_flat
        d_N2 = 4.12e17 * (t_flat** 0.750) / n2_flat
        d_inv = (1.0 / d_O2) + (1.0 / d_O) + (1.0 / d_N2)
        d_mean = 1.0 / d_inv
        
        # Eddy diffusion coefficient
        d_kzz = np.reshape(self.kzz, flat_shape)

        # Set the upper boundary of [H] array to the user input [H]exo
        d_flux = np.zeros(flat_shape)
        d_flux[-1] = np.ravel(self.exo.h_exo) # nH[-1]

        # Set the H flux from user input
        flux = np.reshape(self.phiH, flat_shape)

        # Define the effective scale height
        scal_d = (1/(d_mean+d_kzz))*((d_kzz/Ha) + (d_mean/scal_p) + (((d_mean*(1.0+alpha)) + d_kzz) * (1/t_flat)*(dT/dz)))
        scal_d[0] = scal_d[1]
        
        # Interative solution of diffusion equation to calculate H density profile
        # from 74km to Zexo altitude based on Joshi [2022]
        scal_top = scal_d
        flux_top = flux / (d_mean + d_kzz)
        scal_bot = scal_d
        flux_bot = flux / (d_mean + d_kzz)
        scal_mid = (scal_top + np.roll(scal_bot, 1, axis=0)) / 2.0
        flux_mid = (flux_top + np.roll(flux_bot, 1, axis=0)) / 2.0
        
        for i in range(nz-1, 0, -1):
            # step 1
            dndz_1 = flux_top[i] + d_flux[i] * scal_top[i]
            dmid_1 = d_flux[i] + dndz_1 * dz[i] / 2.0
            # step 2
            dndz_2 = flux_mid[i] + dmid_1 * scal_mid[i]
            dmid_2 = d_flux[i] + dndz_2 * dz[i] / 2.0
            # step 3
            dndz_3 = flux_mid[i] + dmid_2 * scal_mid[i]
            dmid_3 = d_flux[i] + dndz_3 * dz[i]
            # step 4
            dndz_4 = flux_bot[i-1] + dmid_3 * scal_bot[i-1]
            # final summation
            rksum = dndz_1 + 2.0 * dndz_2 + 2.0 * dndz_3 + dndz_4

            d_flux[i - 1] = (d_flux[i] + dz[i] * rksum / 6.0)

        alt_max = self.zmeso  # 85.0
        alt_jnt = self.zjoint  # 112.5
        d_max = self.hmeso

        # Calculate mock-Chapman profile defined in Bishop [2001] for
        # altitudes below zjoint. This profile replaces the H density
        # between 74km and zjoint altitude
        
        if np.any(self.z_flat[0] < alt_jnt):
            iajnt = np.argmin(np.abs(self.z_flat - alt_jnt), axis=0)

            if npt > 1:
                tup = (iajnt, np.indices((npt,))) 
            else:
                tup = iajnt
            zjnt = self.z_flat[tup]
            djnt = d_flux[tup]

            dmock = np.zeros(flat_shape)
            # little iterative relaxation for effective Chapman scale height
            scal_eff = (zjnt - alt_max) / np.log(d_max / djnt)
            for j in range(1, 10):
                xpnt = (zjnt - alt_max) / scal_eff
                chapman = d_max * np.exp(1.0 - xpnt - np.exp(-xpnt))
                scal_eff = scal_eff - np.log(chapman / djnt)

            #   mock-Chapman profile
            xpnt = (self.z_flat - alt_max) / scal_eff
            dmock = d_max * np.exp(1.0 - xpnt - np.exp(-xpnt))

            # impose mock-Chapman profile below the "joint"
            for ii in range(npt):
                d_flux[0:iajnt[ii],ii] = dmock[0:iajnt[ii],ii]
            
            d_flux = np.reshape(d_flux, (nz,) + self.msis.lat.shape)

        # Store H density as a class attribute
        self.h = d_flux

class MLT_SABER:
    """ A class to contain the SABER data """
    
    p = Path("rt_com", "external_files")
    #FIXME: Move files to regular data_files location
    
    def __init__(self):
        """
        Extract MLT parameters on zthermo grid from kzz file.

        Args:

        Returns:
            ndarray: kzz_SCIA

        """
        #TODO: check if these arrays are always fixed
        self.z = np.arange(0, 751)

        phm = Path(self.p, "mlt_kzz.csv").resolve()

        if phm.is_file():
            df = pd.read_csv(phm)
            alt = np.asarray(df.Alt)
            sciamachy = np.asarray(df.SCIAMACHY)
        else:
            alt = np.array(
                [80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 90, 91, 92, 93, 94, 95, 96, 97, 98, 99, 100, 101, 102, 103,
                 104, 105
                 ]
            )
            sciamachy = np.array(
                [102000.00000, 260000.00000, 422000.00000, 271000.00000, 330000.00000, 414000.00000, 431000.00000,
                 454000.00000, 475000.00000, 537000.00000, 715000.00000, 720000.00000, 842000.00000, 846000.00000,
                 905000.00000, 822000.00000, 935000.00000, 928000.00000, 910000.00000, 926000.00000, 962000.00000,
                 913000.00000, 1030000.00000, 1280000.00000, 1330000.00000, 1470000.00000]
            )

        f = interpolate.interp1d(alt, sciamachy * 2, kind='cubic')  # Multiplied by 2 - Swenson 2021

        #TODO: check if these arrays are always fixed
        kzz_scia = np.empty((len(self.z)))
        kzz_scia[:] = np.nan
        kzz_scia[80:100] = f(self.z[80:100])
        z0 = self.z[100]
        kzz_scia[100:] = kzz_scia[99] * np.exp((0.09 * (self.z[100:] - 100)) - (0.14 * (self.z[100:] - z0)))
        self.kzz = kzz_scia


    def get_phiH(self, zthermo):
        """
        The Sounding of the Atmosphere using Broadband Emission Radiometry (SABER) experiment was launched onboard
        the TIMED satellite in December 2001. SABER is designed to provide measurements of the key radiative
        and chemical sources and sinks of energy in the mesosphere and lower thermosphere (MLT). SABER measures
        Earth limb emission in 10 broadband radiometer channels ranging from 1.27 micrometers to 17 micrometers.
        Measurements are made both day and night over the latitude range from 54 deg. S to 87 deg. N with alternating
        hemisphere coverage every 60 days. In this paper we concentrate on retrieved profiles of kinetic
        temperature (T(sub k)) and CO2 volume mixing ratio (vmr), inferred from SABER-observed 15 micrometer
        and 4.3 micrometer limb emissions, respectively. SABER-measured limb radiances are in non-local thermodynamic
        equilibrium (non-LTE) in the MLT region. The complexity of non-LTE radiation transfer combined with the large
        volume of data measured by SABER requires new retrieval approaches and radiative transfer techniques
        to accurately and efficiently retrieve the data products. In this paper we present the salient features
        of the coupled non-LTE T(sub k)/CO2 retrieval algorithm, along with preliminary results.
        https://ntrs.nasa.gov/citations/20030003825

        Extract MLT parameters on zthermo grid from Saber file.

        Args:
            zthermo (ndarray): zthermo grid

        Returns:
            ndarray: phiH

        """
        # Get parameters from Saber files.
        phs = Path(self.p, "mlt_saber.nc").resolve()

        if phs.is_file():
            dsS = xr.open_dataset(phs)
            Z = np.asarray(dsS.Z[0, 0, :])
            H = np.asarray(dsS.H[0, 0, :])
        else:
            # altitude grid in km
            Z = np.array([ 80.79546,  81.17079,  81.54599,  81.92102,  82.29614,  82.67102,
                83.04575,  83.42026,  83.7947 ,  84.16897,  84.54319,  84.91721,
                85.29117,  85.6649 ,  86.03857,  86.41198,  86.78561,  87.15894,
                87.53215,  87.9052 ,  88.2781 ,  88.65093,  89.02361,  89.39615,
                89.7686 ,  90.14093,  90.51315,  90.88517,  91.257  ,  91.62882,
                92.00056,  92.37192,  92.74334,  93.11461,  93.48574,  93.85675,
                94.22773,  94.59847,  94.96913,  95.3397 ,  95.71006,  96.08051,
                96.45082,  96.82075,  97.19068,  97.56046,  97.93019,  98.29961,
                98.66909,  99.03829,  99.40745,  99.77632, 100.     , 100.14529,
               101.     , 102.     , 103.     , 104.     , 105.     , 106.     ,
               107.     , 108.     , 109.     , 110.     , 111.     , 112.     ,
               113.     , 114.     , 115.     , 116.     , 117.     , 118.     ,
               119.     , 120.     , 121.     , 122.     , 123.     , 124.     ,
               125.     , 126.     , 127.     , 128.     , 129.     , 130.     ,
               131.     , 132.     , 133.     , 134.     , 135.     , 136.     ,
               137.     , 138.     , 139.     , 140.     , 141.     , 142.     ,
               143.     , 144.     , 145.     , 146.     , 147.     , 148.     ,
               149.     , 150.     , 151.     , 152.     , 153.     , 154.     ,
               155.     , 156.     , 157.     , 158.     , 159.     , 160.     ,
               161.     , 162.     , 163.     , 164.     , 165.     , 166.     ,
               167.     , 168.     , 169.     , 170.     , 171.     , 172.     ,
               173.     , 174.     , 175.     , 176.     , 177.     , 178.     ,
               179.     , 180.     , 181.     , 182.     , 183.     , 184.     ,
               185.     , 186.     , 187.     , 188.     , 189.     , 190.     ,
               191.     , 192.     , 193.     , 194.     , 195.     , 196.     ,
               197.     , 198.     , 199.     , 200.     , 201.     , 202.     ,
               203.     , 204.     , 205.     , 206.     , 207.     , 208.     ,
               209.     , 210.     , 211.     , 212.     , 213.     , 214.     ,
               215.     , 216.     , 217.     , 218.     , 219.     , 220.     ,
               221.     , 222.     , 223.     , 224.     , 225.     , 226.     ,
               227.     , 228.     , 229.     , 230.     , 231.     , 232.     ,
               233.     , 234.     , 235.     , 236.     , 237.     , 238.     ,
               239.     , 240.     , 241.     , 242.     , 243.     , 244.     ,
               245.     , 246.     , 247.     , 248.     , 249.     , 250.     ,
               251.     , 252.     , 253.     , 254.     , 255.     , 256.     ,
               257.     , 258.     , 259.     , 260.     , 261.     , 262.     ,
               263.     , 264.     , 265.     , 266.     , 267.     , 268.     ,
               269.     , 270.     , 271.     , 272.     , 273.     , 274.     ,
               275.     , 276.     , 277.     , 278.     , 279.     , 280.     ,
               281.     , 282.     , 283.     , 284.     , 285.     , 286.     ,
               287.     , 288.     , 289.     , 290.     , 291.     , 292.     ,
               293.     , 294.     , 295.     , 296.     , 297.     , 298.     ,
               299.     , 300.     ])

            # representative H flux in atoms/cm2/s
            H = np.array([3.37232100e+07, 3.46931384e+07, 3.56627309e+07, 3.66318841e+07,
                       3.76012699e+07, 3.85700354e+07, 3.95384133e+07, 4.05062228e+07,
                       4.14738513e+07, 4.24410405e+07, 4.34081004e+07, 4.43746436e+07,
                       4.53410317e+07, 4.63068254e+07, 4.72724641e+07, 4.82374309e+07,
                       4.92029662e+07, 5.01677262e+07, 5.11321762e+07, 5.20962126e+07,
                       5.30598615e+07, 5.40233294e+07, 5.49864097e+07, 5.59491283e+07,
                       5.69116142e+07, 5.86895251e+07, 6.18059112e+07, 6.49206227e+07,
                       6.80337436e+07, 7.11467806e+07, 7.42591479e+07, 7.73683337e+07,
                       8.04780218e+07, 8.35864541e+07, 8.66937142e+07, 8.97999696e+07,
                       9.29059739e+07, 9.60099687e+07, 9.91132938e+07, 1.02494858e+08,
                       1.05899845e+08, 1.09305659e+08, 1.12710186e+08, 1.16111219e+08,
                       1.19512252e+08, 1.22911907e+08, 1.26311101e+08, 1.29707446e+08,
                       1.33104342e+08, 1.36498664e+08, 1.39892618e+08, 1.43283906e+08,
                       1.45484859e+08, 1.45340358e+08, 1.46334924e+08, 1.47329491e+08,
                       1.48324057e+08, 1.49318624e+08, 1.50313190e+08, 1.51307757e+08,
                       1.52302324e+08, 1.53296890e+08, 1.54291457e+08, 1.55286023e+08,
                       1.59628100e+08, 1.63970177e+08, 1.68312254e+08, 1.72654330e+08,
                       1.76996407e+08, 1.81338484e+08, 1.85680561e+08, 1.90022638e+08,
                       1.94364715e+08, 1.98706791e+08, 2.10110421e+08, 2.21514051e+08,
                       2.32917681e+08, 2.44321310e+08, 2.55724940e+08, 2.67128570e+08,
                       2.78532200e+08, 2.89935829e+08, 3.01339459e+08, 3.12743089e+08,
                       3.22259019e+08, 3.31774948e+08, 3.41290878e+08, 3.50806808e+08,
                       3.60322737e+08, 3.69838667e+08, 3.79354597e+08, 3.88870526e+08,
                       3.98386456e+08, 4.07902386e+08, 4.12919253e+08, 4.17936121e+08,
                       4.22952988e+08, 4.27969855e+08, 4.32986723e+08, 4.38003590e+08,
                       4.43020457e+08, 4.48037325e+08, 4.53054192e+08, 4.58071059e+08,
                       4.60388279e+08, 4.62705499e+08, 4.65022718e+08, 4.67339938e+08,
                       4.69657158e+08, 4.71974378e+08, 4.74291597e+08, 4.76608817e+08,
                       4.78926037e+08, 4.81243256e+08, 4.82339341e+08, 4.83435425e+08,
                       4.84531510e+08, 4.85627594e+08, 4.86723679e+08, 4.87819763e+08,
                       4.88915848e+08, 4.90011932e+08, 4.91108017e+08, 4.92204101e+08,
                       4.92753106e+08, 4.93302110e+08, 4.93851114e+08, 4.94400118e+08,
                       4.94949122e+08, 4.95498126e+08, 4.96047130e+08, 4.96596135e+08,
                       4.97145139e+08, 4.97694143e+08, 4.97991200e+08, 4.98288258e+08,
                       4.98585315e+08, 4.98882373e+08, 4.99179431e+08, 4.99476488e+08,
                       4.99773546e+08, 5.00070603e+08, 5.00367661e+08, 5.00664718e+08,
                       5.00837583e+08, 5.01010448e+08, 5.01183313e+08, 5.01356178e+08,
                       5.01529043e+08, 5.01701909e+08, 5.01874774e+08, 5.02047639e+08,
                       5.02220504e+08, 5.02393369e+08, 5.02451121e+08, 5.02508873e+08,
                       5.02566625e+08, 5.02624377e+08, 5.02682129e+08, 5.02739881e+08,
                       5.02797633e+08, 5.02855385e+08, 5.02913137e+08, 5.02970889e+08,
                       5.03028642e+08, 5.03086394e+08, 5.03144146e+08, 5.03201898e+08,
                       5.03259650e+08, 5.03317402e+08, 5.03375154e+08, 5.03432906e+08,
                       5.03490658e+08, 5.03548410e+08, 5.03606162e+08, 5.03663914e+08,
                       5.03721666e+08, 5.03779418e+08, 5.03837170e+08, 5.03894922e+08,
                       5.03952674e+08, 5.04010427e+08, 5.04068179e+08, 5.04125931e+08,
                       5.04183683e+08, 5.04241435e+08, 5.04299187e+08, 5.04356939e+08,
                       5.04414691e+08, 5.04472443e+08, 5.04530195e+08, 5.04587947e+08,
                       5.04645699e+08, 5.04703451e+08, 5.04761203e+08, 5.04818955e+08,
                       5.04876707e+08, 5.04934460e+08, 5.04992212e+08, 5.05049964e+08,
                       5.05107716e+08, 5.05165468e+08, 5.05223220e+08, 5.05280972e+08,
                       5.05294918e+08, 5.05308864e+08, 5.05322809e+08, 5.05336755e+08,
                       5.05350701e+08, 5.05364647e+08, 5.05378593e+08, 5.05392539e+08,
                       5.05406484e+08, 5.05420430e+08, 5.05434376e+08, 5.05448322e+08,
                       5.05462268e+08, 5.05476213e+08, 5.05490159e+08, 5.05504105e+08,
                       5.05518051e+08, 5.05531997e+08, 5.05545943e+08, 5.05559888e+08,
                       5.05573834e+08, 5.05587780e+08, 5.05601726e+08, 5.05615672e+08,
                       5.05629618e+08, 5.05643563e+08, 5.05657509e+08, 5.05671455e+08,
                       5.05685401e+08, 5.05699347e+08, 5.05713293e+08, 5.05727238e+08,
                       5.05741184e+08, 5.05755130e+08, 5.05769076e+08, 5.05783022e+08,
                       5.05796968e+08, 5.05810913e+08, 5.05824859e+08, 5.05838805e+08,
                       5.05852751e+08, 5.05866697e+08, 5.05880642e+08, 5.05894588e+08,
                       5.05908534e+08, 5.05922480e+08, 5.05936426e+08, 5.05950372e+08,
                       5.05964317e+08, 5.05978263e+08])
        
        f = interpolate.interp1d(Z, H, kind='linear', fill_value='extrapolate')

        phiH = np.apply_along_axis(f, 0, zthermo)            

        return phiH


    def get_k(self, zthermo):
        """
        Extract MLT parameters on zthermo grid from Saber file.

        Args:
            kzz_scia (ndarray): kzz_scia
            zthermo (ndarray): zthermo

        Returns:
            ndarray: k_thermo

        """

        f = interpolate.interp1d(self.z[80:], self.kzz[80:], kind='cubic', fill_value='extrapolate')
        if zthermo.ndim > 1:
            k_thermo = np.apply_along_axis(f, 0, zthermo)            
        else:
            k_thermo = f(zthermo)
        
        return k_thermo
    
def thermospheric_alt_grid(zexos, nzthermo=61):
    """
    Generates thermosphere altitude grid(s) with specified number of points ranging from base (74 km) to z exo(s). Resulting altitude grid(s) will be 
    non-uniformly spaced with denser packing around qpivot (102 km).

    Args:
        zexos (ndarray): Array of altitudes at the exobase.
        nzthermo (int, optional): Number of altitude points to use

    Returns:
        ndarray: zthermo

    """
    z_shape = zexos.shape
    zexos = np.ravel(zexos)
    # Ly-alpha
    # base = 74.0
    # qpivot = 102.0
    # qfctr = 0.118418

    # Ly-beta 
    base = 102.0
    qpivot = 102.0
    qfctr = 0.075128
    quad = 0
    altarr = [base]
    alt = base
    zthermo = np.zeros((len(zexos), nzthermo))
    for i in range(len(zexos)):
        while alt < zexos[i]:
            quad = qfctr * abs(alt - qpivot)
            quad = 1 if quad < 1 else quad
            alt += quad
            altarr.append(alt)
        altarr[-1] = zexos[i]
        idx = np.where(np.asarray(altarr) > 300)[0][0]
        zthermo[i, :idx] = altarr[:idx]
        zthermo[i, idx:] = 10 ** np.linspace(np.log10(altarr[idx]), np.log10(altarr[-1]), nzthermo - idx)
    
    #Reshape
    zthermo = np.swapaxes(zthermo, 0, 1)
    zthermo = np.reshape(zthermo, (nzthermo,) + z_shape)

    return zthermo
    

        
        
