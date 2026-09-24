# !/usr/bin/env python3

from dataclasses import dataclass
import numpy as np
from astropy.time import Time
from astropy import units
from astropy import constants as const
import astropy.coordinates as acoord
import sunpy.coordinates as scoord
from pymsis import msis

from radiative_transfer.rt_common import Exobase
import astropy
R_EARTH = const.R_earth.to_value(astropy.units.km)

# warnings.filterwarnings("ignore")

@dataclass
class SolarActivity:
    """ A class to hold the solar activity parameters necessary to run MSIS """
    date: np.datetime64
    aps: np.ndarray
    f107: float
    f107a: float

class MSIS:
    """
    A class to provide access to pymsis interface and parse outputs.
    """
    
    def __init__(self, solar_activity, lat, lon, do_grid=False):
        """
        Initializes MSIS.

        Args:
            solar_activity (SolarActivity): Solar activity object holding date, aps, f10.7, and f10.7a.
            lat (ndarray): GSE latitude array to calculate MSIS background at.
            lon (ndarray): GSE longitude array to calculate MSIS background at.
            do_grid (bool, optional): If True, lat and lon arrays will be run through meshgrid to create a 2D lat, lon grid. Default is False.

        Returns:
            None.

        """
        
        if do_grid:
            lon, lat = np.meshgrid(lon, lat)
        
        self.solar_activity = solar_activity
        self.lat = np.atleast_1d(lat)
        self.lon = np.atleast_1d(lon) 
        self.msis_version = 2.1
        self.mod_msis = 1 #EXPLAIN: never used
        
        #Make sure date is in datetime64 format
        self.solar_activity.date = np.datetime64(self.solar_activity.date)
        
        #Error check
        if do_grid==False:
            assert self.lat.shape == self.lon.shape, "Latitude and longitude arrays must be same size."
        
        #Error check aps size - force into 1x7 array
        if isinstance(self.aps, int) or isinstance(self.aps, float) or self.aps.size == 1:
            self.aps = self.aps * np.ones((1,7))
        elif self.aps.shape[0] != 1:
            self.aps = np.reshape(self.aps, (1,7))
    
    @property
    def date(self):
        return self.solar_activity.date
    
    @property
    def year(self):
        msis_year = self.date.astype('datetime64[Y]')
        year_str = np.datetime_as_string(msis_year)
        return int(year_str)
    
    @property
    def uday(self):
        #Calculate day of year
        msis_year = self.date.astype('datetime64[Y]')
        msis_doy = self.date.astype('datetime64[D]') - msis_year + 1
        return msis_doy.astype(int)
    
    @property
    def ut(self):
        #EXPLAIN: Calculate msis ut time - where did this come from? Is it used in Fortran??
        msis_hour = self.date.astype('datetime64[h]') - self.date.astype('datetime64[D]')
        msis_minute = self.date.astype('datetime64[m]') - self.date.astype('datetime64[h]')
        return msis_hour.astype(float) + msis_minute.astype(float) * 0.01666667
    
    @property
    def f107(self):
        return self.solar_activity.f107
    
    @property
    def f107a(self):
        return self.solar_activity.f107a
    
    @property
    def aps(self):
        return self.solar_activity.aps
    
    @aps.setter
    def aps(self, aps_val):
        self.solar_activity.aps = aps_val
    
    @property
    def geo_grid(self):
        # gse_to_geo
        GSE_frame = scoord.frames.GeocentricSolarEcliptic
        tan_pts_GSE = acoord.SkyCoord(lon=self.lon * units.deg,
                                      lat=self.lat * units.deg,
                                      distance=R_EARTH * units.km,
                                      obstime=Time(self.date, scale="utc"),
                                      frame=GSE_frame)
        ITRS_frame = acoord.ITRS(obstime=Time(self.date, scale="utc"))
        tan_pts_ITRS = tan_pts_GSE.transform_to(ITRS_frame)
        tan_pts_GEO = acoord.EarthLocation.from_geocentric(tan_pts_ITRS.x,
                                                           tan_pts_ITRS.y,
                                                           tan_pts_ITRS.z).to_geodetic()

        return [tan_pts_GEO.lat.deg, tan_pts_GEO.lon.deg]
    
    @property
    def geo_grid_flat(self):
        
        lat_geo, lon_geo = self.geo_grid
        
        return [np.ravel(lat_geo), np.ravel(lon_geo)]
    
    @property
    def npt(self):
        return self.lat.size

    def set_background(self, z):
        """
        Calculates MSIS background at provided altitudes for each lat, lon point. MSIS return values are parsed
        and saved as properties. Densities are converted from [atoms/m^3] to [atoms/cm^3].

        Args:
            z (ndarray): Altitude array (km).

        Returns:
            nmsis_bkg (ndarray): Raw MSIS return value array. See pymsis for details.

        """
        
        z = np.atleast_1d(z)
        z_shape = z.shape
        nz = z_shape[0]
        
        if len(z_shape) > 1:
            assert np.prod(z_shape) // nz ==self.npt, 'z must either be 1D array or nz x nlat x nlon.'
            z = np.reshape(z, (nz, self.npt))
        elif len(z_shape) == 1:
            nz = z.size
            z = np.expand_dims(z, 1) * np.ones((nz, self.npt))
        else:
            raise ValueError('z must either be 1D array or nz x nlat x nlon.')
        
        #Call msis model
        lats, lons = self.geo_grid_flat
        nmsis_bkg = np.zeros((11, nz, self.npt))
        for ii in range(self.npt):
            nmsis = msis.run(self.date, lons[ii], lats[ii], z[:, ii], self.f107, self.f107a, self.aps,
                             version=self.msis_version)
            nmsis_bkg[:,:,ii] = np.squeeze(nmsis).T
        
        #Reshape everything
        gse_shape = self.lat.shape
        nmsis_bkg = np.reshape(nmsis_bkg, (11, nz) + gse_shape)
        
        #Parse data & convert densities from [atoms/m^3] to [atoms/cm^3]
        self.z = np.reshape(z, (nz,) + gse_shape)
        self.total_dens = nmsis_bkg[0] * 1e-6
        self.N2 = nmsis_bkg[1] * 1e-6
        self.O2 = nmsis_bkg[2] * 1e-6
        self.O = nmsis_bkg[3] * 1e-6
        self.He = nmsis_bkg[4] * 1e-6
        self.H = nmsis_bkg[5] * 1e-6
        self.Ar = nmsis_bkg[6] * 1e-6
        self.N = nmsis_bkg[7] * 1e-6
        self.O_anom = nmsis_bkg[8] * 1e-6
        self.NO = nmsis_bkg[9] * 1e-6
        self.T = nmsis_bkg[10]
        
        return nmsis_bkg
    
    def exobase(self, zfixed=np.arange(200, 800, 1)):
        """
        Calculate MSIS values at the exobase. Returns exobase altitude (km), H density, and temperature.

        Args:
            zfixed (ndarray, optional): Altitude array to use to search for exobase. Defaults to np.arange(200, 800, 1).

        Returns:
            Exobase: Dataclass holding h_exo, z_exo, and t_exo results.

        """
        
        #Calculate background
        self.set_background(zfixed)

        # oxygen constants, cross-section and mass
        sigma_c = 3.5e-15 # cm2, 0.35e-19m2
        M_O = 0.016 # kg mol-1

        # # corresponds to the 6.88e+06 - 6.38e+06m = 0.50e+6 (~500km above Earth surface = exosphere)
        # g = const.G.value * const.M_earth.value / 6.88e+06**2  # 8.421  ms-2
        g = 8.44  # m s-2
        # R = 8.314  # kg m2 mol-1 K-1 s-2
        dens_sum = self.N2 + self.O2 + self.O + self.He + self.H + self.Ar + self.N + self.O_anom + self.NO
        Ha = (const.R.value * self.T * 1e2) / (M_O * g)  # cm
        prod = sigma_c * Ha * dens_sum
        prod = np.reshape(prod, (zfixed.size, self.npt))
        
        idx = np.argmin(np.abs(prod-1), axis=0)
        
        #Extract z_exo, h_exo, and t_exo
        z_exo = zfixed[idx]
        inds = (idx, np.arange(0, self.npt, 1, dtype=int))
        H_flat = np.reshape(self.H, (zfixed.size, self.npt))
        T_flat = np.reshape(self.T, (zfixed.size, self.npt))
        h_exo = H_flat[inds]
        t_exo = T_flat[inds]
        
        #Reshape everything
        z_exo = np.reshape(z_exo, self.lat.shape)
        h_exo = np.reshape(h_exo, self.lat.shape)
        t_exo = np.reshape(t_exo, self.lat.shape)
        
        return Exobase(h_exo, z_exo, t_exo, self.lat, self.lon)
        

if __name__=="__main__":
    
    sa = SolarActivity(np.datetime64('2001-08-15T06:00:00'), 3 * np.ones((1, 7)), 210, 210)
    
    msis_obj = MSIS(sa, 0, 90)
    
    
    exos = msis_obj.exobase()
    
    gselat_wedge = np.array([-7.46304341, -22.30589489, -37.11437768, -51.83217264, -66.27426059,-79.2435324,-79.25898342,-66.29455618,-51.85177357,-37.13214555,
                    -22.32123879, -7.47562649, 7.37761445, 22.21705515, 37.0131819, 51.70176491, 66.07105525, 78.79961508, 78.79975748, 66.07195768, 51.70411914,
                    37.01758453, 22.22396608, 7.38731662])

    gselon_wedge = np.array([84.44887256, 86.68556306, 89.4410583, 93.53490112, 101.61767089, 129.24482077, -141.9249623, -114.22957886, -106.13598766, -102.0388846,
                    -99.28201398, -97.04462436, -94.93672137, -92.65662081, -89.79841002, -85.49346307, -76.94464324, -48.51717889, 35.92240194, 64.34983947, 72.89860143,
                    77.20345208, 80.06150946, 82.34136751])
    
    msis_obj = MSIS(sa, gselat_wedge, gselon_wedge)
    
    
    exos = msis_obj.exobase()

    msis_obj = MSIS(sa, np.arange(-90, 90, 10), np.arange(-180, 180, 10), do_grid=True)

    exos = msis_obj.exobase()
