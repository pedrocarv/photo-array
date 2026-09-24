'''
Class defining photometer observation model
'''
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit, minimize
from scipy.integrate import quad
from importlib import resources

import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, GCRS, ITRS
from astropy.coordinates import SkyCoord
import sunpy.coordinates as scoord
from astroplan import Observer
from astropy.time import Time, TimeDelta

from scip.data_processing.helpers import *

# Physical constants
mh = 1.6738e-24 # Atomic Mass of Hydrogen [g]
k = 1.3806e-16 # Boltzmann constant [erg/K]
c = 2.9979e10 # Speed of light [cm/s]

class Photometer:
    '''
    (FIXME: Only modeling fm2b003004 & 10x10 pixel binning)
    Loads instrument parameters, applies realistic instrument noise, computes expected SNR
    Should store pointing information?
    '''
    def __init__(self, loc=get_loc_info('aeronomy')):
        # Define Instrument parameters
        self.A = np.pi * (2.5)**2 # Area of aperature [cm] (objective diameter = 50mm)
        self.QE = 0.65 # Manufacturer defined QE [e-/phot] NOTE: Eventually want to measure this on-sky
        self.G = 0.27 # manufacturer defined gain [e-/adu]
        self.blocking = 0.69 # CPI-Reported blocking filter transmissivity Near H-alpha and off-band
        self.atms = 0.95 # Atmospheric transmissivity to Halpha (approx)

        datadir = resources.files('scip') / 'calibration/data'

        # Load noise data & fov info 
        on_dframe= np.load(datadir / 'dark_fm2b003004_202510/onband_-6_10.npy')
        off_dframe= np.load(datadir / 'dark_fm2b003004_202510/offband_-6_10.npy')

        m,n = on_dframe.shape
        self.N_pix = m*n

        fov_y = 2.287247518
        fov_x = 2.860205313
        self.omega = (fov_y/m)*(fov_x/n) * (np.pi ** 2) / (180**2) # FIXME? Square degrees to steradians?

        # Define measured instrument noise
        self.d_on = self.G * np.mean(on_dframe) / (30*60) # Per-pixel dark noise [e-/s/pix]
        self.b_on = np.load(datadir / 'read_noise/on_read_noise_-6_10.npy')[5] # Per-pixel read noise [e-/pix]
        self.d_off = self.G * np.mean(off_dframe) / (30*60) # Per-pixel dark noise [e-/s/pix]
        self.b_off = np.load(datadir / 'read_noise/off_read_noise_-6_10.npy')[5] # Per-pixel read noise [e-/pix]
        
        # Load transmissivity & define continuous frequency profile
        df_on = pd.read_csv(datadir / 'transmissivity_curves/fm2b004_normalized_tcurve.csv')
        wl_on = df_on['wavelength'].to_numpy()
        f_t_on = df_on['transmissivity'].to_numpy()

        freq_on = wl_on * (1e-8 / c)

        def func(x, a, x0, sigma):
            return a * np.exp(-(x - x0)**2 / (2 * sigma**2))

        params_on, _ = curve_fit(func, freq_on, f_t_on, p0=[0.35, 6563*(1e-8/c), 1.5 * (1e-8/c)])

        self.f_on = lambda x : func(x, *params_on) # Continuous filter transmissivity (wrt wavelength)

        df_off = pd.read_csv(datadir / 'transmissivity_curves/fm2b003_normalized_tcurve.csv')
        wl_off = df_off['wavelength'].to_numpy()
        f_t_off = df_off['transmissivity'].to_numpy()

        freq_off = wl_off * (1e-8 / c)

        params_off, _ = curve_fit(func, freq_off, f_t_off, p0=[0.35, 6589*(1e-8/c), 1.5 * (1e-8/c)])

        self.f_off = lambda x : func(x, *params_off) # Continuous filter transmissivity (wrt wavelength)

        filter_transmissions = np.load(datadir / 'filter_transmission/total_filter_transmission.npy')
        self.t_ratio = filter_transmissions[3]/filter_transmissions[2]

        self.loc = loc


    def __call__(self, rad, bkg, t_int, t_exo):
        '''
        For now assume that a constant t_int must be some integer number of seconds
        :input:
        rad (ndarray) - Balmer-Alpha radiance profile [R]
        bkg (float) - Background emission brightness [R]
        t_int (int) - integration time [s]
        t_exo (float) - exobase temperature [K]
        '''
        if type(t_int) == int:
            t_int = t_int*np.ones(len(rad))

        S_out = []
        for i in range(len(rad)):
            def doppler(v):
                '''
                Doppler line profile centered at 6562.79 angstrom (H-alpha) wrt LINEAR FREQUENCY
                NOTE: Assumes doppler broadening much larger than natural/pressure/and collisional broadening
                See http://ftp.astro.wisc.edu/~townsend/resource/teaching/astro-310-F08/16-line-profiles-1.pdf
                https://www.slideserve.com/keilah/lecture-9-inhomogeneous-broadening-the-laser-equation-and-threshold-gain
                '''
                v0 = (6562.79 * 1e-8) / c # Line center frequency [Hz]

                dop_width = (v0 / c) * np.sqrt((2*k*t_exo) / mh)

                return rad[i] * (1/(np.sqrt(2*np.pi)*dop_width)) * np.exp(-1 * ((v-v0)**2)/(2 * dop_width**2)) # Area under curve normalized to I 
            
            on_range = [(1e-8/c)*6556, (1e-8/c)*6570] # Assume center at 6563
            off_range = [(1e-8/c)*6582, (1e-8/c)*6596] # Assume center at 6589

            from scipy.integrate import quad
            def on(x):
                return self.f_on(x) * (doppler(x) + (bkg / (on_range[1]-on_range[0])))

            i_on, _ = quad(on, *on_range)

            def off(x):
                return (bkg / (off_range[1] - off_range[0])) * self.f_off(x)
            i_off, _ = quad(off, *off_range)
            
            # Steps 1&2: Defining distribution of signal
            p_on = np.sum((self.blocking * self.atms * 1e6 * self.A * self.omega * i_on) / (4*np.pi)) # Average photon rate [photon/s]
            P_on = self.QE * np.random.poisson(p_on*t_int[i], size=self.N_pix)

            p_off = np.sum((self.blocking * self.atms * 1e6 * self.A * self.omega * i_off) / (4*np.pi)) # Average photon rate [photon/s]
            P_off = self.QE * np.random.poisson(p_off*t_int[i], size=self.N_pix)

            # 3: Define arrays containing other sources of noise
            D_on = np.random.poisson(self.d_on*t_int[i], size=self.N_pix)
            R_on = np.random.normal(0, self.b_on, size=self.N_pix)

            D_off = np.random.poisson(self.d_off*t_int[i], size=self.N_pix)
            R_off = np.random.normal(0, self.b_off, size=self.N_pix)

            # 4: Profit
            S_on = (1/self.G) * np.sum(P_on + D_on + R_on)
            S_off = (1/self.G) * np.sum(P_off + D_off + R_off)

            # Now dark subtract like usual, leaving only on-band signal and some noise
            S_sub = (S_on - (self.d_on/self.G)*t_int[i]*self.N_pix) - self.t_ratio*(S_off - (self.d_off/self.G)*t_int[i]*self.N_pix)
            S_out += [S_sub / t_int[i]] # Dark subtracted (neglect bias... deterministic) & time-averaged... comparable to ER Plots from observations

        return S_out


    def snr(self, rad, t_int, rad_bkg=6):
        '''
        Get photometer SNR based on scene radiance, assumes 6 R Off-Band
        Based on Dawn's formulation of SNR
        '''
        # NOTE: For now, account for 0.95 atmospheric attenuation
        p_on = 0.95*self.blocking*((1e6 * self.A * self.omega * (rad+rad_bkg)) / (4*np.pi)) # Average photon rate [photon/s]
        p_off = 0.95*self.blocking*((1e6 * self.A * self.omega * rad_bkg) / (4*np.pi)) # Average photon rate [photon/s]
        # Mean of poisson <-> Variance of poisson

        # Check box for derivations of these
        snr_on = ((self.N_pix * self.QE * p_on * t_int)/ np.sqrt(self.N_pix * (self.QE**2 * p_on *t_int + self.d_on * t_int + 2*self.b_on**2)))
        snr_off = ((self.N_pix * self.QE * p_off * t_int)/ np.sqrt(self.N_pix * (self.QE**2 * p_off *t_int + self.d_off * t_int + 2*self.b_off**2)))
        snr_sub = ((self.N_pix * self.QE * p_on * t_int) - self.t_ratio*(self.N_pix * self.QE * p_off * t_int)) / np.sqrt(self.QE**2 * p_on *t_int + self.d_on * t_int + 2*self.b_on**2 + (self.t_ratio**2) * (self.N_pix * (self.QE**2 * p_off *t_int + self.d_off * t_int + 2*self.b_off**2)))

        return snr_on, snr_off, snr_sub

    
    def antisolar_pointing(self, dates):
        '''
        Return altitude/azimuth for antisolar pointing for the photometer's location at a given time
        '''
        sol_alt = np.array([get_solalt(dates[i], self.loc, 'deg') for i in range(len(dates))])
        sol_az = np.array([get_solaz(dates[i], self.loc, 'deg') for i in range(len(dates))])
        
        alts = -sol_alt
        azs = (sol_az + 180) % 360 

        return alts, azs


    def sunward_pointing(self, dates):
        '''
        Return alt/az for sunward pointing for a given loc/time
        NOTE: Optimization procedure...beware of funky results
        '''
        alts, azs = [], []

        for d in dates:
            aloc = EarthLocation(lat=self.loc[1]*u.deg, lon=self.loc[0]*u.deg, height=(self.loc[2] * 1e-3)*u.km)
            atime = Time(d, scale='utc')
            gse_frame = scoord.frames.GeocentricSolarEcliptic(obstime=atime)
            
            # Observer location in GSE
            obs_itrs = aloc.get_itrs(obstime=atime) # Location in ITRS
            # Transform to GSE
            obs_gse = obs_itrs.transform_to(gse_frame)
            obs_gse.representation_type = 'cartesian'
            r_obs = obs_gse.cartesian  # ~6371 km offset from Earth center

            # Optimization routine to find sunward pointing
            def shadow_perpendicular_standoff(angles, r_obs_km, atime, aloc, Re=6371):
                alt_deg, az_deg = angles
                altaz_frame = AltAz(obstime=atime, location=aloc)
                
                try:
                    los = SkyCoord(alt=alt_deg*u.deg, az=az_deg*u.deg,
                                distance=1*u.au, frame=altaz_frame)
                    gse = los.transform_to(scoord.GeocentricSolarEcliptic(obstime=atime))
                    v_hat = gse.cartesian.xyz.to(u.km).value
                    v_hat = v_hat / np.linalg.norm(v_hat)
                except:
                    return 1e9

                ry, rz = r_obs_km[1], r_obs_km[2]
                vy, vz = v_hat[1], v_hat[2]

                # GOAL: Minimize magnitude of vector between Obs location and shadow crossing
                # i.e. minimize magnitude s (should be positive and real)
                a = vy**2 + vz**2
                b = 2*(ry*vy + rz*vz)
                c = ry**2 + rz**2 - (Re * u.km)**2

                disc = b**2 - 4*a*c

                if disc < 0:
                    return 9e9 # Must have real magnitude
                
                s1 = (-b + np.sqrt(disc)) / (2*a)
                s2 = (-b - np.sqrt(disc)) / (2*a)

                # Return the smallest crossing magnitude
                if s1 > 0 and s2 > 0:
                    return np.min([s1.value, s2.value])
                else:
                    return np.max([s1.value,s2.value])


            result = minimize(
                shadow_perpendicular_standoff,
                x0=[90, 0], # Start near-zenith looking north
                args=(r_obs.xyz, atime, aloc),
                bounds=[(0, 90), (0, 360)],
                method='L-BFGS-B',
                options={'ftol': 1e-8}
            )

            solar_alt, solar_az = result.x
            alts += [solar_alt]
            azs += [solar_az]
        
        return alts, azs
    
    def standard_pointing(self, pointing, dates):
        '''
        Return pointing for a 'standard pointing' scheme (i.e. a pointing scheme dependent on shadow geometry)
        '''
        alts, azs = [], []
        match pointing:
            case 'zenith':
                alts = [90] * len(dates)
                azs = [180] * len(dates)
            case 'antisolar':
                alts, azs = self.antisolar_pointing(dates)
            case 'sunward':
                alts, azs = self.sunward_pointing(dates)
            case _:
                raise NotImplementedError('Yell at JC')
            
        return alts, azs



    def schedule_obs(self, start, pointing, t_int, profile='full', end=None):
        '''
        Get observation times and angles for a desired night and pointing scheme

        FIXME: Needs further functionality, such as...
        - Ability to accept arbitrary pointing schemes (for observation optimization)
        - Accounting for image download and mount rotation time (For accurate start times)
        - Variable integration time
        - Calibration images

        NOTE: datetime64 should automatically be in UTC

        Args:
            start (np.datetime64) - Start day of observations
            pointing (str) - Pointing Scheme to use
            t_int (int) - Integration time of each observation
            profile (str) - type of profile to schedule for ('full', 'dusk', 'dawn')

        Output: 
            (ndarray) - Array containing obs info for each point, with each line structured as [time, alt, az]
        '''
        obs = Observer(timezone='UTC', longitude=self.loc[0]*u.deg, latitude=self.loc[1]*u.deg, elevation=self.loc[2]*u.m)

        if profile != 'subinterval':
            start = (start.astype('datetime64[D]') + np.timedelta64(12, 'h')).astype('datetime64[s]')
            astart = Time(start)

            dusk = obs.twilight_evening_astronomical(astart, which='next')
            dawn = obs.twilight_morning_astronomical(astart, which='next')
            mid = obs.midnight(astart, which='next')

            dusk = dusk.to_value('datetime64')
            dawn = dawn.to_value('datetime64')

        match profile:
            case 'full':
                current = dusk
                end = dawn
            case 'dusk':
                current = dusk
                end = mid
            case 'dawn':
                current = mid
                end = dawn
            case 'subinterval':
                # Case if you only want to schedule for a subset of the night (i.e. only an hour)
                assert end is not None, 'Must specify end time'
                current = start
                end = end
            case _:
                raise NotImplementedError('Please choose either dawn, dusk, or a full night of observations')


        obs_times = [current]

        # FIXME: The following loop needs to be much more robust
        if pointing != 'custom':
            while current < end:
                current += np.timedelta64(t_int + 5, 's') # Allow 5 seconds for image download
                obs_times += [current]

            alts, azs = self.standard_pointing(pointing, obs_times)
        else: 
            raise NotImplementedError('Yell at JC')
        

        # return list(zip(obs_times, alts, azs))
        return pd.DataFrame(
            {
                'observation_time':obs_times,
                'altitude':alts,
                'azimuth':azs
            }
        )

    
    def rt_los(self, sched):
        '''
        Return vector containing necessary LOS info for lyao_rt for a given observation schedule
        FIXME: Move 2 Scene?
        '''
        dates = sched['observation_time'].to_numpy()
        obs = len(dates)
        alts = sched['altitude'].to_numpy()
        azs = sched['azimuth'].to_numpy()

        szas = [get_sza(d, loc=self.loc, unit='deg') for d in dates]
        sol_azs = [get_solaz(d, loc=self.loc, unit='deg') for d in dates]

        los = []
        for i in range(obs):
            # [distance from earth center (km), Solar (zenith) Angle, Off-Zenith Angle, LOS Azimuth]

            azi = (azs[i] - sol_azs[i]) % 360
            if azi > 180:
                azi = 360 - azi
            azi += 0.180

            los += [[6.371224e3, szas[i], 90 - alts[i] + 0.180, azi]] # FIXME: Check if LOS azimuth determination is correct

        return los
        

        
