'''
Perform sensitivity tests on forward model
JC

NOTE: Since we are still working with a spherically symmetric exosphere model, only need to model one limb profile
'''

import numpy as np
import pandas as pd
from astropy.time import Time
import matplotlib.pyplot as plt

import glide.science.radiative_transfer.rt_common as rt_common
import glide.science.radiative_transfer.rt_msis as rt_msis
import glide.science.radiative_transfer.rt_thermosphere as rt_thermosphere
import glide.science.radiative_transfer.rt_h_dens as rt_h_dens
from glide.common_components.science_pixel_binning import SciencePixelBinning
from glide.common_components.spacecraft import SpaceCraft
from glide.common_components.camera import CameraNFI

def fwd_los(hexo, tpgse_coord, orbit_info, solar_activity, solar_flux, model, tscale=1):
        '''
        Generate a [nr x 1] radiance profile for a given hexo and tpgse coord
        '''
        msis = rt_msis.MSIS(solar_activity, tpgse_coord[0], tpgse_coord[1])
        exo = msis.exobase(tscale=tscale)
        mlt_saber = rt_thermosphere.MLT_SABER()
        thm_obj = rt_thermosphere.Thermosphere(msis, mlt_saber, exo, tscale=tscale)
        
        if model == 'Bishop':
            hdm = rt_h_dens.HDensityBishop(thm_obj)
        elif model == 'Chamberlain':
            hdm = rt_h_dens.HDensityChamberlain(thm_obj)
        else:
            raise ValueError('Select Bishop or Chamberlain')
        
        los_results = rt_common.run_forward_los(hdm, orbit_info, float(hexo), solar_flux)

        return los_results[:,5] # Return radiance from forward model

def geomagnetic_data(msis_hist, current_date, solar_flux):
    """ Get geomagnetic data based on solar flux intensity.

        Args:
            msis_hist (bool): Use historical data.
            current_date (np.datetime64): A current, spacecraft date.
            solar_flux (float): Solar flux.

        Return:
            np.datetime64: msis_date
            np.ndarray: aps
            float: f107
            float: f107a

    """
    aps = 3
    year = 2002  # solar medium conditions, solar flux = 5e11
    f107 = 150
    f107a = 150
    # MSIS historical date - use the same month and day as in date, hour should be either dusk or dawn.
    # For the 36 ground truth generation we should do the following:
    if np.abs((solar_flux - 3e11) / 3e11) < .05:  # math.isclose(solar_flux, 3e11):
        if msis_hist:
            year = 2008  # solar minimum conditions, solar flux = 3e11
        f107 = 80
        f107a = 80
    elif np.abs((solar_flux - 5e11) / 5e11) < .05:  # math.isclose(solar_flux, 5e11):
        if msis_hist:
            year = 2002  # solar medium conditions, solar flux = 5e11
        f107 = 150
        f107a = 150
    elif np.abs((solar_flux - 7e11) / 7e11) < .05:  # math.isclose(solar_flux, 7e11):
        if msis_hist:
            year = 2001  # solar maximum conditions, solar flux = 7e11
        f107 = 210
        f107a = 210

    geomag_date = np.datetime64(str(year))

    msis_date = current_date
    if msis_hist:
        if current_date > np.datetime64('now') and current_date > geomag_date:
            # Pandas to_date is the most reliable way
            dt = pd.to_datetime(current_date)
            gdt = pd.to_datetime(geomag_date)
            dt_new = dt.replace(year=gdt.year)
            msis_date = np.datetime64(dt_new)

    return msis_date, aps, f107, f107a

# For now, focus on varying over hexo or solar flux
hexo = 1e5
sf = 7e11
date_s = '2027-08-09T06:00:00'
model = 'Bishop'
channel = 'NFI'
cam = CameraNFI()
# cam_spec = load_lab_data(cam_spec)

tpgse = [0, 90] # Tan. Pt. in GSE for dawn 

current_date= np.datetime64(date_s)
msis_hist = True
msis_date, aps, f107, f107a = geomagnetic_data(msis_hist, current_date, sf)

solar_activity = rt_msis.SolarActivity(msis_date, aps, f107, f107a)

#generate all tpgses from the binning scheme, then choose the one that I want
nth = 24
scraft = SpaceCraft(Time(current_date, scale="utc"), sensors=[cam])
#Get tangent points
cam = scraft.sensors[channel]
re_pix = scraft.re_to_pixels(camID=channel) #Number of pixels per Re
dr     = .05*re_pix #radial bin size
nr     = 6
rmin   = 1.128
rlim   = (rmin * re_pix, rmin * re_pix + nr*dr) #radial limits
sci_pix_bin = SciencePixelBinning(cam.npix, (nr, nth), rlim = rlim)
gse_lat, gse_lon, sza = scraft.wedge_tan_pts(sci_pix_bin, rlim=(rlim[0], rlim[0]+dr))

#Construct orbit_info, tpgse
ctr = cam.npix/2
r = np.arange(0, nr)
orbit_infos = []
tpgses = []
for ii in range(nth):
    #Get viewing geometry for each radial line
    u = r * np.cos(sci_pix_bin.theta_ctrs[ii]) + ctr
    v = -r * np.sin(sci_pix_bin.theta_ctrs[ii]) + ctr
    uv = np.vstack((u,v))  
    orbit_infos.append(scraft.calc_rt_angles(uv=uv, orbit_info=True))
    tpgses.append([gse_lat[ii], gse_lon[ii]])

orbit_info = orbit_infos[0]
tpgse = tpgses[0]

bins = rmin  + 0.05 * r

plt.clf()
# for ts in [0.85, 0.9, 0.95, 1, 1.05, 1.1, 1.15]:
for ts in [4e4, 6e4, 8e4, 10e4, 12e4]:
    # act = rt_msis.SolarActivity(msis_date, aps, f, f)

    rad = fwd_los(hexo, tpgse, orbit_info, solar_activity, sf, model)

    norm = rad / np.mean(rad)

    plt.plot(bins, rad, linewidth=1.5, label='%.2f'%ts)

plt.legend()
plt.grid()
plt.ylabel('Radiance [10^6 Photons/cm^2/s]')
plt.xlabel('Distance from Earth Surface [Re]')
plt.title('Model Sensitivity to Exobase Temperature Bias')

plt.savefig('rads_temp.png')


