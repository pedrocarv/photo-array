import sys
from pathlib import Path
import numpy as np
import pandas as pd
from astropy.time import Time
import warnings
# warnings.filterwarnings("ignore")

import glide.science.radiative_transfer.rt_common as rt_common
import glide.science.radiative_transfer.rt_msis as rt_msis
from glide.common_components.cam import nadir_nfi_mode, nadir_wfi_mode, CamSpec, CamMode
from glide.common_components.constants import LAMBDA_GRID
from glide.science.radiative_transfer.rt_inversion import inversion_init, optimize_profile
from glide.common_components.science_pixel_binning import SciencePixelBinning
from glide.common_components.spacecraft import SpaceCraft
from glide.common_components.camera import Camera
from glide.validation.instrument import Instrument
from glide.validation.scene import Scene
import glide.calibration.calibration_helpers as ch
from glide.validation.cam import load_lab_data
import glide.calibration.oob as oob_algo
import glide.science.radiative_transfer.rt_thermosphere as rt_thermosphere
import glide.science.radiative_transfer.rt_h_dens as rt_h_dens


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

def science_pixel_bin(scraft, radiance, bin_type="frac",
                      inner_r_boundary=1.078, nrad=6, dth=15, thlim=(-180, 180)):
    """
    Science binning for the HEXO part.

    Args:
        scraft (Spacecraft): A spacecraft object. Make sure a spacecraft is properly initialized with a channel in cam
                             (cam=CameraNFI()) and sensor (sensors=[cam])
        radiance (ndarray): I_exo or noisy image.
        channel (str): The channel to use (NFI or WFI)
        bin_type (str): The type of binning to use on the edges (frac or full)
        inner_r_boundary (float): Inner radial boundary.
        nrad (int): Number of radial bins.
        dth (int): Angular bin size (deg).
        thlim (tup): Angular limits to bin. (inner_azimuth_boundary, outer_azimuth_boundary),
                     Default is -180 to 180 deg.

    Returns:
        ndarray: Binned images. Shape (nr, nth)
        ndarray: Tangent point gse coordinates. Shape  (2, nth)
        ndarray: Solar zenith angles. Shape (nth)
        ndarray: Bin centers. Shape (2, (nr * nth))

    """
    re_pix = scraft.re_to_pixels(scraft.sensor.camID)  # number of pixels per Re

    # Set binning parameters
    dr = 0.05 * re_pix # Radial bin size
    ntheta = int((thlim[1] - thlim[0]) / dth)  # number of Angular (theta) bins, 24

    # Radial limits
    rlim_min = inner_r_boundary * re_pix
    rlim_max = inner_r_boundary * re_pix + nrad * dr
    # Radial limits with dr shift
    #rlim_min = inner_r_boundary * re_pix + (1 * dr)
    #rlim_max = inner_r_boundary * re_pix + nrad * dr + (1 * dr)

    rlim = (rlim_min, rlim_max)  # Radial limits
    sci_pix_bin = SciencePixelBinning(scraft.sensor.npix, shape=(nrad, ntheta), rlim=rlim, bin_type=bin_type,
                                      thlim=thlim)

    image_bin = sci_pix_bin(radiance) # image_bin = np.zeros((6, 24))

    # Calculate tpgse coordinates
    tpgse_lat, tpgse_lon, sza = scraft.wedge_tan_pts(sci_pix_bin, camID=scraft.sensor.camID, rlim=rlim, thlim=thlim)

    sza = np.reshape(sza, (nrad, ntheta))
    tpgse_lat = np.reshape(tpgse_lat, (nrad, ntheta))
    tpgse_lon = np.reshape(tpgse_lon, (nrad, ntheta))
    szas = sza[0, :]
    tp_gses = np.array([tpgse_lat[0, :], tpgse_lon[0, :]])  # (2, ntheta)

    # Calculate bin centers
    bin_ctrs = sci_pix_bin.get_bin_ctrs(rlim=rlim, thlim=thlim)  # 2 x (nrad * ntheta) (2, 144)
    return image_bin, tp_gses, szas, bin_ctrs

def noisy_images(scraft, radiance,
                 iph, oob, stars, moon, blur, distortion,
                 t_op, radiation, activation, lifetime, num_images=1):
    """
    Add noise to the images. No calibration.

    Args:
        scraft (Spacecraft): A spacecraft object.
        radiance (ndarray): I_exo or noisy image.
        iph (bool): IPH mean event rate is generated.
        oob (bool): Out of band mean event rate is generated.
        stars (str or None): Use netCDF file containing stellar spectral data.
        moon (bool): Moon event rate is generated..
        blur (bool): Blurs scene with Voigt kernel using cam_specs.psf_funct.
        distortion (bool): Distorts all incoming photons.
        t_op (int): t_op
        radiation (str): Radiation
        activation (str): Activation baseline
        lifetime (str): BOL beginning of life of a detector due to contamination etc. EOL end of life.
        num_images (int): Number of noisy images. Default = 1

    Returns:
        ndarray: Binned images. Shape (nr, nth)
        obj: Instrument

    """
    # noisy images
    cam_mode = nadir_nfi_mode(filter='LyaN', t_op=t_op)
    cam_spec = CamSpec(cam_mode)
    cam_spec = load_lab_data(cam_spec, radiation=radiation, activation=activation, lifetime=lifetime)
    # cam_spec.set_fw_temp()  # ? AttributeError: 'CamSpec' object has no attribute 'fw_temp'
    cam = Camera(cam_mode, cam_spec)
    scene = Scene(scraft, iph=iph, oob=oob, stars=stars, I_exo=radiance,
                  moon=moon, blur=blur, distortion=distortion)
    # Instrument model
    instr = Instrument(cam, scene)
    nadir_images = instr(num_images=num_images)

    return nadir_images, instr

def create_profile(msis_date, channel, orbit_info, tpgse_coord, hexo, solar_flux, solar_activity, parallel, igeo=0, model='Chamberlain'):
    """ Calculates initial radiance and background files by running the Forward and LOS models.

        Args:
            msis_date (np.datetime64): MSIS date.
            channel (str): Type of data, instrument channel; 1x, NFI, WFI
            orbit_info (ndarray): Orbit info data.
            tpgse_coord (list): tpgse coordinates.
            hexo (float): A Hydrogen density.
            solar_flux (float): Solar flux.
            aps (float): aps.
            f107 (float): F107.
            f107a (float): F107a.
            parallel (bool): LOS code run in parallel with Pool module. Make sure max_os and import los match.
                             Default False

        Return:
            ndarray: Lyao source file with radiance and ratio.

    """

    #Set up H density model
    msis = rt_msis.MSIS(solar_activity, tpgse_coord[0], tpgse_coord[1])
    exo = msis.exobase()
    mlt_saber = rt_thermosphere.MLT_SABER()
    thm_obj = rt_thermosphere.Thermosphere(msis, mlt_saber, exo)
    
    if model == 'Bishop':
        hdm = rt_h_dens.HDensityBishop(thm_obj)
    elif model == 'Chamberlain':
        hdm = rt_h_dens.HDensityChamberlain(thm_obj)
    else:
        raise ValueError('Select Bishop or Chamberlain')

    los_results = rt_common.run_forward_los(hdm, orbit_info, hexo, solar_flux, ncore=8)
    return los_results


if __name__ == "__main__":
    import time
    start_time = time.time()

    # 2027-08-09T06:00:00 1e5 7e11 NFI C
    date_s = '2027-08-09T06:00:00'
    hexo_s = 1e5
    sf_s = 7e11
    channel = 'NFI'
    model = 'C'

    igeo = {"C": 0, "B": 1}
    model_name = {"C": "Chamberlain", "B": "Bishop"}

    date = date_s.replace("-", "").replace(":", "")
    date = date.replace("T", "").replace("Z", "")

    hexo = None
    solar_flux = None
    try:
        hexo = float(hexo_s)
    except ValueError as err:
        print(err)
    try:
        solar_flux = float(sf_s)
    except ValueError as err:
        print(err)

    p = Path(__file__).parent
    pc = Path(p, "rt_com", "results")
    pc.mkdir(mode=0o777, exist_ok=True)

    current_date= np.datetime64(date_s)
    msis_hist = True
    msis_date, aps, f107, f107a = geomagnetic_data(msis_hist, current_date, solar_flux)

    solar_activity = rt_msis.SolarActivity(msis_date, aps, f107, f107a)

    print(current_date, msis_date, hexo, solar_flux, channel, model)
    print(solar_activity)

    # t_op = 60
    t_op = 90
    # Set up SpaceCraft
    cam_mode = None
    if channel == 'NFI':
        cam_mode = nadir_nfi_mode(t_op=t_op) # default t_op=30
    elif channel == 'WFI':
        cam_mode = nadir_wfi_mode() # default t_op=60
    cam_spec = CamSpec(cam_mode)
    cam_spec = load_lab_data(cam_spec)
    # cam_spec.set_fw_temp()  # ? AttributeError: 'CamSpec' object has no attribute 'fw_temp'
    cam = Camera(cam_mode, cam_spec)
    # u,v coordinates for all pixels in focal plane. Coordinates returned: [2 x npix * npix]
    scraft = SpaceCraft(Time(current_date, scale="utc"), sensors=[cam])
    uv = cam.get_pixels()

    # Set up Science Pixel Binning
    re_pix = scraft.re_to_pixels(scraft.sensor.camID)  # number of pixels per Re
    inner_r_boundary=1.078 
    nrad=6 
    dth=15 
    thlim=(-180, 180)
    bin_type='frac'
    # Set binning parameters
    dr = 0.05 * re_pix # Radial bin size
    ntheta = int((thlim[1] - thlim[0]) / dth)  # number of Angular (theta) bins, 24

    # Radial limits
    rlim_min = inner_r_boundary * re_pix
    rlim_max = inner_r_boundary * re_pix + nrad * dr

    rlim = (rlim_min, rlim_max)  # Radial limits
    sci_pix_bin = SciencePixelBinning(scraft.sensor.npix, shape=(nrad, ntheta), rlim=rlim, bin_type=bin_type,
                                      thlim=thlim)
    
    tpgse_lat, tpgse_lon, sza = scraft.wedge_tan_pts(sci_pix_bin,rlim=rlim,thlim=thlim,camID=channel)
    tpgse_lat = np.reshape(tpgse_lat, (nrad, ntheta))
    tpgse_lon = np.reshape(tpgse_lon, (nrad, ntheta))
    
    tp_gses = np.array([tpgse_lat[0, :], tpgse_lon[0, :]])  # (2, ntheta)

    # tpgse_coord_fwd = [0, 90]
    tpgse_coord_fwd = tp_gses[:,0]


    # version = date + "_" + channel + "_" + hexo_s + "_" + sf_s + "_" + model
    # rad_name = "radiance_{0}.txt".format(version)
        
    # Calculate orbit info for the NFI radiance
    # orbit_info = r_obs, sza_sc, los_zenith, los_azimuth, los_tp_r
    orbit_info = scraft.calc_rt_angles(uv, camID=channel,orbit_info=True)

    los_results = create_profile(msis_date, channel, orbit_info, tpgse_coord_fwd,
                                    hexo, solar_flux, solar_activity,
                                    False, igeo[model])
    

    if los_results.all():
        radiance = los_results[:, 5]
        # np.savetxt(Path(pc, rad_name), radiance, delimiter=",")
    else:
        print(los_results)
        exit(1)

    radiance = radiance.reshape((cam.npix, cam.npix))
    print("min:", np.nanmin(radiance), "max:", np.nanmax(radiance))

    import ipdb
    ipdb.set_trace()

    profiles = 1 # number of profiles
    p_offset = 0 # starting profile number - added for saving txt file

    # 0 = only radiance is binned and inverted, >0 = noisy images are calculated and inverted
    res_all_profiles = []
    for profile in range(profiles + 1):
        if profile == 0:
            # Science binning
            # image_bin, tp_gses, szas, bin_ctrs = science_pixel_bin_hexo(scraft, radiance)
            image_bin, tp_gses, szas, bin_ctrs = science_pixel_bin(scraft, radiance)
        else:
            # Noisy image profile
            radiation = "S3"
            activation = "lot"  # baseline
            lifetime = "BOL"  # BOL beginning of life of a detector due to contamination etc. EOL end of life.

            iph = True
            oob = True
            stars = None
            moon = False
            blur = True  # True
            distortion = False

            # Generate noisy images - we don't use stack, number of images = 1
            nadir_image, instr = noisy_images(scraft, radiance, t_op=t_op,
                                               iph=iph, oob=oob, stars=stars, moon=moon, blur=blur,
                                               distortion=distortion, activation=activation, radiation=radiation,
                                               lifetime=lifetime)

            # Out of band (OOB) testing, default = 1 (oob_full false), or use full pre-calculated or calculated file
            oob_calc_full = False
            if oob_calc_full:
                if blur:
                    oob_name = "oob_{0}.txt".format(version + "_blur")
                else:
                    oob_name = "oob_{0}.txt".format(version + "_no_blur")

                oob_file = True
                if oob_file:
                    oob_bkgd = np.loadtxt(Path(pc, oob_name))
                    oob_bkgd = oob_bkgd.reshape(cam.npix, cam.npix)
                else:
                    # calculate oob_bkgd, skipped for now
                    oob_bkgd = oob_full(current_date, radiance,
                                        iph=iph, oob=oob, stars=stars, moon=moon, blur=blur, distortion=distortion,
                                        activation=activation, radiation=radiation, lifetime=lifetime)
                    np.savetxt(Path(pc, oob_name), oob_bkgd.reshape(-1), delimiter=",")
                # print("OOB nans", np.count_nonzero(np.isnan(oob_bkgd)))
                # print("OOB not nans", np.count_nonzero(~np.isnan(oob_bkgd)))

                # Calibrate noisy images without oob and abs_intensity_scale
                nadir_image = instr.calibrate_images(nadir_image, do_distortion_correction=distortion,
                                                      oob_scale=None, abs_intensity_scale=None)

                if blur:
                    # Deblur noisy images
                    psf = cam.spec.psf_func(res=1)
                    nadir_image = ch.deblur(nadir_image[0], psf)

                # Subtract OOB array from an image, oob in digital number/s
                nadir_image_oob = ch.subtract_mean_oob(nadir_image, oob_bkgd,
                                                       instr.cam.mode.t_int, instr.cam.spec.mask_fov)
                # Scale to absolute intensity
                nadir_image = ch.absolute_intensity(nadir_image_oob, instr.cam.spec.lambda_grid,
                                                     instr.cam.spec.R_diffuse, instr.cam.mode.t_int)
            else:
                # Calibrate noisy images
                oob_scale = 1.0
                nadir_image = instr.calibrate_images(nadir_image, do_distortion_correction=distortion,
                                                      oob_scale=oob_scale, abs_intensity_scale=1)

            # Science binning
            image_bin, tp_gses, szas, bin_ctrs = science_pixel_bin(scraft, nadir_image)