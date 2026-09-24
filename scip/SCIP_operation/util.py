# -*- coding: utf-8 -*-
"""
Created on Mon Mar 25 20:26:07 2019

@author: photometer

Utilites for photometer
"""

import matplotlib.pyplot as plt
import numpy as np
import os
#import geomag
#from datetime import date
from datetime import datetime, timezone
#from astropy.coordinates import EarthLocation,SkyCoord
#from astropy.time import Time
#from astropy import units as u
#from astropy.coordinates import AltAz
from astropy.io import fits
from pathlib import Path
#from photutils import DAOStarFinder
from astropy.stats import mad_std
from photutils.aperture import aperture_photometry, CircularAperture
import configparser

from astropy.time import Time, TimeDelta
import astropy.units as u
from astropy.coordinates import SkyCoord, FK4
from astroplan import FixedTarget

def save_as_fits(array, camera, exp_starttime, exptime, binning, b,m=False,s=False, target='test', session=False, comment='None', path='/home/'):
    if camera == 'c1': 
        channel = 'on-band'
        ccd_temperature = b.c1.get_temperature()
    if camera == 'c2': 
        channel = 'off-band'
        ccd_temperature = b.c2.get_temperature()
    hdr = fits.Header()
    hdr['CAMERA'] = camera #c1 or c2
    hdr['CHANNEL'] = channel #on-band or off-band
    hdr['DATE-OBS'] = exp_starttime.strftime("%m/%d/%Y %H:%M:%S %Z")
    if session:
        hdr['OBSERVER'] = session['observer']
        hdr['TELESCOP'] = session['instrument']
        hdr['LATITUDE'] = session['lat']
        hdr['LONGITUD'] = session['lon']
        hdr['RA'] = np.round(m.mount.RightAscension,3)
        hdr['DEC'] = np.round(m.mount.Declination,3)
        hdr['ALT'] = np.round(m.mount.Altitude,3)
        hdr['AZ'] = np.round(m.mount.Azimuth,3)
    if s:
        hdr['AMB-TEMP'] = s.temperature()
        hdr['AMB-HUMI'] = s.humidity()
    hdr['CCD-TEMP'] = np.round(ccd_temperature,3)
    hdr['TARGET'] = target
    hdr['EXP-TIME'] = exptime # exposure time in seconds    
    hdr['BINNING'] = binning
    hdr['COMMENT'] = comment
    hdu = fits.PrimaryHDU(array, header=hdr)
    hdul = fits.HDUList([hdu])
    filename = camera + exp_starttime.strftime("_%m-%d-%Y_%H-%M-%S-%Z") + '.fits'
    if session:
        folder = Path(path)
        if not os.path.exists(folder):
            os.makedirs(folder)
    else:
        folder = Path(path)
        if not os.path.exists(folder):
            os.makedirs(folder)
    hdul.writeto(folder / filename, overwrite=False)
    return

def log(file, string):
    if file == "SCIP":
        filename = "SCIP_log" + datetime.now(timezone.utc).strftime("_%m_%d_%Y") + ".txt"
        filepath = Path(f"C:/Users/photometer/Box Sync/SCIP/Synced_Raw_Data/logs/{filename}")
    if file == 'temp':
        filename = "temp_log" + datetime.now(timezone.utc).strftime("_%m_%d_%Y") + ".txt"
        filepath = Path(f"C:/Users/photometer/Box Sync/SCIP/Synced_Raw_Data/logs/{filename}")
    else:
        filename = file
        filepath = Path(filename)
    with open(filepath, 'a') as logfile:  # Open file to append, create if needed
        print(f"{filepath}")
        time_str = datetime.now(timezone.utc).strftime("%m/%d/%Y %H:%M:%S %Z ")
        log_entry = time_str + string + "\n"
        logfile.write(log_entry)
        print(f"Logged: {log_entry}")
    return None

def identify_stars(image): #TODO clean up
    # Takes a numpy array, locates stars, returns their info, positons, and 
    # apertures to draw on the image.
    image = image.astype(float)
    image -= np.median(image) 
    bkg_sigma = mad_std(image)  
    daofind = DAOStarFinder(fwhm=4., threshold=8*bkg_sigma)  
    sources = daofind(image) 
    if sources != None:
        positions = np.transpose((sources['xcentroid'], sources['ycentroid']))  
        apertures = CircularAperture(positions, r=9.)  
        return positions, apertures
    else:
        return None, None, None

def show_img(nparray):
    # Shows any np array image
    plt.figure(figsize=(20,30))
    plt.subplot(121)
    plt.imshow(nparray, interpolation='none', vmin=np.percentile(nparray, 1), vmax=np.percentile(nparray, 99), cmap='gray')
    
def show_imgs(img1, img2, star_detect = False):
    # Show two images from the cameras
    plt.figure(figsize=(30,20))
    plt.subplot(121)
    plt.imshow(img1, interpolation='none', vmin=np.percentile(img1, 1), vmax=np.percentile(img1, 99), cmap='gray')
    if star_detect == True:
        positions, apertures = identify_stars(img1)
        if positions == None: print("No stars found")
        if apertures != None: apertures.plot(color='red', lw=1.5, alpha=0.5)
    plt.colorbar()
    plt.title("Camera C1")
    plt.subplot(122)
    plt.imshow(img2, interpolation='none', vmin=np.percentile(img2, 1), vmax=np.percentile(img2, 99), cmap='gray')
    if star_detect == True:
        positions, apertures = identify_stars(img2)
        if positions == None: print("No stars found")
        if apertures != None: apertures.plot(color='red', lw=1.5, alpha=0.5)
    plt.colorbar()
    plt.title("Camera C2")
    
def log_session_info(session):
    log('SCIP', str(session))
    for key in session:
        log("SCIP", key + ' ' + session[key])

def crop_im(im, dimy, dimx=None):
    '''
    Crop image to specified dimensions

    :input:
    im (np.ndarray) - Image to crop
    dimy (int) - # # of pixels to crop y dimension to
    dimx (int) - # of pixels to crop x dimension to (if None, square crop with 1st specified dimension)
    :output:
    (ndarray) - cropped image (ydim x xdim)
    '''
    if dimx is None:
        dimx = dimy

    crop_y = im.shape[0] - dimy
    crop_x = im.shape[1] - dimx

    return im[crop_y//2:-1*(int)(np.ceil(crop_y/2)), crop_x//2:-1*(int)(np.ceil(crop_x//2))]

def mark_area(im, dimy, dimx=None):
    '''
    Mark cropped area in image (as specified by passed dims)
    
    :input:
    im (np.ndarray) - Image to crop
    dimy (int) - # # of pixels to crop y dimension to
    dimx (int) - # of pixels to crop x dimension to (if None, square crop with 1st specified dimension)
    :output:
    (matplotlib.pyplot Figure) - image with cropped area marked
    '''

    if dimx is None:
        dimx = dimy

    crop_y = im.shape[0] - dimy
    crop_x = im.shape[1] - dimx

    y_inds = [crop_y//2, im.shape[0]-(int)(np.ceil(crop_y/2))]
    x_inds = [crop_x//2, im.shape[1]-(int)(np.ceil(crop_x/2))]

    plt.plot([x_inds[0], x_inds[0]], [y_inds[0],y_inds[1]],[x_inds[0], x_inds[1]], [y_inds[1],y_inds[1]],[x_inds[1], x_inds[1]], [y_inds[1],y_inds[0]],[x_inds[1], x_inds[0]], [y_inds[0],y_inds[0]], color='red', linewidth=1.5)
    plt.imshow(im, interpolation='none', vmin=np.percentile(im, 1), vmax=np.percentile(im, 99), cmap='gray')
    plt.colorbar()

    return plt.gcf()

def targ_from_name(target_name):
    '''
    Return astroplan Fixed Target from a given star/nebula name
    Some nebular targets are based on a specific region from literature (Nossal '93, Scherb '81, Mierkiewicz '06)

    Known H-alpha point source names: Arcturus, Aldebaran, NGC7662, Meissa (targeting nebular region around Meissa)
    '''
    if target_name == 'NAN':
        ra = '20h56m17s'
        dec='+44d24m03s'
        targ_coord = SkyCoord(ra, dec, frame=FK4(equinox=Time('B1950')))
        target = FixedTarget(coord=targ_coord, name=target_name)
    elif target_name == 'Barnards':
        # FIXME: Weird!
        ra = '00h23m26.73s'
        dec='+06d05m08.00s'
        targ_coord = SkyCoord(ra, dec, frame=FK4(equinox=Time('B1950')))
        target = FixedTarget(coord=targ_coord, name=target_name)
    # TODO: Add the other 2 nebulae
    elif target_name is not None:
        target = FixedTarget.from_name(target_name)

    return target

# Legacy functions that (JC) doesn't find useful
'''
def moves_to_center_stars(x1, y1, x2, y2): #TODO
    # Takes two locations of a stars on the image, and returns delta alt and 
    # delta az needed to center the star in both FOVs as best as possible.
    # Pixel size is 3.04 arcsec =  8.4444e-4 degrees
    
    x = ((x1-x2)**2)**0.5
    y = ((y1-y2)**2)**0.5
    
    
    return(delta_az, delta_alt)
'''
'''
def fieldline_WMM(lat, lon, h, time): # h is height above sea level
    gm = geomag.geomag.GeoMag()
    mag = gm.GeoMag(lat, lon, h, time)
    # mag.dec is the difference between true north and magnetic north
    # west-pointing is negative, east-pointing is positive
    # true north is celestial north within a few arcsec
    # mag.dip is the dip angle, along true noth
    if mag.dec < 0:
        az = 360 + mag.dec
    if mag.dec > 0: 
        az = mag.dec
    # change az to up rather than down along field line
    az = -az
    alt = mag.dip
    return alt, az
    # this only works for normal northern hemisphere (?)
    # I need to get back to this
'''