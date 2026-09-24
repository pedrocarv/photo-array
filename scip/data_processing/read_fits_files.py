# -*- coding: utf-8 -*-
"""
Created on Sun Sep 20 00:13:50 2020

@author: dawnh

Function to read in a folder of fits files and all their associated data.

Input: a path to a folder containing the files

Output: two lists, one for on-band and one for off-band. Each is a list of
dicts, each dict corresponds to one image. Keys are the filename, the image
data, and all the FITS keywords from the files.

Notes: All files in the folder need to have the same FITS keywords. 
"""

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime, timezone

def get_fits_from_folder(directory, sort_by_date=True, open_image=False, pix_sum=False):
    '''
    :input:
    directory (str or pathlib.Path) - path to directory containing files of interest
    [other args are self-explanitory]
    
    :output:
    master_on (list) - info from fits of on-band observations 
    master_off (list) - info from fits of off-band observations
    NOTE: each list item will be a dict containing all info from fits files (this func may not be good for many observations)
    '''
    if type(directory) is str:
        directory = Path(directory)

    master_on = [] # a list of dicts, sorted by datetime, of the data of each file.
    master_off = [] # a list of dicts, sorted by datetime, of the data of each file.
    
    keys = []
    
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        hdul = fits.open(directory / filename, ignore_missing_end=True)
        header = [*hdul[0].header]
        if keys == []:
             keys = header
        if header != keys:
            print(f"Found that {filename} does not match the same keys as previous files. F.")
            hdul.close()
            break
        imagedict = {'filename':file}
        imagedict['directory'] = directory
        if open_image:
            imagedict['image'] = hdul[0].data.astype(float)
        if pix_sum:
            imagedict['pix_sum'] = np.sum(hdul[0].data.astype(float))
        for key in keys: # add all the header file data
            #if key == "CCD_TEMP":
            #    print("warning: underscores in key names deprecated")          
            imagedict[key] = hdul[0].header[key]
        if filename[0:2] == 'c1':    
            master_on.append(imagedict)
        elif filename[0:2] == 'c2':
            master_off.append(imagedict)
        else:
            print('The filename {filename} does not start with "c1" or "c2" so it cannot be sorted into the correct list.')
            hdul.close()
            break
        hdul.close()
        
    print(f'You have opened {len(master_on)} sets of images. Your keywords are: {[*master_on[0]]}')

    if sort_by_date == True:
        # sort lists by date of observation
        master_on = sorted(master_on, key = lambda i: datetime.strptime(i['DATE-OBS'][:-4], "%m/%d/%Y %H:%M:%S"))
        master_off = sorted(master_off, key = lambda i: datetime.strptime(i['DATE-OBS'][:-4], "%m/%d/%Y %H:%M:%S"))
    
    return master_on, master_off

def image_from_filename(directory, filename):
    '''
    :input:
    directory (pathlib.Path or str) - path to directory containing file
    filename (str) - name of fijle to read

    :output:
    image (ndarray) - numpy array containing image
    '''
    if type(directory) is str:
        directory = Path(directory)
        
    hdul = fits.open(directory / filename)
    image = hdul[0].data.astype(float)
    return image

def add_images_to_master_lists(master_on, master_off):
    for master in [master_on, master_off]:
        for imagedict in master:
            imagedict['image'] = image_from_filename(imagedict['directory'], imagedict['filename'])


