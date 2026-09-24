# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 00:27:54 2020

@author: dawnh


Input a folder of fits file observations and calibrate data. Output is
calibration factors and measurements.
"""

from read_fits_files import *
# Open fits files

directory = '/Users/dawnh/Documents/Local Raw Data/09-04-2020_one-night'
master_on, master_off = get_fits_from_folder(directory, sort_by_date=True)

#%% STAR CALIBRATION

# Pick out cal star images
on_calstars_list = []
for i in range(len(master_on)):
    if master_on[i]['TARGET'] != 'airglow': 
        on_calstars_list.append(master_on[i])
        
off_calstars_list = []
for i in range(len(master_off)):
    if master_off[i]['TARGET'] != 'airglow': 
        off_calstars_list.append(master_off[i])

# for each star:
    # find brighestest star (error: no star)
    # get ap phot of star (error: saturation, very dim)
    # get sensitivity using absolute brightness
    # line fit 
    # plot
    # save zenith sensitivity
    # save error on line fit

#%% APPLY BIAS AND DARK CORRECTION



#%% DERIVE J_0 AND J_8446 FOR EACH TIME

# plot and save

