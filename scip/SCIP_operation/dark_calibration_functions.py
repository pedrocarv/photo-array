# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 15:34:02 2023

@author: alexmz2
"""

### Written by Robert Irvin 09/15/2022

from photometer_classes import *
import os
from util import *
import numpy as np
import time

def set_ccd_temp(b, temp, cooldown_wait):
    b.set_cooling(temp)
    c1t, c2t = b.get_CCD_temps()
    had_to_wait = False
    while  c1t > temp+0.1 or c2t > temp+0.1 or c1t < temp-0.1 or c2t < temp-0.1:
        had_to_wait = True
        print("Waiting on:")
        print(temp)
        time.sleep(30)
        c1t, c2t = b.get_CCD_temps()
    if had_to_wait:
        print(f"Reached temp, waiting {cooldown_wait} minutes.")
        time.sleep(60*cooldown_wait)
    print("Reached temp and ready to go.")

# target is the folder where things will be saved
# number of images is how many to take
# pass in arrays of of temps, exposure time, and binning value
def take_images(b, target, exptime, binning, num_images = 1, comment = ""):
    for i in range(num_images):
        print(f'Image #{i+1} of {num_images}, exptime {exptime}, binning {binning}')
        img1, img2, exp_starttime = b.dual_exposure(exptime, False,s=False, show = False, binning = binning, save = 'fits', target = target, session=False, comment=comment)

#wrapper to run all of the calibration tests. Just run this for each camera
#target given is the full target directory
#automatically disconnects after everything finishes.
def dark_calibration_run(target, test_run = False, do_bias = True):
    cooldown_wait = 15 #can be 0, but just in case
    if test_run:
        cooldown_wait = 0
    
    b = Binocular()
    
    if not test_run:
        print("120 seconds to shut lights off and leave the lab!")
        time.sleep(120) #time to shut lights off and get out of slab
    
    if do_bias:
        print("\n \nStarting on biases \n \n")
        
        #start with bias - always 1ms exposure, 0C temp, binning val 1.
        bias_num_images = 500
        if test_run:
            bias_num_images = 5
        
        set_ccd_temp(b, 0, cooldown_wait)
        take_images(b, target + "\\bias_frames", 0.001, 1, num_images = bias_num_images, comment = "bias")
        
    print("\n \nStarting on darks \n \n")
    #then dark sequence of darks, at varying exposures and temps. Set binning val to 1 always.
    dark_temps = [0, 10, 20]
    dark_expo_times = [1, 10, 30, 60, 120, int(3.5*60), 5 * 60]
    dark_binning_vals = [1]
    dark_num_images = 10 #should provide enough spread at each setting to get a good value
    
    if test_run:
        dark_temps = [10, 20]
        dark_expo_times = [10]
        dark_num_images = 2

    for temp in dark_temps:
        set_ccd_temp(b, temp, cooldown_wait)
        for exptime in dark_expo_times:
            for binning in dark_binning_vals:
                temp_target = target + "\\dark_frames\\temp_" + str(temp) + "\\expotime_" + str(exptime)
                take_images(b, temp_target, exptime, binning, num_images = dark_num_images, comment = "dark")

    # disconnect so next camera can be reloaded!
    b.disconnect()