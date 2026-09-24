# -*- coding: utf-8 -*-
"""
Created on Fri Mar 10 15:42:11 2023

@author: alexmz2
"""

from dark_calibration_functions import *

target1 = "PATHS1_images_Mar_23"
target2 = "PATHS2_images_Mar_23"
target3 = "PATHS3_images_Mar_23"

#dark_calibration_run(target, test_run = False, do_bias = True)

dark_calibration_run(target1, test_run = False, do_bias = False)