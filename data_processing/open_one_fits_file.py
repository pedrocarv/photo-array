# -*- coding: utf-8 -*-
"""
Created on Tue Nov  3 12:44:23 2020

@author: dawnh
"""

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime, timezone

filename ='1x1bin_1s_dk1.fit' 
hdul = fits.open(filename)
image = hdul[0].data.astype(float)
# plt.imshow(image, interpolation='none', vmin=np.percentile(image, 1), vmax=np.percentile(image, 99),cmap='gray')
# plt.colorbar()
# print(np.sum(image)/(3379*2703))
print(np.mean(image), np.median(image))
plt.hist(image.flatten(), bins=1000, log=True)