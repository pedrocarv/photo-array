# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 01:29:14 2020

@author: dawnh

Find sensitivity given aperture photometry counts and actaul flux of a star.
"""
import numpy as np

xpix = 3379 # number of pixels in the x-direction
ypix = 2703 # number of pixels in the y-direction
Npix = ypix*xpix # total number of pixels
A = np.pi*(50e-3/2)**2 # aperture area, m^2
xFOV = 171.45/60 # deg
yFOV = 137.15/60 # deg
omega = (xFOV * yFOV) / (180/np.pi)**2 # FOV in sr
omega_p = omega/Npix
t = 30 # IP in s
fon_integral = 3.136206756591953e-10 # in units of m
h = 6.626e-34 # planks constant J-s
c = 299792458 # speed of light m/s
phi = 1e10/(4*np.pi) # photon flux of 1 R [phot/m^2/s/R/4pi]

star = "pollux"
counts = 115469.3
#counts = 116315.9
real_stellar_flux = 8.60e-9 # w/m^2, from CDS portal
lambda_of_real_stellar_flux = 8.78e-7 # m, from CDS portal


# star='alcor'
# counts = 2264.8
# real_stellar_flux = 2.39e-10

# star = 'mizar'
# counts = 21277.9
# real_stellar_flux = 9.97e-10

# star = 'capella'
# counts = 139006.5
# real_stellar_flux = 1.88e-8

# star  = 'arcturus'
# counts = 462712
# real_stellar_flux = 3.62e-8

# efficiency alpha counts per photon
alpha = counts / (t * A * fon_integral * real_stellar_flux*(1/(h*c)))

phot_per_R_sec  = phi*omega*A

cts_per_R_sec = alpha * phot_per_R_sec # counts per Rayleigh-second

print(star)
print(cts_per_R_sec)