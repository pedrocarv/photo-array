# -*- coding: utf-8 -*-
"""
Created on Thu May 28 22:48:51 2020

@author: dawnh
"""
import numpy as np
import matplotlib.pyplot as plt


# INPUTS
target = 'airglow' # 'airglow' or 'star'
I_8446 = # 844.6 intensity in R
I_bg = # background intensity in R
s_flux = # stellar flux in W/m^2/
max_exp = 5 # maximum exposure time in minutes
t = np.linspace(60, 60*max_exp, num = 100)
#t=1

# BINNING --------------------------------------------------------------------

binning = 12

# PHOTOMETER CHARACTERISTICS -------------------------------------------------

xpix = np.floor(3379/binning) # number of pixels in the x-direction
ypix = np.floor(2703/binning) # number of pixels in the y-direction
xlostpix = 3379 - xpix*binning # fraction binned pixels are not read out
ylostpix = 2703 - ypix*binning
xfraction = xpix*binning/(3379)
yfraction = ypix*binning/(2703)
Npix_CCD = 3379*2703*xfraction*yfraction
Npix = ypix*xpix # total number of pixels
#pixsize = 3.69e-6 # length of the side of one pixel, m
A = np.pi*(5/2)**2 # aperture area in cm^2
xFOV = 171.45/60 * (xfraction) # deg. FOV is slightly smaller when binning
yFOV = 137.15/60 * (yfraction) # deg
omega = (xFOV * yFOV) / (180/np.pi)**2 # FOV in sr (0.002?)
omega_p = omega/Npix # FOV in sr of one pixel

gain = 1/0.19 # counts per e-
phi = 1e10/(4*np.pi) # photon flux of 1 R [phot/m^2/s/R/4pi]

# EFFICIENCY FACTORS ALPHA ---------------------------------------------------
QE = 0.35 # approx quantum efficiency at 844.6 nm 
trans_on = 0.50 # peak on-band filter transmission
trans_off = 0.44 # peak off-band filter transmission
trans_dome = 1#0.90 # dome transmissivity estimation
atms = 0.45 # atmospheric transmission
alpha_on = QE*trans_on*trans_dome*atms/gain # counts/phot
alpha_off = QE*trans_off*trans_dome*atms/gain # counts/phot

# FILTER INTEGRALS -----------------------------------------------------------

# Normalized integrated off-band filter transmissivity 
foff_int = 0.2997345826982963

# Normalized integrated on-band filter transmissivity
fon_int = 0.3136206756591936

# Normalized integrated emission line
line_int = 0.04667402571951789

# Normalized integrated on-band multplied with normalized emission line
line_thru_fon = 0.023850237590858202

# if binning == 1:
#     fon_line_sum = 312831.387111227
# elif binning == 4:
#     fon_line_sum = 19565.03587472763   
# elif binning == 8:
#     fon_line_sum = 4885.513809952065  
# elif binning == 12:
#     fon_line_sum = 2173.2956941860702   
# elif binning == 20:
#     fon_line_sum = 781.4992101737781     
# elif binning == 28:
#     fon_line_sum = 402.0650778435002 
# elif binning == 32:
#     fon_line_sum = 305.02870134019577

# I should derive this properly, but from a line fit:
fon_line_sum = np.e**(-1.99940563220494*np.log(binning)+12.652922188575378)


# EXPECTED RADIANCES ---------------------------------------------------------
#I_background = 1 * 10**6 # R to phot/(4*pi*cm^2*s*sr)
Jbb_input = 10 # phot/sec/arcsec^2/nm/m^2
Jbb_input *= 4.25e10 # phot/sec/sr/nm/m^2
Jbb_input *= 1e-4 # phot/sec/sr/nm/cm^2

#Jbb_input = 3.19e6 for 1 R input for on-band filter

I_8446= 30 * 10**6 # R to phot/(4*pi*cm^2*s*sr)
J0_input = I_8446/line_int


# EXPECTED ON-BAND AND OFF-BAND COUNTS ---------------------------------------

N_off = Npix * t * A * omega_p * alpha_off * foff_int * Jbb_input
# print(fon_int * Jbb_input)

N_on = t * A * omega_p * alpha_on * (Npix * Jbb_input * fon_int + J0_input * fon_line_sum)

# NOISE SOURCES --------------------------------------------------------------
'''
Checked for 12/21
'''
if binning == 1:
    read = 6 # read noise, e-/pixel
elif binning == 4:
    read = 7    
elif binning == 8:
    read = 9
elif binning == 12:
    read = 10  
elif binning == 20:
    read = 20      
elif binning == 28:
    read = 23      
elif binning == 32:
    read = 48
else:
    print("Need to characterize that still...")
    read = None      
NR = Npix*read**2
dark = 0.0003 # dark current, e-/pixel/s, -10C
ND = t*Npix_CCD*dark
NPoff = np.sqrt(N_off/gain)
NPon = np.sqrt(N_on/gain)

sigma_Noff = np.sqrt((NPoff)**2 + ND + 2*NR)
sigma_Non = np.sqrt((NPon)**2 + ND + 2*NR)


SNRoff = (N_off/gain)/sigma_Noff
SNRon = (N_on/gain)/sigma_Non


# GET UNCERTAINTY ON BIAS SUM ------------------------------------------------

bias_sum_sigma = 0#np.sqrt(Npix*read**2)


# SOLVING FOR J0 UNCERTAINTY -------------------------------------------------

KJbb = Npix * t * A * omega_p * alpha_off * foff_int
KJ0_A = Npix * t * A * omega_p * alpha_on * fon_int
KJ0_B = t * A * omega_p * alpha_on * fon_line_sum

J0 = N_on/KJ0_B - (KJ0_A * N_off)/(KJ0_B * KJbb) # scale factor

dJ0dNoff = - KJ0_A / (KJ0_B * KJbb) # partial derivative Noff
dJ0dNon = - 1/KJ0_B # partial derivative Non

# variance of J0
sigma_J0 = np.sqrt( (sigma_Noff**2+bias_sum_sigma**2) * gain**2 * dJ0dNoff**2 + (sigma_Non**2+bias_sum_sigma**2)
                   * gain**2 * dJ0dNon**2) 

UJ0 = sigma_J0/J0 # Uncertainty on scale factor J0 
UI8446 = UJ0

fig, ax = plt.subplots(1,2, figsize=(10,4))
ax[1].plot(t, UI8446*100)#, label=f"Scale factor $J_0$")
ax[1].set_xlabel("Integration time (s)")
ax[1].set_ylabel("Relative uncertainty (%)")
ax[1].set_title(f"Relative uncertainty of the 844.6 nm line radiance measurement")#, binning = {binning}")



ax[0].plot(t, SNRon, label="On-band SNR, 30R emission + 13R background")
ax[0].plot(t, SNRoff, label="Off-band SNR, 13R background")
ax[0].legend()
ax[0].set_ylabel("SNR")
ax[0].set_xlabel("Integration time (s)")
ax[0].set_title(f"SNR, on-band and off-band")#, binning = {binning}")


        