# -*- coding: utf-8 -*-
"""
Created on Tue Jun 16 15:42:56 2020

@author: dawnh

Open filter transmission plot data and 8446 triplet plot data, interpolate

"""

import csv
import matplotlib.pyplot as plt
import numpy as np
from scipy import interpolate
from scipy import integrate
from time import time

#%%

'''
------------------------------------------------------------------------------
Open files and create interpolation functions
------------------------------------------------------------------------------
'''
#unit = 1e-9
unit = 1

# Extract triplet data -------------------------------------------------------
with open('slanger-triplet-data.csv') as csvfile:        
    rows=list(csv.reader(csvfile,delimiter=','))
line_lambdas = []
line_intensities = []
for row in rows: 
    line_lambdas.append(float(row[0])*unit)
    line_intensities.append(float(row[1]))
line_lambdas = np.array(line_lambdas[11:46])
line_intensities = np.array(line_intensities[11:46])
line_intensities = np.array(line_intensities)/max(line_intensities)
csvfile.close()


# add zeros to the data to aid interpolation
line_lambdas = np.concatenate((line_lambdas, np.linspace(844.72*unit, 849.0*unit, 500)))
line_intensities = np.concatenate((line_intensities, np.zeros(500)))

line_lambdas = np.concatenate((np.linspace(843.8*unit, 844.6*unit, 100), line_lambdas))
line_intensities = np.concatenate((np.zeros(100), line_intensities))

# FUNCTION FOR EMISSION LINE
interp_line = interpolate.interp1d(line_lambdas, line_intensities, kind = 'cubic')

line_plotting_lambdas = np.linspace(843.8*unit, 849.0*unit, 500)
line_intensity_interp = interp_line(line_plotting_lambdas)

# offband --------------------------------------------------------------------

with open('transmission_plot_digitized_offband.csv') as csvfile:        
    rows=list(csv.reader(csvfile,delimiter=','))
off_lambdas = []
off_trans= []
for row in rows: 
    off_lambdas.append(float(row[0])*unit)
    off_trans.append(float(row[1]))
off_trans = np.array(off_trans)/max(off_trans)
csvfile.close()

# add zeros to the data to aid interpolation
off_lambdas = np.concatenate((off_lambdas, np.linspace(849.2*unit, 849.0*unit, 100)))
off_trans = np.concatenate((off_trans, np.zeros(100)))

off_lambdas = np.concatenate((np.linspace(843.8*unit, 846.8*unit, 100), off_lambdas))
off_trans = np.concatenate((np.zeros(100), off_trans))

# FUNCTION FOR OFFBAND FILTER
interp_off = interpolate.interp1d(off_lambdas, off_trans, kind = 'linear')

offband_plotting_lambdas = np.linspace(843.8*unit, 849.0*unit, 500)
off_trans_interp = interp_off(offband_plotting_lambdas)

# onband ---------------------------------------------------------------------
with open('transmission_plot_digitized_onband.csv') as csvfile:        
    rows=list(csv.reader(csvfile,delimiter=','))
on_lambdas = []
on_trans= []
for row in rows: 
    on_lambdas.append(float(row[0])*unit)
    on_trans.append(float(row[1]))
on_trans = np.array(on_trans)/max(on_trans)
csvfile.close()

# add zeros to the data to aid interpolation
on_lambdas = np.concatenate((on_lambdas, np.linspace(845.9*unit, 850.0*unit, 100)))
on_trans = np.concatenate((on_trans, np.zeros(100)))

on_lambdas = np.concatenate((np.linspace(843.8*unit, 843.5*unit, 100), on_lambdas))
on_trans = np.concatenate((np.zeros(100), on_trans))

#FUNCTION FOR ONBAND FILTER
interp_on = interpolate.interp1d(on_lambdas, on_trans, kind = 'linear')
on_plotting_lambdas = np.linspace(843.8*unit, 849.0*unit, 500)
on_trans_interp = interp_on(on_plotting_lambdas)

#%%

'''
------------------------------------------------------------------------------
Integrals of the emission line, off-band, and on-band independently
------------------------------------------------------------------------------
'''

line_int = integrate.trapz(line_intensity_interp, line_plotting_lambdas)
foff_int = integrate.trapz(off_trans_interp, offband_plotting_lambdas)
fon_int = integrate.trapz(on_trans_interp, on_plotting_lambdas)
fon_line_total = integrate.trapz(line_intensity_interp*on_trans_interp, line_plotting_lambdas)
plt.plot(line_plotting_lambdas, line_intensity_interp*on_trans_interp, label="Emission line through the filter")
plt.plot(on_plotting_lambdas, on_trans_interp, label="On-band filter transmissivity")
plt.plot(line_plotting_lambdas, line_intensity_interp, label="844.6 emission line")
plt.legend()

#%%

# Radiance and transmission plot ---------------------------------------------

Jbb_i = np.ones(500)*0.3
rad = np.add(Jbb_i, line_intensity_interp)
rad_l = np.linspace(843.8, 849.0, 500)

fig, ax = plt.subplots(figsize=(15,5))
ax.plot(rad_l, rad, lw=2, label="OI emission and background", color='firebrick')
ax.plot(offband_plotting_lambdas, off_trans_interp, '--', label=r"Off-band filter transmissivity, $\theta$=0", color='steelblue')
ax.plot(on_plotting_lambdas, on_trans_interp, '--', label=r"On-band filter transmissivity, $\theta_{min}$=0$^{\circ}$", color='seagreen')
on_trans_thetamax = interp_on(on_plotting_lambdas+0.129) #TODO where is this from
ax.plot(on_plotting_lambdas, on_trans_thetamax, '--', label=r"On-band filter transmissivity, $\theta_{max}$=1.84$^{\circ}$", color='limegreen')

ax.set_ylabel("Spectral radiance (arb. units)")
ax.set_xlabel("Wavelength (nm)")
ax2 = ax.secondary_yaxis('right')
ax2.set_ylabel("Normalized filter transmissivity")
plt.legend()
plt.title("Simulated airglow and filter transmissivities")


#%%

# Filter center wavelength and temperature -----------------------------------

# fig = plt.figure(figsize=(10,10))
# onband_l = [844.80, 844.67]
# onband_T = [25.0, 22.0]
# plt.plot(onband_T, onband_l)


with open('center-vs-angle.csv') as csvfile:        
    rows=list(csv.reader(csvfile,delimiter=','))
centers = []
angles= []
for row in rows: 
    angles.append(float(row[0]))
    centers.append(float(row[1]))
csvfile.close()

centers_f = interpolate.interp1d(angles, centers, kind = 'quadratic')

def get_integral(theta):
    c_offset = 844.677 - centers_f(theta) # find center wavelength offset from theta=0
    nms = np.linspace(843.8, 849.0, 500)
    nms_offset = np.linspace(843.8, 849.0, 500) + c_offset
    airglow = interp_line(nms)
    onband_filter = interp_on(nms_offset)
    signal = airglow*onband_filter
    return integrate.trapz(signal, on_plotting_lambdas)

binning = 20

sizex = int(np.round(3379/binning)) #338
FOVx = 171.45/60
scalex = FOVx/sizex

sizey = int(np.round(2703/binning)) #270
FOVy = 137.15/60
scaley = FOVy/sizey

fon_line_sum = 0

t1 = time()
pix = []
offsetx = sizex/2.0 + 0.5
offsety = sizey/2.0 + 0.5
for i in range(sizex):
    for j in range(sizey):
        r = ( ((float(i)-offsetx)*scalex)**2 + ((float(j)-offsety)*scaley)**2)**0.5 # scaled to FOV
        rad = get_integral(r)
        fon_line_sum += rad
        #pix.append([rad, i, j])

t2 = time()
print("fon_line_sum")
print(fon_line_sum)
'''
image = np.zeros((sizex,sizey))

for p in pix:
    image[p[1],p[2]] = p[0]
#pixmax = np.amax(image)
#image = image/pixmax
figimg = plt.figure()
plt.imshow(image)
plt.colorbar()
plt.title("Normalized on-band counts per pixel due to variable angle of incidence")
plt.xlabel("CCD pixels")
plt.ylabel("CCD pixels")
'''