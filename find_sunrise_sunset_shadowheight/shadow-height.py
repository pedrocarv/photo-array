# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 13:42:53 2020

@author: dawnh
"""

import astropy
import astroplan
import astropy.units as u
from datetime import timezone
import datetime
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm

#Areonomy
# lon = -88.159
# lat = 40.167
# elev = 219
# conj_lon = -103.14
# conj_lat = -61.48

#ECEB (found using https://gps-coordinates.org/)
lon = -88.228
lat = 40.115
elev = 224 # NOTE: Probably messed up, need to add height of building as well

#Arecibo
#lon = -66.75
#lat = 18.35
#elev = 317
#conj_lon = -58.31
#conj_lat = -40.41

local = astroplan.Observer(timezone='UTC', longitude=lon*u.deg, latitude=lat*u.deg, elevation=elev*u.m)
# conjugate = astroplan.Observer(timezone='UTC', longitude=conj_lon*u.deg, latitude=conj_lat*u.deg, elevation=0*u.m)
R = 6371

dates = []
mins = np.arange(0, 24*60 - 2, 20)
angles_2d = []
for day in range(30*12-2):
    day_offset = datetime.timedelta(days=1*day)
    angles_1d = []
    for i in range(len(mins)):
        min_offset = datetime.timedelta(minutes=int(mins[i]))
        #now = datetime.datetime.now(timezone.utc) + day_offset
        now = datetime.datetime(2025, 7, 6, hour=0, minute=00, tzinfo = timezone.utc) + day_offset + min_offset
        now_astro = astropy.time.Time(now, format='datetime')  
        
        # For overlap only
        '''
        sun_alt_l = local.sun_altaz(now_astro).alt.radian
        sun_alt_c = conjugate.sun_altaz(now_astro).alt.radian

        if sun_alt_l > 0.0: shadowh_l = 0.0
        else: shadowh_l = R*(1/np.cos( np.abs(sun_alt_l) ) - 1.0)
        if sun_alt_c > 0.0: shadowh_c = 0.0
        else: shadowh_c = R*(1/np.cos( np.abs(sun_alt_c) ) - 1.0)
        
        if shadowh_l > 500 and shadowh_c < 500:
            shadows = 1
        else: shadows = 0
        angles_1d.append(shadows)
        '''

        # For shadow height only, choose local or conj
        
        sun_alt = local.sun_altaz(now_astro).alt.radian
        if sun_alt > 0.0:
            shadowh = 0.0
        if sun_alt < 0.0:
            shadowh = R*(1/np.cos( np.abs(sun_alt) ) - 1.0)
            # if shadowh > 1000: shadowh = 1000 # NOTE: why zero this, should accept shadow alt up to 1e5 km
        angles_1d.append(shadowh)
        

        
    dates.append(now.strftime("%m/%d/%Y"))
    angles_2d.append(angles_1d)
    
    
fig, ax = plt.subplots(figsize=(10,10))

# for shadow height only:    
#ax.contourf(mins/60, dates, angles_2d, np.linspace(0, 1000, 1000, endpoint=True), vmin=0.0, vmax=999, cmap=cm.coolwarm)

# for overlap only:
ax.contourf(mins/60, dates, angles_2d)
    
every_nth = 10
for n, label in enumerate(ax.yaxis.get_ticklabels()):
    if n % every_nth != 0:
        label.set_visible(False)
        
# For shadow height only

m = plt.cm.ScalarMappable(cmap=cm.coolwarm)
m.set_array(angles_2d)
m.set_clim(0., 1000.)
fig.colorbar(m, ax=ax, boundaries=np.linspace(0, 1000, 1000, endpoint=False))


ax.set_xlabel("Hour (UTC)")
ax.set_ylabel("Date (UTC)")
ax.set_title("Local shadow height, Aeronomy, km")

plt.show()