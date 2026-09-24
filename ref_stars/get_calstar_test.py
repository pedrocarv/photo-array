# -*- coding: utf-8 -*-
"""
Created on Wed May 19 00:13:23 2021

@author: dawnh
"""

# def set_calstar(session, NE, m):
#     m.local
#     m.conjugate
#     now = astropy.time.now()
#     sunrise = NE.sunrise()
#     time_range = 
    
#     m.calstar = None
#     m.calstar_ra = None
#     m.calstar_dec = None
#     m.calstar_expiry = None
#     return

import numpy as np
import astropy
from astroplan import Observer, FixedTarget
from astropy.coordinates import SkyCoord
import astropy.units as u
from astroplan import AltitudeConstraint
from astroplan import is_always_observable

lon = -88.159
lat = 40.167
elev = 219
conj_lon = -103.14
conj_lat = -61.48

local = Observer(timezone='UTC', longitude=lon*u.deg, latitude=lat*u.deg, elevation=elev*u.m)
conjugate = Observer(timezone='UTC', longitude=conj_lon*u.deg, latitude=conj_lat*u.deg, elevation=0*u.m)

now = astropy.time.Time.now()
sunrise = now + astropy.time.TimeDelta(11*60*60,format='sec')
endtime = sunrise
time_range = (now, endtime)
stardata = np.load("stardata.npy")
targets = [FixedTarget(coord=SkyCoord(ra=stardata[i,5]*u.deg, dec=stardata[i,6]*u.deg), name = str(stardata[i,0])) for i in range(len(stardata))]
constraints = [AltitudeConstraint(30*u.deg, 90*u.deg)]
observables = np.array([]) 
always_observable = np.array(is_always_observable(constraints, local, targets, time_range=time_range))
observables = stardata[always_observable]
while len(observables) == 0:
    endtime = endtime - astropy.time.TimeDelta(2*60*60,format='sec')
    time_range = (now, endtime)
    try:
        always_observable = np.array(is_always_observable(constraints, local, targets, time_range=time_range))
    except ValueError:
        print("No stars found observable for the next two hours.") 
        break
    observables = stardata[always_observable]
brightest = observables[np.argmin(observables[:,3])]

    
calstar_mag = brightest[3]
calstar_ra = brightest[5]
calstar_dec = brightest[6]
calstar_expiry = endtime


















