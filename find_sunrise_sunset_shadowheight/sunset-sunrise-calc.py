# -*- coding: utf-8 -*-
"""
Created on Thu Feb  6 16:24:07 2020

@author: photometer
"""

import astropy
import astroplan
import astropy.units as u
from datetime import timezone
import datetime

def log(string):
    with open("sunrise-sunset-aeronomy.txt", 'a') as logfile:  # Open file to append, create if needed
        log_entry = string + "\n"
        logfile.write(log_entry)
        print(f"Logged: {log_entry}")
    return None


#lon = -65.30
#lat = 18.32
#elev = 200
#conj_lon = 301.69
#conj_lat = -40.41
lon = -88.159
lat = 40.167
elev = 219
conj_lon = -103.14
conj_lat = -61.48
local = astroplan.Observer(timezone='UTC', longitude=lon*u.deg, latitude=lat*u.deg, elevation=elev*u.m)
conjugate = astroplan.Observer(timezone='UTC', longitude=conj_lon*u.deg, latitude=conj_lat*u.deg, elevation=0*u.m)



for day in range(3):
    day_offset = datetime.timedelta(days=-1*day)
    #now = datetime.datetime(2020, 3, 16, hour=16, minute=00, tzinfo = timezone.utc) + day_offset
    now = datetime.datetime.now(timezone.utc) + day_offset
    now = astropy.time.Time(now, format='datetime')
    
    # If an event is not found in the next 24 hours, add 5 min
    five_min = astropy.time.TimeDelta(5*60,format='sec')
    AST_offset = datetime.timedelta(hours=-6)
    
    # Sunset CONJ
    conjugate_sunset = conjugate.sun_set_time(now, which='next', horizon = -0.8333*u.deg)
    if type(conjugate_sunset.iso) != type('str'):
        # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
        conjugate_sunset = conjugate.sun_set_time(now+five_min, which='next', horizon = -0.8333*u.deg)
    conjugate_sunset = conjugate_sunset.to_datetime(timezone=timezone.utc) + AST_offset
    log(f'Conjugate sunset  : {conjugate_sunset.strftime("%m/%d/%Y %H:%M:%S")}')

    
    # Sunset LOCAL
    sunset = local.sun_set_time(now, which='next', horizon = -0.8333*u.deg)
    if type(sunset.iso) != type('str'):
        # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
        sunset = local.sun_set_time(now+five_min, which='next', horizon = -0.8333*u.deg)
    sunset = sunset.to_datetime(timezone=timezone.utc) + AST_offset
    log(f'Local sunset      : {sunset.strftime("%m/%d/%Y %H:%M:%S")}')
    

    # Sunrise CONJ
    conjugate_sunrise = conjugate.sun_rise_time(now, which='next', horizon = -0.8333*u.deg)
    if type(conjugate_sunrise.iso) != type('str'):
        # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
        conjugate_sunrise = conjugate.sun_rise_time(now+five_min, which='next', horizon = -0.8333*u.deg)
    conjugate_sunrise = conjugate_sunrise.to_datetime(timezone=timezone.utc) + AST_offset
    log(f'Conjugate sunrise : {conjugate_sunrise.strftime("%m/%d/%Y %H:%M:%S")}')
    
    # Sunrise LOCAL
    sunrise = local.sun_rise_time(now, which='next', horizon = -0.8333*u.deg)
    if type(sunrise.iso) != type('str'):
        # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
        sunrise = local.sun_rise_time(now+five_min, which='next', horizon = -0.8333*u.deg)        
    sunrise = sunrise.to_datetime(timezone=timezone.utc) + AST_offset
    log(f'Local sunrise     : {sunrise.strftime("%m/%d/%Y %H:%M:%S")}')
    log('\n')
