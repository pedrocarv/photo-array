# -*- coding: utf-8 -*-
"""
Created on Wed Nov 13 19:06:57 2019

@author: photometer
"""

'''
Sunset: sun is at the local horizon
Dusk: sun at 18 deg below local or geomagnetic conjugate location, whichever is last
Dawn: sun at 18 deg below local or geomagnetic conjugate location, whichever is first
Sunrise: sun is at the local horizon
   _________________________________________________________________
  |            |          |                 |          |            |
  |  daylight  | twilight |      night      | twilight |  daylight  |
  |____________|__________|_________________|__________|____________|
               |          |                 |          |
            sunset       dusk             dawn      sunrise  
            
18 deg below horizon is 'astronomical night'
ALL TIMES USE UTC
'''

import configparser
import numpy as np
import astropy
from astroplan import Observer, FixedTarget, AltitudeConstraint, is_always_observable
from astropy.coordinates import SkyCoord
import astropy.units as u
import time
from datetime import datetime, timezone
from photometer_classes import *
from util import log, log_session_info

FAKETESTINGTIME = False
FAKETIME = datetime(2020, 2, 19, hour=23, minute=50, tzinfo = timezone.utc)

class NextEvents(object):
    '''
    When created, this class calculates the time of the next sunset, sunrise, 
    dawn, and dusk. Then, the user can call NextEvents.sun_status() to find
    out if it is currently daylight, twilight, or night. If the next events 
    need to be re-calculated, sun_status returns "recalculate". It takes some
    processing time to find the times of next events so it's better to only 
    calculate when needed.
    
    '''
    def __init__(self, local, conjugate):
        self.calculate_next(local, conjugate)
        
    def calculate_next(self, local, conjugate):
        
        now = datetime.now(timezone.utc)
        
        # FOR TESTING ONLY.
        if FAKETESTINGTIME:
            now = FAKETIME
        
        # Convert datetime object to astropy time object
        now = astropy.time.Time(now, format='datetime')
        
        # If an event is not found in the next 24 hours, add 5 min. Depending
        # on when this is calcualted, there may be more than 24 hours between
        # the current time and the next sunset. Two consecutive sunsets or
        # sunrises can be nearly 3 min apart in Champaign, IL. At the 
        # conjugate it can be just over 3 min. 
        five_min = astropy.time.TimeDelta(5*60,format='sec')
        
        # Sunset here means local sunset (cannot observe during local daylight)
        sunset = local.sun_set_time(now, which='next', horizon = -0.8333*u.deg)
        if type(sunset.iso) != type('str'):
            # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
            sunset = local.sun_set_time(now+five_min, which='next', horizon = -0.8333*u.deg)
        sunset = sunset.to_datetime(timezone=timezone.utc)
        self.sunset = sunset
        
        # Sunrise here means local sunrise (cannot observe during local daylight)
        sunrise = local.sun_rise_time(now, which='next', horizon = -0.8333*u.deg)
        if type(sunrise.iso) != type('str'):
            # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
            sunrise = local.sun_rise_time(now+five_min, which='next', horizon = -0.8333*u.deg)        
        sunrise = sunrise.to_datetime(timezone=timezone.utc)
        self.sunrise = sunrise
        
        # Dusk is whichever happens LAST - local dusk or conjugate sunset
        local_dusk = local.sun_set_time(now, which='next', horizon = -18*u.deg)
        if type(local_dusk.iso) != type('str'):
            # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
            local_dusk = local.sun_set_time(now+five_min, which='next', horizon = -18*u.deg)
        conjugate_sunset = conjugate.sun_set_time(now, which='next', horizon = -0.8333*u.deg)
        if type(conjugate_sunset.iso) != type('str'):
            # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
            conjugate_sunset = conjugate.sun_set_time(now+five_min, which='next', horizon = -0.8333*u.deg)
        dusk = max(local_dusk, conjugate_sunset)
        dusk = dusk.to_datetime(timezone=timezone.utc)
        self.dusk = dusk
        
        # Dawn is whichever happens FIRST - local dawn or conjugate sunrise
        local_dawn = local.sun_rise_time(now, which='next', horizon = -18*u.deg)
        if type(local_dawn.iso) != type('str'):
            # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
            local_dawn = local.sun_rise_time(now+five_min, which='next', horizon = -18*u.deg)
        conjugate_sunrise = conjugate.sun_rise_time(now, which='next', horizon = -0.8333*u.deg)
        if type(conjugate_sunrise.iso) != type('str'):
            # Bad value - does not occur in the next 24 hours. Add 5 minutes and calculate again.
            conjugate_sunrise = conjugate.sun_rise_time(now+five_min, which='next', horizon = -0.8333*u.deg)
        dawn = min(local_dawn, conjugate_sunrise)
        dawn = dawn.to_datetime(timezone=timezone.utc)
        self.dawn = dawn
        
    def sun_status(self):
        now = datetime.now(timezone.utc) # Get the time when the function is called.
        if FAKETESTINGTIME:
            now = FAKETIME
        sunset_delta = self.sunset-now
        dusk_delta = self.dusk-now
        dawn_delta = self.dawn-now
        sunrise_delta = self.sunrise-now
        deltas = [sunset_delta, dusk_delta, dawn_delta, sunrise_delta]
        smallest_delta = min(d for d in deltas if d.total_seconds()>0) # Nearest one which is in the future.
        
        if smallest_delta == sunset_delta:
            return 'daylight'
        elif smallest_delta == dusk_delta:
            return 'twilight'
        elif smallest_delta == dawn_delta:
            return 'night'
        elif smallest_delta == sunrise_delta:
            return 'twilight'
        else:
            return 'recalculate'

def set_calstar(session, NE, m):
    now = astropy.time.now()
    sunrise = NE.sunrise()
    endtime = sunrise
    time_range = (now, endtime)
    stardata = np.load("stardata.npy")
    
    targets = [FixedTarget(coord=SkyCoord(ra=stardata[i,5]*u.deg, dec=stardata[i,6]*u.deg), name = str(stardata[i,0])) for i in range(len(stardata))]
    constraints = [AltitudeConstraint(30*u.deg, None)]
    observables = np.array([]) 
    always_observable = np.array(is_always_observable(constraints, m.local, targets, time_range=time_range))
    observables = stardata[always_observable]
    while len(observables) == 0:
        endtime = endtime - astropy.time.TimeDelta(2*60*60,format='sec')
        time_range = (now, endtime)
        try:
            always_observable = np.array(is_always_observable(constraints, m.local, targets, time_range=time_range))
        except ValueError:
            log("SCIP", "\n") 
            log("SCIP", "WARNING") 
            log("SCIP", "No calibration stars found observable for the next two hours.") 
            log("SCIP", "This should not happen. More star options, or a differnt method for choosing stars may be needed.") 
            log("SCIP", "\n") 
            m.calstar_mag = None
            m.calstar_ra = None
            m.calstar_dec = None
            m.calstar_expiry = None
            return False
        observables = stardata[always_observable]
    brightest = observables[np.argmin(observables[:,3])]
    
    m.calstar_mag = brightest[3]
    m.calstar_ra = brightest[5]
    m.calstar_dec = brightest[6]
    m.calstar_expiry = endtime
    return True

    
def image_airglow(session, b,m,s):
    '''
    Takes a picture of the zenith with both cameras. 

    Parameters
    ----------
    session : The observation session data
    b : The Binocular object
    m : The CubeProMount object
    s : The TempHumidSensor object

    Returns
    -------
    None.

    '''
    m.goto('altaz', 90.0, 0.0)
    m.track(False)
    exptime = float(session['airglow exposure time'])
    b.dual_exposure(exptime, m,s, show = False, save = 'fits', target = 'airglow', session=session, comment='None')
    return

def image_calstar(session, NE, b,m,s):
    '''
    Takes a picture of the calibration star with both cameras.

    Parameters
    ----------
    session : The observation session data
    NE : The NextEvents object
    b : The Binocular object
    m : The CubeProMount object
    s : The TempHumidSensor object

    Returns
    -------
    None.

    '''
    if m.calstar_mag == None: # If no calibration star has been chosen yet
        set_calstar(session, NE, m) # Choose a new calibration star
        if m.calstar_mag == None: # If STILL no calibration star, it couldn't find one.
            log("SCIP", "No calibration star to image. See above warning.")
            return
    if m.calstar_expiry <= datetime.now(timezone.utc): # If the star has expired
        set_calstar(session, NE, m) # Choose a new calibration star
        if m.calstar_mag == None: # If no calibration star, it couldn't find one.
            log("SCIP", "No calibration star to image. See above warning.")
            return
    m.goto('radec', m.calstar_ra, m.calstar_dec)
    m.track(True)
    exptime = float(session['calibration star exposure time'])
    b.dual_exposure(exptime, m,s, show = False, save = 'fits', target=f'calstar, I mag {m.calstar_mag}', session=session)
    return

def check_schedule(filename):
    '''

    Parameters
    ----------
    filename : The name of the schedule file.

    Returns
    -------
    session : False if there no session to run. Otherwise, session is a list 
        of the parameters for the observation session
    time_to_wait : Number of seconds to wait before re-checking the schedule. 
        Zero if session is not False.

    '''
    schedule = configparser.ConfigParser()
    try:
        schedule.read(filename)
    except:
        log("SCIP", "Schedule file not found.")
    for session in schedule.sections():
        session = schedule[session]
        # Read start and end time and specify UTC
        starttime = datetime.strptime(session['start']+" +0000", "%m/%d/%Y %H:%M:%S %z")
        endtime = datetime.strptime(session['end']+" +0000", "%m/%d/%Y %H:%M:%S %z")
        now = datetime.now(timezone.utc)
        if FAKETESTINGTIME:
            now = FAKETIME
        if starttime > now:
            # Wait either 15 min or until the session is to start.
            waittime = starttime - now
            log("SCIP", f"Next session {session} starts in {waittime}.")
            session = False
            time_to_wait = min(waittime.total_seconds(), 15*60)
            return session, time_to_wait
        if endtime < now:
            # This session is past.
            log("SCIP", f"Session {session} ended at {endtime}.")
        if (starttime < now) and (endtime > now):
            time_to_wait = 0
            return session, time_to_wait
    # if all sessions are looped through, and all are in the past, you get here    
    log("SCIP", "All scheduled sessions in the past. Check again in 15 minutes.")
    session = False
    time_to_wait = 15*60
    return session, time_to_wait
            
def run_session(session, b,m,s):
        # It is time to run a scheduled session.
        # Log the settings of the session for future reference.
        # log_session_info(session)
        # Get the session start time.
        starttime = datetime.strptime(session['start']+" +0000", "%m/%d/%Y %H:%M:%S %z")
        log("SCIP", f"{session} now in progress; start scheduled at {starttime}.")
        # Keep track of the end time of the session.
        endtime = datetime.strptime(session['end']+" +0000", "%m/%d/%Y %H:%M:%S %z")
        # Create the Observers at local and conjugate locations using the 
        # session parameters found in the schedule file
        lon = float(session['lon'])
        lat = float(session['lat'])
        elev = float(session['elev'])
        conj_lon = float(session['conj_lon'])
        conj_lat = float(session['conj_lat'])
        m.local = Observer(timezone='UTC', longitude=lon*u.deg, latitude=lat*u.deg, elevation=elev*u.m)
        m.conjugate = Observer(timezone='UTC', longitude=conj_lon*u.deg, latitude=conj_lat*u.deg, elevation=0*u.m)
        
        # Find out when the next sunset/dusk/dawn/sunrise occur
        NE = NextEvents(m.local, m.conjugate)      

        # While the end time of the session is in the future:
        while endtime > FAKETIME: #datetime.now(timezone.utc): 
            #Use datetime.now(timezone.utc) for actual observations
            status = NE.sun_status() # check sun status
            print(status)
            if status == 'recalculate':
                # if it is time to recalcuate the next events, do so
                log("SCIP", 'Recalculating next sunset/sunrise times')
                NE = NextEvents(m.local, m.conjugate) 
                status = NE.sun_status() 
            if status == 'daylight':
                # Ensure the photometer is pointed away from the sun
                m.point_away_from_sun()
                # Log ambient conditions
                log("SCIP", f'Ambient: {s.temperature()} C, {s.humidity()}% RH')
                # Wait until sunset or wait 15 min, whichever is shorter
                time_until_sunset = NE.sunset - datetime.now(timezone.utc)
                time_to_wait = min(time_until_sunset.total_seconds(), 15*60)
                time_to_wait_minutes = np.round(time_to_wait/60,2)
                log("SCIP", f'Currently daylight, sunset at {NE.sunset}, waiting {time_to_wait_minutes} minutes.')
                time.sleep(time_to_wait)
            if status == 'twilight':
                log("SCIP", 'Starting twilight observations.')
                while NE.sun_status()=='twilight':
                    # image the airglow for the integration time in the schedule
                    image_airglow(session, b,m,s)
                    log("SCIP", 'Airglow imaged')
                    log("SCIP", f'Ambient: {s.temperature()} C, {s.humidity()}% RH')
            if status == 'night':
                log("SCIP", 'Starting night observations.')
                airglow_count = 0
                while NE.sun_status() == 'night':
                    # image the airglow for the integration time in the schedule
                    image_airglow(session, b,m,s)
                    log("SCIP", 'Airglow imaged')
                    log("SCIP", f'Ambient: {s.temperature()} C, {s.humidity()}% RH')
                    airglow_count += 1
                    if airglow_count%session['airglow images per cal star'] == 0: # cal star every 10 airglow images
                        log("SCIP", 'Now would be the time for a cal star')
                        image_calstar(session, NE, b,m,s)
                        log("SCIP", 'Calibration star imaged')
                        log("SCIP", f'Ambient: {s.temperature()} C, {s.humidity()}% RH')
        
        # If you get here, you are out of the while (endtime>now) loop and the
        # the session has ended.
        return 0