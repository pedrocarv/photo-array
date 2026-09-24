'''
helpers.py - collection of helper functions to assist with data processing and analysis

By JC
'''

import numpy as np
from astropy.time import Time
import astroplan
import astropy.units as u
from astropy.coordinates import get_sun
from datetime import datetime

from radiative_transfer.rt_msis import SolarActivity

def fitsdate_to_dt(datestr):
    '''
    Converts datestring from fits file to a numpy datetime64 object

    :input:
    datestr (str) - Date string from fits file (Assumed format: "mm/dd/YYYY HH:MM:SS UTC")

    :output:
    numpy datetime64 object
    '''
    return np.datetime64(f"{datestr[6:10]}-{datestr[:2]}-{datestr[3:5]}T{datestr[11:19]}")


def get_shadowalt(dt, loc):
    '''
    Get earth's shadow altitude (in km) from a given location
    Based on Dawn's code from shadow-height.py

    # TODO: Non-zenith shadow altitude

    :input:
    dt (np.datetime64) - Datetime of observation
    loc (list) - location info of observing site [(decimal) longitude, (decimal) latitude, elevation (m)]
    '''
    local = astroplan.Observer(timezone='UTC', longitude=loc[0]*u.deg, latitude=loc[1]*u.deg, elevation=loc[2]*u.m)
            
    a_time = Time(dt, format='datetime64')

    sun_alt = local.sun_altaz(a_time).alt.radian

    R = 6371 # Earth radius (km) # Note, will add screening altitude of 102 to account for )2 photoabsorption
    theta = np.abs(sun_alt)

    shadowh = 0.0
    if sun_alt < 0.0:
        shadowh = (R)*(1/np.cos(theta) - 1.0) + 102/np.cos(theta)

    return shadowh

def iterative_shadowalt(dt, loc, oza, azi):
    '''
    Compute shadow altitude using the iterative method described in
    Chamberlain's 'Physics of the Aurora and Airglow' Ch 10.
    All variable names that don't have immediate physical interpretation
    follow the naming scheme derived there.

    NOTE: DOES NOT ACCOUNT FOR SCREEN HEIGHT OF O2
    '''
    local = astroplan.Observer(timezone='UTC', longitude=loc[0]*u.deg, latitude=loc[1]*u.deg, elevation=loc[2]*u.m)
            
    a_time = Time(dt, format='datetime64')

    sun_alt = local.sun_altaz(a_time).alt.radian

    R = 6371 # Earth radius (km) # Note, will add screening altitude of 102 to account for )2 photoabsorption
    alpha = np.abs(sun_alt) # FIXME: Is this really solar declination (as described in Chamberlain?)

    if sun_alt < 0.0:
        shadowh = 0.0

    # init with zenith approximation
    shadowh = R * ((1/np.cos(alpha)) - 1)
    beta = np.arccos(1/(shadowh/R + 1))

    for i in range(5): # TODO: Validate with closed form solution FIXME: make while statement that terminates when the difference between successive angles is small
        theta_s = np.arcsin(np.cos(beta) * np.sin(oza))
        gamma = oza - theta_s
        beta = np.arcsin(np.cos(gamma)*np.sin(alpha) - np.sin(gamma)*np.cos(alpha)*np.cos(azi))

    shadowh = R * ((1/np.cos(beta)) - 1)

    return shadowh


def get_solalt(dt, loc, unit='rad'):
    '''
    Get the solar altitude angle (radians) from a given location  

    :input:
    dt (np.datetime64) - Datetime of observation
    loc (list) - location info of observing site [(decimal) longitude, (decimal) latitude, elevation (m)]
    unit (str) - unit of angle to return
    '''
    local = astroplan.Observer(timezone='UTC', longitude=loc[0]*u.deg, latitude=loc[1]*u.deg, elevation=loc[2]*u.m)

    a_time = Time(dt, format='datetime64')

    if unit == 'rad':
        alt = local.sun_altaz(a_time).alt.radian
    else:
        alt = local.sun_altaz(a_time).alt.degree

    return alt

def get_solaz(dt, loc, unit='rad'):
    '''
    Get the solar azimuthal angle (radians) from a given location  

    :input:
    dt (np.datetime64) - Datetime of observation
    loc (list) - location info of observing site [(decimal) longitude, (decimal) latitude, elevation (m)]
    unit (str) - unit of angle to return
    '''
    local = astroplan.Observer(timezone='UTC', longitude=loc[0]*u.deg, latitude=loc[1]*u.deg, elevation=loc[2]*u.m)
            
    a_time = Time(dt, format='datetime64')

    if unit == 'rad':
        az = local.sun_altaz(a_time).az.radian
    else:
        az = local.sun_altaz(a_time).az.degree

    return az
        

def get_sza(dt, loc, unit='rad'):
    '''
    Get the solar zenith angle (radians) from a given location  

    :input:
    dt (np.datetime64) - Datetime of observation
    loc (list) - location info of observing site [(decimal) longitude, (decimal) latitude, elevation (m)]
    unit (str) - unit of angle to return
    '''
    local = astroplan.Observer(timezone='UTC', longitude=loc[0]*u.deg, latitude=loc[1]*u.deg, elevation=loc[2]*u.m)
            
    a_time = Time(dt, format='datetime64')

    #solar zenith angle is the complement of the solar altitude angle -AH
    if unit == 'rad':
        solzen = ((np.pi)/2) - local.sun_altaz(a_time).alt.radian
    else:
        # degrees
        solzen = 90.0 - local.sun_altaz(a_time).alt.degree

    return solzen

def get_solHA(dt, loc, unit='rad', wrap_180=False):
    '''
    Get solar hour angle using Astroplan Observer
    For alternative verison see https://gml.noaa.gov/grad/solcalc/solareqns.PDF

    :input:
    dt (np.datetime64) - Date adn Time of observation
    loc (list) - location info of observing site [(decimal) longitude, (decimal) latitude, elevation (m)]
    unit (str) - Unit to return angle in
    wrap_180 (bool) - wrap HA around 180 (display from -180 to 180 [-pi to pi] instead of 0 to 360 [0 to 2pi])
    '''
    local = astroplan.Observer(timezone='UTC', longitude=loc[0]*u.deg, latitude=loc[1]*u.deg, elevation=loc[2]*u.m)
            
    a_time = Time(dt, format='datetime64')

    sol = get_sun(a_time)

    ang = local.target_hour_angle(a_time, sol)

    if wrap_180:
        ang.wrap_at('180d', inplace=True)

    if unit == 'rad':
        ha = ang.radian
    else:
        ha = ang.degree

    return ha


def get_loc_info(name):
    '''
    Gets the location info (lon, lat, and elevation) for a known observing site
    
    :input:
    name (str) - Name of observing location
    
    :output:
    lon (float) - [decimal] longitude of observing location
    lat (float) - [decimal] latitude of observing location
    elev (float) - elevation of observing location [m]
    '''
    match name:
        case 'eceb':
            lon = -88.228
            lat = 40.115
            elev = 224
        case 'aeronomy':
            lon = -88.159
            lat = 40.167
            elev = 219
        case 'arecibo':
            lon = -66.75
            lat = 18.35
            elev = 317
        case 'pbo':
            lon = -89.6717
            lat = 43.0777
            elev = 362
        case 'millstone':
            lon = -71.495
            lat = 42.61
            elev = 131
        case 'kitt':
            lon = -111.5967
            lat = 31.9583
            elev = 2096
        case 'hao':
            # FIXME: Update with actual site coords
            lon = -105.2457744
            lat = 40.0326591
            elev = 1614
        case _:
            raise Exception("Invalid location, talk to Jackson to get your observatory added")
        
    return lon, lat, elev


def fov_to_pix(fov, m, n):
    '''
    Convert FOV angle (in degrees) to a CCD pixel count

    :input:
    fov (float) - FOV angle (deg.)
    m (int) - pixels in y direction
    n (int) - pixels in x direction
    :output:
    (int) - # of pixels in width that corresponds to that angle
    '''
    # NOTE: Following calculation for deg/pixel is based on the PATHS photometer FOV (~2.3x2.8deg)
    # NOTE: Deg/pix should be about the same in both directions... will verify that here.
    fov_y = 2.287247518
    fov_x = 2.860205313

    deg_per_pix =  fov_y / m

    # Check that degree per pixel is roughly the same in each direction (at least up to 10^-3)
    assert np.abs(deg_per_pix - fov_x/n) / deg_per_pix < 1e-3, 'FOV/pix must be the same up to 10^-3'

    return (int)(fov/deg_per_pix)

def geomagnetic_data(date=np.datetime64("2025-09-24")):
    """
    Written by Michal Ondrejcek

    Get geomagnetic data from
    1) GFZ Potsdam at https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt,
        no F10.7A, an 81-day average
        # Short file description (for a detailed file description, see Kp_ap_Ap_SN_F107_format.txt):
        # 40 header lines, all starting with #
        # ASCII, blank separated and fixed length, missing data indicated by -1.000 for Kp, -1 for ap and SN,
            -1.0 for F10.7
        # YYYY MM DD is date of UT day, days is days since 1932-01-01 00:00 UT to start of UT day, days_m is days
            since 1932-01-01 00:00 UT to midday of UT day
        # BSR is Bartels solar rotation number, dB is day within BSR
        # Kp1 to Kp8 (Kp for the eight eighth of the UT day), ap1 to ap8 (ap for the eight eighth of the UT day),
            Ap, SN, F10.7obs, F10.7adj
        # D indicates if the Kp and SN values are definitive or preliminary. D=0: Kp and SN preliminary; D=1:
            Kp definitive, SN preliminary; D=2 Kp and SN definitive
        #
        #
        # The format for each line is (i stands for integer, f for float):
        #iii ii ii iiiii fffff.f iiii ii ff.fff ff.fff ff.fff ff.fff ff.fff ff.fff ff.fff ff.fff iiii iiii iiii
            iiii iiii iiii iiii iiii  iiii iii ffffff.f ffffff.f i
        # The parameters in each line are:
        #YYY MM DD  days  days_m  Bsr dB    Kp1    Kp2    Kp3    Kp4    Kp5    Kp6    Kp7    Kp8
            ap1  ap2  ap3  ap4  ap5  ap6  ap7  ap8    Ap  SN F10.7obs F10.7adj D
        2025 07 06 34155 34155.5 2617 10  3.667  5.000  3.000  3.333  3.000  1.667  3.000  5.000
            22   48   15   18   15    6   15   48    23  87    118.2    122.2 1
        2025 09 24 34235 34235.5 2620  9  1.667  1.333  1.000  1.667  2.333  2.000  1.333  2.333
            6    5    4    6    9    7    5    9     6 144    184.7    185.8 1
        2025 11 18 34290 34290.5 2622 10  2.667  2.000  2.000  0.667  0.333  0.000  0.333  0.667
            12    7    7    3    2    0    2    3     4  63    120.1    117.4 0

    2) fluxdate at https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt,
        no F10.7A, an 81-day average
        fluxdate    fluxtime    fluxjulian    fluxcarrington  fluxobsflux  fluxadjflux  fluxursi
        ----------  ----------  ------------  --------------  -----------  -----------  ----------
        20250924    170000      2460943.197   2302.56         0177.2       0178.3       0160.5
        20250924    200000      2460943.322   2302.56         0184.7       0185.8       0167.2
        20250924    230000      2460943.447   2302.57         0178.4       0179.5       0161.5
        20251118    200000      2460998.322   2304.58         0120.1       0117.4       0105.6

    3) CelesTrak at https://celestrak.org/SpaceData/SW-Last5Years.txt
        DATATYPE CssiSpaceWeather
        VERSION 1.2
        UPDATED 2025 Dec 08 20:32:34 UTC
        # --------------------------------------------------------------------------------------------------------------------------------
        #                              SPACE WEATHER DATA
        # --------------------------------------------------------------------------------------------------------------------------------
        #
        # See https://celestrak.org/SpaceData/SpaceWx-format.php for format details.
        #
        # FORMAT(I4,I3,I3,I5,I3,8I3,I4,8I4,I4,F4.1,I2,I4,F6.1,I2,5F6.1)
        # -------------------------------------------------------------------------------------------
        #
        # yy mm dd BSRN ND Kp Kp Kp Kp Kp Kp Kp Kp Sum Ap  Ap  Ap  Ap  Ap  Ap  Ap  Ap  Avg Cp C9 ISN
        # --------------------------------------------------------------------------------------------------------------------------------
        -------------------------------------
        Adj     Adj   Adj   Obs   Obs   Obs
        F10.7 Q Ctr81 Lst81 F10.7 Ctr81 Lst81
        -------------------------------------
        #
        NUM_OBSERVED_POINTS 2168
        BEGIN OBSERVED
        2025 09 24 2620  9 17 13 10 17 23 20 13 23 137   6   5   4   6   9   7   5   9   6 0.3 1 144
        185.8 0 151.6 153.6 184.7 150.7 150.1
        2025 11 18 2622 10 23 20 20  3  0  0  0  7  73   9   7   7   2   0   0   0   3   4 0.1 0  62
        117.4 0 147.1 150.3 120.1 150.4 150.5

        Args:
            date (datetime64): A vantage date.

        Return:
            solar_activity: obj

    """
    # datetime object
    dt = date.astype(datetime)

    source = "Potsdam" # Potsdam, Canada, CelesTrak
    keys = None
    date_param = {}
    try:
        import requests

        if source == "Potsdam":
            gfz_url = "https://www-app3.gfz-potsdam.de/kp_index/Kp_ap_Ap_SN_F107_since_1932.txt"

            response = requests.get(gfz_url)
            gfz_data = response.text
            gfz_data = gfz_data.splitlines()

            avg81 = []
            for line in gfz_data:
                if line[0] == "#":
                    if line[1:4] == "YYY":
                        keys = line[1:].split()
                else:
                    data = line.split()
                    if int(data[0]) >= 2020:
                        avg81.append(float(data[-2]))
                    if int(data[0]) == dt.year and int(data[1]) == dt.month and int(data[2]) == dt.day:
                        for i, item in enumerate(data):
                            if i == 4:
                                date_param[keys[i]] = float(item)
                            elif i >= 7 and i <= 14:
                                date_param[keys[i]] = float(item)
                            elif i >= 25 and i <= 26:
                                date_param[keys[i]] = float(item)
                            else:
                                date_param[keys[i]] = int(item)
                        avg81 = np.asarray(avg81)[-81:]
                        break
            date_param["F10.7adj_Lst81"] = np.mean(avg81)
        elif source == "Canada":
            gfz_url = "https://www.spaceweather.gc.ca/solar_flux_data/daily_flux_values/fluxtable.txt"

            response = requests.get(gfz_url)
            gfz_data = response.text
            gfz_data = gfz_data.splitlines()

            avg81 = []
            for line in gfz_data:
                if line.startswith("fluxdate") or line.startswith("------"):
                    continue
                else:
                    data = line.split()
                    # data[1]) == 200000 is 8pm
                    if int(data[0][0:4]) >= 2020 and int(data[1]) == 200000:
                        avg81.append(float(data[-2]))
                    if (int(data[0][0:4]) == dt.year and int(data[0][4:6]) == dt.month and int(data[0][6:8]) == dt.day
                            and int(data[1]) == 200000):
                        for i, item in enumerate(data):
                            if i == 5:
                                date_param["F10.7adj"] = float(item)
                        avg81 = np.asarray(avg81)[-81:]
                        break
            date_param["Ap"] = 4
            date_param["F10.7adj_Lst81"] = np.mean(avg81)
        elif source == "CelesTrak":
            gfz_url = "https://celestrak.org/SpaceData/SW-Last5Years.txt"

            response = requests.get(gfz_url)
            gfz_data = response.text
            gfz_data = gfz_data.splitlines()

            for line in gfz_data:
                if (line.startswith("#") or line.startswith("DATATYPE") or line.startswith("VERSION") or
                        line.startswith("VERSION") or line.startswith("UPDATED") or
                        line.startswith("NUM_OBSERVED_POINTS") or line.startswith("BEGIN OBSERVED")):
                    if line.startswith("# yy mm dd"):
                        keys = line.split()
                        keys = keys[1:]

                        keys[-7] = "F10.7adj"
                        keys[-6] = "Q"
                        keys[-5] = "F10.7adj_Ctr81"
                        keys[-4] = "F10.7adj_Lst81"
                        keys[-3] = "F10.7obs"
                        keys[-2] = "F10.7obs_Ctr81"
                        keys[-1] = "F10.7obs_Lst81"
                else:
                    data = line.split()
                    if int(data[0]) == dt.year and int(data[1]) == dt.month and int(data[2]) == dt.day:
                        for i, item in enumerate(data):
                            if i == 23 or i >= 26:
                                date_param[keys[i]] = float(item)
                            else:
                                date_param[keys[i]] = int(item)
                        break
            date_param["Ap"] = date_param["Avg"]
    except:
        # launch date
        date = np.datetime64("2025-09-24")
        date_param = {"YYY": 2025, "MM": 9, "DD": 24, "days": 34235, "days_m": 34235.5,
                      "Bsr": 2620, "dB": 9,
                      "Kp1": 1.667, "Kp2": 1.333, "Kp3": 1.000, "Kp4": 1.667,
                      "Kp5": 2.333, "Kp6": 2.000, "Kp7": 1.333, "Kp8": 2.333,
                      "ap1": 6, "ap2": 5, "ap3": 4, "ap4": 6, "ap5": 9, "ap6": 7, "ap7": 5, "ap8": 9,
                      "Ap": 6, "SN": 144,
                      "F10.7obs": 184.7, "F10.7adj": 185.8, "D": 1}
        date_param["F10.7adj_Lst81"] = 153.6 # 81 days: 2025-07-06 - 2025-09-24

    msis_date = date
    aps = date_param["Ap"]

    f107 = date_param["F10.7adj"]
    f107a = date_param["F10.7adj_Lst81"]

    aps = aps * np.ones((1, 7))
    solar_activity = SolarActivity(msis_date, aps, f107, f107a)

    return solar_activity