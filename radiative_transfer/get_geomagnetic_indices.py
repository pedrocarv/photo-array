#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
import datetime
import pytplot
from pyspedas.projects.kyoto import dst
import spaceweather as sw
from pymsis import msis
from pymsis import utils
import warnings
warnings.filterwarnings('ignore')

# Define the function to fetch geomagnetic indices -------------------------------------------------------------#
def get_geomagnetic_indices(date):
    '''
    #################################################################################
    
    Function that returns dst, kp, ap, f107, and 107a indices for a user-defined date
    
    Args: 
        date - datetime.datetime object

    Outputs:
         dstd - daily dst index at 1-hour cadence (24 x 1)
         kp - daily kp index at 3-hour cadence (8 x 1)
         aps - daily ap index at 3-hour cadence (8 x 7)
         f107 - daily f107 value (constant)
         f107a - daily 81-day moving average of f107 (constant)
         hour_dst - UT hour for dst index at 1-hour cadence (24 x 1)
         minute_dst - UT minute for dst index at 1-hour cadence (24 x 1)
         hour_kp_aps - UT hour for kp and ap indices at 3-hour cadence (8 x 1)
         minute_kp_aps - UT minute for kp and ap indices at 3-hour cadence (8 x 1)

    Author: Pratik P. Joshi (ppjoshi2@illinois.edu)

    #################################################################################
    '''

    # Extract the day-of-year (doy) associated with the input date
    date_format = "%Y-%m-%d %H:%M:%S"
    date_object = datetime.datetime.strptime(str(date),date_format)
    doy = date_object.timetuple().tm_yday

    # Calculate the next day-of-year (doy) and convert it to a datetime object
    year2=date.year
    if ((year%4>0) & (doy==365)):
        year2=date.year+1
        doy2 = int(1)
    elif ((year%4==0) & (doy==366)):
        year2=date.year+1
        doy2 = int(1)
    else:
        doy2 = doy+1
    date2 = datetime.datetime(year2, 1, 1) + datetime.timedelta(doy2 - 1)
    
    # Get dst index from WDC Kyoto at 1-hour cadence -----------------------------------------------------------------------#
    dst(trange=[str(date),str(date2)])                                    # download the dst database up to the current date
    time, dstd = pytplot.get_data("kyoto_dst")                            # get the hourly dst value for the desired date
    epoch = datetime.datetime(1970, 1, 1, 0, 0)                           # define the 1 Jan 1970 epoch
    seconds_since_epoch = (date - epoch).total_seconds()                  # get the UT time in seconds from 1 Jan 1970 
    hour_dst, minute_dst = [((time-seconds_since_epoch)/3600).astype(int) # get UT hour and minute associated with hourly dst values
                            ,(60*(((time-seconds_since_epoch)/3600)%1)).astype(int)]

    # Get kp index from GFZ GFZ Helmholtz Centre for Geosciences at 3-hour cadence ------------------------------------------#
    df_3h = sw.gfz_3h(update=True)                                         # download the Kp database up to the current date
    kp = df_3h['Kp'].loc[str(date):str(date2)].values                      # get 3-hour Kp value for the desired date
    npdatetime_kp_aps = df_3h.loc[str(date):str(date2)].index.values       # get 3-hour np.datetime64 object for the desired date
    datetime_kp_aps = [datetime.datetime.utcfromtimestamp(t.tolist()/1e9)  # get 3-hour datetime.datetime object for the desired date
                       for t in npdatetime_kp_aps] 
    hour_kp_aps = [i.hour for i in datetime_kp_aps]                        # get hour associated with 3-hour Kp value
    minute_kp_aps = [i.minute for i in datetime_kp_aps]                    # get minutes associated with 3-hour Kp value

    # Get daily f107, f107a at daily cadence from Dominion Radio Astrophysical Observatory and Natural Resources Canada,------#
    # and aps at 3-hour cadence from GFZ Helmholtz Centre for Geosciences
    utils.download_f107_ap()                                               # download the F107 and Ap database up to the current date
    utils._load_f107_ap_data()                                             # load the updated database up to the current date
    f107, f107a, aps = utils.get_f107_ap(datetime_kp_aps)                  # get f107, f107a, and Aps for the desired date

    return dstd, kp, aps, f107[0], f107a[0], hour_dst, minute_dst, hour_kp_aps, minute_kp_aps 
    

# Main Call -------------------------------------------------------------------------------------------------------------------#

## User defined date input as yyyy, mm, dd command line arguments
year = int(sys.argv[1])
month = int(sys.argv[2])
day = int(sys.argv[3])
date = datetime.datetime(year, month, day)
print('User input date = ', date)

## Get geomagnetic indices using a function call
print('====================== Function call =====================================')
dstd, kp, aps, f107, f107a, hour_dst, minute_dst, hour_kp_aps, minute_kp_aps = get_geomagnetic_indices(date)

## Print outputs
print('====================== Outputs =====================================')
print('')
print('---------------------- hourly Dst array ----------------------------')
print('Dst array size = ', dstd.shape)
print('Dst UT hour = ', hour_dst)
print('Dst UT minute = ', minute_dst)
print('Dst = ', dstd)
print('')
print('---------------------- 3-hour Kp array ------------------------------')
print('Kp array size = ', kp.shape)
print('Kp UT hour = ', hour_kp_aps)
print('Kp UT minute = ', minute_kp_aps)
print('Kp = ', kp)
print('')
print('---------------------- 3-hour Ap array ------------------------------')
print('Ap array size = ', aps.shape)
print('Ap UT hour = ', hour_kp_aps)
print('Ap UT minute = ', minute_kp_aps)
print('Ap = ', aps)
print('')
print('---------------------- daily F107 and F107a values ------------------')
print('F107 =', f107)
print('F107a =', f107a)
