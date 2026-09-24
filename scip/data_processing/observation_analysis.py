'''
Script housing functions to analyze nights worth of observations

By JC
'''
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.dates import DateFormatter

from scip.data_processing.helpers import *

def show_mean_on_off_exp(on_band, off_band, start_ind):
    '''
    Basic function to plot and show on band and off band measurements for a given observation set

    :input: 
    on_band (dict of FITS data) - on band images loaded from folder
    off_band (dict of FITS data) - off band images loaded from folder
    start_ind (int) - index of the start of observations (NOTE: often, a test exposure is taken as image 0 for a night)
    '''

    on_er = [np.mean(on_band[i]['image'])/on_band[i]['EXP-TIME'] for i in range(start_ind,len(on_band))]
    off_er = [np.mean(off_band[i]['image'])/off_band[i]['EXP-TIME'] for i in range(start_ind,len(off_band))]

    on_date = [fitsdate_to_dt(on_band[i]['DATE-OBS']) for i in range(start_ind,len(on_band))]
    off_date = [fitsdate_to_dt(off_band[i]['DATE-OBS']) for i in range(start_ind,len(off_band))]
    
    on_salt = [get_shadowalt(fitsdate_to_dt(on_band[i]['DATE-OBS']),'eceb') for i in range(start_ind,len(on_band))]
    off_salt = [get_shadowalt(fitsdate_to_dt(off_band[i]['DATE-OBS']), 'eceb') for i in range(start_ind,len(off_band))]
    

    on_exp = [on_band[i]['EXP-TIME'] for i in range(start_ind,len(on_band))]
    off_exp = [off_band[i]['EXP-TIME'] for i in range(start_ind,len(off_band))]

    plt.figure()
    # plt.plot(on_date, on_er, 'b', label='on-band')
    # plt.plot(off_date, off_er, 'r--', label='off-band')
    plt.plot(on_salt, on_er, 'b', label='on-band')
    plt.plot(off_salt, off_er, 'r--', label='off-band')
    # ax = plt.gca()
    # ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    plt.xticks(rotation=-45, ha='left', size='small')
    plt.grid()
    plt.ylabel('Mean Event Rate [ADU/s]')
    # plt.xlabel('Time of Observation [UTC]')
    plt.xlabel('Shadow Altitude [km]')
    plt.legend()
    plt.title('Ha Observations %s'%(on_band[0]['DATE-OBS'][:10]))

    plt.figure()
    plt.plot(on_date, on_exp, 'b', label='on-band')
    plt.plot(off_date, off_exp, 'r--', label='off-band')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    plt.xticks(rotation=-45, ha='left', size='small')
    plt.grid()
    plt.ylabel('Exposure Time [s]')
    plt.xlabel('Time of Observation [UTC]')
    plt.legend()
    plt.title('Ha Observation exposure times %s'%(on_band[0]['DATE-OBS'][:10]))

    plt.show()

def plot_mean_eventrate(on_band, off_band, start_ind, isolate_onband=False):
    '''
    Plot the mean eventrate for a given set of on/off band images
    NOTE: For now, perform mean and dark subtraction automatically

    :input: 
    on_band (dict of FITS data) - on band images loaded from folder
    off_band (dict of FITS data) - off band images loaded from folder
    start_ind (int) - index of the start of observations (NOTE: often, a test exposure is taken as image 0 for a night)

    :output:
    (matplotlib Figure) - Figure containing generated plot
    '''
    # FIXME: Standard bias and dark frames should be added to the repo, change these paths
    bias_on = np.load('/home/jc/research/photometer/data_analysis/bias_temp_sweeps_03272025_100ims/bias_-5_c1.npy')
    bias_off = np.load('/home/jc/research/photometer/data_analysis/bias_temp_sweeps_03272025_100ims/bias_-5_c2.npy')
    dr_on = np.load('/home/jc/research/photometer/scip-2/scip/data_processing/cal_frames_lisa/dark_rate_-6_c1.npy')
    dr_off = np.load('/home/jc/research/photometer/scip-2/scip/data_processing/cal_frames_lisa/dark_rate_-6_c2.npy')

    #FIXME: Bias and dark subtraction should happen in their own function where you toggle bias/dark/flat/etc.
    on_er = [np.mean(on_band[i]['image']-bias_on-(dr_on*on_band[i]['EXP-TIME']))/on_band[i]['EXP-TIME'] for i in range(start_ind,len(on_band))]
    off_er = [np.mean(off_band[i]['image']-bias_off-(dr_off*off_band[i]['EXP-TIME']))/off_band[i]['EXP-TIME'] for i in range(start_ind,len(off_band))]

    on_salt = [get_shadowalt(fitsdate_to_dt(on_band[i]['DATE-OBS']),'eceb') for i in range(start_ind,len(on_band))]
    off_salt = [get_shadowalt(fitsdate_to_dt(off_band[i]['DATE-OBS']), 'eceb') for i in range(start_ind,len(off_band))]

    plt.figure()
    
    if isolate_onband:
        isolated = np.array(on_er) - np.array(off_er)
        plt.plot(on_salt, isolated, label='on-band - off-band')
    else:
        plt.plot(on_salt, on_er, 'b', label='on-band')
        plt.plot(off_salt, off_er, 'r--', label='off-band')

    plt.xticks(rotation=-45, ha='left', size='small')
    plt.grid()
    plt.ylabel('Mean Event Rate [ADU/s]', fontsize=16)
    plt.xlabel('Shadow Altitude [km]', fontsize=16)
    plt.legend()
    plt.title('Ha Observations %s'%(on_band[0]['DATE-OBS'][:10]), fontsize=20, weight='heavy')
    plt.tight_layout()

    return plt.gcf()

def plot_exptime(fits, start_ind):
    '''
    Plot the exposure time given a collection of fits files from an observation session

    :input:
    fits (list of dict) - list of fits files containing observation information
    start_ind (int) - index of first observation

    :output:
    (matplotlib Figure) - Figure containing plot
    '''


    date = [fitsdate_to_dt(fits[i]['DATE-OBS']) for i in range(start_ind,len(fits))]
    
    exptime = [fits[i]['EXP-TIME'] for i in range(start_ind,len(fits))]


    plt.figure()
    plt.plot(date, exptime, 'b', label='on-band')
    ax = plt.gca()
    ax.xaxis.set_major_formatter(DateFormatter('%H:%M:%S'))
    plt.xticks(rotation=-45, ha='left', size='small')
    plt.grid()
    plt.ylabel('Exposure Time [s]')
    plt.xlabel('Time of Observation [UTC]')
    plt.title('Ha Observation exposure times %s'%(fits[0]['DATE-OBS'][:10]))

    return plt.gcf()