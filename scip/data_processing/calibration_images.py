'''
Functions for generating calibration images
JC
Created: 01.22.2025
'''
import numpy as np
from .read_fits_files import *
from astropy.io import fits
import os
from ccdproc import combine
from scip.SCIP_operation.util import show_imgs
import pathlib


def generate_dark(path, inttime, save_array=False, arr_loc='', sub_Bias = False):
    '''
    Generate dark images from passed observations
    :input:
    path (str) - Path to directory containing images
    inttime (int) - exposure time (minutes)
    save_array (Bool) - save dark image to .npy file
    :output:
    (matplotlib figure) Dark Images
    '''

    # NOTE: m_on == c1; m_off == c2???
    m_on, m_off = get_fits_from_folder(path)
    # # m1_on, m1_off = get_fits_from_folder(path1)

    assert len(m_on) == len(m_off), 'Sensors should have same number of observations!'
    N = len(m_on) # Number of observations in stack

    inttime_s = inttime * 60 # convert exposure time to seconds
    
    # should only have one file per camera
    c1_dark = fits.open(m_on[0]['directory']+'/'+m_on[0]['filename'], ignore_missing_end=True)[0].data.astype(float)
    c2_dark = fits.open(m_off[0]['directory']+'/'+m_off[0]['filename'], ignore_missing_end=True)[0].data.astype(float)

    #FIXME: take look after running negative investigation - AH 4-1-25
    if (sub_Bias): #subtract bias if true
        # load bias frames NOTE: REFACTOR (THESE LINES SUCK!!!)
        b1 = np.load('temp_sweep_frames_100/bias_-7_c1.npy') 
        b2 = np.load('temp_sweep_frames_100/bias_-7_c2.npy')

        c1_dark -= b1
        c2_dark -= b2

    c1_dark /= inttime_s
    c2_dark /= inttime_s 

    # breakpoint()
    fig = show_imgs(c1_dark, c2_dark, stitle='Dark Signal Rate [ADU/s]') # FIXME: switch back after comparison

    if save_array:
        np.save('c1_dark_-7.np', c1_dark)
        np.save('c2_dark_-7.np', c2_dark)

    return fig

def generate_dark_signal(path, save_array=False, arr_loc='', sub_Bias=False):
    #generate dark signal frame (not current)- AH 4/1/25
    # NOTE: m_on == c1; m_off == c2
    m_on, m_off = get_fits_from_folder(path)

    assert len(m_on) == len(m_off), 'Sensors should have same number of observations!'
    N = len(m_on) # Number of observations in stack
    
    # should only have one file per camera
    c1_dark = fits.open(m_on[0]['directory']+'/'+m_on[0]['filename'], ignore_missing_end=True)[0].data.astype(float)
    c2_dark = fits.open(m_off[0]['directory']+'/'+m_off[0]['filename'], ignore_missing_end=True)[0].data.astype(float)

    # NOTE: testing if bias is the cause of negative values
    #FIXME: take look after running negative investigation - AH 4-1-25
    if(sub_Bias):
        b1 = np.load('temp_sweep_frames_100/bias_-7_c1.npy') 
        b2 = np.load('temp_sweep_frames_100/bias_-7_c2.npy')

        c1_dark -= b1
        c2_dark -= b2

    #create plot of dark signal
    fig = show_imgs(c1_dark, c2_dark, stitle='Dark Signal [ADU]') # FIXME: switch back after comparison

    if save_array:

        np.save(f'{arr_loc}_c1', c1_dark)
        np.save(f'{arr_loc}_c2', c2_dark)

    return fig

def generate_bias(path, save_array=False, fname=None):
    '''
    Generate bias images from passed observations
    :input:
    path (str) - Path to directory containing images
    save_array (bool) - save average bias image as .npy file
    fname (pathlib.Path) - filename of bias array to be saved [NOTE: Must include directory of file (if not in current)]
    :output:
    (matplotlib figure) Bias Images
    '''
    m_on, m_off = get_fits_from_folder(path)

    assert len(m_on) == len(m_off), 'Sensors should have same number of observations!'
    N = len(m_on) # Number of observations in stack
    
    # NOTE: image stack too large for mem to save all images to array, must compute average cumulatively
    # FIXME: Yuck!
    c1_bias = np.zeros(fits.open(m_on[0]['directory']+'/'+m_on[0]['filename'], ignore_missing_end=True)[0].data.astype(float).shape)
    c2_bias = np.zeros(fits.open(m_off[0]['directory']+'/'+m_off[0]['filename'], ignore_missing_end=True)[0].data.astype(float).shape)
    
    for i in range(N):
        c1_img = fits.open(m_on[i]['directory']+'/'+m_on[i]['filename'], ignore_missing_end=True)[0].data.astype(float)
        c2_img = fits.open(m_off[i]['directory']+'/'+m_off[i]['filename'], ignore_missing_end=True)[0].data.astype(float)

        c1_bias += c1_img
        c2_bias += c2_img


    c1_bias /= N
    c2_bias /= N 

    fig = show_imgs(c1_bias, c2_bias) # FIXME: switch back after comparison

    if save_array:
        np.save(str(fname) + '_c1', c1_bias)
        np.save(str(fname) + '_c2', c2_bias)

    return fig

def gen_bias(path_c1, path_c2):
    '''
    Generate bias images using ccd proc
    NOTE: Not too useful... maybe if we decide to use median combined bias frames come back to this?
    '''
    fnames = ['', '']

    for file in os.listdir(path_c1):
        filename = os.fsdecode(file)
        fnames[0] += path_c1+'/'+filename+','
    fnames[0] = fnames[0][:-1]

    for file in os.listdir(path_c2):
        filename = os.fsdecode(file)
        fnames[1] += path_c2+'/'+filename+','
    fnames[1] = fnames[1][:-1]

    combined = [combine(l, method='median', unit='adu') for l in fnames]  

    print(type(combined[0]), type(combined[1]))

    fig = show_imgs(np.asarray(combined[0]), np.asarray(combined[1]))
    return fig
