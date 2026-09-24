'''
Folder containing functions to capture calibration images

By JC
'''
from .photometer_classes import *
import time

def take_dark(path, exptime, temp, cooldown_wait=5, comment='dark'):
    '''
    Function that will capture dark exposure (to be used to generate dark-rate image)
    :input:
    path (str) - path to save dark exposures to
    exptime (float) - sensor exposure time [in seconds]
    temp (float) - temperature to set CCD to [degrees C]
    cooldown_wait (float) - time to wait after CCDs reach set temp [in minutes]
    comment (str) - comment for exposure fits file
    TODO: Binning support?
    '''

    b = Binocular()

    b.set_cooling(temp, cooldown_wait)

    print('Starting dark image capture in 2 minutes. Run.')
    time.sleep(2*60)

    _,_,_ = b.dual_exposure(exptime, False,False, show = False, save = 'fits', session=False, comment=comment, path=path)

def take_bias(path, N_imgs,  temp, cooldown_wait=5, comment='bias', exptime=0.001):
    '''
    Function that will capture dark exposure (to be used to generate dark-rate image)
    :input:
    path (str) - path to save dark exposures to
    N_imgs (int) - number of bias images to capture
    temp (float) - temperature to set CCD to [degrees C]
    cooldown_wait (float) - time to wait after CCDs reach set temp [in minutes]
    comment (str) - comment for exposure fits file
    exptime (float) - sensor exposure time [in seconds... should always be 0.001 [1ms]]
    TODO: Binning support? Multiple Trials?
    '''

    b = Binocular()
    
    b.set_cooling(temp, cooldown_wait)
    
    for i in range(N_imgs):
        print('Capturing Bias Frame', i)
        _,_,_ = b.dual_exposure(exptime, False,False, show = False, save = 'fits', session=False, comment=comment, path=path)

    # dis = b.disconnect()
    # print(dis)
    # 3/3/25 (JC) - silly to have automatic disconnect, especially if we want to do multiple BF captures in a row

