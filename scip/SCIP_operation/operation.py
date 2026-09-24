'''
File containing all photometer operation classes (e.g. mount, cams, etc.)

By JC
'''

import AtikSDK
from scip.SCIP_operation.util import save_as_fits, show_imgs, targ_from_name
from scip.data_processing.helpers import get_loc_info

import subprocess
import numpy as np
import matplotlib.pyplot as plt
import time
from pathlib import Path
from datetime import datetime, timezone

from astropy.time import Time
import astropy.units as u
from astroplan import Observer

class Mount():
    '''
    Class containing wrapper functions to control iOptron CubeMount pro
    via the rotctl bash command (https://man.archlinux.org/man/rotctl.1.en)
    If we get a different mount brand, will update this to point different mounts (but i refuse to be proactive)
    TODO: Implement with direct serial comm.
    '''
    
    def __init__(self, dev=Path('/dev/ttyUSB0'), loc=None):
        '''
        Initialization

        :input:
        dev_path (pathlib.Path) - path to mount serial connection FIXME: auto detect ioptron mounts
        TODO: 
            > RA dec to alt az (rotctl only takes alt az as input, will track fixed RADec)
            > easy calstar pointing
        '''
        self.dev = dev
        self.model = '1901'
        self.baud = '9600'
        # Build start of command array compatible with subprocess lib
        self.command_start = ['rotctl', '-r', self.dev, '-m', self.model, '-s', self.baud] 

        self.alt = 0
        self.az = 0
        self.get_pointing() # If no error, then we're good to go, inits alt/az
        
        self.alt_offset = 0
        self.az_offset = 0

        self.obs = None
        
        if loc is not None:
            self.set_location(loc)

        print('mount ready')

    def get_pointing(self, ret=False, radec=False):
        '''
        Get current pointing of photometer

        :input:
        ret (bool) - flag to return values or not FIXME: May be weird to have
        TODO: radec (bool) - return pointing wrt RA/Declination (default alt/az)
        :output:
        alt (float) - altitude of current pointing
        az (float) - azimuth of current pointing
        '''
        proc = subprocess.run(self.command_start + ['p'], capture_output=True)

        self.az, self.alt = [float(var) for var in proc.stdout.decode('utf-8').split('\n')[:-1]]

        if ret:
            if radec:
                raise NotImplementedError('yell at jackson to do this')
            else:
                return self.alt, self.az
    
    def goto(self, coord1, coord2, track=False, radec=False, apply_offset=True):
        '''
        Point at specified coordinates

        :input:
        coord1 (float) - first coord of point (alt/ra)
        coord2 (float) - second coord of point (az/dec)
        track (bool) - continue tracking point (keep radec constant)
        radec (bool) - flag for coords specifying radec
        '''
        if radec:
            # Convert coords back to alt/az
            raise NotImplementedError('yell at jackson to do this')
        
        if apply_offset:
            coord1 += self.alt_offset
            coord2 += self.az_offset
        
        subprocess.run(self.command_start + ['P', str(coord2), str(coord1)], capture_output=True)
        # point.wait() # wait for pointing to finish before continuing

        # FIXME: inefficient(?) way to create block for rotor pointing
        self.get_pointing()
        while np.abs(coord1 - self.alt) > 0.5 or np.abs(coord2 - self.az) > 0.5:
            self.get_pointing()

        if not track:
            self.stop()

        # confirm pointing location
        self.get_pointing()
        print('Staring at %f, %f'%(self.alt, self.az))

    def set_offset(self, target_name):
        '''
        Sets alt/az offset based on current pointing
        NOTE: CAMERA HAS TO BE CENTERED ON PASSED STAR NAME WHEN RAN
        '''
        assert self.obs is not None, 'Need to specify observing location to perform alignment!'

        targ = targ_from_name(target_name)

        self.get_pointing()
        targ_altaz = self.obs.altaz(Time.now(), target=targ)

        self.alt_offset = targ_altaz.alt - self.alt
        self.az_offset = targ_altaz.az - self.az


    def set_location(self, loc):
        '''
        set location based on passed name
        FIXME: Allow coords?
        '''
        loc_info = get_loc_info(loc)
        self.obs = Observer(timezone='UTC', longitude=loc_info[0]*u.deg, latitude=loc_info[1]*u.deg, elevation=loc_info[2]*u.m)

    def stop(self):
        '''
        Stop current tracking
        '''
        subprocess.run(self.command_start + ['S'])


class Binocular():
    '''
    Binocular class rewritten with the AtikSDK python wrapper (flexible for any OS)
    Install lib from the .whl file (Windows, also run the setup .exe)

    Most of these methods are yoinked directly from Dawn Haken's original Binocular class
    (photometer_classes.py)... shoutout Dawn
    '''
    def __init__(self):
        '''
        Initialize (TODO: determine on/off band by serial number)
        '''
        self.c1 = AtikSDK.AtikSDKCamera()
        self.c2 = AtikSDK.AtikSDKCamera()
        
        assert self.c1.is_device_present(0) and self.c2.is_device_present(1), 'Make sure both cameras are connected!'

        self.c1.connect(0)
        self.c2.connect(1)

    def set_cooling(self, setpoint, wait=False):
        '''
        Set CCD temp. in degrees C

        :input:
        setpoint (float) - setpoint in degrees C
        '''
        self.c1.set_cooling(setpoint)
        self.c2.set_cooling(setpoint)

        if wait:
            temp1 = 0
            temp2 = 0

            while temp1 > setpoint+0.2 and temp2 > setpoint+0.2:
                temp1, temp2 = self.get_CCD_temps(ret=True)
                time.sleep(30)

    def get_CCD_temps(self, ret=False):
        '''
        Get current CCD temp and setpoint

        :input:
        ret (bool) - return camera temps
        '''
        c1temp = self.c1.get_temperature()
        c2temp = self.c2.get_temperature()
        c1set = self.c1.cooling_info()[4]
        c2set = self.c2.cooling_info()[4]

        print(f'C1: {c1temp}°C [set to {c1set}°C]\nC2: {c2temp}°C [set to {c2set}°C]')

        if ret:
            return c1temp, c2temp
    
    def dual_exposure(self, exptime, m=None, s=None, show = False, binning = 1, 
                      save = 'fits', target = 'test', session=False, comment='None',path=''):
        '''
        Take exposure with both cameras at once

        :input:
        exptime (float) - length of exposure (in s)
        m (Mount object) - mount object detailing current pointing info
        s (TempSensor object) - External temp sensor (FIXME: not implemented or used)
        show (bool) - plot exposure once completed
        binning (int) - on-chip binning size to use
        save (str) - datatype to save as (None to not save)
        target (str) - imaging target (for fits files)
        session (bool) - observation session?
        comment (str) - comment (for fits files)
        '''
        self.c1.set_binning(binning, binning)
        self.c2.set_binning(binning, binning)
        # log("SCIP", f"C1: {self.c1.state()}; C2: {self.c2.state()}.")
        # log("SCIP", f"Exposing both cameras for {exptime} seconds.")
        exp_starttime = datetime.now(timezone.utc)
        self.c1.start_exposure(exptime)
        self.c2.start_exposure(exptime)
        while not (self.c1.image_ready() and self.c2.image_ready()):
            #log("SCIP", f"C1: {self.c1.state()}; C2: {self.c2.state()}")
            print(self.c1.camera_state(),self.c2.camera_state())
            time.sleep(3)
        img1 = self.c1.get_image()
        img2 = self.c2.get_image()
        
        if save == 'fits':
            save_as_fits(img1, 'c1', exp_starttime, exptime, binning, self,m,s, target=target, session=session, comment=comment, path=path)
            save_as_fits(img2, 'c2', exp_starttime, exptime, binning, self,m,s, target=target, session=session, comment=comment, path=path)
            
        if save == 'numpy':  
            img1_name = "c1" + exp_starttime.strftime("_%d_%b_%Y_%H-%M-%S")
            np.save(f'{path}/{img1_name}', img1)
            img2_name = "c2" + exp_starttime.strftime("_%d_%b_%Y_%H-%M-%S")
            np.save(f'{path}/{img2_name}', img2)
        
        if show == True:
            show_imgs(img1, img2, star_detect = False)
            plt.show()
            
        return img1, img2, exp_starttime   
    
    def warm_up(self):
        '''
        warm up CCDs to prevent CCD shock
        '''
        self.c1.cooler_warmup()
        self.c2.cooler_warmup()

    def is_connected(self):
        '''
        Check if both cameras are connected
        '''
        return self.c1.is_connected(), self.c2.is_connected()
    
    def disconnect(self):
        '''
        (Safely) disconnect the cameras after a session
        '''
        self.warm_up()
        return self.c1.disconnect(), self.c2.disconnect()
    