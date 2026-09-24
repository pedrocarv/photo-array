import win32com.client
import numpy as np
import time
from util import save_as_fits, show_imgs, log
from datetime import datetime, timezone
import serial
from ctypes import cdll, c_void_p, c_bool, c_int, c_uint16, byref, POINTER


class ATIKCamera(object):
    def __init__(self, cam_number, dllpath = 'AtikCameras.dll'):
        # The dll file contains all the functions to communicate with the camera.
        # Anything returned from an atik.ArtemisXxxx function is a ctypes object.
        # Generally using .value to access the value of the returned thing.
        # Ctypes assumes that all functions return int. So if it actually
        # returns something else, specify with restype.
        self.atik = cdll.LoadLibrary(dllpath) # load dll
        self.atik.ArtemisConnect.restype = c_void_p # set return type of ArtemisConnect
        # Handle tells the dll functions which camera to talk to.
        if cam_number == 1:
            self.handle = self.atik.ArtemisConnect(0)
        elif cam_number == 2:
            self.handle = self.atik.ArtemisConnect(1)
        else: log("SCIP", 'Enter camera number 1 or 2.')
        
    def state(self):
        states = {
            0: "Idle",
            1: "Waiting",
            2: "Exposing",
            4: "Downloading",
            5: "Flushing",
            3: "Error",
        }
        current_state = self.atik.ArtemisCameraState(self.handle)
        exp_time_remain = self.atik.ArtemisExposureTimeRemaining(self.handle) #seconds
        download_percent = self.atik.ArtemisDownloadPercent(self.handle)
        
        if current_state == 2:
            return states[current_state] + f', {exp_time_remain} seconds remaining'
        elif current_state == 4:
            return states[current_state] + f', {download_percent}% downloaded'
        else:    
            return states[current_state]
    
    def start_exposure(self, exptime): # exptime in seconds
        print(exptime)
        exptimeMS = int(np.round(exptime,3)*1000) # convert to integer milliseconds
        self.atik.ArtemisStartExposureMS(self.handle, exptimeMS)
    
    def image_ready(self):
        self.atik.ArtemisImageReady.restype = c_bool # set return type
        ready = self.atik.ArtemisImageReady(self.handle)
        return ready

    def get_image(self):
        x_start_position = c_int()
        y_start_position = c_int()
        w = c_int()
        h = c_int()
        xbin = c_int()
        ybin = c_int()
        self.atik.ArtemisGetImageData(self.handle, byref(x_start_position), byref(y_start_position), byref(w),
                    byref(h), byref(xbin), byref(ybin))
    
        self.atik.ArtemisImageBuffer.restype=POINTER(c_uint16*(w.value*h.value))
        out = self.atik.ArtemisImageBuffer(self.handle)
        outcontents = out.contents
        
        # There is probably a faster way to do this but I am scared of memory
        image = []
        for i in outcontents:
            image.append(i)
            
        imagearr = (np.array(image)).reshape(h.value, w.value)
        return imagearr
    
    def set_binning(self, xbin, ybin):
        #xbinmax = 1024
        #ybinmax = 1024
        self.atik.ArtemisBin(self.handle, xbin, ybin)
        
    def get_CCD_temp(self):
        temperature=c_int()
        self.atik.ArtemisTemperatureSensorInfo(self.handle, 1, byref(temperature))
        return temperature.value/100 # temp is returned in 1/100ths of a degree
    
    def set_CCD_temp(self, setpoint_decimal):
        setpoint = int(np.round(setpoint_decimal,2)*100)
        self.atik.ArtemisSetCooling(self.handle, setpoint)
        
    def get_CCD_temp_setpoint(self):
        flags=c_int()
        level=c_int()
        minlvl=c_int()
        maxlvl=c_int()
        setpoint=c_int()
        self.atik.ArtemisCoolingInfo(self.handle, byref(flags), byref(level), 
                                       byref(minlvl), byref(maxlvl), 
                                       byref(setpoint))
        return setpoint.value/100
    
    def warm_up(self):
        self.atik.ArtemisCoolerWarmUp(self.handle)
        return
    
    def is_connected(self):
        self.atik.ArtemisIsConnected.restype = c_bool
        ans = self.atik.ArtemisIsConnected(self.handle)
        return ans
    
    def disconnect(self):
        self.atik.ArtemisDisconnect.restype = c_bool
        ans = self.atik.ArtemisDisconnect(self.handle)
        return ans
    
class Binocular(ATIKCamera):
    def __init__(self):
        self.c1 = ATIKCamera(1)
        self.c2 = ATIKCamera(2)
    
    def set_cooling(self, setpoint):
        self.c1.set_CCD_temp(setpoint)
        self.c2.set_CCD_temp(setpoint)
        
    def get_CCD_temps(self, printit=True):
        c1temp = self.c1.get_CCD_temp()
        c2temp = self.c2.get_CCD_temp()
        c1set = self.c1.get_CCD_temp_setpoint()
        c2set = self.c2.get_CCD_temp_setpoint()
        if printit:
            log("SCIP", f"C1 set to {c1set}, currently {c1temp}. C2 set to {c2set}, currently {c2temp}.")
        return c1temp, c2temp
        
    def dual_exposure(self, exptime, m,s, show = False, binning = 1, 
                      save = 'fits', target = 'test', session=False, comment='None'):
        self.c1.set_binning(binning, binning)
        self.c2.set_binning(binning, binning)
        log("SCIP", f"C1: {self.c1.state()}; C2: {self.c2.state()}.")
        log("SCIP", f"Exposing both cameras for {exptime} seconds.")
        exp_starttime = datetime.now(timezone.utc)
        self.c1.start_exposure(exptime)
        self.c2.start_exposure(exptime)
        while not (self.c1.image_ready() and self.c2.image_ready()):
            #log("SCIP", f"C1: {self.c1.state()}; C2: {self.c2.state()}")
            print(self.c1.state(),self.c2.state())
            time.sleep(3)
        img1 = self.c1.get_image()
        img2 = self.c2.get_image()
        
        if save == 'fits':
            save_as_fits(img1, 'c1', exp_starttime, exptime, binning, self,m,s, target=target, session=session, comment=comment)
            save_as_fits(img2, 'c2', exp_starttime, exptime, binning, self,m,s, target=target, session=session, comment=comment)
            
        if save == 'numpy':  
            img1_name = "c1" + time.strftime("_%d_%b_%Y_%H-%M-%S", exp_starttime)
            np.save(img1_name, img1)
            img2_name = "c2" + time.strftime("_%d_%b_%Y_%H-%M-%S", exp_starttime)
            np.save(img2_name, img2)
        
        if show == True:
            show_imgs(img1, img2, star_detect = False)
            
        return img1, img2, exp_starttime

    def warm_up(self):
        self.c1.warm_up()
        self.c2.warm_up()
        return

    def is_connected(self):
        c1con = self.c1.is_connected()
        c2con = self.c2.is_connected()
        return c1con, c2con

    def disconnect(self):
        self.c1.warm_up()
        self.c2.warm_up()
        c1dis = self.c1.disconnect()
        c2dis = self.c2.disconnect()
        return c1dis, c2dis
    
class CubeProMount(object):
    def __init__(self, driverID = "ASCOM.iOptron2014.Telescope"):
        if driverID == None:
            chooser = win32com.client.Dispatch("ASCOM.Utilities.Chooser")
            chooser.DeviceType = "Telescope"
            driverID = chooser.Choose(None)

        self.driverID = driverID
        mount = self.mount = win32com.client.Dispatch(self.driverID)
        self.Name = mount.Name
        mount.connected = True
        
        self.calstar_mag = None
        self.calstar_ra = None
        self.calstar_dec = None
        self.calstar_expiry = None

        self.local = None
        self.conjugate = None
        
    def stop(self):
        self.mount.AbortSlew()
    
    def track(self, value):
        if value == True:
            self.mount.Tracking = True
        elif value == False:
            self.mount.Tracking = False
         
    def goto(self, coordtype, coord1, coord2):
        # Go to a set of coordinates
        if self.mount.slewing == True:
            log("SCIP","Currently slewing. Try again later.")
        else:
            if coordtype == "altaz":
                alt = coord1
                az = coord2
                # Check if you are already there
                if (abs(self.mount.Altitude - alt) > 0.01) or (abs(self.mount.Azimuth - alt) > 0.01):
                    self.mount.SlewToAltAz(az, alt) #order is az, alt for this function
                    self.mount.Tracking = False
                    log("SCIP", f"Slewed to {alt} alt, {az} az.")
                else:
                    log("SCIP", f"Already pointed at {alt} alt, {az} az.")
                
            elif coordtype == "radec":
                ra = coord1
                dec = coord2
                # Check if you are already there
                if (abs(self.mount.RightAscension - ra) > 0.01) or (abs(self.mount.Declination - dec) > 0.01):
                    #self.mount.TargetDeclination = dec
                    #self.mount.TargetRightAscension = ra
                    #self.mount.SlewToTarget()
                    self.mount.SlewToCoordinates(ra, dec)
                    log("SCIP", f"Slewed to {ra} ra, {dec} dec.")
                log("SCIP", f"Already pointed at {ra} ra, {dec} dec.")
        return
    
    def get_pointing(self, coordtype):
        # Get the current pointing of the mount in the selected coordinate type
        if coordtype == "altaz":
            alt = self.mount.Altitude
            az = self.mount.Azimuth
            return alt,az
        elif coordtype == "radec":
            ra = self.mount.RightAscension
            dec = self.mount.Declination
            return ra,dec
    
    def point_away_from_sun(self):
        # Point the photometer low on the nothern horizon
        self.goto("altaz", 0.0, 0.0)
        self.track(False)    
        return
    
        
# This class communicates with the Arduino Nano thru a COM port and gets the 
# temperature and humidity data from the sensor. Check "Device Manager" to find
# out which COM port it is connected to. 
class TempHumidSensor(object):
    def __init__(self, COMport): 
        # Connect to the port
        self.ser = serial.Serial(COMport, 9600, timeout=1)
        
    def temperature(self):
        # Get the temperature by writing "t" to the Arduino
        t = 't'.encode('utf-8')
        self.ser.write(t)
        return float(self.ser.read(7).decode("utf-8")) # Yeah my code is very readable.
    
    def humidity(self):
        # Get the humidity by writing "h" to the Arduino
        h = 'h'.encode('utf-8')
        self.ser.write(h)
        return float(self.ser.read(7).decode("utf-8"))

    def __del__(self):
        # Close port
        self.ser.close()



