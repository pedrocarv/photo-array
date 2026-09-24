import AtikSDK
import numpy as np
from ctypes import *
from AtikSDK import FastCallback

import time

callback_arr = ""
has_image = False

def fast_callback(handle, x, y, w, h, bx, by, img):
    raw_string = string_at(img, w * h * 2)
    global callback_arr    
    callback_arr = np.fromstring(raw_string, dtype=np.uint16, count=w * h)
    callback_arr = callback_arr.reshape(h, w)
    global has_image 
    has_image = True

try:
    cam = AtikSDK.AtikSDKCamera()

    while(cam.is_device_present(0) == False):
        time.sleep(1)
    
    cam.connect()

    fast_func = FastCallback(fast_callback)

    if(cam.has_fast_mode()):
        cam.set_exposure_speed(AtikSDK.ExposureSpeed.Fast)
        cam.set_fast_callback(fast_func)

        cam.start_fast_exposure(1)

        # exit after receiving 10 images
        num_images = 0
        while(num_images < 10):
            i = 0
            while(has_image == False):
                print(f"Waiting for Image... {i}", end='\r')
                i = i + 1
            
            has_image = False
            num_images = num_images + 1
            print("Image taken:", callback_arr, callback_arr.shape)
        
        cam.stop_exposure()

    cam.disconnect()

    AtikSDK.ArtemisShutdown()

except Exception as e:
    print("Atik Camera error:", str(e))
