import AtikSDK
from AtikSDK import SixteenBitMode
import time

try:
    cam = AtikSDK.AtikSDKCamera()

    while(cam.is_device_present(0) == False):
        time.sleep(1)

    cam.connect()

    #Check if the camera has 16bit mode.
    if(cam.has_16bit_mode()):

        #Get the current 16bit mode setting
        res = cam.get_16bit_mode()
        print(f"Current 16 bit mode: {res}")

        #Set the 16bit mode into CombinedExposure
        cam.set_16bit_mode(SixteenBitMode.CombinedExposure)
        exp_dur = 1

        num_images = 0
        while(num_images < 3):
            arr = cam.take_image(exp_dur)

            try:
                from astropy.io import fits
                hdu = fits.PrimaryHDU(arr)
                filename = 'result%i.fits' %(num_images+1)
                hdu.writeto(filename, overwrite=True)

            except ImportError as e:
                print("Could not import AstroPy to save FITS image:", str(e))
            print("Saved results.fits successfully")
            num_images = num_images + 1
            print("Image taken:", num_images)

    cam.disconnect()

    AtikSDK.ArtemisShutdown()

except Exception as e:
    print("Atik Camera error:", str(e))
