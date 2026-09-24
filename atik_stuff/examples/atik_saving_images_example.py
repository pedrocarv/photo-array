import AtikSDK
import time

try:
    cam = AtikSDK.AtikSDKCamera()

    while(cam.is_device_present(0) == False):
        time.sleep(1)
    
    cam.connect()

    #Take 3 images and save them as FITS Files
    num_images = 0
    while(num_images < 3):
        arr = cam.take_image(1)
        try:
            from astropy.io import fits
            hdu = fits.PrimaryHDU(arr)
            filename = 'result%i.fits' %(num_images+1)
            hdu.writeto(filename, overwrite=True)

        except ImportError as e:
            print("Could not import AstroPy to save FITS image:", str(e))

        num_images = num_images + 1
        print("Image taken:", num_images)    

    #Take 3 images and save them as TIFF Files
    num_images = 0
    while(num_images < 3):
        arr = cam.take_image(1)
        try:
            import tifffile
            filename = 'result%i.tiff' %(num_images+1)
            tifffile.imwrite(filename, arr)
            
        except ImportError as e:
            print("Could not import tifffile to save TIFF image:", str(e))  

        num_images = num_images + 1
        print("Image taken:", num_images)   

    cam.disconnect()

    AtikSDK.ArtemisShutdown()

except Exception as e:
    print("Atik Camera error:", str(e))
