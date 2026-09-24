import AtikSDK
import time

try:
    cam = AtikSDK.AtikSDKCamera()

    while(cam.is_device_present(0) == False):
        time.sleep(1)
    
    cam.connect()

    #Get the max and min gain and offset values
    ranges = cam.get_gain_offset_range()
    print(f"Minimum Gain {ranges[0]} Maximum Gain {ranges[1]} Minimum Offset {ranges[2]} Maximum Offset {ranges[3]}")

    #Get the current gain and offset, then set.
    gainOffset = cam.get_gain_offset()
    print(f"Current Gain: {gainOffset[0]}, Current Offset: {gainOffset[1]}")  
    cam.set_gain_offset(65, 400)

    #Set the exposure speed into Powersave, take a 10 second exposure then save the result out
    cam.set_exposure_speed(AtikSDK.ExposureSpeed.PowerSave)
    arr = cam.take_image(10)
    try:
        from astropy.io import fits
        hdu = fits.PrimaryHDU(arr)
        filename = 'PowerSave_result.fits'
        hdu.writeto(filename, overwrite=True)

    except ImportError as e:
        print("Could not import AstroPy to save FITS image:", str(e))

    #Set the exposure speed into Normal, take a 1 second exposure then save the result out
    cam.set_exposure_speed(AtikSDK.ExposureSpeed.Normal)
    arr = cam.take_image(1)
    try:
        from astropy.io import fits
        hdu = fits.PrimaryHDU(arr)
        filename = 'Normal_result.fits'
        hdu.writeto(filename, overwrite=True)

    except ImportError as e:
        print("Could not import AstroPy to save FITS image:", str(e))

    #Disconnect camera and shutdown
    cam.disconnect()
    AtikSDK.ArtemisShutdown()
    

except Exception as e:
    print("Atik Camera error:", str(e))
