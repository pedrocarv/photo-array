import AtikSDK
import time

try:
    cam = AtikSDK.AtikSDKCamera()

    api_version = cam.get_api_version()

    print("API version is: ", api_version)

    while(cam.is_device_present(0) == False):
        time.sleep(1)

    cam.connect()

    properties = cam.get_properties()
    print("Camera properties:", properties)

    serial = cam.get_serial_str()
    print("Camera serial:", serial)

    fx3, fpga = cam.get_firmware_versions()
    print(f"FX3 version: {fx3}, FPGA version: {fpga}")

    temp_count = cam.get_temperature_sensor_count()
    print("Temperature sensors count:", temp_count)

    if temp_count > 0:
        temp = cam.get_temperature()
        print("Temperature °C:", temp)
    else:
        print("No temperature sensor on camera")

    # Adjust binning, subframing, cooling before taking the image.
    # cam.set_binning(2, 2)
    # cam.set_cooling(-5)
    # cam.set_subframe(60, 60, 640, 640)

    ranges = cam.get_gain_offset_range()
    print(f"Minimum Gain {ranges[0]} Maximum Gain {ranges[1]} Minimum Offset {ranges[2]} Maximum Offset {ranges[3]}")
    cam.set_gain_offset(65, 200)
    
    if cam.has_pad_data():
        cam.set_pad_data(True)

    exposure_secs = 0.01
    arr = cam.take_image(exposure_secs)
        
    print("Image size:", arr.shape)

    # Save TIFF
    try:
        import tifffile

        tifffile.imwrite("result.tiff", arr)
    except ImportError as e:
        print("Could not import tifffile to save TIFF image:", str(e))

    # Save FITS file
    try:
        from astropy.io import fits

        hdu = fits.PrimaryHDU(arr)
        hdu.writeto('result.fits', overwrite=True)
    except ImportError as e:
        print("Could not import AstroPy to save FITS image:", str(e))

    print("Saved results.tiff and results.fits successfully")

except Exception as e:
    print("Atik Camera error:", str(e))
