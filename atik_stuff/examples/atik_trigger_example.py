import AtikSDK
import time

try:
    cam = AtikSDK.AtikSDKCamera()

    while(cam.is_device_present(0) == False):
        time.sleep(1)

    cam.connect()

    # After calling set_triggered_exposure every subsequent call to take_image will
    # arm the trigger in 
    cam.set_triggered_exposure(True)

    # After calling take_image the trigger will be armed
    arr = cam.take_image(2)
    # Once the image is recieved the trigger is de-armed ready for another call to
    # take_image

    try:
        import tifffile

        tifffile.imwrite("result.tiff", arr)
    except ImportError as e:
        print("Could not import tifffile to save TIFF image:", str(e))

    cam.disconnect()
except Exception as e:
    print("Atik Camera error:", str(e))
