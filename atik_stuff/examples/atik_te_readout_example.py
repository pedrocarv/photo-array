import AtikSDK

try:
    cam = AtikSDK.AtikSDKCamera()

    cam.connect()

    if(cam.is_connected()):
        readoutMode = cam.get_te_readout_mode()

        print(f"Readout mode is: {readoutMode}")

        cam.set_te_readout_mode(AtikSDK.ReadoutModeTE.Preview)

        readoutMode = cam.get_te_readout_mode()

        print(f"Readout mode is: {readoutMode}")

        cam.disconnect()

    AtikSDK.ArtemisShutdown()

except Exception as e:
    print("Atik Camera error:", str(e))