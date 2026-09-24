import AtikSDK
import random
import time

try:
    cam = AtikSDK.AtikSDKCamera()

    while(cam.is_device_present(0) == False):
        time.sleep(1)

    cam.connect(0)
 
    if(cam.is_connected()):
        info = cam.internal_filterwheel_info()

        num_of_pos = info[0]
        print(f"Number of positions: {num_of_pos}")

        current_Pos = info[2]
        print(f"Current cam position: {current_Pos}")
        
        #Move the filter wheel to 3 random positions 
        numOfMoves = 0
        while(numOfMoves < 3):
            randPos = random.randrange(0, num_of_pos - 1)
            info = cam.internal_filterwheel_info()
            if(randPos != info[2]):
                cam.move_internal_filterwheel(randPos)
                while True:
                    time.sleep(.5)
                    info = cam.internal_filterwheel_info()
                    moving = info[1]
                    if not moving:
                        break
                print(f"New Position is: {info[2]}")
                numOfMoves += 1              

    cam.disconnect()

    AtikSDK.ArtemisShutdown()

except Exception as e:
    print("Atik Filter Wheel error:", str(e))
