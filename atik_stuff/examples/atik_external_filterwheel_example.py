import AtikSDK
import random
import time

try:
    efw = AtikSDK.AtikSDKCamera()

    while(efw.is_efw_present(0) == False):
        time.sleep(1)

    efw.connect_efw(0)
 
    if(efw.is_efw_connected() == True):
        details = efw.efw_device_details()
        print(f"EFW Type {details[0]}, EFW Serial Number {details[1]}")

        num_of_pos = efw.efw_num_positons()
        print(f"Number of EFW positions: {num_of_pos}")

        current_Pos = efw.get_current_efw_positon()
        print(f"Current EFW position: {current_Pos}")
        
        #Move the filter wheel to 3 random positions 
        numOfMoves = 0
        while(numOfMoves < 3):
            randPos = random.randrange(0, num_of_pos - 1)
            if(randPos != efw.get_current_efw_positon()):
                efw.set_efw_position(randPos)
                while(efw.check_efw_moving() == True):
                    time.sleep(.5)
                print(f"New Position is: {efw.get_current_efw_positon()}")
                numOfMoves += 1              

    efw.disconnect()

    AtikSDK.ArtemisShutdown()

except Exception as e:
    print("Atik Filter Wheel error:", str(e))
