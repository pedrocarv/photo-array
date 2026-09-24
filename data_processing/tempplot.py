# -*- coding: utf-8 -*-
"""
Created on Thu Feb 27 22:52:36 2020

@author: Dawn
"""
from datetime import datetime
import matplotlib.pyplot as plt

#%% Plot where temperature only is recorded
'''
dates = []
temps = []
file = open('SCIP_temp_log_start_2_18_2020.txt', 'r') 
for line in file:
    date = datetime.strptime( line[0:24], "%m/%d/%Y %H:%M:%S %Z " )
    dates.append(date)
    temp = float(line[24:])
    temps.append(temp)
    
plt.plot(dates,temps,',')
'''
#%% Plot for temerature and humidity, temperature first
dates = []
temps = []
humids = []
file = open('SCIP_temp_log_test_8_3_2020.txt', 'r') 
for line in file:
    date = datetime.strptime( line[0:24], "%m/%d/%Y %H:%M:%S %Z " )
    dates.append(date)
    stringvals = line[24:]
    i = 1
    for char in stringvals:
        if char == ' ':
            space_index = i
            break
        i += 1
    temp = float(line[24:24+i])
    temps.append(temp)
    humid = float(line[24+i:])
    humids.append(humid)

fig, ax = plt.subplots(2,1, figsize=(20,15))    
ax[0].plot(dates,temps, label="Sensor temperature")
ax[0].set_title("Dome box temperature, with reflective cover")
ax[0].axhline(25,color='r', label="Max operational filter temp")
#ax[0].axhline(22.222,color='g', label="Room temp (72 F)")
sunset = datetime.strptime("08/04/2020 00:03:54 UTC ", "%m/%d/%Y %H:%M:%S %Z " )
sunrise = datetime.strptime("08/04/2020 09:53:17 UTC ", "%m/%d/%Y %H:%M:%S %Z " )
ax[0].axvline(sunset, linestyle='--', linewidth=2, color = '#454545', label = "Sunset")
ax[0].axvline(sunrise, linestyle='--', linewidth=2, color = '#999999', label = "Sunrise")
ax[0].legend()
ax[0].set_xlabel("Date and time (UTC)")
ax[0].set_ylabel("Temperature (C)")

ax[1].plot(dates,humids, label="Humidity")
ax[1].set_title("Dome box humidity")
ax[1].axvline(sunset, linestyle='--', linewidth=2, color = '#454545', label = "Sunset")
ax[1].axvline(sunrise, linestyle='--', linewidth=2, color = '#999999', label = "Sunrise")
ax[1].legend()
ax[1].set_xlabel("Date and time (UTC)")
ax[1].set_ylabel("% Humidity")