# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 13:59:29 2019

@author: Dawn
"""
import pandas as pd
import math
import csv
import matplotlib.pyplot as plt
from datetime import datetime as dt
from datetime import timedelta

xl = pd.read_excel("Arecibo_8446_data_notes_combined.xls")

with open('Ap.txt') as csvfile:        
    ap=list(csv.reader(csvfile,delimiter=' '))
  
allap = []
alldates = []
    
for j in range(len(ap)-17):
    timeobject = dt.strptime( ap[j+17][0]+" "+ap[j+17][1], "%Y-%m-%d %H:%M")
    if timeobject > dt(1996, 1, 1, 0, 0):
        allap.append(float(ap[j+17][2]))
        alldates.append(timeobject)
    



dataap = []
datadates = []

isrdata = []
isrdates = []

winddata = []
winddates = []

photdata = []
photdates = []

# Add AP index to a copy of the master Excel file
for i in range(720):
    year = xl['year'][i]
    if math.isnan(year) != True:
        if (1995 < year and 2010 > year):
            date = str(year)+'-'+xl['textdate'][i][:2]+'-'+xl['textdate'][i][3:]
            
            # Get daily average AP's for dates with data
            for j in range(len(ap)-17):
                #print(type(ap[j+17][0]))
                if ap[j+17][0] == date:
                    ap_sum = 0
                    ap_ave = 0
                    for k in range(8):
                        ap_sum += float(ap[j+17+k][2])
                        ap_ave = ap_sum/8
                    xl['Ap index'][i] = ap_ave
                    datadates.append(dt.strptime( ap[j+17][0]+" 12:00", "%Y-%m-%d %H:%M"))
                    dataap.append(ap_ave)
                    break
                
            # ISR data?  
            if xl['ISR?'][i] == 'ISR':
                isrdates.append(dt.strptime( date+" 12:00", "%Y-%m-%d %H:%M"))
                isrdata.append(-5)
                
            # O8446 data?
            if xl['PHOT#1'][i] == 'OI 8446':
                photdates.append(dt.strptime( date+" 12:00", "%Y-%m-%d %H:%M"))
                photdata.append(-25)
                
            # Wind data?
            if xl['WINDS?'][i] == 'OI 6300 winds':
                winddates.append(dt.strptime( date+" 12:00", "%Y-%m-%d %H:%M"))
                winddata.append(-15)
            
            
#xl.to_excel("output.xlsx")
plt.figure(figsize=(20,5))
plt.plot(alldates, allap, linewidth=1, color='cadetblue')                   
plt.plot(datadates, dataap, 'o', color='teal', label="AP values for dates with Arecibo data")
plt.plot(isrdates, isrdata, 'o', color='limegreen', label="Dates with ISR data")
plt.plot(photdates, photdata, 'o', color='orange', label="Dates with OI 8446 data")
plt.plot(winddates, winddata, 'o', color='mediumpurple', label="Dates with wind data")

allgood = []
for x in datadates:
    if (x in isrdates) and (x in photdates) and (x in winddates):
        allgood.append(x)

plt.plot(allgood, [-35]*len(allgood), '*', color='red', label="Dates with all data")

plt.hlines(30, alldates[0], alldates[len(alldates)-1])

plt.legend()
               
