# -*- coding: utf-8 -*-
"""
Created on Mon Jan  4 17:46:47 2021

@author: dawnh
"""

import numpy as np
import matplotlib.pyplot as plt
from read_fits_files import *

#% Read noise from bias images

directory = 'C:/Users/photometer/Box Sync/SCIP/Synced_Raw_Data/binning_v_read_noise'
master_on, master_off = get_fits_from_folder(directory, sort_by_date=False, open_image=True, pix_sum=False)
fig,ax = plt.subplots(2,1)

#xaxis_len = 211
exptime = 0.001

masters = master_on, master_off
for m_index, master in enumerate(masters):
    #master = [item for item in master if item['EXP-TIME']==exptime]
    #master = [item for item in master if item['NAXIS1']==xaxis_len]
    if m_index == 0:
        camera='On-band'
    else: 
        camera='Off-band' 
    rnoise_list = []  
    total_noise_list = []
    binning = []
    for item in master:
        image=item['image']*0.19
        print(np.max(item['image']))
        rnoise = np.sqrt(np.mean((image- np.mean(image))**2))
        total_noise = image.size*rnoise**2
        total_noise_list.append(total_noise)
        rnoise_list.append(rnoise)
        binning.append(item['BINNING'])
    ax[m_index].plot(binning, rnoise_list, 'o-')
    ax[m_index].set_ylabel("Read noise (e-)")
    
    #ax[m_index].plot(binning, total_noise_list, 'o-')
    ax[m_index].set_xlabel('Binning (X and Y)')
    #ax[m_index].set_ylabel("Total read noise (e$^-$)")
    ax[m_index].set_title(camera)
    #ax[m_index].set_yscale("log")

    
"""
#%% Plot sum or mean vs approximate binning value

directory = 'C:/Users/dawnh/Box Sync/SCIP/Synced_Raw_Data/bias_darks_vs_potential_best_bins'
master_on, master_off = get_fits_from_folder(directory, sort_by_date=False, open_image=False, pix_sum=True)
fig,ax = plt.subplots(2,1)

masters = master_on, master_off
comment = 'variable binning, inside'
for m_index, master in enumerate(masters):
    #master = [item for item in master if item['COMMENT']==comment]
    master = [item for item in master if item['EXP-TIME']==300]
    #master = [item for item in master if item['NAXIS1']!=3379]
    temps = []
    sums = []
    bins = []
    for item in master:
        #temps.append(item['CCD_TEMP'])
        #sums.append(np.max(item['image'].flatten()))
        sums.append(item['pix_sum'])
        bins.append(3379/item['NAXIS1'])
    if m_index == 0:
        camera='On-band'
    else: 
        camera='Off-band' 
    ax[m_index].hist(sums, bins=100)
    #ax[m_index].plot(bins,sums)
    ax[m_index].set_xlabel('Bin number')
    ax[m_index].set_ylabel("CCD sum")
    ax[m_index].set_title(camera + ", T = 3C")
    #ax[m_index].set_xscale('log')
   
#%% Compare same set CCD temperature, indoors vs outdoors test

directory = 'C:/Users/dawnh/Box Sync/SCIP/Synced_Raw_Data/bias'
master_on, master_off = get_fits_from_folder(directory, sort_by_date=False, open_image=True, pix_sum=True)

fig,ax = plt.subplots(2,1)
hist=True

masters = master_on, master_off
comments = ['inside binning', 'outside binning']
for m_index, master in enumerate(masters):
    for comment in comments:
        masterc = [item for item in master if item['COMMENT']==comment]
        masterc = [item for item in masterc if item['EXP-TIME']==0.001]
        temps = []
        sums = []
        for item in masterc:
            temps.append(item['CCD_TEMP'])
            sums.append(item['pix_sum'])
        if m_index == 0:
            camera='On-band'
        else: 
            camera='Off-band'
        if comment=='inside binning':    
            marker='x'
            color='blue'
        elif comment=='outside binning':
            marker='.'
            color='orange'
        if not hist:
            ax[m_index].plot(temps, sums, marker, label = comment)
        if hist:    
            ax[m_index].hist(sums, color=color, bins=200)
        print(f'{comment} {camera} {np.mean(sums)}')
    if not hist:
        ax[m_index].set_title(f'{camera} indoors vs. outdoors')
        ax[m_index].legend()
        ax[m_index].set_xlabel('Temperature')
        ax[m_index].set_ylabel('Sum over CCD')
    if hist:
        ax[m_index].set_title(f'{camera} indoors vs. outdoors')
        ax[m_index].legend()
        ax[m_index].set_xlabel('Sum over CCD')
        ax[m_index].set_ylabel('Number of images')
        
        

#%% Bias analysis - x and y slices

#directory = 'C:/Users/dawnh/Documents/Local Raw Data/bias'
directory = 'C:/Users/dawnh/Box Sync/SCIP/Synced_Raw_Data/bias'
master_on, master_off = get_fits_from_folder(directory, sort_by_date=False, open_image=False, pix_sum=False)

binning = True
imagetype = 'bias'
comment = 'outside binning'
min_images_per_temp = 10
max_images_per_temp = 10
fig,ax = plt.subplots(2,2)

flag_xy_means = True
flag_median_sums = False

if flag_xy_means:
    fig,ax = plt.subplots(2,2)
if flag_median_sums:
    fig,ax = plt.subplots(2,1)

masters = master_on, master_off
for m_index, master in enumerate(masters):
    # Get images with desired target:
    if imagetype == 'bias':
        master = [item for item in master if item['EXP-TIME']==0.001]
    else:
        print('Not set up for this yet bruh.')
        break
    # Get images with desired binning:
    if binning == True:
        master = [item for item in master if item['NAXIS1']==337 and item['NAXIS2']==422]
    elif binning == False:
        master = [item for item in master if item['NAXIS1']==2703 and item['NAXIS2']==3379]
    # Pick out only images with this comment
    master = [item for item in master if item['COMMENT']==comment]   
    # Get a list of all tempertures in the set
    temps=[]
    [temps.append(np.round(master[index]['CCD_TEMP'], decimals=1)) for index in range(len(master)) if np.round(master[index]['CCD_TEMP'], decimals=1) not in temps]
    temps.sort()
    print(temps)
    # Create median images for each temp
    median_sums = []
    median_sums_temps = []
    for temp in temps:
        # Get all the filenames of the images at this temperature
        filenames_at_temp = [master[index]['filename'] for index in range(len(master)) if np.round(master[index]['CCD_TEMP'], decimals=1) == temp]
        # Get the images at this temp
        images_at_temp = []
        for f in filenames_at_temp:
            images_at_temp.append(image_from_filename(directory, f))
        # Trim according to max/min number of images per median
        print(f"T = {temp}C has {len(images_at_temp)} images.")
        if len(images_at_temp) >= min_images_per_temp:
            
            if len(images_at_temp) > max_images_per_temp:
                images_at_temp = images_at_temp[0:max_images_per_temp]
                
            image_median = np.median( np.stack(images_at_temp), axis = 0 )
            
            if flag_xy_means:
                x_means = np.mean(image_median, axis=0)
                y_means = np.mean(image_median, axis=1)
                ax[m_index,0].plot(x_means, label=f"{temp} C")
                ax[m_index,1].plot(y_means, label=f"{temp} C")
            
            if flag_median_sums:
                median_sums_temps.append(temp)
                median_sums.append(np.sum(image_median))
                     
    if m_index == 0:
        camera='On-band'
    else: camera='Off-band'
    
    if flag_xy_means:
        ax[m_index,0].set_title(f"{camera}, median image")
        ax[m_index,0].legend()
        ax[m_index,0].set_xlabel("X-axis coordinate")
        ax[m_index,0].set_ylabel("Mean value of row")
        ax[m_index,1].set_title(f"{camera}, median image")
        ax[m_index,1].legend()
        ax[m_index,1].set_xlabel("Y-axis coordinate")
        ax[m_index,1].set_ylabel("Mean value of column")        
        
    if flag_median_sums:
        ax[m_index].plot(median_sums_temps, median_sums, 'o')
        ax[m_index].set_title(f"{camera} sum of CCD")
        ax[m_index].set_xlabel('Temperature (C)')
        ax[m_index].set_ylabel('Sum of median bias image')
        for i in range(len(median_sums_temps)):
            logfilename = camera + "-sums-v-temp.csv"
            logdata(logfilename, str(median_sums_temps[i])+','+str(median_sums[i]))

    

#%% Show some images
            
directory = 'C:/Users/dawnh/Box Sync/SCIP/Synced_Raw_Data/binning_v_read_noise'
master_on, master_off = get_fits_from_folder(directory, sort_by_date=False, open_image=True, pix_sum=False)

# comment = 'inside binning'

# master_on = [item for item in master_on if item['COMMENT']==comment]
# master_off = [item for item in master_off if item['COMMENT']==comment]

num_images = len(master_on)
#num_images = 3
for index in range(num_images):
    fig, ax = plt.subplots(1,2)
    img1 = master_on[index]['image']
    img2 = master_off[index]['image']
    im1 = ax[0].imshow(img1, interpolation='none', cmap='gray')#, vmin=np.percentile(img1, 1), vmax=np.percentile(img1, 99))
    im2 = ax[1].imshow(img2, interpolation='none', cmap='gray')#, vmin=np.percentile(img2, 1), vmax=np.percentile(img2, 99))
    ax[0].set_title(f"On-band {np.round(3379/master_on[index]['NAXIS1'])} image")
    ax[1].set_title(f"Off-band {np.round(3379/master_off[index]['NAXIS1'])} image")
    fig.colorbar(im1, ax=ax[0])
    fig.colorbar(im1, ax=ax[1])
"""
