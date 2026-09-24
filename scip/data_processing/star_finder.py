# -*- coding: utf-8 -*-
"""
Created on Sun Nov 22 22:56:41 2020

@author: dawnh
"""

import os
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
from datetime import datetime, timezone
from matplotlib.patches import Circle
from matplotlib.patheffects import withStroke
import scipy.signal
from scipy.stats import multivariate_normal
from scipy import signal
from skimage.measure import label, regionprops
from skimage.morphology import remove_small_objects
import time
import pywt



# get star postions and brightnesses

def get_star_vals(stars, camera):
    if stars == '[]':
        return None
    split_stars = stars[2:-2].split('], [')
    star_numbers = []
    if camera == 'sim':
        for star in split_stars:
            vals = star.split(',')
            #print(vals)
            numbers = [float(vals[1])/10, float(vals[0])/10, float(vals[2])]
            star_numbers.append(numbers)
    else:
        for star in split_stars:
            vals = star.split(',')
            #print(vals)
            numbers = [float(i) for i in vals]
            star_numbers.append(numbers)
    return star_numbers




# open images

def get_fitslist(folder):
    directory = folder
    fitslist = [] # a list of dicts, sorted by datetime, of the data of each file.
    keys = []
    
    for file in os.listdir(directory):
        filename = os.fsdecode(file)
        hdul = fits.open(directory+'/'+filename)
        header = [*hdul[0].header]
        if keys == []:
             keys = header
        if header != keys:
            print(f"Found tft {filename} does not match the same keys at previous files. F.")
            hdul.close()
            break
        imagedict = {'FILENAME':file}
        imagedict['IMAGE'] = hdul[0].data.astype(float)
        for key in keys: # add all the header file data
            imagedict[key] = hdul[0].header[key]
        fitslist.append(imagedict)
        hdul.close()
        
    print(f'Your keywords for the files in {folder} are: {[*fitslist[0]]}')
    
    return(fitslist)




# draw circles around real stars

def circle_real_stars(fitslist, index, potential_centroids=None, image = None):
    if type(image) == type(None):
        image = fitslist[index]['IMAGE']
    #plt.hist(image.flatten(), list(np.arange(0,2500,10)), log=True)
    stars = get_star_vals(fitslist[index]['STARS'], fitslist[index]['CAMERA'])
    fig,ax = plt.subplots()
    plt.imshow(image, interpolation='None', vmin=np.percentile(image, 1), vmax=np.percentile(image, 99), cmap='gray')
    #plt.imshow(image, interpolation='None', cmap='gray')
    plt.title(fitslist[index]['FILENAME'])
    if stars == None:
        return
    # if fitslist[index]['CAMERA'] == 'sim':
    #     for coords in stars:
    #         x=coords[1]/10
    #         y=coords[0]/10
    #         radius=30
    #         print(x,y)
    #         circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=1,
    #                         edgecolor='red', facecolor=(0, 0, 0, .0125))
    #         ax.add_artist(circle)
    # else:
    for coords in stars:
        x=coords[0]
        y=coords[1]
        radius=45
        circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=3,
                        edgecolor='limegreen', facecolor=(0, 0, 0, .0125))
        ax.add_artist(circle)
    if potential_centroids != None:
        for coords in potential_centroids:
            x=coords[0]
            y=coords[1]
            radius=30
            circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=2,
                            edgecolor='red', facecolor=(0, 0, 0, .0125))
            ax.add_artist(circle)




#Functions for pre-processing, finding stars

def subtract_median(image):
    return image-np.median(image)
    
def subtract_median_filtered_image(image):
    median_filt_image = scipy.ndimage.median_filter(image, size=10)
    return image - median_filt_image

def background_estimate(image, index):
    k=4
    original_std = np.std(image)
    image_noisetrim = image[ abs(image) < 3*np.std(image) ]
    while(abs(np.std(image_noisetrim) - original_std) > 1e-10):
        original_std = np.std(image_noisetrim)
        image_noisetrim = image_noisetrim[ abs(image_noisetrim) < 3*np.std(image_noisetrim) ]
    fitslist[index]['NOISESTD'] = np.std(image_noisetrim)
    
def matched_filter(image, kernel):
    return signal.convolve2d(image, kernel, mode='same')


def threshold(multiplier, image, index, draw=False):
    threshold = multiplier*fitslist[index]['NOISESTD']
    potential_trutharray = image > threshold
    if draw:
        potential_coords = np.where(potential_trutharray)
        # circle the points
        fig,ax = plt.subplots()
        #plt.imshow(image, interpolation='None', vmin=np.percentile(image, 1), vmax=np.percentile(image, 99), cmap='gray')
        plt.imshow(potential_trutharray, interpolation='None', cmap='gray')
        #print(f'{len(potential_coords[0])} pixels above threshold found.')
        # for i in range(len(potential_coords[0])):
        #     x=potential_coords[1][i]
        #     y=potential_coords[0][i]
        #     radius=30
        #     circle = Circle((x, y), radius, clip_on=False, zorder=10, 
        #                     linewidth=1, edgecolor='yellow', facecolor=(0, 0, 0, 0))
        #     ax.add_artist(circle)
    return potential_trutharray

def find_single_regions(potential_trutharray):
    regions = label(potential_trutharray, background = 0)
    multipixel_regions = remove_small_objects(regions, min_size = 2)
    single_pixel_regions = [a-b for a,b in zip(regions, multipixel_regions)]
    return single_pixel_regions

def find_regions(potential_trutharray, draw=False):
    regions = label(potential_trutharray, background = 0)
    multipixel_regions = remove_small_objects(regions, min_size = 2)
    #multipixel_regions = regions
    potential_source_regions = regionprops(multipixel_regions) 
    if draw:
        fig,ax=plt.subplots()
        plt.imshow(multipixel_regions, cmap='gray') 
        plt.colorbar
    potential_centroids = []
    for i in range(len(potential_source_regions)):
        x=potential_source_regions[i]['centroid'][1]
        y=potential_source_regions[i]['centroid'][0]
        potential_centroids.append([x,y])
        if draw:
            radius=30
            circle = Circle((x, y), radius, clip_on=False, zorder=10, 
                            linewidth=2, edgecolor='red', 
                            facecolor=(0, 0, 0, .0125))
            ax.add_artist(circle)
    return potential_source_regions, potential_centroids

def neighborhood_check(potential_centroids, potential_trutharray, disk_coords, image):
    neighborhood_means = []
    for centroid in potential_centroids:
        # find coords around the centroid
        neighborhood_coords = np.add(np.array(disk_coords), np.array(centroid)).astype(int)
        pixel_vals = []
        for c in neighborhood_coords:
            # if it is not a thresholded pixel
            try:
                if not potential_trutharray[c[1],c[0]]:
                    # add it's pixel value to the list
                    pixel_vals.append(image[c[1],c[0]])
            except IndexError:
                pass
        #print(pixel_vals)
        neighborhood_means.append(np.mean(pixel_vals))
    return neighborhood_means
        



# Compare found stars to actual stars

def compare_results(stars, potential_centroids, filename):
    showit = False
    true_positives = []
    false_positives = []
    false_negatives = []
    
    max_offset = 20
    for coords in potential_centroids:
        found_match = False
        matches = []
        if stars != None:
            for star in stars:
                offset = ((coords[0]-star[0])**2 + (coords[1]-star[1])**2 )**0.5
                if offset < max_offset:
                    found_match = True
                    star.append('matched')
                    star.append('matched')
                    matches.append([coords, offset])
            if found_match == True:
                true_positives.append([coords, len(matches), matches])
            if found_match == False:
                false_positives.append([coords])
    if stars != None:
        for star in stars:
            try:
                test = star[3]
            except IndexError:
                #print(filename)
                false_negatives.append([star])
                #showit=True

    return true_positives, false_positives, false_negatives, showit

#true_positives = [(actualx, actualy, brightness), num_of_matches, [(match1x, match1y, offset), (match1x, match1y, offset)] ]
#false_positives = [foundx, foundy, min_offset]
#false_negatives = [actualx, actualy, brightness]





# aperture phot

def ap_phot(image, potential_centroids, disk_coords, annulus_coords):
    star_counts = []
    for centroid in potential_centroids:
        
        neighborhood_coords = np.add(np.array(disk_coords), np.array(centroid)).astype(int)
        background_coords = np.add(np.array(annulus_coords), np.array(centroid)).astype(int)
        
        neighborhood_vals = []
        for c in neighborhood_coords:
            # if it is not a thresholded pixel
            try:
                    neighborhood_vals.append(image[c[1],c[0]])
            except IndexError:
                pass
        neighborhood_sum = np.sum(neighborhood_vals)
        neighborhood_pixel_count = len(neighborhood_vals)
        
        background_vals = []
        for c in background_coords:
            # if it is not a thresholded pixel
            try:
                    star_val = image[c[1],c[0]]
                    if star_val >= 2**16:
                        print('fuck')
                    background_vals.append(star_val)
            except IndexError:
                pass
        background_mean = np.mean(background_vals)
        
        star_count = neighborhood_sum - neighborhood_pixel_count*background_mean
        star_counts.append(star_count)

    return(star_counts)
        

# draw circles that match photometric circles

def circle_phot_stars(fitslist, index, ri, ro, star_counts, potential_centroids=None):
    image = fitslist[index]['IMAGE']
    #plt.hist(image.flatten(), list(np.arange(0,2500,10)), log=True)
    stars = get_star_vals(fitslist[index]['STARS'], fitslist[index]['CAMERA'])
    fig,ax = plt.subplots()
    plt.imshow(image, interpolation='None', vmin=np.percentile(image, 1), vmax=np.percentile(image, 99), cmap='gray') #vmin=np.percentile(image, 1), vmax=np.percentile(image, 99),
    #plt.imshow(image, interpolation='None', cmap='gray')
    plt.title(fitslist[index]['FILENAME'])
    if stars == None:
        return
    # if fitslist[index]['CAMERA'] == 'sim':
    #     for coords in stars:
    #         x=coords[1]/10
    #         y=coords[0]/10
    #         radius=30
    #         print(x,y)
    #         circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=1,
    #                         edgecolor='red', facecolor=(0, 0, 0, .0125))
    #         ax.add_artist(circle)
    # else:
    for coords in stars:
        x=coords[0]
        y=coords[1]
        radius=50
        circle = Circle((x, y), radius, clip_on=False, zorder=10, linewidth=3,
                        edgecolor='limegreen', facecolor=(0, 0, 0, .0125))
        ax.add_artist(circle)
    if potential_centroids != None:
        for i, coords in enumerate(potential_centroids):
            x=coords[0]
            y=coords[1]
            circle = Circle((x, y), ri, clip_on=False, zorder=10, linewidth=2,
                            edgecolor='red', facecolor=(0, 0, 0, .0125))
            ax.add_artist(circle)
            circle = Circle((x, y), ro, clip_on=False, zorder=10, linewidth=2,
                            edgecolor='yellow', facecolor=(0, 0, 0, .0125))
            ax.add_artist(circle)
            ax.annotate(f'counts = {star_counts[i]}', xy=(x+ri, y), fontsize = 13, bbox={'facecolor': 'white', 'alpha': 1, 'pad': 10})
            fig.texts.append(ax.texts.pop())
#%% Get data


"""
TODO: Once how/when to use the above functions, delete code below (will be kept for reference for now)
folder = 'fits_stars'
#folder='simulation_output_friday_9pm'
fitslist = get_fitslist(folder)

#%% View histograms
# fig = plt.figure(figsize=(12, 3))
# num_images = 5
# for i in range(num_images):
#     ax = fig.add_subplot(5, 1, i+1)
#     image = fitslist[i]['IMAGE']
#     ax.hist(image.flatten(), list(np.arange(0,5000,20)), log=True)
#     ax.set_title(fitslist[i]['FILENAME'])
#     print(fitslist[i])

#%% run ap_phot
 
num_images = len(fitslist)
#num_images = 1

disk_coords = []
annulus_coords = []
ri = 20
ro = 30
multiplier=4

for x in range(-ro, ro, 1):
    for y in range(-ro, ro, 1):
        if (x**2 + y**2)**0.5 <= ro:
            if (x**2 + y**2)**0.5 <= ri:
                disk_coords.append([x,y])
            else:
                annulus_coords.append([x,y])

for index in range(num_images):
    image = fitslist[index]['IMAGE']
    image = subtract_median(image)
    background_estimate(image, index)
    potential_trutharray = threshold(multiplier, image, index, draw=False)  
    potential_source_regions, potential_centroids = find_regions(potential_trutharray, draw=False)
    neighborhood_means = neighborhood_check(potential_centroids, potential_trutharray, disk_coords, image)
    good_neighbor_centroids = []
    for i,m in enumerate(neighborhood_means):
        if abs(m) > 8 or potential_source_regions[i]['area'] > 5:
            good_neighbor_centroids.append(potential_centroids[i])    
    star_counts = ap_phot(fitslist[index]['IMAGE'], good_neighbor_centroids, disk_coords, annulus_coords)
    #print(f'{star_counts}')
    print(fitslist[index]['COMMENT'])
    #print(fitslist[index]['NOISESTD'])
    circle_phot_stars(fitslist, index, ri, ro, star_counts, good_neighbor_centroids)

#%% Run star finding test
'''
num_images = len(fitslist)
#num_images = 9
list_true_positives = []
list_false_positives = []
list_false_negatives = []
list_times = []
list_true_positives_ncheck = []
list_false_positives_ncheck = []
list_false_negatives_ncheck = []
list_times_ncheck = []
#multipliers = np.arange(3.2, 4.7, step=0.1)
multipliers = [4]

disk_coords = []
annulus_coords = []
ri = 5
ro = 30

for x in range(-ro, ro, 1):
    for y in range(-ro, ro, 1):
        if (x**2 + y**2)**0.5 <= ro:
            if (x**2 + y**2)**0.5 <= ri:
                disk_coords.append([x,y])
            else:
                annulus_coords.append([x,y])

options = [True]
for option in options:
    for multiplier in multipliers:
        num_true_positives = 0
        num_false_positives = 0
        num_false_negatives = 0
        times = []
        total_stars = 0
        for index in range(num_images):
            index+=0
            
            t1 = time.time()
            
            image = fitslist[index]['IMAGE']
            #total_stars += float(fitslist[index]['NUMSTARS'])
            
            image = subtract_median(image)
            #image = subtract_median_filtered_image(image)
            background_estimate(image, index)
            
            potential_trutharray = threshold(multiplier, image, index, draw=False)
            
            # Optional remove single pixels and match filter
            if False:
                #get singles
                single_pixel_regions = find_single_regions(potential_trutharray)
                #set singles to zero
                #image = [a*(not b) for a,b in zip(image.ravel(), np.array(single_pixel_regions).ravel())]
                image = np.multiply(image, np.logical_not(np.array(single_pixel_regions)))
                #image = image.reshape((3379,2703))
                x, y = np.mgrid[-4:4:1, -4:4:1]
                pos = np.dstack((x, y))
                rv = multivariate_normal([0, 0], [[2, 0], [0, 2]])
                kernel = rv.pdf(pos)
                image = matched_filter(image, kernel)
                #threshold again
                potential_trutharray = threshold(2, image, index, draw=False)
            
            potential_source_regions, potential_centroids = find_regions(potential_trutharray, draw=False)
            
            # Optional check neighborhood
            if option:
                neighborhood_means = neighborhood_check(potential_centroids, potential_trutharray, disk_coords, image)
                
                good_neighbor_centroids = []
                for i,m in enumerate(neighborhood_means):
                    if abs(m) > 8 or potential_source_regions[i]['area'] > 5:
                        good_neighbor_centroids.append(potential_centroids[i])
                #print(neighborhood_means)
                #print(potential_centroids)
                #print(good_neighbor_centroids)
            else:
                good_neighbor_centroids = potential_centroids
            
            t2 = time.time()
            times.append(t2-t1)
            stars = get_star_vals(fitslist[index]['STARS'], fitslist[index]['CAMERA'])
            true_positives, false_positives, false_negatives, showit = compare_results(stars, good_neighbor_centroids, index)
            #if showit:
            circle_real_stars(fitslist,index, good_neighbor_centroids, image)
            num_true_positives += len(true_positives)
            num_false_positives += len(false_positives)
            num_false_negatives += len(false_negatives)
        print(f'Total real stars {total_stars}')
        print(f'True positives: {num_true_positives}')
        print(f'False positives: {num_false_positives}')
        print(f'False negatives: {num_false_negatives}')
        print(f'runtime for 1 image is {np.mean(times)}')
        
        print(f'Threshold noise multiplier is {multiplier}')
        
        if option:
            list_true_positives_ncheck.append(num_true_positives)
            list_false_positives_ncheck.append(num_false_positives)
            list_false_negatives_ncheck.append(num_false_negatives)
            list_times_ncheck.append(np.mean(times))
        else:
            list_true_positives.append(num_true_positives)
            list_false_positives.append(num_false_positives)
            list_false_negatives.append(num_false_negatives)
            list_times.append(np.mean(times))
'''
"""