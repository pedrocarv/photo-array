# -*- coding: utf-8 -*-
"""
Created on Mon May 10 13:45:57 2021

@author: dawnh
"""
import matplotlib.pyplot as plt
import numpy as np

# Data source: UBVRIJKLMNH Photoelectric Catalogue (Morel+ 1978)
# https://vizier.u-strasbg.fr/viz-bin/VizieR?-source=II/7A

# Header:
# JP11 LID n_LID I ref RA DE

lines = []
file = open('asu (1).tsv', 'r') 
for line in file:
    if line[0] == "#":
        continue
    linesplit = line.split(";")
    if len(linesplit) == 7:
        lines.append(line.split(';'))
# remove header
lines = lines[3:]
# remove stars without Johnson I magnitude         
lines = [a for a in lines if a[3] != '      ']
# remove variable stars
lines = [a for a in lines if a[2] == ' ']
# limiting magnitude
lines = [a for a in lines if float(a[3]) < 1]
RA = []
DEC = []
# Format values as numbers
for line in lines:
    line[0] = float(line[0]) #index
    line[1] = float(line[1]) # identifier
    line[2] = 0
    line[3] = float(line[3]) # empty
    line[4] = float(line[4]) # I magnitude
    line[5] = float(line[5]) # RA
    RA.append(line[5])
    line[6] = float(line[6][:-3])
    DEC.append(line[6]) # DEC
stardata = np.array(lines,dtype=float)    
np.save("stardata.npy", stardata)
plt.scatter(RA, DEC)