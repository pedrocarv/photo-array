# Generate albedo netcdf from given albedo txt file

import xarray as xr
import numpy as np
import csv

# params
r = 32 # grid radial resolution
sza = 32 # grid azimuthal resolution
h_exo = '5e4' # H density used to generate albedo
season = 'Fall' # Season used to generate albedo
re = True # elevation has + 1Re?

a_path = 'rt_com/results/albedo_20250930_5e4_7e11_32_RE.txt' # Path to albedo used
nc_path = f'albedo_r{r}_sza{sza}_{h_exo}_{season}.nc' # location & filename of albedo netcdf

f = open(a_path, 'r')
reader = csv.reader(f, delimiter=',')
alb = []
for row in reader:
    alb.append([float(i) for i in row])

alb = np.array(alb)

# R and SZA grids NOTE: Change these if using different resolution
zgrid = np.array([74.5, 76.0, 78.5, 82.0, 86.5, 92.0, 98.5, 106.0, 116.5, 133.5,
                          161.0, 205.5, 277.5, 394.0, 582.5, 887.5, 1381.0, 2179.5, 3471.5, 5562.0,
                          8945.0, 14417.0, 23273.0, 37601.0, 60785.0, 98299.0, 158999.0, 257199.0, 416109.0, 673219.0,
                          1089229.0, 1762329.0])

if re:
    zgrid += 6371

sza_grid = np.linspace(0, 180, num=32)

ds = xr.Dataset(
    data_vars=dict(
        ALBEDO = (['r','sza'], alb),
        R = (['r'],zgrid),
        SZA = (['sza'],sza_grid)
    )
)

ds.to_netcdf(nc_path)