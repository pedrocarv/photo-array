'''
LOS visualization
for a given observation location, time, and oza
TODO: given LOS specification as...
 dist from earth center | SZA | OZA | Az

By JC
'''

import plotly.graph_objects as go
import numpy as np
from astropy.coordinates import get_sun
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import AltAz, EarthLocation, GCRS, ITRS
from astropy.coordinates import SkyCoord
import sunpy.coordinates as scoord

from scip.data_processing.helpers import *
from osse.Photometer import Photometer


# Input Params
# date = '2026-01-19T00:00:00' #YYYY-MM-DDTHH:MM:SS
# date = '2026-03-20T00:00:00' # Spring Equinox
# date = '2026-09-23T03:00:00'
date = '2025-10-23T00:30:00'
# date = '2026-06-21T03:00:00' # Summer Solstice
# date = '2026-12-21T00:00:00' # Winter Solstice
dur = 5 # Duration [hr] of observation session
points = 100 # Number of LOSs to plot (over 100 gets laggy)


# lon, lat, elev = get_loc_info('aeronomy') # Or manually define latitude longitude and elevation


# NOTE: OZA is defined below (can require knowledge of location and date e.g. for antisolar pointing)

# Constants et al.
Re = 6371 # (Approx) earth radius in km
dates = [np.datetime64(date) + np.timedelta64(int(i * (dur*3600/points)), 's') for i in range(points)] * 3 # More than 100 points gets laggy
# dates = [np.datetime64(date)]
pointing = ['zenith'] * points + ['antisolar'] * points + ['solar'] * points


photom = Photometer() # Defaults to the aeronomy field site
loc = photom.loc

# locs = [(0,0,0)] * points

# dates = []
# locs = []
# dates = [np.datetime64(date) + np.timedelta64(3,'h')] * 3
# locs = [(get_loc_info('aeronomy')), (get_loc_info('kitt')), (get_loc_info('millstone'))]
# locs = [(0,0,0), (0,40,0), (0,-40,0)]
# dates = [np.datetime64(date)] * 3

# Add cylindrical Earth shadow
theta, phi = np.mgrid[0:np.pi:100j, 0:2*np.pi:100j]

x_shdw, ang = np.mgrid[-3*Re:0:100j, 0:2*np.pi:100j]
# NOTE: Including Screen Height
y_shdw = (Re+105) * np.cos(ang)
z_shdw = (Re+105) * np.sin(ang)

fig = go.Figure()

fig.add_trace(go.Surface(x=Re*np.sin(theta)*np.cos(phi),
                        y=Re*np.sin(theta)*np.sin(phi),
                        z=Re*np.cos(theta)*np.ones(np.shape(phi)),
                        colorscale=['blue','blue']))


fig.add_trace(go.Surface(x=x_shdw, y=y_shdw, z=z_shdw, colorscale=['gray','gray'],opacity=0.65))

# Add earth rotational axis and equator
# NOTE: Assumes that the GSE coordinate system doesn't vary too much over the course of a night
sol = get_sun(Time(date))
r_ax = np.linspace(-1.5*Re, 1.5*Re, 100)
tilt = 23.44 * (np.pi/180)
r_az = (sol.ra.deg-90) * (np.pi/180)

x_ax = r_ax*np.sin(tilt)*np.cos(r_az)
y_ax = r_ax*np.sin(tilt)*np.sin(r_az)
z_ax = r_ax*np.cos(tilt)

# fig.add_scatter3d(x=x_ax,y=y_ax,z=z_ax,mode='lines',line=dict(color='orange',width=5))

# Trace north and south pole, get GSE coords
r_ax = np.linspace(0, 0.5*Re, 10)
t0 = Time(dates[0], scale='utc')

# North
alocs = [EarthLocation(lat=90*u.deg, lon=0*u.deg, height=r*u.km) for r in r_ax] 

north_itrs = [aloc.get_itrs(obstime=t0) for aloc in alocs]

# Transform to GSE
north_gse = []
gse_frame = scoord.frames.GeocentricSolarEcliptic(obstime=t0)
for i in range(len(r_ax)):
    pole = north_itrs[i].transform_to(gse_frame)
    pole.representation_type = 'cartesian'
    north_gse += [pole.cartesian.xyz]  # ~6371 km offset from Earth center

north_gse = np.array(north_gse)
fig.add_scatter3d(x=north_gse[:,0],y=north_gse[:,1],z=north_gse[:,2],mode='lines',line=dict(color='orange',width=5))

# South
alocs = [EarthLocation(lat=-90*u.deg, lon=0*u.deg, height=r*u.km) for r in r_ax] 

south_itrs = [aloc.get_itrs(obstime=t0) for aloc in alocs]

# Transform to GSE
south_gse = []
gse_frame = scoord.frames.GeocentricSolarEcliptic(obstime=t0)
for i in range(len(r_ax)):
    pole = south_itrs[i].transform_to(gse_frame)
    pole.representation_type = 'cartesian'
    south_gse += [pole.cartesian.xyz]  # ~6371 km offset from Earth center

south_gse = np.array(south_gse)
fig.add_scatter3d(x=south_gse[:,0],y=south_gse[:,1],z=south_gse[:,2],mode='lines',line=dict(color='orange',width=5))

eq_az = np.linspace(-180, 180, 100)
# Trace equator, get GSE coords
alocs = [EarthLocation(lat=0*u.deg, lon=lon*u.deg, height=0*u.km) for lon in eq_az] 
t0 = Time(dates[0], scale='utc')

eq_itrs = [aloc.get_itrs(obstime=t0) for aloc in alocs]

# Transform to GSE
eq_gse = []
gse_frame = scoord.frames.GeocentricSolarEcliptic(obstime=t0)
for i in range(len(eq_az)):
    eq = eq_itrs[i].transform_to(gse_frame)
    eq.representation_type = 'cartesian'
    eq_gse += [eq.cartesian.xyz]  # ~6371 km offset from Earth center

eq_gse = np.array(eq_gse)

fig.add_scatter3d(x=eq_gse[:,0],y=eq_gse[:,1],z=eq_gse[:,2],mode='lines',line=dict(color='gray',width=10,dash='dash'))

# Define LOS angles of interest
# ozas = ((lat - sol.dec.deg) * (np.pi/180) - np.pi/2) * np.ones(len(dates)) # Off-Zenith Angle

# sol_alt = np.array([get_solalt(dates[i], unit='deg') for i in range(len(dates))])
# szas = np.array([get_sza(dates[i]) for i in range(len(dates))])
# sol_az = np.array([get_solaz(dates[i], unit='deg') for i in range(len(dates))])
# has = np.array([get_solHA(dates[i], unit='deg') for i in range(len(dates))])
# # ozas = -1*(sol_alt+np.pi/2)
# ozas = 0* sol_az
# # ozas = np.array([0,60,60])

# # Compute Antisolar pointing
# alts = -sol_alt
# azs = (sol_az + 180) % 360 




# Plot LOSs
ps = []

for i in range(len(dates)):
    lon, lat, elev = photom.loc
    h = get_solHA(dates[i], [lon, lat, elev], 'deg')
    
    # Determine radius of LOS
    r = np.linspace(0,Re,100) # Plot out to 1 Re of line

    '''
    TODO:
    Find observation location as a rotation from GSE... save xyz coords for cartesian translation in GSE
    Rotate Local observation angles to be in GSE:
        +z will be rotated by zenith angle... 
        +x will be rotated by hour angle???
    Translate using saved observation location

    ... or just use astropy...
    '''
    aloc = EarthLocation(lat=lat*u.deg, lon=lon*u.deg, height=(elev * 1e-3)*u.km)
    atime = Time(dates[i], scale='utc')
    gse_frame = scoord.frames.GeocentricSolarEcliptic(obstime=atime)
    
    # Observer position in ITRS (ECEF equivalent)
    obs_itrs = aloc.get_itrs(obstime=atime)

    # Transform to GSE
    obs_gse = obs_itrs.transform_to(gse_frame)
    obs_gse.representation_type = 'cartesian'

    r_obs = obs_gse.cartesian  # ~6371 km offset from Earth center

    # Define line of sight in Alt-Az
    altaz_frame = AltAz(obstime=atime, location=aloc)

    if pointing[i] == 'zenith':
        solar_alt = 90
        solar_az = 180
        color='red'
    elif pointing[i] == 'antisolar':
        solar_alt, solar_az = photom.antisolar_pointing([dates[i]])
        solar_alt = solar_alt[0]
        solar_az = solar_az[0]
        color='green'
    else:
        solar_alt, solar_az = photom.sunward_pointing([dates[i]])
        solar_alt = solar_alt[0]
        solar_az = solar_az[0]
        color='magenta'


    los = SkyCoord(alt=solar_alt*u.deg, az=solar_az*u.deg, distance=1.495979e+8*u.km, frame=altaz_frame)

    # Convert to GSE (GeocentricMeanEcliptic = GSE when obstime provided)
    gse = los.transform_to(gse_frame)

    # Cartesian GSE components
    gse.representation_type = 'cartesian'
    # print(gse.cartesian)  # x, y, z unit vector in GSE

    los_astro = r_obs.xyz[:, np.newaxis] + gse.cartesian.xyz[:,np.newaxis] / np.linalg.norm(gse.cartesian.xyz) * r * u.km

    fig.add_scatter3d(x=los_astro[0],y=los_astro[1],z=los_astro[2],mode='lines',line=dict(color=color,width=10))

'''
END TO FIX
'''


fig.update_layout(title_text='LOS Visualization', showlegend=False)

fig.show()

fig.write_html('aeronomy_all.html')
