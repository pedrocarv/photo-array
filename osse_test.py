'''
Script to test the entire OSSE workflow.

JC
'''
#%% imports
import numpy as np
from time import time

from osse.Photometer import Photometer
from osse.Scene import Scene
from scip.data_processing.helpers import geomagnetic_data
from radiative_transfer.rt_inversion import optimize_profile

#%% init
obs_date = np.datetime64('now') - np.timedelta64(1, 'D') # [UTC]
t_int = 3 * 60 # [s]
t_int2 = 1 * 30 # [s]
phot = Photometer() # Photometer @ Aeronomy

# TODO Exosphere initialization (Scene needs a major overhaul)
hexo = 5e4
sflux = 9e9
bkg = 12 # Background BB Signal [R]
solar_activity = geomagnetic_data(obs_date)
exo = Scene(hexo=hexo, solar_flux=sflux, solar_activity=solar_activity)
exo2 = Scene(hexo=hexo, solar_flux=sflux, solar_activity=solar_activity, twoD=True)

#%% sched
point = 'zenith'
kind = 'dusk' # Don't schedule full night
schedule = phot.schedule_obs(obs_date, point, t_int, profile=kind)
schedule2 = phot.schedule_obs(obs_date, point, t_int2, profile=kind)

#%% scene_rad
rt_los = phot.rt_los(schedule)
rt_los2 = phot.rt_los(schedule2)


t0 = time()
rad = exo(rt_los)
tf  =time()
print('Forward time: ', tf-t0)
# rad2 = exo(rt_los2)
rad2 = exo2(rt_los)

#%% noise
inst = phot(rad, bkg, t_int, exo.exo_obj.t_exo[0])
inst2 = phot(rad2, bkg, t_int, exo.exo_obj.t_exo[0])

#%% retr
# print('Retrieving Spherical')
t0 = time()
retr = optimize_profile([0,90], inst, np.array(rt_los), sflux, solar_activity, 'Bishop')
tf  =time()
print('Retrieval time: ', tf-t0)
# print('Retrieving Assymetric')
#retr2 = optimize_profile([0,90], inst2, np.array(rt_los), sflux, solar_activity, 'Bishop')

# print('Spherical Exosphere: Hexo=%f; SFlux = %fe9'%(retr[0], retr[1] / 1e9))
# print('Assymetric Exosphere: Hexo=%f; SFlux = %fe9'%(retr2[0], retr2[1] / 1e9))

#%% plot
import matplotlib.pyplot as plt

plot_rads = True
if plot_rads:
    # Plot unnormalized radiance
    # plt.plot(schedule['observation_time'].to_numpy(), rad, label='1D')
    # plt.plot(schedule['observation_time'].to_numpy(), rad2, label='2D')
    # plt.ylabel('Radiance [R]')

    # Plot normalized radiance
    # plt.plot(schedule['observation_time'].to_numpy(), rad/np.max(rad), label='1D')
    # plt.plot(schedule['observation_time'].to_numpy(), rad2/np.max(rad2), label='2D')
    plt.plot(schedule['observation_time'].to_numpy(), inst/np.max(inst), label='1D')
    plt.plot(schedule['observation_time'].to_numpy(), inst2/np.max(inst2), label='2D')
    plt.ylabel('Max-Normalized Radiance')

    # Plot instrument profile
    # plt.plot(schedule['observation_time'].to_numpy(), inst, label='1D')
    # plt.plot(schedule['observation_time'].to_numpy(), inst2, label='2D')
    # plt.ylabel('Measurement [ADU/s]')

    plt.legend()
    plt.grid()
    plt.show()
