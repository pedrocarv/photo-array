'''
All plotting scripts for osse results and scheduling visualization
'''

import matplotlib.pyplot as plt

def plot_sched_alt_az(sched):
    '''
    Plot of Altitude and Azimuth from passed list of schedule tuples (obstime, alt, az)
    '''
    times = [s[0] for s in sched]
    alts = [s[1] for s in sched]
    azs = [s[2] for s in sched]