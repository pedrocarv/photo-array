'''
File containing all calibration steps (e.g. bias, dark, EPS removal, etc.)
'''

import numpy as np
import matplotlib.pyplot as plt
import xarray as xr

from scip.data_processing.read_fits_files import *
from scip.data_processing.observation_analysis import *

def inst_cal(on_band, off_band, start_ind=0, end_ind=None, binoc=2, loc=None, fit_bkg=False, out_path=None):
    '''
    Perform instrument calibration on passed on- & off-band measurements

    Versioning:
        v0.0- Uniform Dark, Bias, and Hot Pixel removal
    
    Args:
        on_band (dict) - Dictionary containing all on-band information
        off_band (dict) - Dictionary containing all off-band information
        start_ind (int) - index to process observations from
        end_ind (int) - index to end processing at
        binoc (int) - binocular pair being processed (0 = 001&002, 2 = 003&004, 4 = 005&006)
        fit_bkg (bool) - whether to apply linear fit to background signal (Counters wackiness of of fm2b003's chip cooler)
    '''
    
    if end_ind is None:
        end_ind = len(on_band) - 1

    # Get observation info
    s_alts = np.array([get_shadowalt(fitsdate_to_dt(on_band[i]['DATE-OBS']), loc=loc) for i in range(start_ind,end_ind+1)])

    times = np.array([fitsdate_to_dt(on_band[i]['DATE-OBS']) for i in range(start_ind,end_ind+1)])

    on_exp = np.array([on_band[i]['EXP-TIME'] for i in range(start_ind, end_ind+1)])
    off_exp = np.array([off_band[i]['EXP-TIME'] for i in range(start_ind, end_ind+1)])
    
    
    # Load calibration frames
    # FIXME: Don't hardcode the calibration directory
    temps = -0.5 * np.arange(14)
    path = Path('/home/jc/research/photometer/scip-2/scip/calibration/data/bias_fm2b005006_20251119')
    bias_on_frames = []
    bias_off_frames = []

    for t in temps:
        bias_on_frames += [np.load(path / f'onband_-6_10.npy')]
        bias_off_frames += [np.load(path / f'offband_-6_10.npy')]
    
    on_temps = np.array([on_band[i]['CCD-TEMP'] for i in range(len(on_band))])
    off_temps = np.array([off_band[i]['CCD-TEMP'] for i in range(len(off_band))])
    on_b_ind, off_b_ind = np.abs(np.round(on_temps * 2)).astype(int), np.abs(np.round(off_temps * 2)).astype(int)

    # Variable Dark Subtraction
    on_ims, off_ims = [], []
    dpath = Path('/home/jc/research/photometer/scip-2/scip/calibration/data/dark_fm2b003004_20260125')
    temps = -1 * np.arange(7)
    for b in temps:
        on_ims += [np.load(dpath/f'onband_{b}_10.npy')]
        off_ims += [np.load(dpath/f'offband_{b}_10.npy')]

    t_int = 30*60 # 30 minute integration time for each calibration image

    dark_on = np.array([np.mean(im)/t_int for im in on_ims])
    dark_off = np.array([np.mean(im)/t_int for im in off_ims])

    on_d_ind, off_d_ind = np.abs(np.round(on_temps)).astype(int), np.abs(np.round(off_temps)).astype(int)
    # If below operating temp, set dark frame to operating temp
    on_d_ind[on_d_ind > 6] = 6
    off_d_ind[off_d_ind > 6] = 6

    stencil =  np.ones(on_band[0]['image'].shape)

    filter_transmissions = np.load('/home/jc/research/photometer/scip-2/scip/calibration/data/filter_transmission/total_filter_transmission.npy')
    t_ratio = filter_transmissions[binoc+1]/filter_transmissions[binoc]


    on_adu = np.array([np.sum(extreme_pixel_removal(on_band[i]['image']- bias_on_frames[on_b_ind[i]] - dark_on[on_d_ind[i]]*stencil*on_band[i]['EXP-TIME'],passes=2))/on_band[i]['EXP-TIME']  for i in range(start_ind, end_ind+1)])
    off_adu = np.array([np.sum(extreme_pixel_removal(off_band[i]['image'] - bias_off_frames[off_b_ind[i]] - dark_off[off_d_ind[i]]*stencil*off_band[i]['EXP-TIME'],passes=2)/off_band[i]['EXP-TIME']) for i in range(start_ind, end_ind+1)])

    
    if fit_bkg:
        off_fit = np.polyfit(s_alts, off_adu, 1) # Linear fit
        off_adu = off_fit[0] * s_alts + off_fit[1]

    '''
    FIXME: Clean l8r
    # Remove indicies with galactic bkg
    rem_wham = True
    if rem_wham:
        eloc = EarthLocation(lon=loc[0]*u.deg, lat=loc[1]*u.deg, height=loc[2]*u.m)
            # Pull galatic coordinates for zenith observations
        lons = []
        lats = []
        lons_deg = []
        lats_deg = []
        mid = np.datetime64('2026-01-23T00:00:00')
        time_since = []
        for dt in times:
            time = Time(dt)
            zen_coords = AltAz(az=180*u.deg, alt=90*u.deg, obstime=time, location=eloc)
            zen_g = SkyCoord(zen_coords).transform_to('galactic')
            lons += [-zen_g.l.wrap_at('180d').rad]
            lats += [zen_g.b.rad]
            lons_deg += [-zen_g.l.wrap_at('180d').deg]
            lats_deg += [zen_g.b.deg]
            time_since += [(dt - mid).astype('int') / 3600]

        # Load wham file to pandas dataframe
        wham_path = '/home/jc/research/photometer/scip-2/data/wham-ss-DR1-v161116-170912-int.txt'
        wham_df = pd.read_csv(wham_path, sep='\s+', skipinitialspace=True, skiprows=45)

        # Pull relevant data to numpy arrays
        inten = wham_df['INTEN'].to_numpy()
        gal_coords = wham_df[['GAL-LON', 'GAL-LAT']].to_numpy()
        gal_coords = SkyCoord(gal_coords[:, 0], gal_coords[:,1], frame='galactic', unit = u.deg) # Easier to work with astropy
        l = -gal_coords.l.wrap_at('180d') # Galactic Longitude NOTE: Reversed!
        b = gal_coords.b # Galactic latitude

        # Create map from points (NOTE: WHAM has 1deg resolution)
        # Boundaries of the map will be at whole degrees (i.e. -180, -179)
        # thus, the value at [0,0] will contain the average intensity measured between lats -180:-179 and lons -90:-89
        # value at [1,0] will contain data from lats -180:-179 and lons -89:-88
        # matrix should have dimension of (180, 360)
        # ==> To accomplish this mapping, we round down and add either 180 or 90
        g_lon = (np.pi/180) * np.linspace(-180, 180, 362)[1:-1]
        g_lat = (np.pi/180) * np.linspace(-90, 90, 182)[1:-1]
        L, B = np.meshgrid(g_lon, g_lat) # Meshgrids in radians

        l_ind = (np.floor(np.array(l)) + 180).astype(int)
        b_ind = (np.floor(np.array(b)) + 90).astype(int)

        ha_map = np.zeros((180,360))
        ha_map[b_ind, l_ind] = inten

        # indicies of observations
        lon_ind = (np.floor(np.array(lons_deg)) + 180).astype(int)
        lat_ind = (np.floor(np.array(lats_deg)) + 90).astype(int)

        obs_bkg = [ha_map[lat_ind[i], lon_ind[i]] for i in range(end_ind+1 - start_ind)]
    '''
    # Redefine here so they save to netcdf properly
    on_temps = np.array([on_band[i]['CCD-TEMP'] for i in range(start_ind, end_ind+1)])
    off_temps = np.array([off_band[i]['CCD-TEMP'] for i in range(start_ind, end_ind+1)])



    obs_adu = on_adu - t_ratio * off_adu

    if out_path is not None:
        ds = xr.Dataset(
            {
                'obs_adu':  ('time', obs_adu),
                'on_adu':   ('time', on_adu),
                'off_adu':  ('time', off_adu),
                'on_exp':   ('time', on_exp),
                'off_exp':  ('time', off_exp),
                's_alts':   ('time', s_alts),
                'on_temps': ('time', on_temps),
                'off_temps':('time', off_temps),
            },
            coords={'time': times},
        )
        ds.to_netcdf(out_path)

    return obs_adu

def extreme_pixel_removal(im, tile_size=9, thresh=2, passes=1):
    '''
    Remove extreme pixel values (e.g. hot pixels, cosmic ray hits) from passed image
    
    :input:
    im (ndarray) - Observation image
    tile_size (int) - tile size to scan over
    thresh (int) - # of stds above mean pixel has to be to remove
    '''
    h,w = im.shape
    for p in range(passes):
        for i in range(h):
            for j in range(w):
                h_ind = i*tile_size
                w_ind = j*tile_size
                tile = im[h_ind:h_ind+tile_size, w_ind:w_ind+tile_size]
                hot = np.abs(tile) > np.median(tile) + thresh*np.std(tile)
                im[h_ind:h_ind+tile_size, w_ind:w_ind+tile_size][hot] = np.median(tile[~hot])

    return im
 
# FIXME: Remove if you're nt doing anything crazy
def dark_removal():
    pass

def bias_removal():
    pass

def crop_im(im, dim1, dim2=None, shape='rect', center=(0,0)):
    '''
    Crop/zero image to specified dimensions

    :input:
    im (np.ndarray) - Image to crop
    dim1 (int) - # # of pixels to crop y dimension to / radius of circular crop
    dim2 (int) - # of pixels to crop x dimension to (if None, square crop with 1st specified dimension)
    shape (str) - shape of area to crop out (rect or circ)
    center (tuple:int,int) - center index of crop area (x,y)
    :output:
    (ndarray) - cropped image (ydim x xdim)
    '''
    m,n = im.shape

    if shape == 'rect':
        if dim2 is None:
            dim2 = dim1

        crop_y = m - dim1
        crop_x = n - dim2

        return im[crop_y//2:-1*(int)(np.ceil(crop_y/2)), crop_x//2:-1*(int)(np.ceil(crop_x//2))]
    
    elif shape == 'circle':
        # Circular crop
        x = np.arange(-1*n//2,n//2).reshape(n,1)
        y = np.arange(m//2, -1*m//2,-1).reshape(m,1)

        mask = ((y-center[1])**2) + ((x-center[0])**2).T <= dim1**2

        return im * mask

    
    else:
        raise Exception('Invalid crop shape')

def mark_area(im, dimy, dimx=None):
    '''
    Mark cropped area in image (as specified by passed dims)
    
    :input:
    im (np.ndarray) - Image to crop
    dimy (int) - # # of pixels to crop y dimension to
    dimx (int) - # of pixels to crop x dimension to (if None, square crop with 1st specified dimension)
    :output:
    (matplotlib.pyplot Figure) - image with cropped area marked
    '''

    if dimx is None:
        dimx = dimy

    crop_y = im.shape[0] - dimy
    crop_x = im.shape[1] - dimx

    y_inds = [crop_y//2, im.shape[0]-(int)(np.ceil(crop_y/2))]
    x_inds = [crop_x//2, im.shape[1]-(int)(np.ceil(crop_x/2))]

    plt.plot([x_inds[0], x_inds[0]], [y_inds[0],y_inds[1]],[x_inds[0], x_inds[1]], [y_inds[1],y_inds[1]],[x_inds[1], x_inds[1]], [y_inds[1],y_inds[0]],[x_inds[1], x_inds[0]], [y_inds[0],y_inds[0]], color='red', linewidth=1.5)
    plt.imshow(im, interpolation='none', vmin=np.percentile(im, 1), vmax=np.percentile(im, 99))
    plt.colorbar()

    return plt.gcf()

