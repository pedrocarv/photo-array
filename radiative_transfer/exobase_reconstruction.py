import numpy as np
from scipy.interpolate import RectSphereBivariateSpline
from sklearn.metrics.pairwise import haversine_distances
import glide.science.radiative_transfer.rt_msis as rt_msis

import warnings
warnings.filterwarnings("ignore")

class GlobalExobase:
    """
        ---------------------------------------------------------------------------------
        This class implements the global H exobase reconstruction algorithm using the
        retrieved Hexo contraints at the limb

        author: Pratik Joshi (ppjoshi2@illinois.edu)
        ----------------------------------------------------------------------------------
    """
    
    def __init__(self, solar_activity, tpgse, hexo):
        """
            Args:
                solar_activity (SolarActivity) = SolarActivity object containing geomagnetic indices
                tpgse (ndarray) = Latitude, longitude tangent points for all wedges
                hexo (ndarray) = Exobase H density for all wedges
        """
        self.solar_activity = solar_activity
        self.tpgse = tpgse
        
        #Set exobase for each wedge
        msis_obj = rt_msis.MSIS(self.solar_activity, tpgse[0], tpgse[1])
        self.wedge_exo = msis_obj.exobase()
        self.h_exo = hexo
    
    @property
    def nwedges(self):
        return self.wedge_exo.npt
        
            
    def calc_global_hexo(self, gselat, gselon, power=.5):
    
        """
        Reconstruct the global exobase H density

        Args:
            gselat (ndarray): GSE latitude grid (degrees)
            gselon (ndarray): GSE longitude grid (degrees)

        Returns:
            merged_hexo_2d (ndarray): reconstructed global exobase H density in cm^-3
        """
        
        global_msis = rt_msis.MSIS(self.solar_activity, gselat, gselon, do_grid=True)
        global_exo = global_msis.exobase()
        
        ngselat = gselat.size
        ngselon = gselon.size
        
        # Calculate the a parameter -------------------------------------------------------------------------------------------
        
        apar = global_exo.h_exo * (global_exo.t_exo**2.5)

        ## Calculate b parameter -----------------------------------------------------------------------------------------------
        
        apar_wedge = self.wedge_exo.h_exo * (self.wedge_exo.t_exo**2.5)
        bpar = (np.log10(apar_wedge) - np.log10(self.h_exo)) / np.log10(self.wedge_exo.t_exo)

        # Calculate 2D Hexo for each wedge using a and b parameters ------------------------------------------------------------
        hexowedge = np.zeros((self.nwedges, ngselat, ngselon))
        for w in range(self.nwedges):
            hexowedge[w,:,:] = apar / (global_exo.t_exo**bpar[w])
            
        ## Merge 2D Hexo from 24 wedges using inverse distance weighted average function -----------------------------------------     
        points = np.stack([global_exo.lat, global_exo.lon], axis=0)
        weights_arr = np.zeros((self.nwedges, ngselat, ngselon))
        values_arr = np.zeros((self.nwedges, ngselat, ngselon))

        for w in range(self.nwedges):
            weights = distance_weights_haversine(points, self.tpgse[:,w], power=power)
            weights_arr[w] = weights
            values_arr[w] = hexowedge[w]
            
        weighted_sum = np.sum(weights_arr*values_arr,axis=0)
        total_weight = np.sum(weights_arr,axis=0)
        merged_hexo = weighted_sum/total_weight
        
        return merged_hexo

def distance_weights_haversine(points, target_point, power=.5):
    """
    Calculates the distance-dependent weights at a target point using haversine distance.
    Returns 1/(distances ** power).

    Args:
        points (ndarray): Lat, lon array (2 x nlat x nlon).
        target_point (ndarray): Lat, lon point to compute distances to.
        power (float, optional): Exponent used in distance weighting. Defaults to .5.

    Returns:
        weights (ndarray): Weights array.

    """
    
    lat = np.ravel(points[0]) * np.pi / 180
    lon = np.ravel(points[1]) * np.pi / 180
    
    distances = haversine_distances(np.vstack((lat, lon)).T, np.reshape(target_point * np.pi/180, (1,2)))
    distances = np.reshape(distances, np.shape(points[0]))
    weights = 1 / (distances ** power)

    return weights

## ******************************************* Main call ********************************************************

if __name__=="__main__":

    ## User defined date input as yyyy, mm, dd command line arguments
    # year = int(sys.argv[1])
    # month = int(sys.argv[2])
    # day = int(sys.argv[3])
    # uthh = int(sys.argv[4])
    # utmm = int(sys.argv[5])
    gselat = np.linspace(-90,90,19) # matches Evan's gselon grid of total density above 3 Re
    gselon = np.linspace(-180,180,37) # matches Evan's gselon grid of total density above 3 Re
    # ngselat = len(gselat) # number of gselats
    # ngselon = len(gselon) # number of gselons

    # Calculate the date array in Numpy 64 ISO format required for background atmosphere calculations
    # py_datetime = datetime.datetime(year, month, day, uthh, utmm)
    # date = np.datetime64(py_datetime)
    date = np.datetime64('2001-08-15T06:00:00')

    # Get geomagnetic indices for the given date
    # dstd, kp, aps, f107, f107a, datetime_dst, datetime_kp_aps = geomagnetic_indices(date)
    # msis_version=2.1 #msis model version
    aps = 3 * np.ones((1,7))
    f107 = 210
    f107a = 210
        
    # Get limb retrieval parameters for the given date
   
    gselat_wedge = np.array([-7.46304341, -22.30589489, -37.11437768, -51.83217264, -66.27426059,-79.2435324,-79.25898342,-66.29455618,-51.85177357,-37.13214555,
                    -22.32123879, -7.47562649, 7.37761445, 22.21705515, 37.0131819, 51.70176491, 66.07105525, 78.79961508, 78.79975748, 66.07195768, 51.70411914,
                    37.01758453, 22.22396608, 7.38731662])
    
    gselon_wedge = np.array([84.44887256, 86.68556306, 89.4410583, 93.53490112, 101.61767089, 129.24482077, -141.9249623, -114.22957886, -106.13598766, -102.0388846,
                    -99.28201398, -97.04462436, -94.93672137, -92.65662081, -89.79841002, -85.49346307, -76.94464324, -48.51717889, 35.92240194, 64.34983947, 72.89860143,
                    77.20345208, 80.06150946, 82.34136751])
    tpgse = np.vstack([gselat_wedge, gselon_wedge])
    
    
    zexo_wedge = np.array([547, 533, 520, 510, 480, 469, 466, 464, 453, 451, 444, 443, 455, 462, 474, 477, 491, 510, 532, 545, 550, 546, 553, 553])
    
    texo_wedge = np.array([1321.48181152, 1265.87341309, 1196.96386719, 1136.0090332, 1001.85577393, 973.27355957, 986.8560791, 993.5925293, 950.19830322, 927.1920166, 922.88067627,
                  937.47210693, 1028.64807129, 1076.20336914, 1160.22436523, 1185.20007324, 1253.38220215, 1310.4597168, 1349.66809082, 1369.36621094, 1385.24731445,
                  1381.37060547, 1396.24609375, 1378.13671875])
            
    hexo_wedge = np.array([63763.999560, 70733.895450, 84025.523680, 92527.837405, 99394.099165, 97187.454365, 106043.778128, 127885.783959, 157089.213202, 142404.978842, 116335.945338,
                   109535.317567,91441.796965, 82280.031307, 66458.638213, 61133.820673, 46641.369965, 38085.103243, 35193.667080, 36999.178064, 46300.186325, 51788.372057,
                   53870.261131, 56506.995493])
    
    # solar flux in ph/s/cm2/A
    solar_flux = 7e11
            
    # Create object of the GlobalHexo class
    solar_activity = rt_msis.SolarActivity(date, aps, f107, f107a)
    obj = GlobalExobase(solar_activity, tpgse, hexo_wedge)

    # Get global Hexo by call to class function
    global_hexo = obj.calc_global_hexo(gselat, gselon)
    
    # Print the scalars and array size of all parameters
    print('date =',date)
    print('ap =',aps)
    print('f107 =',f107)
    print('f107a =',f107a)
    print('number of wedges = ', obj.nwedges)
    print('gselat grid size = ', gselat.shape[0])
    print('gselon grid size = ', gselon.shape[0])
    print('reconstructed global hexo size =', global_hexo.shape)
    
    # # Save to NETCDF
    # p = Path(__file__).parent
    # print(p)
    # ph = Path(p, "rt_com", "input_data")
    # RE=6371 #km
    # da = xr.DataArray(np.asarray((obj.global_zexo+RE)/RE))
    # ds = da.to_dataset(name ='REXO')
    # ds['REXO']= ds.REXO.rename({'dim_0': 'el','dim_1':'az'})
    # ds = ds.assign(TEXO = (['el','az'],obj.global_texo))
    # ds = ds.assign(HEXO = (['el','az'],obj.global_hexo))
    # ds = ds.assign(GSELAT = (['el'],gselat))
    # ds = ds.assign(GSELON = (['az'],gselon))
    # ds = ds.assign(SF = obj.solar_flux)
    # ds = ds.assign(DATE = obj.date)
    # ds.to_netcdf(Path(ph,'global_hexo.nc'),mode='w')

