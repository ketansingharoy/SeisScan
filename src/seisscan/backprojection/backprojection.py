import importlib
import copy
import time
import itertools
import pickle
import numpy as np
from scipy import signal
from scipy import interpolate
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib import dates as mdates

from obspy import Stream, UTCDateTime
from obspy.geodetics.base import gps2dist_azimuth, degrees2kilometers, kilometers2degrees
from obspy.taup import TauPyModel
from obspy.signal.trigger import recursive_sta_lta

import utm

import dask
from dask.distributed import Client as dask_Client
from dask.distributed import LocalCluster as dask_LocalCluster

from .sliding_window_calculation import sliding_window_view, sliding_window_calculation
from .brightness import Brightness4




def prepare_traveltime_lookup_table(min_dist_k, max_dist_k, step_dist_k,
                                    min_dep_k, max_dep_k, step_dep_k,
                                    phase, model_name='iasp91', dask_client=None):
    """Computes seismic wave travel time table.
    
    It computes seismic waves travel times for a given epicentral distance range, depth range and 1-D earth model. These travel times are used by the backprojection program.
    
    Parameters
    ----------
    min_dist_k: float
        Minimum epicentral distance (km).
    max_dist_k: float
        Maximum epicentral distance (km).
    step_dist_k: float
        Epicentral distance increment (km).
    min_dep_k: float
        Minimum depth distance (km).
    max_dep_k: float
        Maximum depth distance (km).
    step_dep_k: float
        Depth increment (km).
    phase: str
        Seismic phase name. "P" or "S".
    model_name: str
        1-D earth model name. Default is 'iasp91'.
    dask_client: dask.Client
        A dask client for parallel processing.
        
    Returns
    -------
    dist_k_r1: numpy.ndarray of float
        An array epicentral distances (km).
    dep_k_r1: numpy.ndarray of float
        An array deths (km).
    tt_r2: numpy.ndarray of float
        An travel times (seconds).
    computation_time: float
        Total time (seconds) taken by the cmputation.
    """

    #--- model
    model = TauPyModel(model=model_name)
    
    phase_list = [phase.lower(), phase.upper()]

    #--- distance and depth array
    dist_k_r1 = np.arange(min_dist_k, max_dist_k+step_dist_k, step_dist_k)    # distance array in km
    tmp_dist_k_r1 = np.linspace(min_dist_k, max_dist_k, 100)
    dep_k_r1 = np.arange(min_dep_k, max_dep_k+step_dep_k, step_dep_k)         # distance array in km

    #--- initialize 2d array of travel times. ndist x ndep
    ndist = dist_k_r1.size
    ndep = dep_k_r1.size

    tt_r2 = np.zeros((ndist, ndep))

    #--- function to get travel time for a depth (km) and an array of distances (deg)
    def get_tt(tmp_dist_k_r1, dist_k_r1, dep_k, phase_list, model):
        
        tmp_dist_d_r1 = kilometers2degrees(tmp_dist_k_r1)     # distance array in degree
        
        tmp_tt_list = []
        for dist_d in tmp_dist_d_r1:
            arrivals = model.get_travel_times(source_depth_in_km=dep_k, distance_in_degree=dist_d, phase_list=phase_list)
            tt = arrivals[0].time
            tmp_tt_list.append(tt)
        tmp_tt_r1 = np.array(tmp_tt_list)
        
        f = interpolate.interp1d(tmp_dist_k_r1, tmp_tt_r1, kind='cubic')
        tt_r1 = f(dist_k_r1)
        
        return tt_r1

    #--- calculate travel time
    time_start = time.time()
    
    if dask_client == None:
        for i, dep_k in enumerate(dep_k_r1):
            tt_r1 = get_tt(tmp_dist_k_r1, dist_k_r1, dep_k, phase_list, model)
            tt_r2[:,i] = tt_r1
            
    else:
        delayed_functions = []

        for dep_k in dep_k_r1:
            delayed_function = dask.delayed(get_tt)(tmp_dist_k_r1, dist_k_r1, dep_k, phase_list, model)
            delayed_functions.append(delayed_function)

        tt_r1_list = dask.compute(*delayed_functions)

        for i in range(dep_k_r1.size):
            tt_r2[:,i] = tt_r1_list[i]

    result_dict = {
        'distances': dist_k_r1,
        'depths': dep_k_r1,
        'traveltimes': tt_r2
    }
    
    time_end = time.time()
    computation_time = time_end - time_start
    computation_time
    
    return dist_k_r1, dep_k_r1, tt_r2, computation_time



def make_source_grid(st, min_lon=None, max_lon=None, min_lat=None, max_lat=None, min_dep=0.0, max_dep=10.0,
                     step_x=0.5, step_y=0.5, step_z=0.5):
    
    stlo_list = [tr.stats.sac.stlo for tr in st]
    stla_list = [tr.stats.sac.stla for tr in st]
    
    if min_lon == None:
        min_lon = min(stlo_list)
    if max_lon == None:
        max_lon = max(stlo_list)
    if min_lat == None:
        min_lat = min(stla_list)
    if max_lat == None:
        max_lat = max(stla_list)
        
    min_x, min_y, ZONE_NUMBER, ZONE_LETTER = utm.from_latlon(min_lat, min_lon)
    max_x, max_y, _, _ = utm.from_latlon(max_lat, max_lon, force_zone_number=ZONE_NUMBER, force_zone_letter=ZONE_LETTER)
    
    min_x /= 1000
    max_x /= 1000
    min_y /= 1000
    max_y /= 1000
    
    ref_x, ref_y = copy.deepcopy(min_x), copy.deepcopy(min_y)
    
    min_x -= ref_x
    max_x -= ref_x
    min_y -= ref_y
    max_y -= ref_y
    
    x_r1 = np.arange(min_x, max_x+step_x, step_x)
    y_r1 = np.arange(min_y, max_y+step_y, step_y)
    z_r1 = np.arange(min_dep, max_dep+step_z, step_z)
    
    _, lon_r1 = utm.to_latlon(   (x_r1 + ref_x)*1000,                np.ones(x_r1.shape)*ref_y*1000,    ZONE_NUMBER,    ZONE_LETTER)
    lat_r1, _ = utm.to_latlon(   np.ones(y_r1.shape)*ref_x*1000,     (y_r1 + ref_y)*1000,               ZONE_NUMBER,    ZONE_LETTER)
    
    return ZONE_NUMBER, ZONE_LETTER, ref_x, ref_y, x_r1, y_r1, lon_r1, lat_r1, z_r1


def get_station_info(st, ref_x, ref_y, ZONE_NUMBER, ZONE_LETTER):

    net_sta_list = [tr.stats.network+'.'+tr.stats.station for tr in st]
    net_sta_list = sorted(list(set(net_sta_list)))

    stlo_list = []
    stla_list = []
    stel_list = []

    for net_sta in net_sta_list:
        net, sta = net_sta.split('.')
        stats = st.select(network=net, station=sta)[0].stats

        stlo_list.append(stats.sac.stlo)
        stla_list.append(stats.sac.stla)
        stel_list.append(stats.sac.stel)

    stlo_r1 = np.array(stlo_list)
    stla_r1 = np.array(stla_list)
    stel_r1 = np.array(stel_list)

    stx_r1, sty_r1, _, _ = utm.from_latlon(stla_r1, stlo_r1, force_zone_number=ZONE_NUMBER, force_zone_letter=ZONE_LETTER)

    stx_r1 /= 1000
    sty_r1 /= 1000

    stx_r1 = stx_r1 - ref_x
    sty_r1 = sty_r1 - ref_y

    return net_sta_list, stlo_r1, stla_r1, stel_r1, stx_r1, sty_r1


def distance(x1,y1,x2,y2):
    '''
    Computes the distance using Pythagoras theorem
    '''
    dist = ((x2-x1)**2 + (y2-y1)**2)**0.5
    return dist


def get_travel_times(mod_dist_r1, mod_dep_r1, mod_tt_r2, x_r1, y_r1, z_r1, stx_r1, sty_r1):
    
    '''
    ttimes: shape=(ny,nx,nz,nsta)
    '''
    
    mod_points = (mod_dist_r1, mod_dep_r1)
    
    nx = x_r1.shape[0]
    ny = y_r1.shape[0]
    nz = z_r1.shape[0]
    nsta = stx_r1.shape[0]
    
    x_r2, y_r2 = np.meshgrid(x_r1, y_r1)
    dists = np.zeros((ny, nx, nsta))
    
    for i in range(nsta):
        dists[:,:,i] = distance(x_r2, y_r2, stx_r1[i], sty_r1[i])
        
    dists_ravel = dists.ravel()
    
    f_1d = np.empty(mod_dep_r1.size, dtype='object')
    
    for i in range(0, mod_dep_r1.size):
        f_1d[i] = interpolate.interp1d(mod_dist_r1, mod_tt_r2[:,i], kind='nearest')
        
    ttimes = np.zeros((ny, nx, nz, nsta))
    
    for j in range(0, nz):
        idx = np.where(z_r1[j]==mod_dep_r1)[0]
        
        if len(idx) == 1:
            i = idx[0]
            ttimes[:,:,j,:] = f_1d[i](dists)
            
        else:
            points = [[d, z_r1[j]] for d in dists_ravel]
            ttimes[:,:,j,:] = interpolate.interpn(mod_points, mod_tt_r2, points).reshape(dists.shape)
        
        # if j == 1:
        #     break # idx
    
    return ttimes


def stream_to_matrix(st, component, net_sta_list):

    npts = min([tr.stats.npts for tr in st])
    nsta = len(net_sta_list)
    data = np.zeros((nsta, npts))

    for i in range(nsta):
        net_sta = net_sta_list[i]
        net, sta = net_sta.split('.')
        st_sel = st.select(network=net, station=sta, component=component)

        if len(st_sel) > 0:
            tr = st_sel[0]
            data[i] = tr.data[:npts]

    return data

def get_data_with_tt(iy, ix, iz, nt0, fs, data_dict, phase_components, ttimes):

    ny, nx, nz, nsta = ttimes.shape
    
    data = []
    
    for i in range(nsta):
        tt = ttimes[iy,ix,iz,i]
        
        i1 = round(tt * fs)
        i2 = i1 + nt0
        
        for component in phase_components:
            data_ = data_dict[component][i,i1:i2]
            data.append(data_)
    
    data = np.array(data)
    
    return data

def calculate_brightness(iy, ix, iz, nt0, fs, data_dict,
                         p_components, s_components, ttimes_p, ttimes_s,
                         nperseg, noverlap=None, pos='start'):

    calculation='mean'

    data = []
    
    if len(p_components) > 0:
        p_data = get_data_with_tt(iy, ix, iz, nt0, fs, data_dict, p_components, ttimes_p)
        data.append(p_data)
    
    if len(s_components) > 0:
        s_data = get_data_with_tt(iy, ix, iz, nt0, fs, data_dict, s_components, ttimes_s)
        data.append(s_data)
    
    data = np.vstack(data)
    data = data.mean(axis=0)
    
    _, b = sliding_window_calculation(data, nperseg, noverlap=noverlap, fs=fs, calculation=calculation, pos=pos)

    return b


def calculate_brightness_serial(idx_array, nt0, fs, data_dict,
                                p_components, s_components, ttimes_p, ttimes_s,
                                nperseg, noverlap=None, pos='start'):

    ny, nx, nz, nsta = ttimes_p.shape

    b_list = []
    
    for idx in idx_array:
        iy, ix, iz = np.unravel_index(idx, (ny,nx,nz))
        
        b = calculate_brightness(iy, ix, iz, nt0, fs, data_dict, p_components, s_components,
                                    ttimes_p, ttimes_s, nperseg, noverlap=noverlap, pos=pos)
        
        b_list.append(b)

    return b_list

def calculate_brightness_parallel(dask_client, nt0, fs, nt0_bp, data_dict,
                                    p_components, s_components, ttimes_p, ttimes_s,
                                    nperseg, noverlap=None, pos=None):

    ny, nx, nz, nsta = ttimes_p.shape
    # nt0 = len(t0_r1)
    # nt0_bp = len(bp_t0_r1)
    
    source_size = nx * ny * nz
    idx_array_main = np.arange(source_size)
    
    n_workers = len(dask_client.scheduler_info()['workers'])
    
    idx_array_list = np.array_split(idx_array_main, n_workers)
    
    delayed_functions = []
    
    for idx_array in idx_array_list:
        
        delayed_function = dask.delayed(calculate_brightness_serial)(idx_array, nt0, fs, data_dict,
                                        p_components, s_components, ttimes_p, ttimes_s,
                                        nperseg, noverlap=noverlap, pos=pos)
    
        delayed_functions.append(delayed_function)
    
    out = dask.compute(*delayed_functions)
    b_r1_list = list(itertools.chain.from_iterable(out))
    B_r4= np.reshape(b_r1_list, (ny, nx, nz, nt0_bp))

    return B_r4

def do_bp(st_dls, dask_client, w,
         min_lon, max_lon, min_lat, max_lat, min_dep, max_dep, step_x, step_y, step_z,
         mod_dist_r1, mod_dep_r1, mod_ttp_r2, mod_tts_r2,
         o=None,
         pos='start',
         p_components=['Z'],
         s_components=['N', 'E'],
         p_cor_dict=None,
         s_cor_dict=None):
    """Perform backprojection to calculate 4-D brightness volume.
    
    Parameters
    ----------
    st_dls: ObsPy.Stream
        Waveform stream.
        For example a differentiated local similarity stream.
    dask_client: dask.Client
        A dask client for parallel processing.
    w: float
        Window size in seconds.
    min_lon: float
        Minimum longitude in degrees.
    max_lon: float
        Maximum longitude in degrees.
    min_lat: float
        Minimum latitude in degrees.
    max_lat: float
        Maximum latitude in degrees.
    min_dep: float
        Minimum depth in kilometers.
    max_dep: float
        Maximum depth in kilometers.
    step_x: float
        Grid size (kilometers) along X-direction.
    step_y: float
        Grid size (kilometers) along Y-direction.
    step_z: float
        Grid size (kilometers) along Z-direction.
    mod_dist_r1: numpy.ndarray
        Model epicentral distance (kilometers) in 1-D.
    mod_dep_r1: numpy.ndarray
        Model depth (kilometers) in 1-D.
    mod_ttp_r2: numpy.ndarray
        Model travel time (seconds) in 2-D for P-wave.
    mod_tts_r2: numpy.ndarray
        Model travel time (seconds) in 2-D for S-wave.
    o: float
        Overlap fraction for successive window.
        Default is None.
    pos: str
        Position of the brightness value in each window.
        Possible values are 'start', 'mid' or 'end'.
        Default value is 'start'.
    p_components: list of str
        A list of components for P-wave brightness.
        Default is ['Z'].
    s_components: list of str
        A list of components for S-wave brightness.
        Default is ['N', 'E'].
    p_cor_dict: dict
        A dictionary of station correction for P-wave.
        Default is None.
    s_cor_dict: dict
        A dictionary of station correction for S-wave.
        Default is None.
        
    Returns
    -------
    brightness4: seisscan.Brightness4
    """
    
    
    #--- get sampling rate and interval
    fs = st_dls[0].stats.sampling_rate
    dt = st_dls[0].stats.delta
    
    #--- define nperseg and noverlap
    nperseg = round(w * fs) + 1
    
    if o == None:
        noverlap = nperseg - 1
    else:
        noverlap = int(w * fs * o) + 1
    
    #--- npts
    # npts = min([tr.stats.npts for tr in st_dls])
    
    
    #--- prepare source grid
    ZONE_NUMBER, ZONE_LETTER, ref_x, ref_y, x_r1, y_r1, lon_r1, lat_r1, z_r1 = make_source_grid(st_dls,
                                                                    min_lon=min_lon, max_lon=max_lon,
                                                                    min_lat=min_lat, max_lat=max_lat,
                                                                    min_dep=min_dep, max_dep=max_dep,
                                                                    step_x=step_x, step_y=step_y, step_z=step_z)
    
    
    #--- get station information
    #--- relative to source coordinates
    net_sta_list, stlo_r1, stla_r1, stel_r1, stx_r1, sty_r1 = get_station_info(st_dls, ref_x, ref_y, ZONE_NUMBER, ZONE_LETTER)
    
    
    #--- make travel time array
    ttimes_p = get_travel_times(mod_dist_r1, mod_dep_r1, mod_ttp_r2, x_r1, y_r1, z_r1, stx_r1, sty_r1)
    ttimes_s = get_travel_times(mod_dist_r1, mod_dep_r1, mod_tts_r2, x_r1, y_r1, z_r1, stx_r1, sty_r1)
    
    
    #--- station correction
    if p_cor_dict == None:
        p_cor_dict = {}
        for net_sta in net_sta_list:
            p_cor_dict[net_sta] = 0.0
    
    if s_cor_dict == None:
        s_cor_dict = {}
        for net_sta in net_sta_list:
            s_cor_dict[net_sta] = 0.0
            
    
    #--- modify travel time array based on station correction
    for i, net_sta in enumerate(net_sta_list):
        ttimes_p[:,:,:,i] += p_cor_dict[net_sta]
        ttimes_s[:,:,:,i] += s_cor_dict[net_sta]
        
    
    #--- starttime and endtime & cut
    starttime = max([tr.stats.starttime for tr in st_dls])
    endtime = min([tr.stats.endtime for tr in st_dls])
    
    st_ = st_dls.copy()
    st_.trim(starttime, endtime)
    
    #--- maximum of modified travel time
    max_tt = max([ttimes_p.max(), ttimes_s.max()])
    
    #--- total duraion of origin times
    total_duration = endtime - starttime - max_tt
    total_duration = total_duration - 5             # reduce by 5 seconds
    
    #--- origin time
    t0_r1 = np.arange(0, total_duration, dt)
    nt0 = t0_r1.shape[0]
    
    #--- components
    components = list(set([tr.stats.channel[-1].upper() for tr in st_dls]))
    
    
    #--- prepare data
    data_dict = {}
    for component in components:
        data_dict[component] = stream_to_matrix(st_, component, net_sta_list)
        
    
    #--- bp origin time
    idx, _ = sliding_window_view(t0_r1, nperseg, noverlap=noverlap, pos=pos)
    bp_t0_r1 = idx * dt
    nt0_bp = len(bp_t0_r1)
    
    #--- calculate brightness parallel
    B_r4 = calculate_brightness_parallel(dask_client, nt0, fs, nt0_bp, data_dict,
                                        p_components, s_components, ttimes_p, ttimes_s,
                                        nperseg, noverlap=noverlap, pos=pos)
    
    #--- prepare brightness object
    brightness4 = Brightness4(B_r4, starttime, bp_t0_r1, lon_r1, lat_r1, z_r1)

    return brightness4


