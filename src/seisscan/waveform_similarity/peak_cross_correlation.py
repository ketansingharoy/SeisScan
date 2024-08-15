from obspy import Stream
from obspy.core.util import AttribDict
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees, degrees2kilometers

import dask
from dask.distributed import Client as dask_Client
from dask.distributed import LocalCluster as dask_LocalCluster

# from .correlation import sliding_max_lagged_ncc
from seisscan.waveform_similarity.correlation import sliding_max_lagged_ncc

def do_pcc_(st, channel, reference, secondary, w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0):
    '''Create peak cross-correlation function between a reference and a secondary station.
    
    Parameters
    ----------
    st: ObsPy.Stream
        Waveform stream.
    channel: str
        Channel code.
    reference: str
        A reference station code
    secondary: str
        A secondary station code.
    w: float
        Window length in seconds.
        Default is 1.0 second.
    dt: float
        Sampling interval (seconds) of the returned characteristic function.
        Default is 0.1 second.
    max_lag: float
        Maximum allowed lag in seconds.
        Default is 0.1 second.
    pos: str
        Position of the max normalized cross-correlation value in each window.
        Possible values are 'start', 'mid' or 'end'.
        Default value is 'end'.
    method: int
        One of the integer (0, 1, 2, 3). It determines type of normalized cross-correlation.
        If method = 0, it returns C.
        If method = 1, it returns |C|.
        If method = 2, it returns C*C.
        If method = 3, it returns C*|C|.
        Default method is 0.
        
    Returns
    -------
    st_r: ObsPy.Stream
        Waveform stream for reference stations.
    st_s: ObsPy.Stream
        Waveform stream for secondary stations.
    st_pcc: ObsPy.Stream
        Peak cross-correlation stream.
    st_dpcc: ObsPy.Stream
        Differentiated peak cross-correlation stream.
    '''

    #--- select stream
    st_sel = Stream()
    st_sel += st.select(station=reference, channel=channel).copy()
    st_sel += st.select(station=secondary, channel=channel).copy()
    st = st_sel
        
    #--- Get reference and secondary stream
    st_r = st.select(station=reference, channel=channel).copy()
    st_s = st.select(station=secondary, channel=channel).copy()
        
    #--- get important header and data
    delta = st_r[0].stats.delta
    x = st_r[0].data
    y = st_s[0].data

    #--- fix data based on minumum length
    n = min([len(x), len(y)])
    x = x[:n]
    y = y[:n]

    # #--- adjust dt
    dt = int(dt / delta) * delta

    #--- prepare parameters for sliding_max_lagged_ncc
    nperseg = int(w / delta) + 1
    stride = int(dt / delta)
    max_ilag = int(max_lag / delta)

    #--- signal shift in seconds due to cross-correlation
    if pos == 'start':
        shift = 0
    if pos == 'mid':
        shift = w / 2
    if pos == 'end':
        shift = w
        
    #--- Calculate master-neighbor distance, azimuth and backazimuth
    stlo_r = st_r[0].stats.sac.stlo
    stla_r = st_r[0].stats.sac.stla
    stlo_s = st_s[0].stats.sac.stlo
    stla_s = st_s[0].stats.sac.stla

    dist_meter, az, baz = gps2dist_azimuth(stla_r, stlo_r, stla_s, stlo_s)
    dist_km = dist_meter / 1000
        
    #--- perform sliding maximum normalized cross-correlation coefficient
    pcc = sliding_max_lagged_ncc(x, y, nperseg=nperseg, stride=stride, max_ilag=max_ilag, method=method)

    #--- prepare pcc stream
    st_pcc = st_s.copy()

    try:
        del st_pcc[0].stats.response
    except:
        pass

    st_pcc[0].data = pcc
    st_pcc[0].stats.delta = dt
    st_pcc[0].stats.starttime += shift

    st_pcc[0].stats.pcc = AttribDict()
    st_pcc[0].stats.pcc.reference = reference
    st_pcc[0].stats.pcc.secondary = secondary
    st_pcc[0].stats.pcc.distance = dist_meter
    st_pcc[0].stats.pcc.az = az
    st_pcc[0].stats.pcc.baz = baz
    st_pcc[0].stats.pcc.w = w

    #--- compute dpcc stream
    st_dpcc = st_pcc.copy()
    st_dpcc.differentiate(method='gradient')
    
    st_dpcc[0].stats.dpcc = st_dpcc[0].stats.pcc
    del st_dpcc[0].stats.pcc

    return st_r, st_s, st_pcc, st_dpcc


def do_pcc(st, channel, reference, secondaries=[], w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0, dask_client=None):
    '''Create peak cross-correlation function between a reference and a list of secondary stations.
    
    Parameters
    ----------
    st: ObsPy.Stream
        Waveform stream.
    channel: str
        Channel code.
    reference: str
        Reference station code
    secondaries: list
        A list of secondary station codes.
    w: float
        Window length in seconds.
        Default is 1.0 second.
    dt: float
        Sampling interval (seconds) of the returned characteristic function.
        Default is 0.1 second.
    max_lag: float
        Maximum allowed lag in seconds.
        Default is 0.1 second.
    pos: str
        Position of the max normalized cross-correlation value in each window.
        Possible values are 'start', 'mid' or 'end'.
        Default value is 'end'.
    method: int
        Default method is 0.
        
    Returns
    -------
    st_r: ObsPy.Stream
        Waveform stream for reference stations.
    st_s: ObsPy.Stream
        Waveform stream for secondary stations.
    st_pcc: ObsPy.Stream
        Peak cross-correlation stream.
    st_dpcc: ObsPy.Stream
        Differentiated peak cross-correlation stream.
        
    Example
    -------
    >>> import seisscan as ss
    >>>
    >>> event_dict, st, inventory, subnetworks, model_name = ss.read_example()
    >>>
    >>> subnetwork = subnetworks[0]
    >>> reference = subnetwork["reference"]
    >>> secondaries = subnetwork["secondaries"]
    >>>
    >>> st_r, st_s, st_pcc, st_dpcc = ss.do_pcc(st, "DPZ", reference, subnetworks=subnetworks)
    '''
    #--- select stream
    st_sel = Stream()
    st_sel += st.select(station=reference, channel=channel).copy()
    for secondary in secondaries:
        st_sel += st.select(station=secondary, channel=channel).copy()
    st = st_sel
        
    #--- Compute pcc and dpcc between master and each neighbor
    st_s = Stream()
    st_pcc = Stream()
    st_dpcc = Stream()
    
    #--- serial processing
    if dask_client == None:
        
        for secondary in secondaries:
            st_r, st_s_, st_pcc_, st_dpcc_ = do_pcc_(st, channel, reference, secondary=secondary, w=w, dt=dt, max_lag=max_lag, pos=pos, method=method)
            st_s += st_s_
            st_pcc += st_pcc_
            st_dpcc += st_dpcc_
            
    #--- parallel processing   
    else:
        delayed_functions = []
        
        for secondary in secondaries:
            delayed_function = dask.delayed(do_pcc_)(st, channel, reference, secondary=secondary,
                                                    w=w, dt=dt, max_lag=max_lag, pos=pos, method=method)
            delayed_functions.append(delayed_function)
            
        out = dask.compute(*delayed_functions)
        
        for out_ in out:
            st_r, st_s_, st_pcc_, st_dpcc_ = out_
            st_s += st_s_
            st_pcc += st_pcc_
            st_dpcc += st_dpcc_
        
        
    return st_r, st_s, st_pcc, st_dpcc