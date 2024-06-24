from obspy import Stream
from obspy.core.util import AttribDict
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees, degrees2kilometers

from utils import fix_data
from correlation import sliding_max_lagged_ncc

def do_pcc_(st, component, reference, secondary, w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0, fix=True):
    '''
    Description:
            Perform sliding window max normalized cross-correlation given an obspy stream.
            
        Input:
            st: obspy stream
            component: Component of seismogram
            reference: Reference station code
            secondary: Secondary station code
            w: Window length in seconds. Default is 1.0 s.
            dt: Sampling rate (seconds) of the output 1-D max normalized cross-correlation. Default is 0.1 s.
            max_lag: Maximum allowed lag in seconds. Default is 0.1.
            pos: Position of the max normalized cross-correlation value in each window.
                Possible values are 'start', 'mid' or 'end'.
                Default value is 'end'
            method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                    If method = 0, it returns C
                    If method = 1, it returns |C|
                    If method = 2, it returns C*C
                    If method = 3, it returns C*|C|
                    Default method is 0.
            
        Output:
            st_r: Reference waveform stream
            st_s: Secondary waveform stream
            st_pcc: Stream of sliding window max normalized cross-correlation coefficient.
            st_dpcc: Differentiated st_pcc
            
        Note:
            Both st_pcc[0].stats.pcc and st_pcc[0].stats.dpcc contain the following information
                Reference: Reference station code
                neighbor: Secondary station code
                distance: distance between the Reference and Secondary stations in meters
                az: Reference Secondary azimuth in degrees
                baz: master Secondary backazimuth in degrees, or Secondary Reference azimuth in degrees.
                w: window size in seconds for the sliding window max normalized cross-correlation
    '''

    #--- select stream
    st_sel = Stream()
    st_sel += st.select(station=reference, component=component).copy()
    st_sel += st.select(station=secondary, component=component).copy()
    st = st_sel

    #--- fix the data
    if fix:
        st = fix_data.fix_stream(st, min_npts=100)
        st = fix_data.fix_times(st)
        
        sampling_rate = min([tr.stats.sampling_rate for tr in st])
        st.interpolate(sampling_rate=sampling_rate)
        
        st.detrend().merge(method=1, interpolation_samples=0, fill_value='interpolate')
        
    #--- Get reference and secondary stream
    st_r = st.select(station=reference, component=component).copy()
    st_s = st.select(station=secondary, component=component).copy()
        
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


def do_pcc(st, component, reference, secondaries=[], w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0, fix=True):
    '''
    Description:
            Perform sliding window max normalized cross-correlation given an obspy stream.
            
        Input:
            st: obspy stream
            component: Component of seismogram
            reference: Reference station code
            secondaries: List of Secondary station code
            w: Window length in seconds. Default is 1.0 s.
            dt: Sampling rate (seconds) of the output 1-D max normalized cross-correlation. Default is 0.1 s.
            max_lag: Maximum allowed lag in seconds. Default is 0.1.
            pos: Position of the max normalized cross-correlation value in each window.
                Possible values are 'start', 'mid' or 'end'.
                Default value is 'end'
            method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                    If method = 0, it returns C
                    If method = 1, it returns |C|
                    If method = 2, it returns C*C
                    If method = 3, it returns C*|C|
                    Default method is 0.
            
        Output:
            st_r: Reference waveform stream
            st_s: Secondary waveform stream
            st_pcc: Stream of sliding window max normalized cross-correlation coefficient.
            st_dpcc: Differentiated st_pcc
            
        Note:
            Both st_pcc[0].stats.pcc and st_pcc[0].stats.dpcc contain the following information
                Reference: Reference station code
                neighbor: Secondary station code
                distance: distance between the Reference and Secondary stations in meters
                az: Reference Secondary azimuth in degrees
                baz: master Secondary backazimuth in degrees, or Secondary Reference azimuth in degrees.
                w: window size in seconds for the sliding window max normalized cross-correlation
    '''
    #--- select stream
    st_sel = Stream()
    st_sel += st.select(station=reference, component=component).copy()
    for secondary in secondaries:
        st_sel += st.select(station=secondary, component=component).copy()
    st = st_sel
    
    #--- fix the data
    if fix:
        st = fix_data.fix_stream(st, min_npts=100)
        st = fix_data.fix_times(st)
        
        sampling_rate = min([tr.stats.sampling_rate for tr in st])
        st.interpolate(sampling_rate=sampling_rate)
        
        st.detrend().merge(method=1, interpolation_samples=0, fill_value='interpolate')
        
    #--- Compute pcc and dpcc between master and each neighbor
    st_s = Stream()
    st_pcc = Stream()
    st_dpcc = Stream()
    
    for secondary in secondaries:
        st_r, st_s_, st_pcc_, st_dpcc_ = do_pcc_(st, component, reference, secondary=secondary, w=w, dt=dt, max_lag=max_lag, pos=pos, method=method, fix=False)
        st_s += st_s_
        st_pcc += st_pcc_
        st_dpcc += st_dpcc_
        
    return st_r, st_s, st_pcc, st_dpcc