import numpy as np
from obspy import Stream

from utils import fix_data
from peak_cross_correlation import do_pcc

def do_ls_(st, component, reference, secondaries=[], w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0, fix=True):
    '''
    Description:
        Perform local similarity given an obspy stream.
            
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
        st_ls: Local similarity stream
        st_dls: Differentiated local similarity stream
        
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
        
    #--- compute pcc and dpcc
    st_r, st_s, st_pcc, st_dpcc = do_pcc(st, component, reference, secondaries=secondaries, w=w, dt=dt, max_lag=max_lag, pos=pos, method=method, fix=False)

    #--- Compute ls
    n = min([tr.stats.npts for tr in st_pcc])
    ls = np.array([tr.data[:n] for tr in st_pcc]).mean(axis=0)

    #--- prepare ls stream
    st_ls = st_pcc.select(station=secondaries[0]).copy()

    st_ls[0].data = ls
    st_ls[0].stats.station = reference
    st_ls[0].stats.sac = st_r[0].stats.sac
    st_ls[0].stats.ls = st_ls[0].stats.pcc
    st_ls[0].stats.ls.secondaries = secondaries
    del st_ls[0].stats.pcc
    del st_ls[0].stats.ls.secondary
    del st_ls[0].stats.ls.distance
    del st_ls[0].stats.ls.az
    del st_ls[0].stats.ls.baz

    #--- compute dls stream
    st_dls = st_ls.copy()
    st_dls.differentiate(method='gradient');

    st_dls[0].stats.dls = st_dls[0].stats.ls
    del st_dls[0].stats.ls
    
    return st_r, st_s, st_pcc, st_dpcc, st_ls, st_dls



def do_ls(st, component, rs_list=[], w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0, fix=True):
    '''
    Description:
        Perform local similarity given an obspy stream.
            
    Input:
        st: obspy stream
        component: Component of seismogram
        rs_list: list of reference_secondaries dictionary
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
        st_ls: Local similarity stream
        st_dls: Differentiated local similarity stream
        
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

    for reference_secondaries in rs_list:
        reference = reference_secondaries['reference']
        secondaries = reference_secondaries['secondaries']
        
        st_sel += st.select(station=reference).copy()
        for secondary in secondaries:
            st_sel += st.select(station=secondary).copy()

    st = st_sel

    #--- fix the data
    if fix:
        st = fix_data.fix_stream(st, min_npts=100)
        st = fix_data.fix_times(st)
        
        sampling_rate = min([tr.stats.sampling_rate for tr in st])
        st.interpolate(sampling_rate=sampling_rate)
        
        st.detrend().merge(method=1, interpolation_samples=0, fill_value='interpolate')
        
    #--- compute pcc, dpcc,ls and dls
    st_r = Stream()
    st_s = Stream()
    st_pcc = Stream()
    st_dpcc = Stream()
    st_ls = Stream()
    st_dls = Stream()

    for reference_secondaries in rs_list:
        reference = reference_secondaries['reference']
        secondaries = reference_secondaries['secondaries']
        
        st_r_, st_s_, st_pcc_, st_dpcc_, st_ls_, st_dls_ = do_ls_(st, component, reference, secondaries=secondaries, w=w, dt=dt, max_lag=max_lag, pos=pos, method=method, fix=False)
        
        st_r += st_r_
        st_s += st_s_
        st_pcc += st_pcc_
        st_dpcc += st_dpcc_
        st_ls += st_ls_
        st_dls += st_dls_
        
        
    return st_r, st_s, st_pcc, st_dpcc, st_ls, st_dls


