import numpy as np
from obspy import Stream

def _split_trace_on_nan(trace):
    """
    Splits an ObsPy Trace object into multiple Trace objects at NaN values.
    
    Parameters:
    trace (obspy.core.trace.Trace): The Trace object to split.
    
    Returns:
    obspy.core.stream.Stream: A Stream object containing the split Trace objects.
    """
    data = trace.data
    nans = np.isnan(data)
    
    # Find indices where the data transitions from NaN to non-NaN and vice versa
    segments = np.where(np.diff(nans.astype(int)))[0] + 1
    
    # If the entire trace is NaN, return an empty stream
    if len(segments) == 0 and nans.all():
        return Stream()
    
    # If there are no NaNs, return the trace as-is in a Stream
    if len(segments) == 0:
        return Stream(trace)
    
    # Split the trace at these segment indices
    traces = []
    start_idx = 0
    for idx in segments:
        if not nans[start_idx]:  # Only consider non-NaN segments
            new_trace = trace.copy()
            new_trace.data = data[start_idx:idx]
            new_trace.stats.starttime += start_idx * trace.stats.delta
            traces.append(new_trace)
        start_idx = idx
    
    # Handle the final segment
    if not nans[start_idx]:
        new_trace = trace.copy()
        new_trace.data = data[start_idx:]
        new_trace.stats.starttime += start_idx * trace.stats.delta
        traces.append(new_trace)
    
    return Stream(traces)

def _split_stream_on_nan(st):
    '''
    Description: Splits ObsPy Trace object from stream object into multiple Trace objects at NaN values.
    '''
    
    st_tmp = Stream()
    for tr in st:
        st_tmp += _split_trace_on_nan(tr)
    st = st_tmp
    
    return st

def _remove_small_traces_from_stream(st, min_npts=100):
    '''
    Description:
        Fix stream. It removes the traces with very small npts.
        It trims the stream so that starttime or endtime is same for all the traces.
    '''
    
    tr_list = [tr for tr in st if tr.stats.npts >= min_npts]
    st = Stream(tr_list)
    
    return st

def fix_stream(st, min_npts=100):
    '''
    Description:
        Fix obspy stream.
        Splits ObsPy Trace object from stream object into multiple Trace objects at NaN values.
        It then removes the traces with very small npts.
    '''
    
    st = _split_stream_on_nan(st)
    st = _remove_small_traces_from_stream(st, min_npts=min_npts)
    
    return st

def fix_times(st):
    
    #--- fix for start and end time
    t1 = min([tr.stats.starttime for tr in st])
    t2 = max([tr.stats.endtime for tr in st])

    st.detrend(type='linear')
    st.merge(method=1, fill_value=0)
    st.trim(t1, t2, pad=True, fill_value=0);
    
    return st