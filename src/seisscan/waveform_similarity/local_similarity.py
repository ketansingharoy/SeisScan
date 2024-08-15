import numpy as np
from obspy import Stream

import dask
from dask.distributed import Client as dask_Client
from dask.distributed import LocalCluster as dask_LocalCluster

# from .peak_cross_correlation import do_pcc_
from seisscan.waveform_similarity.peak_cross_correlation import do_pcc_


def do_ls(st, channel, subnetworks=[], w=1.0, dt=0.1, max_lag=0.1, pos='end', method=0, dask_client=None):
    '''Create characteristic function (Local Similarity).
    
    Parameters
    ----------
    st: ObsPy.Stream
        Waveform stream.
    channel: str
        Channel code.
    subnetworks: list
        station clusters
    w: float
        A win length in seconds.
        Default is 1.0 seconds.
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
        One of the integer (0, 1, 2, 3). It determines type of normalized cross-correlation (cor).
        If method = 0, it returns cor.
        If method = 1, it returns |cor|.
        If method = 2, it returns cor*cor.
        If method = 3, it returns cor*|cor|.
        Default method is 0.
    dask_client: dask.Client
        A dask client for parallel processing.
        Default is None.
        
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
    st_ls: ObsPy.Stream
        Local similarity stream.
    st_dls ObsPy.Stream
        Differentiated local similarity stream.
        
    Example
    -------
    >>> import seisscan as ss
    >>>
    >>> event_dict, st, inventory, subnetworks, model_name = ss.read_example()
    >>>
    >>> st_r, st_s, st_pcc, st_dpcc, st_ls, st_dls = do_ls(st, "DPZ", subnetworks=subnetworks)
    '''

    #--- select stream
    st_sel = Stream()

    for subnetwork in subnetworks:
        reference = subnetwork['reference']
        secondaries = subnetwork['secondaries']
        
        st_sel += st.select(station=reference, channel=channel).copy()
        for secondary in secondaries:
            st_sel += st.select(station=secondary, channel=channel).copy()

    st = st_sel
    
    #--- compute pcc, dpcc,ls and dls
    st_r = Stream()
    st_s = Stream()
    st_pcc = Stream()
    st_dpcc = Stream()
    st_ls = Stream()

    #--- serial processing
    if dask_client == None:

        for subnetwork in subnetworks:
            reference = subnetwork['reference']
            secondaries = subnetwork['secondaries']
            
            for secondary in secondaries:
                _, st_s_, st_pcc_, st_dpcc_ = do_pcc_(st, channel, reference, secondary=secondary,
                                                    w=w, dt=dt, max_lag=max_lag, pos=pos, method=method)
                st_s += st_s_
                st_pcc += st_pcc_
                st_dpcc += st_dpcc_
                
    #--- parallel processing 
    else:
        delayed_functions = []

        for subnetwork in subnetworks:
            reference = subnetwork['reference']
            secondaries = subnetwork['secondaries']

            for secondary in secondaries:
                delayed_function = dask.delayed(do_pcc_)(st, channel, reference, secondary=secondary,
                                                        w=w, dt=dt, max_lag=max_lag, pos=pos, method=method)
                delayed_functions.append(delayed_function)
                
        out = dask.compute(*delayed_functions)

        for out_ in out:
            _, st_s_, st_pcc_, st_dpcc_ = out_
            st_s += st_s_
            st_pcc += st_pcc_
            st_dpcc += st_dpcc_
            
    #--- for reference stream
    for subnetwork in subnetworks:
        reference = subnetwork['reference']
        st_r += st.select(station=reference, channel=channel).copy()
        
    #--- Compute ls
    npts_ = min([tr.stats.npts for tr in st_pcc])

    for subnetwork in subnetworks:
        reference = subnetwork['reference']
        secondaries = subnetwork['secondaries']
        
        ls_ = np.array([st_pcc.select(station=secondary, channel=channel)[0].data[:npts_] for secondary in secondaries]).mean(axis=0)

        #--- prepare ls stream
        st_ls_ = st_pcc.select(station=secondaries[0], channel=channel).copy()

        st_ls_[0].data = ls_
        st_ls_[0].stats.station = reference
        st_ls_[0].stats.sac = st_r.select(station=reference, channel=channel)[0].stats.sac
        st_ls_[0].stats.ls = st_ls_[0].stats.pcc
        st_ls_[0].stats.ls.secondaries = secondaries
        del st_ls_[0].stats.pcc
        del st_ls_[0].stats.ls.secondary
        del st_ls_[0].stats.ls.distance
        del st_ls_[0].stats.ls.az
        del st_ls_[0].stats.ls.baz

        st_ls += st_ls_
        
    #--- compute dls stream
    st_dls = st_ls.copy()
    st_dls.differentiate(method='gradient');

    for tr in st_dls:
        tr.stats.dls = tr.stats.ls
        del tr.stats.ls
        
        
    return st_r, st_s, st_pcc, st_dpcc, st_ls, st_dls


