import pickle
import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
from matplotlib import colormaps as cm
from obspy.geodetics.base import gps2dist_azimuth
from obspy.signal.invsim import cosine_taper

def sliding_coherence(st, station1, station2, component, win_len, ovl_frac=None, window=None, delf=5.0, eps=10):
    '''
    win_len = 1.0
    ovl_frac = None
    ovl_frac = 0.75
    taper_frac = 1.0
    delf = 5.0 # Hz
    eps = 10
    
    window = ('boxcar')
    window = ('hamming')
    window = ('hann')
    window = ('hann')
    window = ('blackman')
    window = ('bartlett')
    window = ('flattop')
    window = ('tukey', M, alpha)
    window = ('kaiser', 14.0) # second param (beta): 0:boxcar, 5:hamming, 6:hann, 8.6:blackman
    window = ('gaussian', 25.0) # second param (std)
    '''
    
    st1 = st.select(station=station1, component=component)
    st2 = st.select(station=station2, component=component)
    
    fs = st1[0].stats.sampling_rate
    dt = st1[0].stats.delta

    starttime = max([tr.stats.starttime for tr in st1 + st2])
    npts = min([tr.stats.npts for tr in st1 + st2])
    x = st1[0].data[:npts]
    y = st2[0].data[:npts]
    times = st1[0].times(type='relative', reftime=starttime)[:npts]

    nperseg = int(win_len * fs)

    if ovl_frac == None:
        stride = 1
    else:
        stride = int(win_len * (1-ovl_frac) * fs)
        
    x_v = np.lib.stride_tricks.sliding_window_view(x,nperseg)[::stride,:]
    y_v = np.lib.stride_tricks.sliding_window_view(y,nperseg)[::stride,:]
    times_v = np.lib.stride_tricks.sliding_window_view(times,nperseg)[::stride,:]

    x_v = signal.detrend(x_v, type='linear', axis=1)
    y_v = signal.detrend(y_v, type='linear', axis=1)
    
    if window == None:
        window = ('boxcar')
        
    taper = signal.get_window(window, nperseg)
    
    for i in range(x_v.shape[0]):
        x_v[i,:] *= taper
        y_v[i,:] *= taper
        
    nfft = nperseg

    X_v = np.fft.rfft(x_v, n=nfft, axis=1)
    Y_v = np.fft.rfft(y_v, n=nfft, axis=1)
    Xcnj_v = np.conjugate(X_v)
    Ycnj_v = np.conjugate(Y_v)
    f = np.fft.rfftfreq(nfft, d=dt)
    df = f[1] - f[0]

    ############
    
    times_out = times_v[:,0]
    # times_out = times_out - times_out[1]
    times_out = times_out + win_len
    # times_out = times_out + win_len / 2
    
    freq_out = np.arange(delf/2+df, f.max()-delf/2+df, df)
    fbins_i = [(int((fc-delf/2)/df), int(fc/df), int((fc+delf/2)/df)) for fc in freq_out]
    lags = np.arange(-eps, eps+1, 1) * dt

    nf = len(fbins_i)
    nt = times_v.shape[0]
    C = np.zeros((nf, nt))

    for f_i in range(nf):

        i1, i0, i2 = fbins_i[f_i]
        # i2 += 1

        f_tile = np.tile(f[i1:i2], (nt,1))

        A = X_v[:,i1:i2]
        Acnj = Xcnj_v[:,i1:i2]
        B = Y_v[:,i1:i2]
        Bcnj = Ycnj_v[:,i1:i2]

        numerator_bkp = Acnj * B
        denominator = np.abs(np.mean(Acnj*A, axis=1) * np.mean(Bcnj*B, axis=1))

        ratio_list = []
        for l_i in range(len(lags)):
            
            lag = lags[l_i]

            numerator = numerator_bkp *  np.exp(2.0 * np.pi * 1j * f_tile * lag)
            numerator = np.abs(numerator.mean(axis=1)) ** 2

            ratio_list.append(numerator / denominator)
            
        ratio_arr = np.array(ratio_list)

        C[f_i,:] = ratio_arr.max(axis=0)
        
    return starttime, times_out, freq_out, C
###########################