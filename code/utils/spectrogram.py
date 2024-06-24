import numpy as np
from scipy import signal
from matplotlib import pyplot as plt
from matplotlib import dates as mdates

def spectrogram(st, station=None, channel=None, component=None, network=None, location=None, seed_id=None,
                w=1.0, dt=0.1, nfft=None, taper='hann', vmin=None, vmax=None, cmap='viridis', ttype='absolute', reftime=None, dB=True, handle=False):
    '''
    network = None
    station = 'SM056'
    location = None
    channel = 'BDF'
    component = None
    seed_id = None
    w = 1.0 # window length in seconds
    dt = 0.1 # time interval in seconds
    taper=('tukey', 0.05)
    taper='hann'
    vmin_p = 1.0
    vmax_p = 99.0
    ttype = 'relative' # 'relative' or 'absolute'
    reftime = None
    '''

    #--- select stream
    st_ = st.select(network=network, station=station, location=location, channel=channel, component=component, id=seed_id).copy()

    #--- get important headers
    fs = st_[0].stats.sampling_rate
    delta = 1 / fs

    #--- fix dt
    dt = int(dt / delta) * delta

    #--- prepare parameters for sliding calculation
    nperseg = int(w / delta) + 1
    noverlap = int((w - dt) / delta) + 1
    
    if nfft == None:
        nfft = nperseg

    #--- loop over traces to calculate segmented spectrograms
    starttime = min([tr_.stats.starttime for tr_ in st_])

    #--- fix reftime
    if reftime == None:
        reftime = starttime

    starttime_list = []
    freq_list, time_list, Sxx_list = [], [], []

    for tr_ in st_:
        starttime_ = tr_.stats.starttime
        data_ = tr_.data
        
        freq_, time_, Sxx_ = signal.spectrogram(data_, fs=fs, window=taper, nperseg=nperseg, noverlap=noverlap, detrend='linear', scaling='density')
        time_ = starttime_ - starttime + time_
        
        if dB == True:
            Sxx_ = 10 * np.log10(Sxx_)
        
        freq_list.append(freq_)
        time_list.append(time_)
        Sxx_list.append(Sxx_)
        
    #--- vmin and vmax for colormap
    Sxx_tmp = np.hstack(Sxx_list)
    # vmin = np.percentile(Sxx_tmp, vmin_p)
    # vmax = np.percentile(Sxx_tmp, vmax_p)
    if vmin == None:
        vmin = Sxx_tmp.min()
        
    if vmax == None:
        vmax = Sxx_tmp.max()

    #--- modify time
    if ttype == 'relative':
        for i_ in range(len(time_list)):
            time_list[i_] = starttime - reftime + time_list[i_]
            
    if ttype == 'absolute':
        for i_ in range(len(time_list)):
            time_list[i_] = mdates.date2num(starttime) + time_list[i_] / 86400
            

    ############# PLOT ############
    #--- figure and axes
    fig = plt.figure(figsize=(12,5))
    fig.subplots_adjust(hspace=0)
    gs = fig.add_gridspec(6, 20)

    ax1 = fig.add_subplot(gs[0,:19])
    ax2 = fig.add_subplot(gs[1:,:19], sharex=ax1)
    cax = fig.add_subplot(gs[1:,19:])

    #--- seismogram
    for tr_ in st_:
        if ttype == 'absolute':
            t_ = tr_.times(type='matplotlib', reftime=reftime)
        if ttype == 'relative':
            t_ = tr_.times(type='relative', reftime=reftime)
            
        d_ = tr_.data
        
        ax1.plot(t_, d_, color='k', lw=0.75)

    #--- spectrogram
    for freq_, time_, Sxx_ in zip(freq_list, time_list, Sxx_list):
        s = ax2.pcolormesh(time_, freq_, Sxx_, vmin=vmin, vmax=vmax, cmap=cmap)
        # break
        
    fig.colorbar(s, ax=ax2, cax=cax, label='PSD')

    #--- title
    if ttype == 'relative':
        reftime_str = reftime.strftime('%Y-%m-%d %H:%M:%S')
        title = f'Ref. time: {reftime_str}  |  Station: {st_[0].id}'
        ax1.set_title(title)
        
    if ttype == 'absolute':
        title = f'Station: {st_[0].id}'
        ax1.set_title(title)

    #--- xlabel
    if ttype == 'relative':
        ax2.set_xlabel('Time (s)')
        
    #--- ylabel
    ax1.set_ylabel('Seis.')
    ax2.set_ylabel('Frequency (Hz)')
    
    #--- hide xticklabels from ax1
    [label.set_visible(False) for label in ax1.get_xticklabels()]

    #--- yticks
    ax1.set_yticks([])

    #--- auto xlim
    ax1.autoscale(enable=True, axis='x', tight=True)

    #--- format absolute time
    if ttype == 'absolute':
        locator = mdates.AutoDateLocator()
        formatter = mdates.ConciseDateFormatter(locator)
        ax2.xaxis.set_major_locator(locator)
        ax2.xaxis.set_major_formatter(formatter)
        
        
    #---- return handle
    if handle:
        return fig