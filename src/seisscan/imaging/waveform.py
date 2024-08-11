import numpy as np
from obspy.geodetics.base import gps2dist_azimuth, kilometers2degrees, degrees2kilometers
from obspy.taup import taup_create, TauPyModel
from matplotlib import pyplot as plt

def compute_dist_tt(dmin, dmax, evdp, phase, model_name='iasp91', dist_unit='km'):
    '''
    Description:
        Given minimum and maximum distance, compute travel times for a gven phase.
        
    Input:
        dmin: minumum distance
        dmax: maximum distance
        evdp: Event depth (km).
        phase: Seismic phase. 'p' or 's'
        model_name: Earth model. Default is 'iasp91'
        dist_unit: Unit of dmin and dmax. 'meter', 'km' or 'deg'. Default is 'km'
        
    Output:
        dist_array: array of distances
        tt_array: array of travel times
    '''

    dist_array = np.linspace(dmin, dmax, 100)

    if dist_unit == 'meter':
        dist_d_array = kilometers2degrees(dist_array / 1000)
    if dist_unit == 'km':
        dist_d_array = kilometers2degrees(dist_array)
    if dist_unit == 'deg':
        dist_d_array = dist_array
        
    #--- Compute travel time
    taup_model = TauPyModel(model=model_name)
    tt_array = []
    
    for dist_d in dist_d_array:
        arrivals = taup_model.get_travel_times(source_depth_in_km=evdp, distance_in_degree=dist_d, phase_list=[phase.lower(), phase.upper()])
        
        if len(arrivals) > 0:
            tt = arrivals[0].time
        else:
            tt = np.nan
            
        tt_array.append(tt)
            
    tt_array = np.array(tt_array)
    
    return dist_array, tt_array


def prs(st, evt0, evlo, evla, evdp, scale=0.5, ax=None, width=None, height=None,
        xmin=None, xmax=None, phases=['p', 's'], model_name='iasp91',
        w_color='k', w_lw=0.75, phase_colors=['b', 'g'], phase_lws=[0.75, 0.75],
        xlabel_fs=12, ylabel_fs=12, title_fs=12, xticklabels_fs=10, yticklabels_fs=10,
        legend=True, handle=False):
    
    """
    Plot a record section of waveform stream.
    
    Parameters
    ---------
    st: ObsPy.Stream
        Waveform stream.
    evt0: ObsPy.UTCDateTime
        Event origin time.
    evlo: float
        Event longitude.
    evla: float
        Event latitude.
    evdp: float
        Event depth (km).
    scale: float
        A scale factor for waveform amplitude.
    model_name: str
        An earth model name.
        Default is "iasp91"
    phases: list of str
        A list of seismic phases.
        Default is ["p", "s"]
    phase_colors: list of str
        A list of colors.
        Default is ["b", "g"]
    ax: matplotlib.axes
        Matplotlib axes object.
    handle: bool
        If True, returns a matplotlib.Figure object.
        
    Returns
    -------
    fig: matplotlib.Figure
        A matplotlib.Figure if parameter handle=True.
    """

    #--- Create figure object
    if ax == None:
        fig = plt.figure()
        ax = fig.add_subplot()
    else:
        fig = ax.figure
        
    if width != None:
        fig.set_figwidth(width)
        
    if height != None:
        fig.set_figheight(height)
        
    # ax = fig.add_subplot()

    #--- Prepare lists of station code, longitude and latitude
    sta_list = [tr.stats.station for tr in st]
    stlo_list = [tr.stats.sac.stlo for tr in st]
    stla_list = [tr.stats.sac.stla for tr in st]

    #--- Calculate epicentral distances
    dist_km_list = []
    for stlo, stla in zip(stlo_list, stla_list):
        dist_meter, _, _ = gps2dist_azimuth(evla, evlo, stla, stlo)
        dist_km = dist_meter / 1000
        dist_km_list.append(dist_km)
        
    #--- Plot waveforms
    for i_ in range(len(sta_list)):
        sta = sta_list[i_]
        dist_km = dist_km_list[i_]
        
        if i_ == 0:
            label = 'Waveform'
        else:
            label = None
        
        for tr in st.select(station=sta):
            times = tr.times(type='relative', reftime=evt0)
            data = tr.data
            data_modified = data / np.abs(data).max() * scale + dist_km
            ax.plot(times, data_modified, color=w_color, lw=w_lw, label=label)
            
    #--- Plot travel times
    dmin, dmax = ax.get_ylim()

    for i_ in range(len(phases)):
        phase = phases[i_]
        phase_color = phase_colors[i_]
        phase_lw = phase_lws[i_]

        dist_array, tt_array = compute_dist_tt(dmin, dmax, evdp, phase, model_name=model_name, dist_unit='km')
        ax.plot(tt_array, dist_array, color=phase_color, lw=phase_lw, label=phase.upper())
        
    #--- legend
    if legend:
        ax.legend(loc=1, framealpha=1.0)

    #--- Format xlabel
    ax.set_xlabel('Time (s)', fontsize=xlabel_fs)

    #--- Format ylabel
    ax.set_ylabel('Distance (km)', fontsize=ylabel_fs)

    #--- Format xticklabels
    ax.tick_params(axis='x', labelsize=xticklabels_fs)

    #--- Format yticklabels
    ax.tick_params(axis='y', labelsize=yticklabels_fs)

    #--- Format title
    evt0_str = evt0.strftime('%Y-%m-%d %H:%M:%S')
    title = f'Origin time: {evt0_str}'
    ax.set_title(title, fontsize=title_fs)

    #--- Auto scale yaxis
    ax.autoscale(enable=True, axis='y', tight=True)

    #--- limit xaxis
    if xmin is None:
        xmin, _ = ax.get_xlim()
        
    if xmax is None:
        _, xmax = ax.get_xlim()
    ax.set_xlim(xmin, xmax)
    
    if handle:
        return fig