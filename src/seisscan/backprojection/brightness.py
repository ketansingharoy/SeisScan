import copy
import pickle
import numpy as np
import scipy
from scipy import signal
import pandas as pd
from obspy import UTCDateTime
from matplotlib import dates as mdates
from matplotlib import pyplot as plt


def _load_stack(filename):
    stack = pickle.load(open(filename, 'rb'))
    return stack

def load_stack(files, duration):

    if isinstance(files, str):
        stack = _load_stack(filename=files)

    if isinstance(files, list):

        stack_list = [_load_stack(filename=file) for file in files]
        
        starttime = min([stack_.starttime for stack_ in stack_list])
        
        bp_t0_list = []
        evlo_list = []
        evla_list = []
        evdp_list = []
        b_list = []
        
        for stack_ in stack_list:
            
            idx = np.where(stack_.bp_t0_r1 < duration)[0]
            times = stack_.get_times(reftime=starttime)
            
            bp_t0_ = times[idx]
            evlo_ = stack_.evlo_r1[idx]
            evla_ = stack_.evla_r1[idx]
            evdp_ = stack_.evdp_r1[idx]
            b_ = stack_.b_r1[idx]
        
            bp_t0_list.append(bp_t0_)
            evlo_list.append(evlo_)
            evla_list.append(evla_)
            evdp_list.append(evdp_)
            b_list.append(b_)
        
        bp_t0_r1 = np.hstack(bp_t0_list)
        evlo_r1 = np.hstack(evlo_list)
        evla_r1 = np.hstack(evla_list)
        evdp_r1 = np.hstack(evdp_list)
        b_r1 = np.hstack(b_list)
        
        stack = Stack(b_r1, starttime, bp_t0_r1, evlo_r1, evla_r1, evdp_r1)

    return stack

def get_catalog(stack, bth, spacing=5.0):

    dt = stack.bp_t0_r1[2] - stack.bp_t0_r1[1]
    spacing_npts = round(spacing / dt)
    
    idx, _ = signal.find_peaks(stack.b_r1, height=bth, distance=spacing_npts)
    
    starttime = stack.starttime
    
    bp_t0_r1 = stack.bp_t0_r1[idx]
    evlo_r1 = stack.evlo_r1[idx]
    evla_r1 = stack.evla_r1[idx]
    evdp_r1 = stack.evdp_r1[idx]
    b_r1 = stack.b_r1[idx]
    
    evt0_list = [(starttime + bp_t0).datetime for bp_t0 in bp_t0_r1]
    
    #---------
    
    lsbp_cat_dict = {
        'origin':evt0_list,
        'longitude':evlo_r1,
        'latitude':evla_r1,
        'depth':evdp_r1,
        'brightness':b_r1
    }
    
    df_lsbp = pd.DataFrame(lsbp_cat_dict)
    
    return df_lsbp


class Stack():
    """A class for holding 1-D backprojected stack.
    
    This class provides methods to get backprojected solution, brightness stack, plot brightness slice etc.
    
    Attributes
    ----------
    b_r1: numpy.ndarray of floats
        1-D numpy array of brightness values.
    starttime: ObsPy.UTCDateTime
        Starttime of the stack
    bp_t0_r1: numpy.ndarray of floats
        1-D numpy array origin times (seconds).
    evlo_r1: numpy.ndarray of floats
        1-D numpy array event longitudes (degree)
    evla_r1: numpy.ndarray of floats
        1-D numpy array event latitudes (degree)
    evdp_r1: numpy.ndarray of floats
        1-D numpy array event depths (km)
        
    Methods
    -------
    copy:
        Returns a deepcopy of the 1-D stack.
    trim:
        Trims (cuts) the 1-D stack given a starttime and endtime
    get_solution:
        Returns event longitude, latitude, depth and brightness given an origin time.
    get_times:
        Returns a 1-D numpy array of times.
    plot:
        Plots backprojected stack on a matploib axes.
    write:
        Writes the stack in a pickle file.
    """

    def __init__(self, b_r1, starttime, bp_t0_r1, evlo_r1, evla_r1, evdp_r1):
        '''Initialize the Stack class.
        
        Parameters
        ----------
        b_r1: numpy.ndarray of floats
             1-D numpy array of brightness values.
        starttime: ObsPy.UTCDateTime
            Starttime of the stack
        bp_t0_r1: numpy.ndarray of floats
            1-D numpy array origin times (seconds).
        evlo_r1: numpy.ndarray of floats
            1-D numpy array event longitudes (degree)
        evla_r1: numpy.ndarray of floats
            1-D numpy array event latitudes (degree)
        evdp_r1: numpy.ndarray of floats
            1-D numpy array event depths (km)
        '''

        self.b_r1 = b_r1
        self.starttime = starttime
        self.bp_t0_r1 = bp_t0_r1
        self.evlo_r1 = evlo_r1
        self.evla_r1 = evla_r1
        self.evdp_r1 = evdp_r1

    def copy(self):
        """Returns a deepcopy of the backprojected stack.
        
        Returns:
        ClassName
            A deepcopy of the current instance.
        """
        return copy.deepcopy(self)

    def trim(self, t1, t2):
        """Trims the backprojected stack given a starttime and endtime.
        
        Parameters
        ----------
        t1: ObsPy.UTCDateTime
            A starttime.
        t2: ObsPy.UTCDateTime
            An endtime.
        """
        
        t1_ = t1 - self.starttime
        t2_ = t2 - self.starttime
        
        idx = np.where((self.bp_t0_r1 >= t1_) & (self.bp_t0_r1 <= t2_))[0]
        
        self.starttime = self.starttime + self.bp_t0_r1[idx[0]]
        self.b_r1 = self.b_r1[idx]
        self.bp_t0_r1 = self.bp_t0_r1[idx] - self.bp_t0_r1[idx[0]]
        self.evlo_r1 = self.evlo_r1[idx]
        self.evla_r1 = self.evla_r1[idx]
        self.evdp_r1 = self.evdp_r1[idx]

    def get_solution(self, t=None):
        """Computes and returns a solution for a given origin time.
        
        Parameters
        ----------
        t: None or float or ObsPy.UTCDateTime
            A time that represents an event origin time. If None, it returns a solution corresponding to the maximum brightness. If float, it returns a solution corresponding to the time (t) relative to the starttime of the stack. If ObsPy.UTCDateTime, it returns a solution corresponding to the time. Default is None.
            
        Returns:
        evt0: ObsPy.UTCDateTime
            Represents origin time of the event.
        evlo: float
            Represents longitude (deg) of the event.
        evla: float
            Represents latitude (deg) of the event.
        evdp: float
            Represents depth (km) of the event.
        b: float
            Represents brightness of the event.
        """

        if t == None:
            idx = np.argmax(self.b_r1)
            
        elif isinstance(t, float):
            idx = np.argmin(np.abs(self.bp_t0_r1 - t))

        elif isinstance(t, UTCDateTime):
            dif = t - self.starttime
            idx = np.argmin(np.abs(self.bp_t0_r1 - dif))

        evt0 = self.starttime + self.bp_t0_r1[idx]
        evlo = self.evlo_r1[idx]
        evla = self.evla_r1[idx]
        evdp = self.evdp_r1[idx]
        b = self.b_r1[idx]

        return evt0, evlo, evla, evdp, b

    def get_times(self, reftime=None):
        """Computes and returns a numpy array of relative times given a reference time.
        
        Parameters
        ----------
        reftime: ObsPy.UTCDateTime
            A reference time. Default is None. If None, it is the starttime of the stack.
            
        Returns
        -------
        times: numpy.ndarray of floats
            An array of times relative to the reference time.
        """

        if reftime == None:
            starttime_mpl = mdates.date2num(self.starttime)
            times = self.bp_t0_r1 / 86400 + starttime_mpl
        else:
            dif = self.starttime - reftime
            times = self.bp_t0_r1 + dif

        return times

    def plot(self, ax=None, ttype='relative', reftime=None, handle=False, **kwarg):
        """Plots the backprojected stack
        
        Parameters
        ----------
        ax: matplotlib.Axes
            A matplotlib.Axes object. Default None.
        ttype: str
            Time type. Options are 'relative' and 'absolute'. Default is 'relative'.
        reftime: ObsPy.UTCDateTime
            A reference time. Default is None.
        handle: bool
            If True, returns a matplotlib.Figure object. Default is False
        
        Returns
        -------
        handle: matplotlib.Figure
            S matplotlib.Figure object
        """

        if ax == None:
            fig = plt.figure(figsize=(10,3))
            ax = fig.add_subplot()
        else:
            fig = ax.figure

        if ttype == 'relative':
            if reftime == None:
                reftime = self.starttime
            times = self.get_times(reftime=reftime)

        if ttype == 'absolute':
            times = self.get_times()

        b_r1 = self.b_r1

        ax.plot(times, b_r1, **kwarg)

        if ttype == 'absolute':
            locator = mdates.AutoDateLocator()
            formatter = mdates.ConciseDateFormatter(locator)
            ax.xaxis.set_major_locator(locator)
            ax.xaxis.set_major_formatter(formatter)

        if handle == True:
            return fig

    def write(self, filename='brightness_stack.p'):
        """Writes the backprojected stack in a pickle file.
        
        Parameters
        ----------
        filename: str
            A filename to write the backprojected stack. Default is 'brightness_stack.p'
        """
        pickle.dump(self, open(filename, 'wb'))




class Brightness4():
    """
    A class for holding 4-D brightness values.
    
    This class provides methods to get backprojected solution, brightness stack, plot brightness slice etc.
    
    Attributes
    ----------
    B_r4: numpy.ndarray
        4-D brightness values.
    starttime: ObsPy.UTCDateTime
        Starttime of the brightness.
    bp_t0_r1: numpy.ndarray
        1-D array of origin times.
    lon_r1: numpy.ndarray
        1-D array of longitudes (degrees).
    lat_r1: numpy.ndarray
        1-D array of latitudes (degrees).
    z_r1: numpy.ndarray
        1-D array of depths (kilometers).
        
    Methods
    -------
    copy():
        Returns a deep copy of the 4-D brightness4 object.
    get_stack():
        Returns stack time series of the 4-D brightness volume.
    get_solution():
        Returns best solution of the 4-D brightness volume.
    get_islice():
        Returns a 2-D space slice of the 4-D brightness volume given a time index.
    get_slice():
        Returns a 2-D space slice of the 4-D brightness volume given a time.
    plot_slice():
        Plots a 2-D space slice of the 4-D brightness volume given a time.
    gaussian_filter():
        Smooths the 4-D brightness volume.
    write():
        Writes a pickle file for the 4-D brightness object.
    """

    def __init__(self, B_r4, starttime, bp_t0_r1, lon_r1, lat_r1, z_r1):
        """Initializes the Brightness4 class with brightness values, time and space.
        
        Parameters
        ----------
        B_r4: numpy.ndarray
            A 4-D array containing brightness values.
        starttime: ObsPy.UTCDateTime
            Starttime of the brightness volume.
        bp_t0_r1:
            1-D array of origin times.
        lon_r1: numpy.ndarray
            1-D array of longitudes (degrees).
        lat_r1: numpy.ndarray
            1-D array of latitudes (degrees).
        z_r1: numpy.ndarray
            1-D array of depths (kilometers).
        """
        
        self.B_r4 = B_r4
        self.starttime = starttime
        self.bp_t0_r1 = bp_t0_r1
        self.lon_r1 = lon_r1
        self.lat_r1 = lat_r1
        self.z_r1 = z_r1

    def copy(self):
        """
        Returns a deep copy of the 4-D brightness4 object.
        
        Returns
        -------
        ClassName
            A deepcopy of the current instance.
        """
        return copy.deepcopy(self)

    def get_stack(self):
        """
        Returns tack time series of the 4-D brightness volume.
        
        Returns
        -------
        stack: seisscan.Stack
            Backprojected stack time series.
        """
        
        bp_nt = len(self.bp_t0_r1)
        nx = len(self.lon_r1)
        ny = len(self.lat_r1)
        nz = len(self.z_r1)

        b_r1 = []
        evt0_r1 = []
        evlo_r1 = []
        evla_r1 = []
        evdp_r1 = []

        for i in range(bp_nt):
            idx = self.B_r4[:,:,:,i].argmax()
            iy, ix, iz = np.unravel_index(idx, (ny,nx,nz))
            b_ = self.B_r4[:,:,:,i].max()
            evlo_ = self.lon_r1[ix]
            evla_ = self.lat_r1[iy]
            evdp_ = self.z_r1[iz]

            b_r1.append(b_)
            evlo_r1.append(evlo_)
            evla_r1.append(evla_)
            evdp_r1.append(evdp_)

        b_r1 = np.array(b_r1)
        evlo_r1 = np.array(evlo_r1)
        evla_r1 = np.array(evla_r1)
        evdp_r1 = np.array(evdp_r1)

        stack = Stack(b_r1, self.starttime, self.bp_t0_r1, evlo_r1, evla_r1, evdp_r1)

        return stack

    def get_solution(self, t=None):
        """
        Returns best solution of the 4-D brightness volume for a given origin time.
        
        Parameters
        ----------
        t: ObsPy.UTCDateTime
            An origin time.
            Default is None.

        Returns
        -------
        evt0: ObsPy.UTCDateTime
            Event origin time
        evlo: float
            Event longitude (degrees)
        evla: float
            Event latitude (degrees)
        evdp: float
            Event depth (kilometers)
        """

        if t == None:
            ii_y, ii_x, ii_z, ii_t = np.unravel_index(self.B_r4.argmax(), self.B_r4.shape)
            
        elif isinstance(t, float):
            ii_t = np.argmin(np.abs(self.bp_t0_r1 - t))
            B_r3 = self.B_r4[:,:,:,ii_t]
            ii_y, ii_x, ii_z = np.unravel_index(B_r3.argmax(), B_r3.shape)

        elif isinstance(t, UTCDateTime):
            dif = t - self.starttime
            ii_t = np.argmin(np.abs(self.bp_t0_r1 - dif))
            B_r3 = self.B_r4[:,:,:,ii_t]
            ii_y, ii_x, ii_z = np.unravel_index(B_r3.argmax(), B_r3.shape)

        evt0 = self.starttime + self.bp_t0_r1[ii_t]
        evlo = self.lon_r1[ii_x]
        evla = self.lat_r1[ii_y]
        evdp = self.z_r1[ii_z]

        return evt0, evlo, evla, evdp

    def get_islice(self, plane='xy', ii_t=None):
        """
        Returns a 2-D space slice of the 4-D brightness volume given a time index.
        
        Parameters
        ----------
        plane: str
            'xy', 'yz' or 'xz' plane
            Default is 'xy'.
        ii_t: int
            time index
            
        Returns
        -------
        slice_plane: numpy.ndarray
            A 2-D slice of the 4-D brightness volume.
        """

        if ii_t == None:
            _, _, _, ii_t = np.unravel_index(self.B_r4.argmax(), self.B_r4.shape)

        B_r3 = self.B_r4[:,:,:,ii_t]

        if plane =='xy' or plane == 'yx':
            _, _, ii_z = np.unravel_index(B_r3.argmax(), B_r3.shape)
            slice_plane = B_r3[:,:,ii_z]

        if plane =='yz' or plane == 'zy':
            _, ii_x, _ = np.unravel_index(B_r3.argmax(), B_r3.shape)
            slice_plane = B_r3[:,ii_x,:].T

        if plane =='xz' or plane == 'zx':
            ii_y, _, _ = np.unravel_index(B_r3.argmax(), B_r3.shape)
            slice_plane = B_r3[ii_y,:,:].T

        return slice_plane

    def get_slice(self, plane='xy', t=None):
        """
        Returns a 2-D space slice of the 4-D brightness volume given a time.
        
        Parameters
        ----------
        plane: str
            'xy', 'yz' or 'xz' plane
            Default is 'xy'.
        t: ObsPy.UTCDateTime
            A time.
            Default is None.
            
        Returns
        -------
        slice_plane: numpy.ndarray
            A 2-D slice of the 4-D brightness volume.
        """

        if t == None:
            ii_t = None

        elif isinstance(t, float):
            ii_t = np.argmin(np.abs(self.bp_t0_r1 - t))

        elif isinstance(t, UTCDateTime):
            dif = t - self.starttime
            ii_t = np.argmin(np.abs(self.bp_t0_r1 - dif))

        slice_plane = self.get_islice(plane=plane, ii_t=ii_t)

        return slice_plane

    def plot_slice(self, ax, plane='xy', t=None):
        """
        Plots a 2-D space slice of the 4-D brightness volume given a time.
        
        Parameters
        ----------
        ax: matplotlib.Axes
            A matplotlib axes
        plane: str
            'xy', 'yz' or 'xz' plane
            Default is 'xy'.
        t: ObsPy.UTCDateTime
            A time.
            Default is None.
        
        Returns
        -------
        s: matplotlib.collections.QuadMesh
            A matplotlib quadmesh.
        """

        slice_plane = self.get_slice(plane=plane, t=t)

        if plane =='xy' or plane == 'yx':
            s = ax.pcolormesh(self.lon_r1, self.lat_r1, slice_plane)

        if plane =='yz' or plane == 'zy':
            s = ax.pcolormesh(self.lat_r1, self.z_r1, slice_plane)

        if plane =='xz' or plane == 'zx':
            s = ax.pcolormesh(self.lon_r1, self.z_r1, slice_plane)

        return s

    def gaussian_filter(self, sigma, radius=None, axes=None):
        """
        Smooths the 4-D brightness volume.
        
        Parameters
        ----------
        sigma: float
            Standard deviation for Gaussian kernel.
        radius: int
            Radius of the Gaussian kernel.
            Default is None.
        axes: tuple of int or None, optional
            Filter along the axes.
            Default is None.
        """
        
        self.B_r4 = scipy.ndimage.gaussian_filter(self.B_r4, sigma, radius=radius, axes=axes)

    def write(self, filename='brightness4.p'):
        """
        Writes a pickle file for the 4-D brightness object.
        
        Parameters
        ----------
        filename: str
            A pickle filename to save the 4-D brightness object.
        """
        pickle.dump(self, open(filename, 'wb'))
        