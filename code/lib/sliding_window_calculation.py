import numpy as np


def sliding_window_view(x, nperseg, noverlap=None, pos='mid'):
    '''
    Description:
        Form a 2D numpy.ndarray sliding window view from an 1D numpy.ndarray.
        This helps to perform moving window calculation.

    Input:
        x: 1d numpy.ndarray. (float)
        nperseg: number of elements in each sliding window. (int)
        noverlap: number of elements overlapped by two consecutive sliding windows. Default is None. (int)
        pos: Position of index. 'Start', 'mid' or 'end' for each sliding window. Default is 'mid'.
            if nperseg is odd, mid is nperseg//2
            if nperseg is even, mid is nperseg//2 - 1

    Output:
        idx: 1d numpy.ndarray of indices. (int)
        x_v: 2D numpy.ndarray sliding window view. (float)
            Total number of row is len(idx).
            Total number of column is nperseg.

    Usage:
        np.random.seed(0)
        x = np.random.normal(0, 1, 30)
        nperseg = 5
        noverlap = 3
        pos = 'mid
        
        idx, x_v = sliding_window_view(x, nperseg, noverlap=noverlap, pos=pos)
    '''

    if noverlap == None:
        noverlap = nperseg - 1

    if pos == None:
        pos = 'mid'

    #--- stride
    stride = nperseg - noverlap

    #--- Prepare output index
    i = np.arange(len(x), dtype='int')
    i_v = np.lib.stride_tricks.sliding_window_view(i, nperseg)
    i_v = i_v[::stride,:]

    if pos == 'start':
        idx = i_v[:,0]

    if pos == 'mid':
        if nperseg % 2 == 1:
            idx = i_v[:,nperseg//2]
        if nperseg % 2 == 0:
            idx = i_v[:,nperseg//2 - 1]

    if pos == 'end':
        idx = i_v[:,-1]

    #--- prepar data view
    x_v = np.lib.stride_tricks.sliding_window_view(x, nperseg)
    x_v = x_v[::stride,:]

    return idx, x_v



def sliding_window_calculation(x, nperseg, noverlap=None, fs=None, calculation = 'mean', pos=None, smallest_number = 1e-10):
    '''
    Description:
        Perform moving window calculation on a 1D numpy array.

    Input:
        x: 1d numpy.ndarray. (float)
        nperseg: number of elements in each sliding window. (int)
        noverlap: number of elements overlapped by two consecutive sliding windows. Default is None. (int)
        calculation: 'sum', 'mean', 'gmean' or 'hmean'
        pos: 'Start', 'mid' or 'end'
        smallest_number: smallest number close to zero.

    Output:
        idx: 1d numpy.ndarray of indices. (int)
        y: Output 1d numpy.ndarray. (float)

    Usage:
        np.random.seed(0)
        x = np.random.normal(0, 1, 30)
        nperseg = 5
        noverlap = 3
        calculation = 'mean'
        pos = 'mid

        idx, y = moving_window_calculation(x, nperseg, noverlap=noverlap, calculation=calculation, pos=pos)
    '''

    idx, x_v = sliding_window_view(x, nperseg, noverlap=noverlap, pos=pos)

    if fs != None:
        dt = 1 / fs
        t = idx * dt

    if calculation == 'sum':
        y = np.sum(x_v, axis=1)

    if calculation == 'mean':
        y = np.mean(x_v, axis=1)

    if calculation == 'gmean':
        y = np.prod(x_v, axis=1) ** (1/nperseg)

    if calculation == 'hmean':
        x_v = np.where(x_v < smallest_number, smallest_number, x_v)
        y = nperseg / np.sum(1 / x_v, axis=1)

    if fs == None:
        return idx, y
    else:
        return t, y

