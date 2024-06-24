import numpy as np

def normalize_data(x):
    x_norm = (x - x.mean()) / x.std()
    return x_norm


def ncc(x, y, method=0):
    '''
    Description:
        Computes normalized correlation coefficient (NCC).
        Given two 1D ndarray of equal length (n) such as x and y NCC is as follows,
        
                          sum[   ( x - mean(x) )  *  ( y - mean(y) )   ]
            C = ------------------------------------------------------------------------------
                [  sum(  ( x - mean(x) )**2  )    *    sum(   ( x - mean(x) )**2   )  ]**(0.5)
            
            or,
            
            C = (1/n)  *  sum[  (x - mean(x)) / std(x)   *   (y - mean(y)) / std(y)  ]
            
            or,
            
            C = (1/n)  *  sum[  x_norm   *   y_norm  ]
        
    Input:
        x: 1D numpy ndarray of length n. n is a positive integer.
        y: 1D numpy ndarray of length n. n is a positive integer.
        method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                If method = 0, it returns C       (Rogers and Nicewander, 1987)
                If method = 1, it returns |C|     (used by Li et al., 2018)
                If method = 2, it returns C*C
                If method = 3, it returns C*|C|   (Gibbons, 2022)
                Default method is 0.
                
    Output:
        out: One of C, |C|, C*C and C*|C|
    '''
    
    #--- length of input data
    n = len(x)

    #--- normalize input data
    x_norm = normalize_data(x)
    y_norm = normalize_data(y)

    #--- computation
    c = (1/n) * np.sum(x_norm * y_norm)

    #--- different methods
    if method == 0:
        out = c
        
    if method == 1:
        out = np.abs(c)
        
    if method == 2:
        out = c * c
        
    if method == 3:
        out = c * np.abs(c)
        
    return out


def lagged_ncc(x, y, method=0):
    '''
    Description:
        Computes normalized correlation coefficient (NCC) at different lags.
        
    Input:
        x: 1D numpy ndarray of length n. n is a positive integer.
        y: 1D numpy ndarray of length n. n is a positive integer.
        method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                If method = 0, it returns C
                If method = 1, it returns |C|
                If method = 2, it returns C*C
                If method = 3, it returns C*|C|
                Default method is 0.
                
    Output:
        lags: 1D numpy ndarray of all possible lags
        out: 1D numpy ndarray of one of lagged C, |C|, C*C and C*|C|
    '''

    #--- length of input data
    n = len(x)

    #--- normalize input data
    x_norm = normalize_data(x)
    y_norm = normalize_data(y)

    #--- perform cross-correlation
    cc = (1/n) * np.correlate(x_norm, y_norm, mode='full')

    #--- different methods
    if method == 0:
        out = cc
        
    if method == 1:
        out = np.abs(cc)
        
    if method == 2:
        out = cc * cc
        
    if method == 3:
        out = cc * np.abs(cc)

    #--- Generate all lag [-n+1, .... -2, -1, 0, 1, 2, ..., n-1]
    lags = np.arange(-n+1, n)
    
    return lags, out


def max_lagged_ncc(x, y, max_ilag=0, method=0):
    '''
    Description:
        Computes maximum normalized correlation coefficient (NCC) at within certain lag.
        
    Input:
        x: 1D numpy ndarray of length n. n is a positive integer.
        y: 1D numpy ndarray of length n. n is a positive integer.
        max_ilag: A non-negative integer. An array of NCC will be claculated
                   for lags [-max_ilag, -max_ilag+1, ..., -2, -1, 0, 1, 2, ..., max_ilag-1, max_ilag]
                  Default max_ilag is 0. It is simply NCC (one of C, |C|, C*C, C*|C|).
                  If max_ilag = -1, it computes NCC for all possible lags.
        method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                If method = 0, it returns C
                If method = 1, it returns |C|
                If method = 2, it returns C*C
                If method = 3, it returns C*|C|
                Default method is 0.
                
    Output:
        max_cc: maximum NCC within max_ilag
    '''
    
    #--- length of input data
    n = len(x)
    
    #--- check max_ilag
    if max_ilag == -1:
        max_ilag = n - 1
    
    #--- perform lagged NCC
    _, cc = lagged_ncc(x, y, method=method)
    
    #--- zeroth lag
    ilag_0 = n - 1
    
    #--- get max cc within certain lag (max_ilag)
    max_cc = cc[ilag_0-max_ilag:ilag_0+max_ilag+1].max()
    
    return max_cc


def sliding_max_lagged_ncc(x, y, nperseg=100, stride=1, max_ilag=0, method=0):
    '''
    Description:
        Computes sliding maximum normalized correlation coefficient (NCC) at within certain lag.
        
    Input:
        x: 1D numpy ndarray of length n. n is a positive integer.
        y: 1D numpy ndarray of length n. n is a positive integer.
        nperseg: number of samples in each sliding window. Default is 100.
        stride: number of points shift between successive windows (int). Default is 1.
        max_ilag: A non-negative integer. An array of NCC will be claculated
                   for lags [-max_ilag, -max_ilag+1, ..., -2, -1, 0, 1, 2, ..., max_ilag-1, max_ilag]
                  Default max_ilag is 0. It is simply NCC (one of C, |C|, C*C, C*|C|).
                  If max_ilag = -1, it computes NCC for all possible lags.
        method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                If method = 0, it returns C
                If method = 1, it returns |C|
                If method = 2, it returns C*C
                If method = 3, it returns C*|C|
                Default method is 0.
                
    Output:
        max_cc: maximum NCC within max_ilag
    '''

    #--- create view for x and y
    x_v = np.lib.stride_tricks.sliding_window_view(x, nperseg)[::stride,:]
    y_v = np.lib.stride_tricks.sliding_window_view(y, nperseg)[::stride,:]

    #--- perform sliding max_lagged_ncc
    max_cc_list = []

    for x_, y_ in zip(x_v, y_v):
        max_cc = max_lagged_ncc(x_, y_, max_ilag=max_ilag, method=method)
        max_cc_list.append(max_cc)

    max_cc_array = np.array(max_cc_list)
    
    return max_cc_array


def sliding_template_ncc(x, y, stride=1, method=0):
    '''
    Description:
        Computes sliding normalized correlation coefficient (NCC) using a template.
        
    Input:
        x: 1D numpy ndarray of length n. n is a positive integer.
        y: 1D numpy ndarray of length nperseg. nperseg is a positive integer. This is a template.
        stride: number of points shift between successive windows (int). Default is 1.
        method: One of the integer (0, 1, 2, 3). It determines type of NCC (C) to return.
                If method = 0, it returns C
                If method = 1, it returns |C|
                If method = 2, it returns C*C
                If method = 3, it returns C*|C|
                Default method is 0.
                
    Output:
        cc_array: 1D numpy ndarray of NCC
    '''
    
    #--- define nperseg
    nperseg = len(y)

    #--- create view for x
    x_v = np.lib.stride_tricks.sliding_window_view(x, nperseg)[::stride,:]

    #--- perform sliding ncc
    cc_list = []

    for x_ in x_v:
        cc = ncc(x_, y, method=method)
        cc_list.append(cc)

    cc_array = np.array(cc_list)
    
    return cc_array