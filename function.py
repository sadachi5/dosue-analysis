import numpy as np
import csv
import pandas as pd
import struct

# Constants
v_c = 220.e+3 # [m/sec] speed of solar system
v_E = v_c # [m/sec] speed of earth --> TODO: Should be calculated for the measurement data
c = 299792458. # [m/sec] speed of light from wikipedia
k_B = 1.380649e-23 # [J/K] boltzmann constant
rbw = 3.e+2 # [Hz]
binwidth = 2.e+3 # [Hz]
TLN2 = 77 # [K]
nrun = 12 # nRun in search measurement

check_freq = np.array([])


# Read file

def dat_to_array(path, doRebin=False, rebinmethod=0, rbw=rbw, binwidth=binwidth, binary=False):
    if binary:
        f = open(path, 'rb')
        freq = []
        dBm = []
        while True:
            read = f.read(16)
            if len(read) == 0:
                break
            v1, v2, = struct.unpack('dd', read)
            freq.append(v1)
            dBm.append(v2)
            pass
        freq = np.array(freq)
        dBm = np.array(dBm)
    else:
        f = open(path, "r")
        data = f.read().split("\n")[3:-1]
        freq = np.array([float(val.split(" ")[0]) for val in data])
        dBm = np.array([float(val.split(" ")[1]) for val in data])
        pass
    W = dBm_to_W(dBm)
    f.close()

    Werr = None
    if doRebin:
        freq, W, Werr = rebin_func_consider_rbw(
                    freq, W, rebin=binwidth, rbw=rbw, method=rebinmethod, verbose=0)
    else:
        Werr = np.full(len(W), 0.)
        pass

    return freq, W, Werr


def csv_to_array(path):
    result = {}
    df = pd.read_csv(path)
    for col in df.columns:
        result[str(col)] =  np.array(df[col])
    
    return result

def cpy_to_dict(path):
    return np.load(path, allow_pickle=True).item()

def get_one_freq_data(
    start, nRun=12,
    datadir='/data/ms2840a/dosue-j/signal_data/2023-03',
    dataprefix='scan_FFT_',
    datasuffix='GHz_span2.50MHz_rbw300Hz_2.0sec_1counts_12runs',
    binary_data=True,
    doRebin=True, rebinmethod=0, binwidth=binwidth,
    onlyAverage=True, cutEdges=True
    ):

    path = f"{datadir}/{dataprefix}{start}{datasuffix}"

    f_list = []
    W_list = []
    Werr_list = []
    for n in range(nRun):
        filesuffix = f'_{n}' if nRun > 1 else  ''
        _path = f'{path}{filesuffix}.dat'
        f, W, Werr = dat_to_array(_path, doRebin=doRebin, rebinmethod=rebinmethod, rbw=rbw, binwidth=binwidth, binary=binary_data)
        if cutEdges: f, W, Werr = cut_data(f, W, Werr)
        f_list.append(f)
        W_list.append(W)
        Werr_list.append(Werr)
        if n > 1:
            if np.sum((f != f_list[0])) > 0:
                print(f"Error! There are different frequency range in {nRun} data files.")
                return -1
            pass
        pass
    f_list = np.array(f_list)
    W_list = np.array(W_list)
    Werr_list = np.array(Werr_list)
    
    f_ave = f_list[0]
    W_ave = np.average(W_list, axis=0)
    Werr_ave = np.average(Werr_list, axis=0)
    
    if onlyAverage:
        return f_ave, W_ave, Werr_ave
    else:
        return {'average':[f_ave, W_ave, Werr_ave], 'raw':[f_list, W_list, Werr_list]}


def read_data(
    freq_min, freq_max, nRun=12, 
    datadir='/data/ms2840a/dosue-j/signal_data/2023-03',
    dataprefix='scan_FFT_',
    datasuffix='GHz_span2.50MHz_rbw300Hz_2.0sec_1counts_12runs',
    binary_data=True,
    doRebin=True, rebinmethod=0, binwidth=binwidth,
    onlyAverage=True, cutEdges=True, flatten=True
):
    all_freq_ave = []
    all_W_ave = []
    all_Werr_ave = []
    all_freq_list = []
    all_W_list = []
    all_Werr_list = []
    
    i_start100MHz, i_end100MHz, i_100MHz = get_loop_index_each100MHz(freq_min, freq_max)
    for i in range(i_start100MHz, i_end100MHz, i_100MHz):
        i_start2MHz, i_end2MHz, i_2MHz = get_loop_index_each2MHz(i)
        for j in range(i_start2MHz, i_end2MHz, i_2MHz):
            if (j-i_start2MHz) % (10*i_2MHz) == 0: print(f'freq = {j*1e-9:.6f} GHz')
            start_Hz = j
            center_Hz = start_Hz + i_2MHz/2
            start_str, start_100MHz_str, is_add_data\
                = get_file_freq(center_Hz, verbose=0)
            
            data = get_one_freq_data(start_str, nRun=nRun,
                                           datadir=datadir, dataprefix=dataprefix, datasuffix=datasuffix,
                                           binary_data=binary_data,
                                           doRebin=doRebin, rebinmethod=rebinmethod, binwidth=binwidth,
                                           onlyAverage=onlyAverage, cutEdges=cutEdges)
            if onlyAverage:
                _f_ave, _W_ave, _Werr_ave = data[0], data[1], data[2]
                all_freq_ave.append(_f_ave)
                all_W_ave.append(_W_ave)
                all_Werr_ave.append(_Werr_ave)
            else:
                data_ave = data['average']
                data_raw = data['raw']
                _f_ave, _W_ave, _Werr_ave = data_ave[0], data_ave[1], data_ave[2]
                _f_list, _W_list, _Werr_list = data_raw[0], data_raw[1], data_raw[2]
                all_freq_ave.append(_f_ave)
                all_W_ave.append(_W_ave)
                all_Werr_ave.append(_Werr_ave)
                all_freq_list.append(_f_list)
                all_W_list.append(_W_list)
                all_Werr_list.append(_Werr_list)
                pass
            pass
        pass
    all_freq_ave = np.array(all_freq_ave)
    all_W_ave = np.array(all_W_ave)
    all_Werr_ave = np.array(all_Werr_ave)
    if flatten:
        all_freq_ave = all_freq_ave.flatten()
        all_W_ave = all_W_ave.flatten()
        all_Werr_ave = all_Werr_ave.flatten()
        pass
    
    if onlyAverage:
        return all_freq_ave, all_W_ave, all_Werr_ave
    else:
        return {'average':[all_freq_ave, all_W_ave, all_Werr_ave], 'raw':[all_freq_list, all_W_list, all_Werr_list]}


def original_signal_to_array(path, n_data=nrun, verbose=0):
    signal = csv_to_array(path)
    freq = signal['freq']
    keys_W = [ f'W_{i}' for i in range(n_data) ]
    W_array = np.array([ signal[key] for key in keys_W ])
    W = np.average(W_array, axis=0) # mean of (W_0, W_1,... W_{n_data})
    if verbose > 0:
        print(f'data keys = {signal.keys()}')
        print(f'power keys = {keys_W}')
        print(f'freq (size={len(freq)}) = ',freq)
        print(f'W (nominal data, size={len(W)}) = ',W)
        pass
    return freq, W, W_array


# Y-factor method

def yfactor_analysis(freq, Wamb, WLN2, Wamb_err, WLN2_err, Tamb, rbw=rbw):
    # ２点を通る直線　y = ax + b
    a = (Wamb - WLN2) / (Tamb - TLN2)
    b = Wamb - a * Tamb
    if (Wamb_err is not None) and (WLN2_err is not None):
        a_err = np.sqrt( np.power(Wamb_err, 2.) + np.power(WLN2_err, 2.) ) / (Tamb - TLN2)
        b_err = np.sqrt( np.power(Wamb_err, 2.) + np.power(a_err*Tamb, 2.) )
    else:
        a_err = np.full(len(a), 0.)
        b_err = np.full(len(b), 0.)
        pass
    #print_list(a, 'a')
    #print_list(b, 'b')
    #print_list(a_err, 'a_err')
    #print_list(b_err, 'b_err')

    gain = a/k_B/rbw
    Trx = b/a
    gain_err = a_err/k_B/rbw
    Trx_err = np.sqrt( np.power(b_err/a, 2.) + np.power(a_err*b/(a*a), 2.) )

    return gain, Trx, gain_err, Trx_err


# Data select for fit

def cut_data(x, y, yerr=None):
    freq = []
    W = []
    Werr = []
    for i, (a, b) in enumerate(zip(x, y)):
        if a >= x[0] + 250.e+3 and a < x[0] + 250.e+3 + 2.e+6:
            freq.append(a)
            W.append(b)
            if yerr is not None:
                Werr.append(yerr[i])
                pass
            pass
        pass
        
    return np.array(freq), np.array(W), np.array(Werr)


# Rebin or average

def rebin_func(freq, data, rebin=binwidth):
    rebin_freq = []
    rebin_data = []
    rebin_data_std = []
    freq_0 = freq[0]
    save_data = []
    save_freq = []
    for x, y in zip(freq, data):
        if x < freq_0 + rebin:
            save_data.append(y)
            save_freq.append(x)
        else:
            rebin_data.append(np.mean(np.array(save_data)))
            rebin_freq.append(np.mean(np.array(save_freq)))
            rebin_data_std.append(np.std(np.array(save_data))/len(save_data)**0.5)
            freq_0 += rebin
            save_data = [y]
            save_freq = [x]

    return np.array(rebin_freq), np.array(rebin_data), np.array(rebin_data_std)

def rebin_func_consider_rbw(freq, data, rebin=binwidth, rbw=rbw, method=0, verbose=0):
    '''
        freq: original freq array [Hz]
        data: original power array [W]
        rebin: rebinning width [Hz] = 2000 Hz
        rbw: RBW of taken data [Hz] = 300 Hz
        verbose: verbosity level (int)
    '''
    rebin_freq = []
    rebin_data = []
    rebin_data_err = []

    freq = np.array(freq)
    freq1 = np.array(freq - rbw/2.) # lower edge of each data
    freq2 = np.array(freq + rbw/2.) # upper edge of each data
    data = np.array(data)
    if verbose > 0:
        print(f'freq  = {freq}')
        print(f'freq1 = {freq1}')
        print(f'freq2 = {freq2}')
        print(f'dfreq = {freq2-freq1}')
        print(f'data = {data}')
        pass

    # Create rebinned freq
    rebin_freq = np.arange(freq[0]+rebin/2, freq[-1]+rebin/2, rebin)
    if verbose > 0: print_list(rebin_freq, 'rebin_freq')

    # Loop over new freq
    for i, rebin_x in enumerate(rebin_freq):
        rebin_x1 = rebin_x - rebin/2. # lower edeg
        rebin_x2 = rebin_x + rebin/2. # upper edge

        if   method == 0: 
            in_range = np.where((freq2 > rebin_x1) & (freq1 < rebin_x2)) # Check lower edge and upper edge (consider bin edges of the original binning)
        elif method == 1:
            in_range = np.where((freq >= rebin_x1) & (freq < rebin_x2)) # Check lower edge and upper edge (consider only bin centers of the original binning)
            pass
        if verbose > 0: print_list(in_range, 'in_range')

        # Retrieve original data in the range
        x = freq[in_range]
        x1 = freq1[in_range]
        x2 = freq2[in_range]
        y = data[in_range]
        n = len(x)
        width = x[1]-x[0]
        if verbose > 0: print_list(x, 'x')

        # Weight considering width in range
        weight = np.full(n, 1.)
        if method == 0:
            # Considering lower edge
            weight[0] = (x2[0] - rebin_x1)/rbw
            # Considering upper edge
            weight[-1] = (rebin_x2 - x1[-1])/rbw
            pass
        
        total_y = np.sum(y*weight)
        total_weight = np.sum(weight)
        total_width = total_weight *rbw
        _rebin_data = total_y * rebin/total_width
        rebin_data.append( _rebin_data )
        rebin_data_err.append( np.sqrt( np.sum( np.power(y*rebin/rbw - _rebin_data, 2.)*weight )/total_weight ) )
        pass

    return np.array(rebin_freq), np.array(rebin_data), np.array(rebin_data_err)

def average_list(data, naverage=100):

    ndata = len(data)
    npoints = int(ndata/naverage)
    
    data_ave = []
    data_err = []
    
    for i in range(npoints):
        data_subset = data[i*naverage:(i+1)*naverage]
        average = np.mean(data_subset)
        average_err = np.std(data_subset)/np.sqrt(naverage) #  = 1/N * sqrt( sum((y-mean)^2))  ( std = sqrt( sum((y-mean)^2) / N) )
        data_ave.append(average)
        data_err.append(average_err)
        pass
    
    return np.array(data_ave), np.array(data_err)

def rebin_hist(bins, hist, rebin, histerr=None):
    nbins = len(hist)
    nbins_new = int(nbins/rebin)
    hist_new = np.zeros( nbins_new )
    bins_new = np.zeros( nbins_new + 1 )

    if histerr is not None: 
        histerr = np.array(histerr)
        histerr_new = np.zeros( nbins_new )
    else:
        histerr_new = None
        pass

    for i in range(nbins_new):
        i_start = i*rebin
        i_end = i_start + rebin
        # bins
        bin_lowedge = bins[i_start]
        bins_new[i] = bin_lowedge
        # hist
        hist_subset = hist[i_start:i_end]
        _y = np.sum(hist_subset)
        hist_new[i] = _y
        # histerr
        if histerr is not None:
            histerr_subset = histerr[i_start:i_end]
            _yerr = np.sqrt(np.sum(np.power( histerr_subset, 2. )))
            histerr_new[i] = _yerr
            pass
        # Add the most right bin edge at last bin
        if i == nbins_new - 1:
            bins_new[i+1] = bins[i_end]
            pass
        pass
    
    return bins_new, hist_new, histerr_new 


# Misc.

def W_to_dBm(W, Werr=None):
    dBm = np.log10(W*1.0e+3)*10.
    if Werr is None:
        return dBm
    else:
        dBm_up = np.log10((W+Werr)*1.0e+3)*10.
        dBm_down = np.log10((W-Werr)*1.0e+3)*10.
        dBm_err = np.array( [dBm_up - dBm, dBm - dBm_down] )
        return dBm, dBm_err 

def dBm_to_W(dBm):
    return np.power(10., dBm*0.1) * 1.0e-3

def ratio_to_dB(ratio, ratio_err=None):
    dB = np.log10(ratio)*10.
    if ratio_err is None:
        return dB
    else:
        dB_up = np.log10(ratio+ratio_err)*10.
        dB_down = np.log10(ratio-ratio_err)*10.
        dB_err = np.array( [dB_up - dB, dB - dB_down] )
        return dB, dB_err

def dB_to_ratio(dB):
    return np.power(10., dB*0.1)

def round_GHz_upto_kHz(GHz):
    '''
        input: GHz = float [GHz]
        output:  string [GHz] up to kHz (ex. 20.000000 GHz)
    '''
    if isinstance(GHz, str): GHz = (float)(GHz) 
    GHz_upto_kHz = round(GHz, 6)
    return f'{GHz_upto_kHz:.06f}'

def get_file_freq(freq0, verbose=0):
    '''
    Arguments:
        freq0: frequency to be used in analysis [Hz]
    
    Return:
        start_str, start_100MHz_str, is_add_data

        start_str: GHz frequency string for 2MHz span file name [ex. 17.999750]
        start_100MHz_str: GHz frequency string for 100MHz span file name [ex. 19.1]
        is_add_data: if this frequency is in additional data or not
    '''
    freq0_MHz = freq0*1e-6
    
    is_add_data = False
    for _check_freq in check_freq:
        # Check if freq0 is in 2MHz measured data span in additional datas
        if _check_freq <= freq0_MHz and freq0_MHz < _check_freq + 2.:
            is_add_data = True
            pass
        pass
    if verbose>0: 
        print('This frequency has an additional data? -->', is_add_data)
        pass 
    
    # Get start freq of 100MHz span
    start_100MHz = (freq0//100e+6)*100 # MHz (100MHz 毎にする)
    start_100MHz *= 1e-3 # GHz
    start_100MHz = round(start_100MHz, 6) # GHz で小数点以下6桁(kHz)までにする
    start_100MHz_str = f'{start_100MHz:.1f}'
    if verbose>0: print(f'start_100MHz = {start_100MHz} GHz')
    
    # Get start freq of 2MHz span
    start_2MHz = (freq0//2e+6)*2 # MHz (2MHz 毎にする)
    start_2MHz *= 1e-3 # GHz
    start_2MHz = round(start_2MHz, 6) # GHz で小数点以下6桁(kHz)までにする
    if verbose>0: print(f'start_2MHz = {start_2MHz} GHz')
    
    # Get start freq of data span (start_2MHz - 0.25 MHz)
    start = start_2MHz - 0.00025
    start_str = round_GHz_upto_kHz(start)
    if verbose>0: print(f'start_str = {start_str}')

    return start_str, start_100MHz_str, is_add_data

def get_loop_index_each100MHz(freq_min, freq_max):
    # freq_min: GHz with 1 decimal place / ex) 10.0
    # freq_max: GHz with 1 decimal place / ex) 10.1
    # return: (int)i_start [Hz], (int)i_end [Hz], (int)interval [Hz]
    i_start = int(freq_min*1e+9) # [Hz]
    i_end = int(freq_max*1e+9) # [Hz]
    interval = int(1e+8) # 100MHz [Hz]
    return i_start, i_end, interval

def get_loop_index_each2MHz(index_100MHz, interval=int(2e+6)):
    # index_100MHz: [Hz] output of get_loop_index_each100MHz()
    # return: (int)i_start [Hz], (int)i_end [Hz], (int)interval [Hz]
    start0_Hz = index_100MHz # [Hz] --> [Hz]
    start_GHz_str, start_100MHz_str, is_add_data = get_file_freq(start0_Hz)

    i_start = int(float(start_GHz_str)*1e+9) # [GHz] --> [Hz]
    i_end = int(i_start + int(1e+8)) # 100MHz span [Hz]
    return i_start, i_end, interval


def listup_low_plocal(freq, p_local, lowP = 1.e-5):
    freq_lowP = []
    p_local_lowP = []
    for i, p in enumerate(p_local):
        if p < lowP:
            freq_lowP.append(freq[i])
            p_local_lowP.append(p)
            print(
                str(freq[i]/1e9) + " & " + 
                str(round(p*1e6,2))
                )
            pass
        pass
    print(f'# of fit with p_local<1e-5 = {len(freq_lowP)}')
    print()
    
    # Remove adjacent results
    print('###############################################')
    print('Remove adjacent results')
    diff_freq_lowP = np.diff(freq_lowP)
    close_lowP = (diff_freq_lowP < 2.1e+3) # Hz
    _plocal_close = []
    _index_close = []
    remove_index = []
    for i, isclose in enumerate(close_lowP):
        if isclose:
            _plocal_close.append(p_local_lowP[i])
            _index_close.append(i)
        else:
            if len(_plocal_close)>0:
                # check lowest p_local
                _plocal_close.append(p_local_lowP[i])
                _index_close.append(i)
                _plocal_min = min(_plocal_close)
                _index_min = _index_close[_plocal_close.index(_plocal_min)]
                # append removing index datas
                for _ind in _index_close:
                    if _ind != _index_min:
                        remove_index.append(_ind)
                        pass
                # clear 
                _plocal_close = []
                _index_close = []
                pass
            pass
        pass
    
    i = len(close_lowP)
    if len(_plocal_close)>0:
        # check lowest p_local
        _plocal_close.append(p_local_lowP[i])
        _index_close.append(i)
        _plocal_min = min(_plocal_close)
        _index_min = _index_close[_plocal_close.index(_plocal_min)]
        # append removing index datas
        for _ind in _index_close:
            if _ind != _index_min:
                remove_index.append(_ind)
                pass
            pass
        pass
    
    freq_lowP = [ f for i, f in enumerate(freq_lowP) if i not in remove_index ]
    p_local_lowP = [ f for i, f in enumerate(p_local_lowP) if i not in remove_index ]
    
    # print
    for i in range(len(freq_lowP)):
        print(
            str(freq_lowP[i]/1e9) + " & " + 
            str(round(p_local_lowP[i]*1e6,2))
            )
        pass
            
    
    freq_lowP = np.array(freq_lowP)
    p_local_lowP = np.array(p_local_lowP)
    print(f'# of fit with p_local<1e-5 = {len(freq_lowP)}')

    return freq_lowP, p_local_lowP

import inspect
def get_var_name(var, back_vars=None):
    name = ''
    if back_vars == None:
        back_vars = inspect.currentframe().f_back.f_globals
        back_vars.update( inspect.currentframe().f_back.f_globals )
        pass;
    for k, v in back_vars.items():
        if id(v) == id(var):
            name=k
            pass
        pass
    return name

def print_list(var, varname=''):
    if isinstance(var, np.ndarray): size = var.shape
    else: size = len(var)
    if varname != '':
        print(f'{varname} (size={size}) = {var}')
    else:
        back_vars = inspect.currentframe().f_back.f_globals
        back_vars.update( inspect.currentframe().f_back.f_locals )
        varname = get_var_name(var, back_vars)
        print(f'{varname} (size={size}) = {var}')
        pass
    pass

def isNoneAny_array(array):
    isNone = [ _x is None for _x in array ]
    return np.any(isNone)
