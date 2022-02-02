#!/bin/env python
import os, argparse
import pathlib

import numpy as np
np.set_printoptions(threshold=10)
from matplotlib import pyplot as plt

g_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:olive','tab:cyan','tab:gray','red','royalblue','turquoise','darkolivegr    een', 'magenta', 'blue', 'green']*5;

class DATA:
    def __init__(self, powDBm=None, amp=None, freq=None, freq_start=None, freq_span=None, MesTimePoint=None, rebin=1, rbw=None):
        if freq is not None: self._freq = freq
        else               : self._freq = np.linspace(freq_start,freq_start+freq_span,len(powDBm))
        if powDBm is not None: 
            self._powDBm = np.array(powDBm)
            self._amp    = np.power(10., self._powDBm*0.1)*1.e-3
        else                 : 
            self._amp    = np.array(amp)
            self._powDBm = np.log10(self._amp*1.e+3)*10.
            pass
        self._amp_err      = None
        self._powDBm_err_p = None
        self._powDBm_err_m = None
        self._rebin = rebin
        if rebin>1:
            new_nbin = (int)(len(self._freq)/rebin)
            ampEachRebins  = self._amp[:new_nbin*rebin] # ignore several bins in the end of array, which are the remainders after dividing by 'rebin'
            ampEachRebins  = ampEachRebins.reshape(rebin, new_nbin) # [[a_1,a_2,...,a_rebin],[a_rebin+1,a_rebin+2,...,a_rebin+rebin],...]
            freqEachRebins = self._freq[:new_nbin*rebin] # ignore several bins in the end of array, which are the remainders after dividing by 'rebin'
            freqEachRebins = freqEachRebins.reshape(rebin, new_nbin) # [[a_1,a_2,...,a_rebin],[a_rebin+1,a_rebin+2,...,a_rebin+rebin],...]
            freqRebin     = np.mean(freqEachRebins, axis=1)
            ampRebin      = np.mean(ampEachRebins, axis=1)
            ampRebin_err  = np.std(ampEachRebins, axis=1)
            # Convert the freq, amp, powDBm
            self._freq       = freqRebin
            self._amp        = ampRebin
            self._amp_err    = ampRebin_err
            self._powDBm       = np.log10(self._amp*1.e+3)*10.
            self._powDBm_err_p = np.abs( np.log10((self._amp+self._amp_err)*1.e+3)*10. - self._powDBm )
            self._powDBm_err_m = np.abs( np.log10((self._amp-self._amp_err)*1.e+3)*10. - self._powDBm )
            pass

        self._rbw    = rbw
        self._time   = MesTimePoint # measurement start unix-time
        self._npoints= len(self._freq)
        if len(self._freq) != len(self._powDBm):
            print(f'DATA:__init__(): Error! The sizes of power and freq. array are NOT the same!')
            pass

    @property
    def powDBm(self):
        return self._powDBm

    @property
    def powDBm_err(self):
        if self._powDBm_err_p is None or self._powDBm_err_m is None: return None
        return [self._powDBm_err_m, self._powDBm_err_p]

    @property
    def amp(self):
        return self._amp 

    @property
    def amp_err(self):
        return self._amp_err

    @property
    def freq(self):
        return self._freq

    @property
    def freq_binwidth(self):
        return np.diff(self._freq)

    @property
    def rbw(self):
        return self._rbw

    @property
    def npoints(self):
        return self._npoints

    @property
    def time(self):
        return self._time

def get_rbw(filepath, startstr, index=0, scale=1.):
    rbw = None
    with open(filepath) as f:
        for line in f:
            if line.startswith(startstr):
                rbw = (float)(line.strip().split()[index])*scale
                break
            pass
        pass
    if rbw is None:
        print(f'get_rbw(): Error! Could not get the RBW value from {filepath}.')
        pass
    return rbw

def read_data(filepath, rebin=1, verbose=0):

    print(f'read_data(): input filepath = {filepath}')
    csv = np.loadtxt(filepath, delimiter=' ', comments='#', dtype=float)
    rbw = get_rbw(filepath, '#RBW', index=2, scale=1.)
    data = DATA(freq=csv[:,0],powDBm=csv[:,1],rbw=rbw,rebin=rebin)
    if verbose>1:
        print('read_data(): Data RBW       [Hz] =',data.rbw)
        print('read_data(): Data Frequency [Hz] =',data.freq)
        print('read_data(): Data Power    [dBm] =',data.powDBm)
        print('read_data(): Data Power      [W] =',data.amp)
        pass

    return data

def plot(datas, outdir='./', outname='aho', doAve=False, verbose=0):

    # expand outdir for '~' (home directory)
    outdir = pathlib.Path(outdir).expanduser()

    # make directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        pass

    fig, ax = plt.subplots(2,1)
    fig.set_size_inches(6,8)
    for i, data in enumerate(datas):
        if data.amp_err is None:
            ax[0].plot(data.freq*1.e-9, data.powDBm, color=g_colors[i], linewidth=1)
            ax[1].plot(data.freq*1.e-9, data.amp   , color=g_colors[i], linewidth=1)
        else :
            ax[0].errorbar(data.freq*1.e-9, data.powDBm, yerr=data.powDBm_err, color=g_colors[i], linewidth=1)
            ax[1].errorbar(data.freq*1.e-9, data.amp   , yerr=data.amp_err   , color=g_colors[i], linewidth=1)
            pass
        #npoints1MHz = np.sum(data.freq<data.freq[0]+100.e+3)
        #nep = np.std(data.amp[0:npoints1MHz])
        #if verbose>0: print(f'NET first {npoints1MHz} points / {data.npoints} = {nep} W')
        pass
    if doAve: 
        data_ave = DATA(freq=datas[0].freq, amp=np.mean([data.amp for data in datas], axis=0), rbw=datas[0].rbw)
        if data.amp_err is None:
            ax[0].plot(data_ave.freq*1.e-9, data_ave.powDBm, color='red', linewidth=2)
            ax[1].plot(data_ave.freq*1.e-9, data_ave.amp   , color='red', linewidth=2)
        else :
            ax[0].errorbar(data.freq*1.e-9, data.powDBm, yerr=data.powDBm_err, color=g_colors[i], linewidth=1)
            ax[1].errorbar(data.freq*1.e-9, data.amp   , yerr=data.amp_err   , color=g_colors[i], linewidth=1)
            pass
        pass
    ax[0].set_xlabel('Frequency [GHz]')
    ax[0].set_ylabel('Power [dBm]')
    ax[0].grid()
    ax[1].set_xlabel('Frequency [GHz]')
    ax[1].set_ylabel('Power [W]')
    ax[1].grid()
    fig.tight_layout()

    fig.savefig(f'{outdir}/{outname}.pdf')

    return fig, ax


if __name__=='__main__':
    outname='aho2'
    #input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/2021-10-21/data'
    #input_files= ['plate8_10500.0MHz_RBW1000.0kHz_100times_1.0sec_20211021185051.dat']
    '''
    input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/scan3/2021-10-26/data'
    input_files= [
            'scan_FFT_20.0GHz_span1.00MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            'scan_FFT_20.0GHz_span2.50MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            'scan_FFT_20.0GHz_span5.00MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            #'scan_FFT_20.0GHz_span10.00MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            ]
    '''

    '''
    input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/2021-10-27/data'
    input_files= [
            'FFT_2.5MHz_10sec_12times_0.dat',
            'FFT_2.5MHz_10sec_12times_1.dat',
            'FFT_2.5MHz_10sec_12times_2.dat',
            'FFT_2.5MHz_10sec_12times_3.dat',
            'FFT_2.5MHz_10sec_12times_4.dat',
            'FFT_2.5MHz_10sec_12times_5.dat',
            'FFT_2.5MHz_10sec_12times_6.dat',
            'FFT_2.5MHz_10sec_12times_7.dat',
            'FFT_2.5MHz_10sec_12times_8.dat',
            'FFT_2.5MHz_10sec_12times_9.dat',
            ]
    '''

    '''
    input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/signal_data_test/2021-10-28/data'
    input_files= [
                'test_FFT_26.399750GHz_span2.50MHz_rbw0.3kHz_1.0sec_1counts_2runs_0.dat',
                'test_FFT_26.399750GHz_span2.50MHz_rbw0.3kHz_1.0sec_1counts_2runs_1.dat',
                ]
    '''

    input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/2021-10-28/data'
    input_files= [
                'SWEEP_plate_test.dat',
                ]

    #datas = [ read_data(f'{input_dir}/{infile}', rebin=100, verbose=2) for infile in input_files ]
    datas = [ read_data(f'{input_dir}/{infile}', rebin=1, verbose=2) for infile in input_files ]

    plot(datas, 'output/', outname, doAve=True, verbose=2)

    pass
