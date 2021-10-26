#!/bin/env python
import os, argparse
import pathlib

import numpy as np
np.set_printoptions(threshold=10)
from matplotlib import pyplot as plt

g_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:olive','tab:cyan','tab:gray','red','royalblue','turquoise','darkolivegr    een', 'magenta', 'blue', 'green']*5;

class DATA:
    def __init__(self, powDBm, freq=None, freq_start=None, freq_span=None, MesTimePoint=None, rbw=None):
        if freq is not None: self._freq = freq
        else               : self._freq = np.linspace(freq_start,freq_start+freq_span,len(powDBm))
        self._powDBm = np.array(powDBm)
        self._rbw    = rbw
        self._time   = MesTimePoint # measurement start unix-time
        self._npoints= len(freq)
        if len(freq) != len(powDBm):
            print(f'DATA:__init__(): Error! The sizes of power and freq. array are NOT the same!')
            pass

    @property
    def powDBm(self):
        return self._powDBm

    @property
    def amp(self):
        return np.power(10., self._powDBm*0.1)*1.e-3

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

def read_data(filepath, verbose=0):

    print(f'read_data(): input filepath = {filepath}')
    csv = np.loadtxt(filepath, delimiter=' ', comments='#', dtype=float)
    rbw = get_rbw(filepath, '#RBW', index=2, scale=1.)
    data = DATA(freq=csv[:,0],powDBm=csv[:,1],rbw=rbw)
    if verbose>1:
        print('read_data(): Data RBW       [Hz] =',data.rbw)
        print('read_data(): Data Frequency [Hz] =',data.freq)
        print('read_data(): Data Power    [dBm] =',data.powDBm)
        print('read_data(): Data Power      [W] =',data.amp)
        pass

    return data

def plot(datas, outdir='./', outname='aho', verbose=0):

    # expand outdir for '~' (home directory)
    outdir = pathlib.Path(outdir).expanduser()

    # make directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        pass

    fig, ax = plt.subplots(2,1)
    fig.set_size_inches(6,8)
    for i, data in enumerate(datas):
        ax[0].plot(data.freq*1.e-9, data.powDBm, color=g_colors[i], linewidth=1)
        ax[1].plot(data.freq*1.e-9, data.amp   , color=g_colors[i], linewidth=1)
        npoints1MHz = np.sum(data.freq<data.freq[0]+100.e+3)
        nep = np.std(data.amp[0:npoints1MHz])
        if verbose>0: print(f'NET first {npoints1MHz} points / {data.npoints} = {nep} W')
        pass
    ax[0].set_xlabel('Frequency [GHz]')
    ax[0].set_ylabel('Power [dBm]')
    ax[0].grid()
    ax[1].set_xlabel('Frequency [GHz]')
    ax[1].set_ylabel('Power [W]')
    ax[1].grid()
    fig.tight_layout()

    fig.savefig(f'{outdir}/{outname}.pdf')

    return 0


if __name__=='__main__':
    #input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/2021-10-21/data'
    #input_files= ['plate8_10500.0MHz_RBW1000.0kHz_100times_1.0sec_20211021185051.dat']
    input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/scan3/2021-10-26/data'
    input_files= [
            'scan_FFT_20.0GHz_span1.00MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            'scan_FFT_20.0GHz_span2.50MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            'scan_FFT_20.0GHz_span5.00MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            #'scan_FFT_20.0GHz_span10.00MHz_rbw0.3kHz_1.0sec_10counts_1runs.dat',
            ]
    

    datas = [ read_data(f'{input_dir}/{infile}', verbose=2) for infile in input_files ]

    plot(datas, 'output/', 'aho', verbose=2)

    pass
