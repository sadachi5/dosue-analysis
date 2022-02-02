#!/bin/env python
import os, argparse
import pathlib

import numpy as np
np.set_printoptions(threshold=10)
from matplotlib import pyplot as plt

from read_data import read_data

g_colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:olive','tab:cyan','tab:gray','red','royalblue','turquoise','darkolivegr    een', 'magenta', 'blue', 'green']*5;

# Constants
c = 299792458. # [m/sec] speed of light from wikipedia
k_B = 1.380649e-23 # [J/K] boltzmann constant


class YFactor:
    def __init__(self, filepath_77K, filepath_300K, filepath_plate, rebin=1, verbose=0):
        self.data_77K   = read_data(filepath_77K, rebin=rebin)
        self.data_300K  = read_data(filepath_300K, rebin=rebin)
        self.data_plate = read_data(filepath_plate, rebin=rebin)
        self.datas  = [self.data_77K, self.data_300K, self.data_plate]
        self.labels = ['77K', '300K', 'Plate']
        self.verbose = verbose

        self.Trx_str = '$T_{\mathrm{rx}}$'
        self.Tsys_str = '$T_{\mathrm{sys}}$'
        self.Tplate_str = '$T_{\mathrm{plate}}$'
        pass

    def plot_raw(self):
        fig, ax = plt.subplots(2,1)
        fig.set_size_inches(6,8)
        for i, data in enumerate(self.datas):
            ax[0].plot(data.freq*1.e-9, data.powDBm, color=g_colors[i], linewidth=1, alpha=0.5, label=self.labels[i])
            ax[1].plot(data.freq*1.e-9, data.amp   , color=g_colors[i], linewidth=1, alpha=0.5, label=self.labels[i])
            pass
        ax[0].set_xlabel('Frequency [GHz]')
        ax[0].set_ylabel('Power [dBm]')
        ax[0].grid()
        ax[0].legend()
        ax[1].set_xlabel('Frequency [GHz]')
        ax[1].set_ylabel('Power [W]')
        ax[1].legend()
        ax[1].grid()
        fig.tight_layout()
        return fig, ax


    def calculate(self, doPlot=False):
        # y = a*x+b
        p1 = (77., self.data_77K.amp)
        p2 = (300., self.data_300K.amp)
        self.a = (p2[1]-p1[1])/(p2[0]-p1[0])
        self.b = p1[1]-self.a*p1[0]

        # Pout = gain*conv_P*(Tin + Trx)
        #      = a*T + b
        self.rbw = self.data_plate.rbw
        conv_P = k_B * self.rbw # P = k_B*T*delta_nu(RBW)
        self.gain = self.a/conv_P
        self.gainDB = np.log10(self.gain)*10.
        self.Trx    = self.b/(self.gain*conv_P)
        self.Tplate = (self.data_plate.amp - self.b)/(self.a)
        self.Tsys   = self.Trx + self.Tplate

        self.Psys   = k_B*self.Tsys*self.rbw
        self.NEP    = np.sqrt(2.)*self.Psys/np.sqrt(self.rbw)

        if doPlot:
            fig, ax = plt.subplots(3,3)
            fig.set_size_inches(12,12)

            _ax = ax[0][0]
            _ax.plot(self.data_plate.freq*1.e-9, self.a, color=g_colors[0], linewidth=1, alpha=1., label='a')
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel('a')
            _ax.grid()
            _ax.legend()

            _ax = ax[0][1]
            _ax.plot(self.data_plate.freq*1.e-9, self.gainDB, color=g_colors[0], linewidth=1, alpha=1., label='Gain [dB]')
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel('Gain [dB]')
            _ax.grid()
            _ax.legend()
            
            _ax = ax[0][2]
            _ax.plot(self.data_plate.freq*1.e-9, self.b, color=g_colors[0], linewidth=1, alpha=1., label='b')
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel('b (intercept)')
            _ax.grid()
            _ax.legend()

            _ax = ax[1][0]
            _ax.plot(self.data_plate.freq*1.e-9, self.Trx, color=g_colors[0], linewidth=1, alpha=1., label=self.Trx_str)
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel(self.Trx_str)
            _ax.grid()
            _ax.legend()

            _ax = ax[1][1]
            _ax.plot(self.data_plate.freq*1.e-9, self.Tplate, color=g_colors[0], linewidth=1, alpha=1., label=self.Tplate_str)
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel(self.Tplate_str)
            _ax.grid()
            _ax.legend()

            _ax = ax[1][2]
            _ax.plot(self.data_plate.freq*1.e-9, self.Tsys, color=g_colors[0], linewidth=1, alpha=1., label=self.Tsys_str)
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel(self.Tsys_str)
            _ax.grid()
            _ax.legend()

            _ax = ax[2][0]
            _ax.plot(self.data_plate.freq*1.e-9, self.Psys, color=g_colors[0], linewidth=1, alpha=1., label='$P_{sys}$')
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel('$P_{sys}$')
            _ax.grid()
            _ax.legend()

            _ax = ax[2][1]
            _ax.plot(self.data_plate.freq*1.e-9, self.NEP, color=g_colors[0], linewidth=1, alpha=1., label='$NEP$ [W/$\sqrt{Hz}$]')
            _ax.set_xlabel('Frequency [GHz]')
            _ax.set_ylabel('$NEP$ [W/$\sqrt{Hz}$]')
            _ax.grid()
            _ax.legend()



            fig.tight_layout()
            return fig, ax
        else :
            return 0


def main(input_dir, input_files, outdir='./', outname='aho', rebin=1, verbose=0):
    # expand outdir for '~' (home directory)
    outdir = pathlib.Path(outdir).expanduser()

    # make directory
    if not os.path.exists(outdir):
        os.makedirs(outdir)
        pass

    filepath_77K   = f'{input_dir}/{input_files[0]}'
    filepath_300K  = f'{input_dir}/{input_files[1]}'
    filepath_plate = f'{input_dir}/{input_files[2]}'
    yfactor = YFactor(filepath_77K, filepath_300K, filepath_plate, rebin=rebin)
    fig, ax = yfactor.plot_raw()
    fig.savefig(f'{outdir}/{outname}_raw.pdf')

    fig, ax = yfactor.calculate(doPlot=True)
    fig.savefig(f'{outdir}/{outname}_result.pdf')

    return 0


if __name__=='__main__':

    input_dir = '/Users/shadachi/Experiment/DOSUE/data/ms2840a/2021-10-21/data'
    input_files= [
            '77K8_10500.0MHz_RBW1000.0kHz_100times_1.0sec_20211021185011.dat', # 77K
            '300K8_10500.0MHz_RBW1000.0kHz_100times_1.0sec_20211021184923.dat', # 300K
            'plate8_10500.0MHz_RBW1000.0kHz_100times_1.0sec_20211021185051.dat', # plate
            ]

    main(input_dir, input_files, outdir='output', outname='yfactor', rebin=100, verbose=2)
    pass
