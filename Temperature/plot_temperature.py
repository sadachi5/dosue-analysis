#!/bin/env python3
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import utils

# set colorful lines
cmap = plt.get_cmap('jet')

def read_datafile(filename, 
        columns=['unittime', 'date', 'time', 'name1', 'temp1', 'raw1'],
        delimiter=' ', comment='#',
    ):
    df = pd.read_csv(filename, delimiter=delimiter, comment=comment)
    if len(columns) == len(df.columns):
        df.columns = columns
        pass
    return df
    

def plot(channels, data_configs, xcolumn='time', 
        logx=False, logy=False, 
        outdir='./figure', outname='aho'):
    utils.makedir(outdir)
    # check data_configs
    for i, _dc in enumerate(data_configs):
        if not ('label' in _dc.keys()):
            data_configs[i]['label']=f'{i}'
            pass
        if not ('lw' in _dc.keys()):
            data_configs[i]['lw']=2
            pass
        if not ('ls' in _dc.keys()):
            data_configs[i]['ls']='-'
            pass
        pass

    data_list = []
    for _dc in data_configs:
        filenames = utils.list_datafiles(
            starttime=_dc['start'], stoptime=_dc['stop'], 
            data_dir=_dc['dir'], add_dir=True)

        _df_list = []
        for _filename in filenames:
            print(_filename)
            _df = read_datafile(_filename, columns=_dc['columns'], delimiter=' ', comment='#')
            _df_list.append(_df)
            pass

        _df_all = pd.concat(_df_list, join='inner')

        # cut time range
        _start = utils.time2unixTime(_dc['start'])
        _stop = utils.time2unixTime(_dc['stop'])
        print(_start, _stop)
        print(_df_all['unixtime'])
        _select = (_df_all['unixtime'] >= _start) & \
                  (_df_all['unixtime'] <= _stop)
        print(_select)
        _df_all = _df_all[_select]
        time = (_df_all['unixtime'] - _start)
        _df_all['time'] = time
        print(_df_all.keys())

        data_list.append(_df_all)
        pass

    # initialize figure
    fig = plt.figure(figsize=(16, 4))
    fig.tight_layout()
    plt.rcParams["font.size"] = 16 # Nominal font size
    plt.subplots_adjust()
    ax = fig.add_subplot(1,1,1)

    for i, _df in enumerate(data_list):
        _dc = data_configs[i]
        for j, _c in enumerate(channels):
            _column = _c['column']
            _dc_label = _dc['label']
            _c_label = _c['label']
            _label = f"{_dc_label}: {_c_label}"
            _lw = _dc['lw']
            _ls = _dc['ls']
            ax.plot(_df[xcolumn], _df[_column], label=_label, lw=_lw, ls=_ls)
            pass
        pass
    if xcolumn == 'time':
        ax.set_xlabel('Time [sec]')
        pass
    ax.set_ylabel('Temperature [K]')
    ax.legend(fontsize=10).get_frame().set_alpha(0)
    if logx: ax.set_xscale('log')
    if logy: ax.set_yscale('log')

    fig.set_tight_layout(True)
    fig.savefig(f'{outdir}/{outname}.png')
    return True


if __name__=='__main__':

    default_data_columns = [
         'unixtime', 'date', 'time', 
         'name1', 'temp1', 'raw1',
         'name2', 'temp2', 'raw2',
         'name3', 'temp3', 'raw3',
         'name4', 'temp4', 'raw4',
         'name5', 'temp5', 'raw5',
         'name6', 'temp6', 'raw6',
         'name7', 'temp7', 'raw7',
         'name8', 'temp8', 'raw8',
         ]

    outname='temperature_compare_20220806'
    channels = [
        {'column':'temp1', 'label':'CH1 2nd stage'},
        {'column':'temp4', 'label':'CH4 2nd stage HEMT'},
        {'column':'temp2', 'label':'CH2 1st stage'},
        {'column':'temp3', 'label':'CH3 Wave-guide @ 40K'},
        ]
    data_configs = [
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/08/02 16:23', 'stop':'2022/08/03 12:00', 
            'label':'2nd cooling', 'ls':'--', 'lw':2},
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/08/06 10:07', 'stop':'2022/08/12', 
            'label':'3rd cooling', 'ls':'-', 'lw':2},
        ]

    '''
    outname='temperature_compare_20220802'
    channels = [
        {'column':'temp1', 'label':'CH1 2nd stage'},
        {'column':'temp3', 'label':'CH3 1st stage'},
        {'column':'temp2', 'label':'CH2 40K shield bottom'},
        {'column':'temp4', 'label':'CH4 40K shield top'},
        ]
    data_configs = [
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/07/30 08:24', 'stop':'2022/07/31 08:00', 
            'label':'1st cooling', 'ls':'--', 'lw':2},
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/08/02 16:23', 'stop':'2022/08/10', 
            'label':'2nd cooling', 'ls':'-', 'lw':2},
        ]
    '''

    plot(channels=channels, data_configs=data_configs, 
            logy=True,
            outdir='./figure', outname=outname)

    pass
