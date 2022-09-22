#!/bin/env python3
import numpy as np
import argparse
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
        logx=False, logy=False, xtype='sec', 
        figsize = (16, 4), grid=True, 
        ymin=None, ymax=None,
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
    fig = plt.figure(figsize=figsize)
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
            _x = _df[xcolumn]
            if   xtype=='hour': _x = _x/3600.
            elif xtype=='min': _x = _x/60.
            ax.plot(_x, _df[_column], label=_label, lw=_lw, ls=_ls)
            pass
        pass
    if xcolumn == 'time':
        if   xtype=='hour': ax.set_xlabel('Time [hour]')
        elif xtype=='min' : ax.set_xlabel('Time [min]')
        else:               ax.set_xlabel('Time [sec]')
        pass
    ax.set_ylabel('Temperature [K]')
    ax.legend(fontsize=10).get_frame().set_alpha(0)
    ax.grid(grid)
    if logx: ax.set_xscale('log')
    if logy: ax.set_yscale('log')
    if ymin is not None: ax.set_ylim(ymin=ymin)
    if ymax is not None: ax.set_ylim(ymax=ymax)

    fig.set_tight_layout(True)
    fig.savefig(f'{outdir}/{outname}.png')
    return True


if __name__=='__main__':

    xcolumn = 'time'
    xtype = 'sec'
    figsize = (16, 4)
    grid = True
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

    '''
    outname='temperature_single'
    xtype = 'hour'
    figsize = (12,4)
    channels = [
        {'column':'temp1', 'label':'4K stage'},
        ]
    data_configs = [
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/08/15 16:45', 'stop':'2022/08/16 10:00', 
            'label':'', 'ls':'-', 'lw':2},
        ]
    '''

    outname='temperature_compare'
    channels = [
        {'column':'temp1', 'label':'CH1 2nd stage'},
        {'column':'temp2', 'label':'CH2 1st stage'},
        {'column':'temp3', 'label':'CH3 Wave-guide @ 40K'},
        {'column':'temp4', 'label':'CH4 40K shield top'},
        ]
    data_configs = [
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/08/15 16:45', 'stop':'2022/08/17', 
            'label':'4th cooling (Wt 4K MLI)', 'ls':'--', 'lw':2},
        {'dir':'/data/lakeshore218', 'columns':default_data_columns,
            'start':'2022/08/24 08:48', 'stop':'2022/08/31 12:00', 
            'label':'5th cooling (No 4K MLI)', 'ls':'-', 'lw':2},
        ]

    '''
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

    parser = argparse.ArgumentParser()
    parser.add_argument('--channel_columns', nargs='*', default=None, 
                            help=f'column names for channels (ex: "temp1 temp2 temp3 temp4")')
    parser.add_argument('--channel_labels', nargs='*', default=None,
                            help=f'labels for channels (ex: "CH1 CH2 CH3 CH4")')
    parser.add_argument('--data_dir', default='/data/lakeshore218',
                            help=f'temprature data directory (default: "/data/lakeshore218")')
    parser.add_argument('--data_starts', nargs='*', default=None,
                            help=f'start datetime list (ex: "\'2022/08/15 16:15\' \'\2022/08/24 08:48\'")')
    parser.add_argument('--data_stops', nargs='*', default=None,
                            help=f'stop datetime list (ex: "\'2022/08/17 12:00\' \'\2022/08/26 12:00\'")')
    parser.add_argument('--data_labels', nargs='*', default=None,
                            help=f'label list (ex: "\'4th cooling (Wt 4K MLI)\' \'5th cooling (No 4K MLI)\'")')
    parser.add_argument('--data_linestyles', nargs='*', default=["--", "-"],
                            choices=['', 'solid', 'dashed', 'dashdot', 'dotted'],
                            help=f'line-style list (default: "\'--\' \'-\'")')
    parser.add_argument('--data_linewidth', default=2, type=float, help=f'Line width (default: 2)')
    parser.add_argument('--xcolumn', default=xcolumn, help=f'column name of time data (default: {xcolumn})')
    parser.add_argument('--xtype', default=xtype, choices=['sec', 'min', 'hour'], help=f'time type (default: {xtype})')
    parser.add_argument('--grid', action='store_true', default=False if grid else True, help=f'show grid (default: {grid})')
    parser.add_argument('--logx', action='store_true', default=False, help=f'log x-axis (default: False)')
    parser.add_argument('--logy', action='store_true', default=False, help=f'log y-axis (default: False)')
    parser.add_argument('--ymin', default=None, type=float, help=f'Y-axis min (default: None)')
    parser.add_argument('--ymax', default=None, type=float, help=f'Y-axis max (default: None)')
    parser.add_argument('--figsizeW', default=figsize[0], type=float, help=f'Width of figure (default: {figsize[0]})')
    parser.add_argument('--figsizeH', default=figsize[1], type=float, help=f'Height of figure (default: {figsize[1]})')
    parser.add_argument('--outname', default=outname, help=f'output filename (default: {outname})')
    parser.add_argument('-v', '--verbose', dest='verbose', default=0, type=int, help=f'Print out verbosity (default: 0)')
    args = parser.parse_args()

    if args.channel_columns is not None\
        and args.channel_labels is not None:
        channels = []
        for _column, _label in zip(args.channel_columns, args.channel_labels):
            channels.append({'column':_column, 'label':_label})
            pass
        pass

    if args.data_starts is not None\
        and args.data_stops is not None\
        and args.data_labels is not None:
        data_configs = []
        for _start, _stop, _label, _ls \
            in zip(args.data_starts, args.data_stops, 
                    args.data_labels, args.data_linestyles):
            data_configs.append({
                'dir':args.data_dir,
                'columns':default_data_columns,
                'start':_start, 'stop':_stop,
                'label':_label, 'ls':_ls,
                'lw':args.data_linewidth
                })
            pass
        pass

    xcolumn = args.xcolumn
    xtype = args.xtype
    grid = args.grid
    logx = args.logx
    logy = args.logy
    ymin = args.ymin
    ymax = args.ymax
    figsize = (args.figsizeW, args.figsizeH)
    outname = args.outname

    plot(channels, data_configs, xcolumn=xcolumn, 
        logx=logx, logy=logy, xtype=xtype, 
        ymin=ymin, ymax=ymax,
        figsize=figsize, grid=grid, 
        outdir='./figure', outname=outname);
    pass
