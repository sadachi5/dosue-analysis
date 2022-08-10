#!/bin/env python3
import os
import sys
from datetime import datetime as datetime

# Input types
# - time0 = '2022/08/01' --> set to 00:00
# - time0 = '2022/08/01 00:00'
# - time0 = datetime(2022, 8, 1, 0, 0)
def time2unixTime(time0='2022/08/01 00:00'):
    if isinstance(time0, str):
        split_time0 = time0.split(' ')
        if len(split_time0) == 2:
            split_HHMMSS = split_time0[1].split(':')
            if len(split_HHMMSS) == 1:
                time0 = datetime.strptime(time0, '%Y/%m/%d %H')
            elif len(split_HHMMSS) == 2:
                time0 = datetime.strptime(time0, '%Y/%m/%d %H:%M')
            elif len(split_HHMMSS) == 3:
                time0 = datetime.strptime(time0, '%Y/%m/%d %H:%M:%S')
                pass
            pass
        elif len(split_time0) == 1:
            time0 = datetime.strptime(time0, '%Y/%m/%d')
            pass
        pass

    if isinstance(time0, int) or isinstance(time0, float):
        print(f'Warning!: The input time looks unix-time. (time0={time0})')
        print(f'          --> return the input time!')
        return time0

    if isinstance(time0, datetime):
        return datetime.timestamp(time0)
    else:
        print(f"Error!: (time0={time0})"\
               "        The input time did not match to any expected types: \n"\
               "         - time0 = '2022/08/01' \n"\
               "         - time0 = '2022/08/01 00:00' \n"\
               "         - time0 = '2022/08/01 00:00:00' \n"\
               "         - time0 = datetime(2022, 8, 1, 0, 0)")
        return False


# input filename should be 'data_YYYY-MM-DD.dat'
def list_datafiles(starttime='2022/07/28 00:00', stoptime='2023/01/01 00:00', 
        data_dir='/data/lakeshore218', add_dir=True):
    start_unixtime = time2unixTime(starttime)
    stop_unixtime = time2unixTime(stoptime)
    filenames = os.listdir(data_dir)
    select_filenames = []
    for _filename in filenames:
        _file_date = _filename.split('_')[1].split('.')[0]
        _file_start = datetime.timestamp(datetime.strptime(_file_date, '%Y-%m-%d'))
        _file_stop = _file_start + 60*60*24
        if start_unixtime <= _file_start and _file_start <= stop_unixtime:
            # start edge file
            select_filenames.append(_filename)
        elif start_unixtime <= _file_start and _file_stop <= stop_unixtime:
            # in-range file
            select_filenames.append(_filename)
        elif start_unixtime <= _file_stop and _file_stop <= stop_unixtime:
            # stop edge file
            select_filenames.append(_filename)
            pass
        pass

    # sort
    select_filenames = sorted(select_filenames)

    # add directory path
    if add_dir:
        select_filenames = [ f'{data_dir}/{_filename}' for _filename in select_filenames ] 
        pass

    return select_filenames


def makedir(dirname='./aho'):
    if os.path.isdir(dirname):
        print(f'Warning! Directory of "{dirname}" already exists!')
        print(f'         --> No action!')
        return True
    else:
        print(f'Making directory: {dirname}')
        os.makedirs(dirname)
        return True


if __name__=='__main__':
    print(time2unixTime('2022/08/01 00:00'))
    print(time2unixTime('2022/08/01 00:00:01'))
    print(time2unixTime('2022/08/01 01:00'))
    print(time2unixTime(datetime(2022, 8, 1, 0, 0)))

    print(list_datafiles(starttime='2022/07/28 00:00', stoptime='2023/01/01 00:00', 
        data_dir='/data/lakeshore218', add_dir=True))
    pass
