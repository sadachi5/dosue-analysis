import csv
import numpy as np
from matplotlib import pyplot as plt

# start_freq, stop_freq, npoints are only used in OneColumn type
def read_csv(filename, csvType='Anritsu', start_freq=None, stop_freq=None, npoints=None, verbose=0):
    
    freq = [] # frequency list [GHz]
    power = [] # power list  [mW]
    
    if verbose > 0:
        print(f'Read csv file: {filename}')
        pass
    f = open(filename, 'r');
    if csvType=='TwoColumn':
        fin = list( csv.reader(f, delimiter=' ') )
    else:
        fin = list(csv.reader(f))
    if verbose > 1: 
        print(fin)  #リストの中身を出力
        pass
    isData = False
    
    if csvType=='Anritsu': # Anritsu : NOTE: only for RMS detection
        
        start_freq = 0
        stop_freq = 0
        npoints = 0
        for line in fin:
            if len(line)==0 : continue
            first = line[0].strip()
            # Search for frequency range
            if first == 'Trace-A':
                start_freq = int(line[1])
                stop_freq  = int(line[2])
                continue
            # Search for npoints
            if first == 'RMS':
                npoints = int(line[1])
                continue
            # Search for data starting point (Anritsu: Wave Data)
            if first.startswith('Wave Data'):
                isData = True
                continue
            # Get data
            if isData:
                power.append(10 ** (float(line[0])*0.1)) # dBm --> mW
                pass
            pass
        freq = np.linspace(start_freq,stop_freq,npoints) * 1.e-9 # Hz --> GHz
            
    elif csvType=='Keysight' : # Keysight
        
        for line in fin:
            if len(line)==0 : continue
            first = line[0].strip()
            # Search for data starting point (Keysight: DATA)
            if verbose > 1:
                print(f'first = {first}')
                pass
            if first == 'DATA':
                if verbose > 0:
                    print(f'Find DATA!')
                    pass
                isData = True
                continue
            # Get data
            if isData:
                freq.append( float(line[0]) * 1.e-9 ) # Hz --> GHz
                power.append(10 ** (float(line[1])*0.1)) # dBm --> mW
                pass
            pass
        
    elif csvType=='TwoColumn' : # Hz, dBm
        
        for line in fin:
            if len(line)==0 : continue
            first = line[0].strip()
            #print(f'first = {first}')
            if first[0]=='#':
                # skip line
                continue
            # Get data
            freq.append( float(line[0]) * 1.e-9 ) # Hz --> GHz
            power.append(10 ** (float(line[1])*0.1)) # dBm --> mW
            pass
        
    elif csvType=='OneColumn' : # dBm
        
        for line in fin:
            if len(line)==0 : continue
            first = line[0].strip()
            #print(f'first = {first}')
            if first[0]=='#':
                # skip line
                continue
            # Get data
            power.append(10 ** (float(line[0])*0.1)) # dBm --> mW
            pass
        if (start_freq is None) or (stop_freq is None) or (npoints is None):
            print('Error! There is no arguments for frequency information (start_freq, stop_freq, npoints).')
            print('Error! Please specify them!')
            return None
        freq = np.linspace(start_freq,stop_freq,npoints) * 1.e-9 # Hz --> GHz
        
        pass
    
    return np.array(freq), np.array(power)
                
def freq_average(data, n_average=100):

    n_data = len(data)
    n_points = int(n_data/n_average)
    
    data_ave = []
    data_err = []
    
    for i in range(n_points):
        data_subset = data[i*n_average:(i+1)*n_average]
        average = np.median(data_subset)
        average_err = np.std(data_subset)/np.sqrt(n_average) #  = 1/N * sqrt( sum((y-mean)^2))  ( std = sqrt( sum((y-mean)^2) / N) )
        data_ave.append(average)
        data_err.append(average_err)
        pass
    
    return np.array(data_ave), np.array(data_err)

def plot_power(
    is_dBm, freq_list, power_list, label_list, color_list, 
    suffix='', outdir='aho', 
    power_err_list=None, 
    logy=False, show=True,
    markersize=0.5, linestyle='', linewidth=0.5,
    figsize=(6.4, 4.8),
    xmin=None, xmax=None, ymin=None, ymax=None):
    # dBm or mW
    unit = 'dBm' if is_dBm else 'mW'
    
    if is_dBm:
        # Nominal value
        new_power_list = [ mW2dBm(power) for power in power_list ]
        power_list = new_power_list
        # Error (up, down)
        new_power_err_list = None
        if power_err_list is not None:
            new_power_err_list = []
            for power, err in zip(power_list, power_err_list):
                nom  = mW2dBm(power)
                up   = mW2dBm(power + err)
                down = mW2dBm(power - err)
                up_err = up - nom
                down_err = nom - down
                new_power_err_list.append([up_err,down_err])
                pass
            power_err_list = new_power_err_list
            pass
        pass
    # Plot
    plt.figure(figsize=figsize)
    if power_err_list is not None:
        # With error bar
        for freq, power, err, label, color in zip(freq_list, power_list, power_err_list, label_list, color_list):
            plt.errorbar(
                freq, power, yerr=err, label=label, 
                color=color, marker='o', markersize=markersize, linestyle=linestyle, linewidth=linewidth, 
                ecolor=color, capsize=3, fmt='o')
            pass
    else:
        # Without error bar
        for freq, power, label, color in zip(freq_list, power_list, label_list, color_list):
            plt.plot(freq, power, label=label, color=color, marker='o', markersize=markersize, linestyle=linestyle, linewidth=linewidth)
            pass
    plt.xlabel('Frequency [GHz]',fontsize=16) #x軸の名前
    plt.ylabel(f'Power [{unit}]',fontsize=16) #y軸の名前
    if (not is_dBm) and (not logy): plt.ylim(ymin=0.)
    if xmin is not None: plt.xlim(left=xmin)
    if xmax is not None: plt.xlim(right=xmax)
    if ymin is not None: plt.ylim(ymin=ymin)
    if ymax is not None: plt.ylim(ymax=ymax)
    plt.xticks(fontsize=12)
    plt.yticks(fontsize=12)
    plt.grid(True)
    plt.title(f'Power[{unit}] {suffix}',fontsize=16)
    if logy:
        plt.yscale('log')
        pass
    plt.legend()
    savefile = f'{outdir}/power_{unit}_{suffix}.png'
    print(f'save to {savefile}')
    plt.savefig(savefile)
    if show:
        plt.show()
        pass
    plt.close()

def mW2dBm(mW):
    return to_dB(mW)

def to_dB(G):
    if G is None:
        return None
    else:
        return np.log10(G)*10.
