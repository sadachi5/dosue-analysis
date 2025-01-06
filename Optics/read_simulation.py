#!/bin/env python3
# read_simulation.py
import numpy as np
import lmfit
import csv
from matplotlib import pyplot as plt
from scipy.interpolate import griddata
from matplotlib import ticker, cm, colors
import function as func
from ellipsoid_mirror import *


#####################
# Data modification #
#####################

def read_CST(filepath, datatype='', verbose=0):
    data = {}
    keys = []
    n_skipline = 0
    if datatype == '3D_E-field':
        keys = ['x', 'y', 'z', 'ExRe', 'ExIm', 'EyRe', 'EyIm', 'EzRe', 'EzIm']
        data = { key: [] for key in keys }
        n_skipline = 2
        split = None
    elif datatype == '3D_twovalue':
        keys = ['x', 'y', 'z', 'val0', 'val1']
        data = { key: [] for key in keys }
        n_skipline = 2
        split = None
    elif datatype == 'farfield1D':
        keys = ['angle', 'val']
        data = { key: [] for key in keys }
        split = '\t'
    elif datatype == '1D':
        keys = ['x', 'y']
        data = { key: [] for key in keys }
        split = ' '
        n_skipline = 2
    else:
        print(f'Error!! There is no datatype for "{datatype}".')
        return -1
        pass
    if verbose > -1: print(f'input file = {filepath} (type={datatype})')
    if verbose > 0: print(f'keys = {keys}')

    with open(filepath) as f:
        fin = list( csv.reader(f, delimiter=' ', skipinitialspace=True) )
        for n, line in enumerate(fin):
            if split is not None:
                new_line = []
                for _l in line:
                    new_line+=_l.split(split)
                    pass
                line = new_line
                pass
            if n < n_skipline:
                if verbose > 1: print(f'skip line {n}')
                continue
            if len(line) == 0:
                continue
            if line[0][0] == '#':
                if verbose > 1: print(f'skip line starting with #')
                continue
            if verbose > 1: print('line:', line)
            for k, key in enumerate(keys):
                data[key].append( float( line[k] ) )
                if verbose > 3: print(f'append data [{key}] = {data[key][-1]}')
                pass
            pass
        pass

    # Change to array
    for key in keys:
        data[key] = np.array( data[key] )
        pass

    return data


def read_CST_complex_field(
    indir='input/CST/215GHz_v8-4-215GHz/215GHz',
    prefix='E-Field_',
    suffix='_W300_215GHz.dat',
    grid=( np.linspace(-300, 300, 600+1), np.linspace(-300, 300, 600+1) ),
    levels=np.power(10, np.linspace(-30, 30, 21)/10),
    logz=True, cmap='jet', colorbarformat='%.1e',
    doPlot=False,
    verbose=0 ):
    
    infiles = {
        'U-Re': prefix + 'U_Re' + suffix,
        'U-Im': prefix + 'U_Im' + suffix,
        'V-Re': prefix + 'V_Re' + suffix,
        'V-Im': prefix + 'V_Im' + suffix,
        'W-Re': prefix + 'W_Re' + suffix,
        'W-Im': prefix + 'W_Im' + suffix,
    }
    
    data_dict = {}
    for _k, _f in infiles.items():
        fig, ax, _data_grid = plotCST2D(
            f'{indir}/{_f}', xkey='x', ykey='y', newgrid=grid, levels=levels, title=_k,
            colorbarformat=colorbarformat, cmap=cmap, doPlot=doPlot, logz=logz, verbose=verbose )
        _data_grid['x'] = _data_grid['x']*mm
        _data_grid['y'] = _data_grid['y']*mm
        data_dict[_k] = _data_grid
        pass
    
    return data_dict


def get_index_data(data, indices):
    newdata = {}
    for _c in data.keys():
        if isinstance(data[_c], list) or isinstance(data[_c], np.ndarray):
            newdata[_c] = data[_c][indices]
        else:
            newdata[_c] = data[_c]
            pass
        pass
    return newdata


def getSlice(data, slicekey='x', slice_val=0., verbose=0):
    if verbose > 0:
        print(slicekey, slice_val)
        pass
    val = data[slicekey]
    diff = np.abs(val - slice_val)
    min_diff = min(diff)+0.01
    slice_index = np.where( diff < min_diff )[0]
    if verbose > 0:
        print(f'getSlice min_diff = {min_diff}')
        pass
    data_slice = get_index_data(data, slice_index)
    return data_slice


def slice_file(infile, indir, slicekey='x', slice_val=0.):
    _infile_path = f'{indir}/{infile}'
    _data = read_CST(_infile_path, '3D_twovalue' )
    _data_slice = getSlice(_data, slicekey=slicekey, slice_val=slice_val)
    return _data_slice


#######################
# Beamsize Estimation #
#######################

def gauss_fit(x, y, x_fit_range=None, xlim=None, title='', verbose=0):

    if x_fit_range is not None:
        _index = np.where( (x >x_fit_range[0] ) & (x < x_fit_range[1]) )
        _x = x[_index]
        _y = y[_index]
    else:
        _x = x
        _y = y
        pass

    model = lmfit.models.ConstantModel() +  lmfit.models.GaussianModel()

    params = model.make_params()
    params['c'].set(value=0., min=0., max=1000, vary=False)
    params['amplitude'].set(value=max(y), min=0., max=None, vary=True)
    params['sigma'].set(value=20., min=0., max=1000., vary=True)
    params['center'].set(value=x[np.where(y==max(y))][0], min=-1000., max=1000, vary=True)

    result = model.fit(x=_x, data=_y, weights=None, params=params, method='leastsq')

    if verbose > 0:
        print(params)
        print(result.fit_report())
        print(result.ci_report())
        pass

    fig = result.plot(data_kws={'markersize': 5})
    fig.suptitle(title)
    fig.set_figwidth(8)
    fig.set_figheight(4)
    axes = fig.get_axes()

    xlim = [min(x), max(x)] if xlim is None else xlim
    axes[1].scatter(x, y, lw=1, c='tab:blue', s=10)
    axes[0].set_xlim(xlim)
    axes[1].set_xlim(xlim)
    axes[0].grid(True)
    axes[1].grid(True)

    fig.tight_layout()

    return result


# Convert: gaussian sigma in power ( P = P0 exp(-x^2/2sigma^2) ) to  beam size w ( P = P0 exp(-x^2/w^2) )
def sigma2beamsize(sigma):
    return np.sqrt(2.) * sigma


def get_fwhm_from_half(x, y):
    x_half = func.getX(x, y, 0.5*max(y))
    n_x_half = len(x_half)
    if n_x_half != 2:
        print('Warning! There are several candidates for the x with y=0.5*ymax.')
        print(f'x candidates = {x_half}')
        print('Use two x closer to 0.')
        i_order = np.argsort( np.abs(x_half) )
        new_x_half = [0., 0.]
        new_x_half[0] = x_half[i_order[0]]
        new_x_half[1] = x_half[i_order[1]]
        x_half = new_x_half
        print(f'x_half = {x_half}')
        pass
    fwhm = np.abs(x_half[1] - x_half[0])
    return fwhm


# Fit (z, w) with a gaussian beam propagation
# unit: [m]
def fit_beamsize(z, w, wavelength, w0_init=None):
    def residual(pars, x, data=None, eps=None):
        # unpack parameters: extract .value attribute for each parameter
        parvals = pars.valuesdict()
        w0 = parvals['w0']
        z0 = parvals['z0']

        model = calc_beamwaist(z=x+z0, w0=w0, wavelength=wavelength)

        if data is None:
            return model
        if eps is None:
            return model - data
        return (model-data) / eps

    params = lmfit.Parameters()
    if w0_init is None: w0_init = 1.*mm
    params.add('w0', value=w0_init, vary=True, min=0.*mm, max=100*mm)
    params.add('z0', value=0.*mm, vary=True, min=-100.*mm, max=100*mm)
    result = lmfit.minimize(residual, params, args=(z, w))

    print(f'w0 = {result.params["w0"].value/mm} +- {result.params["w0"].stderr/mm} mm')
    print(f'z0 = {result.params["z0"].value/mm} +- {result.params["z0"].stderr/mm} mm')
    print(f'redchi = {result.redchi}')
    print(f'success = {result.success}')
    w0 = result.params['w0'].value
    z0 = result.params['z0'].value

    fig, ax = plt.subplots(1, 1, figsize=(12, 4))
    _z = np.arange(0, np.max(z)+50*mm, 1*mm)
    ax.scatter(z/mm, w/mm, label='Data')
    ax.plot(_z/mm, calc_beamwaist(z=_z+z0, w0=w0, wavelength=wavelength)/mm, label='Fit', c='tab:red')
    ax.legend(frameon='False')

    return result


###########
# Plotter #
###########


def plot2D_Eyz(data, logz=False, figsize0=8., colorbar_frac = 0.2, levels=None):
    figsize = (figsize0*(1.+colorbar_frac), figsize0)
    x_list = np.unique(data['x'])
    y_list = np.unique(data['y'])
    z_list = np.unique(data['z'])
    func.print_list(x_list)
    func.print_list(y_list)
    func.print_list(z_list)

    print(y_list.tolist())
    print(z_list.tolist())

    y_grid, z_grid = np.meshgrid(y_list, z_list)

    EyRe = griddata((data['y'], data['z']), data['EyRe'], (y_grid, z_grid))
    EzRe = griddata((data['y'], data['z']), data['EzRe'], (y_grid, z_grid))
    EyIm = griddata((data['y'], data['z']), data['EyIm'], (y_grid, z_grid))
    EzIm = griddata((data['y'], data['z']), data['EzIm'], (y_grid, z_grid))

    Ey_Mag = np.sqrt( np.power(EyRe, 2.) + np.power(EyIm, 2.) )
    Ez_Mag = np.sqrt( np.power(EzRe, 2.) + np.power(EzIm, 2.) )
    E_Mag = np.sqrt( np.power(Ey_Mag, 2.) + np.power(Ez_Mag, 2.) )

    fig, ax = plt.subplots(1,1,figsize=figsize)
    
    locator=ticker.LogLocator() if logz else ticker.MaxNLocator()
    cs = ax.contourf(y_grid, z_grid, E_Mag, levels=levels, cmap='jet', locator=locator)
    cbar = fig.colorbar(cs, ax=ax, format=colorbarFormat, fraction=colorbar_frac)
    cbar.ax.tick_params(labelsize=12)
    ax.set_aspect('equal')
    plt.gca().set_aspect('equal')
    fig.tight_layout()
    
    return fig, ax


def plot2D_key(data, keys=['y', 'z', 'val0'], figsize0=8., colorbar_frac = 0.2, logz=False, levels=None, newgrid=None, cmap='jet', alpha=1, 
              xlabel='', ylabel='', zlabel='', title='', colorbarformat='%.1e', grid=True, contourtype=None, nconv2D=None, doPlot=True, equal_aspect=True, verbose=0):
    if equal_aspect:
        figsize = (figsize0*(1.+colorbar_frac), figsize0)
    else:
        figsize = (figsize0*1.5, figsize0)
        pass
    x = data[keys[0]]
    y = data[keys[1]]
    if newgrid is None:
        x_list = np.unique(x)
        y_list = np.unique(y)
    else:
        x_list = newgrid[0]
        y_list = newgrid[1]
        pass
    func.print_list(x_list)
    func.print_list(y_list)

    x_grid, y_grid = np.meshgrid(x_list, y_list)

    z_grid = griddata((x, y), data[keys[2]], (x_grid, y_grid))
    #print('z_grid shape', z_grid.shape)
    if nconv2D is not None:
        from scipy import signal
        nconv2D_array = np.ones(nconv2D)/(nconv2D[0]*nconv2D[1])
        z_grid = signal.convolve2d(z_grid, nconv2D_array, 'same')
        pass

    fig, ax = None, None
    if doPlot:
        fig, ax = plt.subplots(1,1,figsize=figsize)
        fig.suptitle(title)

        locator=ticker.LogLocator() if logz else ticker.MaxNLocator()
        if contourtype is None:
            cs = ax.contourf(x_grid, y_grid, z_grid, levels=levels, cmap=cmap, alpha=alpha, locator=locator)
        elif contourtype=='line':
            cs = ax.contour(x_grid, y_grid, z_grid, levels=levels, cmap=cmap, alpha=alpha, locator=locator)
        elif contourtype=='fill-line':
            ax.contour(x_grid, y_grid, z_grid, levels=levels, colors='black', alpha=alpha, locator=locator)
            cs = ax.contourf(x_grid, y_grid, z_grid, levels=levels, cmap=cmap, locator=locator)
            pass
        cbar = fig.colorbar(cs, ax=ax, format=colorbarformat, fraction=colorbar_frac, label=zlabel)
        cbar.ax.tick_params(labelsize=12)
        ax.set_xlabel(xlabel)
        ax.set_ylabel(xlabel)
        ax.grid(grid)
        if equal_aspect:
            ax.set_aspect('equal')
            plt.gca().set_aspect('equal')
            pass
        fig.tight_layout()
        pass
    
    return fig, ax, {'x':x_grid.flatten(), 'y':y_grid.flatten(), 'z':z_grid.flatten()}


def plotCST2D(infile_path, xkey='z', ykey='y', newgrid=None, 
              levels=np.power(10, np.linspace(10, 25, 15)/10), cmap='jet', alpha=1, colorbarformat='%.1e', title='', logz=True, doPlot=True, equal_aspect=True, verbose=0 ):
    _raw_data = read_CST(f'{infile_path}', '3D_twovalue', verbose=0)
    #if newgrid is None:
    #    newgrid = ( np.linspace(-100, 500, 600+1), np.linspace(-500, 500, 500+1) )
    #    pass
    fig, ax, data_grid = plot2D_key(_raw_data, keys=[xkey, ykey, 'val0'], logz=logz,  # file column (x, y, z, val0) --> (x', y', z'(val0))
                     levels=levels, newgrid=newgrid, colorbarformat=colorbarformat, title=title, cmap=cmap, alpha=alpha, doPlot=doPlot, equal_aspect=equal_aspect, verbose=verbose)
    data_grid['z'][np.isnan(data_grid['z'])] = 0.
    return fig, ax, data_grid


######################
# Complex Beam Calculation #
######################

def calc_mag2(data_dict):
    mag2 = data_dict['U-Re']['z']**2. + data_dict['U-Im']['z']**2. \
            + data_dict['V-Re']['z']**2. + data_dict['V-Im']['z']**2. \
            + data_dict['W-Re']['z']**2. + data_dict['W-Im']['z']**2.
    return mag2

def integral2D(x, y, z, rmax = None, doPlot=False, verbose=0):
    x_uniq = np.unique(x)
    y_uniq = np.unique(y)
    if verbose > 1:
        print(x_uniq)
        print(y_uniq)
        pass
    dx = abs(x_uniq[1] - x_uniq[0])
    dy = abs(y_uniq[1] - y_uniq[0])
    if verbose > 0:
        print(f'dx = {dx}, dy = {dy}')
        pass
    
    if rmax is not None:
        in_range = ( np.power(x, 2.) + np.power(y, 2.) <= np.power(rmax, 2.) )
        x = x[in_range]
        y = y[in_range]
        z = z[in_range]
        pass
    if doPlot:
        plot2D_key({'x':x/mm, 'y':y/mm, 'z':z}, keys=['x', 'y', 'z'])
        pass
    
    integral = np.sum(z) * dx * dy
    return integral

# Output: integral(real)^2 + integral(imag)^2
def calc_complex_integral_mag2(real_griddata, imag_griddata, rmax=None, verbose=0):
    integral_real = integral2D(real_griddata['x'], real_griddata['y'], real_griddata['z'], rmax=rmax)
    integral_imag = integral2D(imag_griddata['x'], imag_griddata['y'], imag_griddata['z'], rmax=rmax)
    integral_mag2 = integral_real**2. + integral_imag**2.
    if verbose > 0:
        print(f'integral real (r<{rmax}) = {integral_real}')
        print(f'integral imag (r<{rmax}) = {integral_imag}')
        print(f'integral real^2 + integral imag^2 (r<{rmax}) = {integral_mag2}')
        pass
    return integral_mag2

