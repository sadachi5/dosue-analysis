#!/bin/env python3 
import os
import sys
from matplotlib import pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
from scipy.constants import c

mm = 1e-3 # m
cm = 1e-2 # m
GHz = 1e+9 # Hz


def calc_wavelength(freq):
    return c/freq

def calc_z_c(w, wavelength):
    return np.pi * w * w / wavelength

def calc_R(d, w, wavelength):
    z_c = calc_z_c(w, wavelength)
    return d + z_c**2 / d

# m base
# d1/z_c should be larger than w_out/w_in=m.
def calc_R1R2d2(
        freq = 200.*GHz, # frequency
        w_in=5*mm,       # input beam waist
        w_out=100*mm,    # output beam waist
        d1=1000*mm,       # input distance
        f0fix = False, # must not use (wrong calculation)
        verbose=1,
        ):

    wavelength = calc_wavelength(freq)
    m = w_out/w_in # magnification
    alpha = 1.-1./(m*m)
    z_c = np.pi*(w_in**2.)/wavelength
    f0 = np.pi * w_in * w_out /wavelength

    if verbose > 0:
        print(f'wavelength = {wavelength/mm:.1f} mm')
        print(f'w_in = {w_in/mm:.1f} mm')
        print(f'w_out = {w_out/mm:.1f} mm')
        print(f'm = {m} ')
        print(f'd1/z_c = {d1/z_c:.1f}')
        print(f'z_c = {z_c/mm:.1f} mm')
        print(f'f0 = {f0/mm:.1f} mm')
        pass

    if f0fix:
        d1 = f0
        d2 = f0
        d2_p = d2_m = d2
        print(f'd1 = {d1/cm:.1f} cm')
        print(f'd2 = {d2/cm:.1f} cm')
        pass
    
    # calculate R1
    R1 = calc_R(d1, w_in, wavelength)

    # calculate f
    if w_in == w_out :
        f = z_c * (1.+(d1/z_c)**2.)/(2.*d1/z_c)
        f_p = f_m = f
    else:
        f2 = np.sqrt( 1. - alpha*(1. + (z_c/d1)**2.) )
        f1_p = 1+f2
        f1_m = 1-f2
        f_p = z_c * d1/(z_c*alpha) * f1_p
        f_m = z_c * d1/(z_c*alpha) * f1_m
        pass

    # calculate R2
    R2_p = R1*f_p/(R1 - f_p)
    R2_m = R1*f_m/(R1 - f_m)

    if verbose > 0:
        print(f'R1 = {R1/cm:.1f} cm')
        print(f'd1 = {d1/cm:.1f} cm')
        print(f'f0 = {f0/cm:.1f} cm')
        print(f'f_p = {f_p/cm:.1f} cm')
        print(f'f_m = {f_m/cm:.1f} cm')
        print(f'R2_p = {R2_p/cm:.1f} cm')
        print(f'R2_m = {R2_m/cm:.1f} cm')
        pass

    # calculate d2
    if not f0fix:
        d2_p = f_p + m * np.sqrt(f_p**2. - f0**2.)
        d2_m = f_m - m * np.sqrt(f_m**2. - f0**2.)

        d2_R2p_p = ( R2_p + np.sqrt( R2_p**2. - 4.*(m**4.)*(z_c**2.) ) )/2.
        d2_R2m_p = ( R2_m + np.sqrt( R2_m**2. - 4.*(m**4.)*(z_c**2.) ) )/2.

        d2_R2p_m = ( R2_p - np.sqrt( R2_p**2. - 4.*(m**4.)*(z_c**2.) ) )/2.
        d2_R2m_m = ( R2_m - np.sqrt( R2_m**2. - 4.*(m**4.)*(z_c**2.) ) )/2.
 
        #_d2_pp = 1/2.*(R2_p + np.sqrt( R2_p**2 - 4.*m**4*z_c**2) )
        #_d2_pm = 1/2.*(R2_p - np.sqrt( R2_p**2 - 4.*m**4*z_c**2) )
        #_d2_mp = 1/2.*(R2_m + np.sqrt( R2_m**2 - 4.*m**4*z_c**2) )
        #_d2_mm = 1/2.*(R2_m - np.sqrt( R2_m**2 - 4.*m**4*z_c**2) )

        if verbose > 0:
            print(f'd2_p = {d2_p/cm:.1f} cm')
            print(f'd2_m = {d2_m/cm:.1f} cm')

            print(f'd2_R2p_p from R2_p = {d2_R2p_p/cm:.1f} cm')
            print(f'd2_R2m_p from R2_m = {d2_R2m_p/cm:.1f} cm')
            print(f'd2_R2p_m from R2_p = {d2_R2p_m/cm:.1f} cm')
            print(f'd2_R2m_m from R2_m = {d2_R2m_m/cm:.1f} cm')
 
            #print(f'_d2_pp = {_d2_pp/cm:.1f} cm')
            #print(f'_d2_pm = {_d2_pm/cm:.1f} cm')
            #print(f'_d2_mp = {_d2_mp/cm:.1f} cm')
            #print(f'_d2_mm = {_d2_mm/cm:.1f} cm')
            pass

        # calculate R2 from d2, w_out
        R2_from_d2p = calc_R(d2_p, w_out, wavelength)
        R2_from_d2m = calc_R(d2_m, w_out, wavelength)
        if verbose > 0:
            print(f'R2 calculated from d2_p = {R2_from_d2p/mm} mm')
            print(f'R2 calculated from d2_m = {R2_from_d2m/mm} mm')
            pass
        pass



    return (d1, d2_p, d2_p), (d1, d2_m, d2_m)

def scan(
        w_in_min=3*mm, w_in_max=4*mm, w_in_scan=1*mm,
        d1_min=62*cm, d1_max=200*cm, d1_scan=1*cm,
        w_out = 100*mm
        ):

    w_in_list = np.arange(w_in_min, w_in_max+w_in_scan, w_in_scan)
    d1_list = np.arange(d1_min, d1_max+d1_scan, d1_scan)

    X, Y = np.meshgrid(w_in_list, d1_list)
    d1 = calc_R1R2d2(w_in=X, d1=Y, w_out=w_out, verbose=0)[0][0]
    d2_p = calc_R1R2d2(w_in=X, d1=Y, w_out=w_out, verbose=0)[0][2]
    d2_m = calc_R1R2d2(w_in=X, d1=Y, w_out=w_out, verbose=0)[1][2]

    #for w_in, d1 in xy:
    #    ret = calc_R1R2d2(w_in=w_in, d1=d1, verbose=0)
    #    ret_pos = ret[0]
    #    R1 = ret_pos[0]
    #    R2_p = ret_pos[1]
    #    d2_p = ret_pos[2]
    #    pass

    Xcm = X/cm
    Ycm = Y/cm
    d1cm = d1/cm
    d2_p_cm = d2_p/cm
    d2_m_cm = d2_m/cm

    plt.pcolormesh(Xcm, Ycm, d2_p, cmap='hsv', vmin=0, vmax=100)
    plt.xlabel('$w_{\mathrm{in}}$ [cm]')
    plt.ylabel('$d_{1}$ [cm]')
    plt.title('$d2_p$ [m]')
    plt.colorbar(orientation='vertical')
    plt.savefig('d2_p.pdf')
    plt.close()

    plt.pcolor(Xcm, Ycm, d2_m, cmap='hsv', vmin=0, vmax=100)
    plt.xlabel('$w_{\mathrm{in}}$ [cm]')
    plt.ylabel('$d_{1}$ [cm]')
    plt.title('$d2_m$ [m]')
    plt.colorbar(orientation='vertical')
    plt.savefig('d2_m.pdf')
    plt.close()

    plt.plot(Y, d2_p, label = 'd2_p')
    plt.plot(Y, d2_m, label = 'd2_m')
    plt.ylim(0,100)
    plt.grid(True)
    plt.savefig('d2_1D.pdf')
    plt.close()
    return 0


def calc_R(z, w, wavelength):
    z_c = calc_z_c(w, wavelength)
    return z + z_c * z_c / z

def calc_waist(z, w0, wavelength):
    z_c = calc_z_c(w, wavelength)
    return w0 * np.sqrt( 1. + np.power(z/z_c, 2.) )




if __name__ == '__main__':

    # Old
    #calc_R1R2d2(w_in=3*mm, w_out=2.99*mm, d1=10*cm, freq=100*GHz, verbose=1)
    #calc_R1R2d2(w_in=5*mm, w_out=5*mm, d1=50*mm, freq=200*GHz, verbose=1)
    #calc_R1R2d2(w_in=5*mm, w_out=100*mm, d1=1.2*m, freq=200*GHz, verbose=1)
    #calc_R1R2d2(w_in=2*mm, w_out=100*mm, d1=0.43, freq=200*GHz, verbose=1)
    #scan()

    # 2024/04/01
    #print(beam_waist(w_in=2.7*mm, freq=200*GHz, focus=100))
    #print(focus_from_waist(w_in=2.7*mm, w_out=200*mm, freq=200*GHz))
    print('### w_out = 10 mm #########################################')
    calc_R1R2d2(w_in=5*mm, w_out=10*mm, d1=10*cm, verbose=1)
    print('### w_out = 100 mm #########################################')
    calc_R1R2d2(w_in=5*mm, w_out=100*mm, d1=200*cm, verbose=1)






