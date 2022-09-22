#!/bin/env python3 
import os
import sys
from matplotlib import pyplot as plt
import numpy as np
from scipy.constants import c

mm = 1e-3 # m
cm = 1e-2 # m
GHz = 1e+9 # Hz

# m base
def calc_R1R2d2(
        freq = 200.*GHz, # frequency
        w_in=5*mm,       # input beam waist
        w_out=150*mm,    # output beam waist
        d1=500*mm,       # input distance
        verbose=1,
        ):

    wavelength = c/freq 
    m = w_out/w_in # magnitude
    z_c = np.pi*(w_in**2.)/wavelength

    # calculate R1
    R1 = d1 + z_c**2./d1

    # calculate f
    alpha = 1.-1./m
    if m==1. :
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
    R2_p = R1*f_p/(f_p - R1)
    R2_m = R1*f_m/(f_m - R1)

    # calculate d2
    f0 = np.pi * w_in * w_out /wavelength
    if verbose > 0:
        print(f'f0 = {f0}')
        print(f'f_p = {f_p}')
        print(f'f_m = {f_m}')
        pass
    d2_p = f_p + m * np.sqrt(f_p**2. - f0**2.)
    d2_m = f_m - m * np.sqrt(f_m**2. - f0**2.)

    if verbose > 0:
        print(f'R1 = {R1/mm} mm')
        print(f'R2_p = {R2_p/mm} mm')
        print(f'R2_m = {R2_m/mm} mm')
        print(f'd2_p = {d2_p/mm} mm')
        print(f'd2_m = {d2_m/mm} mm')
        pass

    return (R1, R2_p, d2_p), (R1, R2_m, d2_m)

def scan(
        w_in_min=3*mm, w_in_max=10*mm, w_in_scan=1*mm,
        d1_min=50*cm, d1_max=200*cm, d1_scan=1*cm,
        ):

    w_in_list = np.arange(w_in_min, w_in_max+w_in_scan, w_in_scan)
    d1_list = np.arange(d1_min, d1_max+d1_scan, d1_scan)

    X, Y = np.meshgrid(w_in_list, d1_list)
    R1 = calc_R1R2d2(w_in=X, d1=Y, verbose=0)[0][0]
    R2_p = calc_R1R2d2(w_in=X, d1=Y, verbose=0)[0][1]
    d2_p = calc_R1R2d2(w_in=X, d1=Y, verbose=0)[0][2]

    #for w_in, d1 in xy:
    #    ret = calc_R1R2d2(w_in=w_in, d1=d1, verbose=0)
    #    ret_pos = ret[0]
    #    R1 = ret_pos[0]
    #    R2_p = ret_pos[1]
    #    d2_p = ret_pos[2]
    #    pass

    Xcm = X/cm
    Ycm = Y/cm
    R1cm = R1/cm
    R2_p_cm = R2_p/cm
    d2_p_cm = d2_p/cm

    plt.pcolor(Xcm, Ycm, R1cm, cmap='hsv')
    plt.xlabel('$w_{\mathrm{in}}$ [cm]')
    plt.ylabel('$d_{1}$ [cm]')
    plt.title('$R1$ [cm]')
    plt.colorbar(orientation='vertical')
    plt.savefig('R1.pdf')
    plt.close()

    plt.pcolor(Xcm, Ycm, R2_p_cm, cmap='hsv')
    plt.xlabel('$w_{\mathrm{in}}$ [cm]')
    plt.ylabel('$d_{1}$ [cm]')
    plt.title('$R2$ [cm]')
    plt.colorbar(orientation='vertical')
    plt.savefig('R2_p.pdf')
    plt.close()

    plt.pcolor(Xcm, Ycm, d2_p_cm, cmap='hsv')
    plt.xlabel('$w_{\mathrm{in}}$ [cm]')
    plt.ylabel('$d_{1}$ [cm]')
    plt.title('$d2$ [cm]')
    plt.colorbar(orientation='vertical')
    plt.savefig('d2_p.pdf')
    plt.close()
    return 0

if __name__ == '__main__':
    #calc_R1R2d2(w_in=8*mm, d1=195*cm, verbose=1)
    calc_R1R2d2(w_in=3*mm, w_out=2.99*mm, d1=10*cm, freq=100*GHz, verbose=1)
    #scan()






