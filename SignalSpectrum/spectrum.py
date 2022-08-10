#!/bin/python
import argparse;

import numpy as np;
np.set_printoptions(threshold=10);
from scipy import special;
from matplotlib import pyplot as plt;

from minuitfit import *;
from lmfit.models import GaussianModel;

colors = ['tab:blue','tab:orange','tab:green','tab:red','tab:purple','tab:brown','tab:pink','tab:olive','tab:cyan','tab:gray','red','royalblue','turquoise','darkolivegr    een', 'magenta', 'blue', 'green']*5;

# Constants
v_c = 220.e+3 # [m/sec] speed of solar system
v_E = v_c # [m/sec] speed of earth
c = 299792458. # [m/sec] speed of light from wikipedia
k_B = 1.380649e-23 # [J/K] boltzmann constant


def cummulative_velocity(v):
    # C * (exp_p - exp_m) + 1/2(erf_p + erf_m)
    C = v_c/(2.*np.sqrt(np.pi)*v_E) 
    exp_p = np.exp( -1. * np.power((v+v_E)/v_c, 2.) );
    exp_m = np.exp( -1. * np.power((v-v_E)/v_c, 2.) );
    erf_p = special.erf((v+v_E)/v_c);
    erf_m = special.erf((v-v_E)/v_c);

    #print(C, exp_p-exp_m, erf_p, erf_m);
    f = C*(exp_p-exp_m) + 1./2. * (erf_p + erf_m);
    return f;

def freq_to_nu(freq, freq_0):
    if np.isscalar(freq): 
        v = c * np.sqrt( 1. - np.power(freq_0/freq, 2.)) if freq>freq_0 else 0.;
    else : # list or array
        ok = (freq>freq_0)
        v  = np.full(len(freq), 0.);
        v[ok] = c * np.sqrt( 1. - np.power(freq_0/freq[ok], 2.));
        pass;
    return v;

def integral_binwidth_velocity(freq, freq_0, binwidth):
    v_0 = freq_to_nu(freq_0,freq_0);
    v_p = np.where(freq+binwidth/2.>freq_0, freq_to_nu(freq+binwidth/2., freq_0), v_0);
    v_m = np.where(freq-binwidth/2.>freq_0, freq_to_nu(freq-binwidth/2., freq_0), v_0);
    integral = cummulative_velocity(v_p) - cummulative_velocity(v_m);
    return integral;

def peak_spectrum(freq, freq_0, P_DP, b, a=0., binwidth=5.e+3):
    integral = integral_binwidth_velocity(freq, freq_0, binwidth);
    peak = P_DP * integral;
    power = peak + a*(freq-freq_0) + b;
    return power;

# P_DP [W] / A_eff [m^2] / rho_CDM (density) [GeV/cm^3]
def calc_chi(P_DP, A_eff, rho_CDM=0.39, alpha=1./np.sqrt(3.)):
    chi = 4.5e-14 * np.sqrt( P_DP/1.e-23  * 1./A_eff * 0.3/rho_CDM ) * np.sqrt(2./3.)/alpha ;
    return chi;

def calc_P_DP(chi, A_eff, rho_CDM=0.39, alpha=1./np.sqrt(3.)): 
    C    = A_eff * rho_CDM/0.3 * alpha*alpha/(2./3.);
    P_DP = 1.e+5/(4.5*4.5) * np.power(chi, 2.) * C;
    return P_DP;

# T [K] / d_nu [Hz]
def calc_noise(T, d_nu, t):
    return np.sqrt(2.) * k_B*T*d_nu / np.sqrt(d_nu * t); # J*Hz = J/sec = W

def calc_rebin(rebin, y, iserr=False):
    if rebin>1:
        new_nbin = (int)(len(y)/rebin)
        # ignore several bins in the end of array, which are the remainders after dividing by 'rebin'
        yEachRebins = y[:new_nbin*rebin]; 
        yEachRebins = yEachRebins.reshape(new_nbin, rebin);
        if not iserr :
            yRebin = np.mean(yEachRebins, axis=1);
        else :
            #print(f'rebin = {rebin}');
            #print(f'yEachRebins = {yEachRebins}');
            yRebin = np.sqrt(np.sum(yEachRebins*yEachRebins, axis=1))/rebin;# sum{(err_i)^2}/N
            #print(f'yRebin = {yRebin}');
            pass;
        return yRebin;
    else:
        return y;

def create_spectrum(
        chi  = 2.e-10, # kumodes = 2.e-10
        T_noise = 10, # [K] System temperature
        gain = 1.e+6, # 1e+6 = +60dB
        time = 60., # [sec] 1min
        freq_0 = 20.e+9, # [Hz] peak position
        freq_binwidth = 5.e+3, # [Hz] frequency bin width
        rebin = 1,
        A_eff= (6.*1.e-2/2.)**2.*np.pi, # [m^2]
        doPlot=True, outdir='', outname='spectrum.pdf', verbosity=1, freq_half_scale = 1e-6*10.):

    # Prepare frequency array
    nfreq_half = (int)((freq_0*freq_half_scale)/freq_binwidth); # determine the fit range
    nfreq      = nfreq_half * 2 + 1; # +1 is center bin
    freq_min_center  = freq_0 - nfreq_half*freq_binwidth; # determine the fit range
    freq_max_center  = freq_0 + nfreq_half*freq_binwidth; # determine the fit range
    freq_center = np.arange(freq_min_center, freq_max_center+freq_binwidth, freq_binwidth);
    if verbosity>1 :
        print(f'nfreq = {nfreq}, nfreq_half={nfreq_half}');
        print(f'freq_min_center = {freq_min_center}');
        print(f'freq_max_center = {freq_max_center}');
        print(f'freq_center = {freq_center}');
        pass;
    freq = np.arange(freq_min_center-freq_binwidth/2., freq_max_center+freq_binwidth/2.+freq_binwidth, freq_binwidth);
    # Hz --> GHz
    freq_conv = 1.e-9;
    freq           *= freq_conv;
    freq_center    *= freq_conv;
    freq_0         *= freq_conv;
    freq_min_center *= freq_conv;
    freq_max_center *= freq_conv;
    freq_binwidth   *= freq_conv;

    # Calculate P_DP
    P_DP = gain * calc_P_DP(chi, A_eff);
    d_nu = freq_binwidth/freq_conv;
    noise   = calc_noise(T_noise, d_nu=d_nu, t=time) * gain; # noise fluctuation propotional to T_noise after amplifications
    noise_floor = k_B * T_noise * d_nu * gain; # measure in spectrum analyzer after amplifications
    if verbosity>0 :
        print(f'chi     = {chi}');
        print(f'A_eff   = {A_eff*1.e+4} [cm^2]');
        print( 'gain    = {:+f} [dB]'.format(np.log10(gain)*10.));
        print(f'T_noise = {T_noise} [K]');
        print(f'd_nu    = {freq_binwidth/freq_conv*1.e-3} [kHz]');
        print(f'P_DP    = {P_DP} [W] (After *gain)');
        print(f'Noise   = {noise} [W] (After *gain)');
        print(f'Noise floor  = {noise_floor} [W] (After *gain)');
        pass;

    #v_noise = np.random.normal(0., np.sqrt(noise), nfreq);
    #white_noise = v_noise*v_noise;
    white_noise = np.random.normal(0., noise, nfreq);
    v_noise = np.sqrt(abs(white_noise));
    integral_velocity = integral_binwidth_velocity(freq_center/freq_conv, freq_0/freq_conv, binwidth=freq_binwidth/freq_conv);

    if verbosity>1: print(f'integral_velocity = {integral_velocity}');
    y = peak_spectrum(freq_center/freq_conv, freq_0=freq_0/freq_conv, P_DP=P_DP, a=0., b=noise_floor, binwidth=freq_binwidth/freq_conv) + white_noise;
    if verbosity>0:
        print('freq_center', freq_center);
        print('power', y);
        pass;

    # Create error array from std of white_noise
    std_white_noise = np.std(white_noise);
    y_err = np.full(nfreq, std_white_noise);

    # Rebinning
    if rebin>1:
        # freq_center
        freq_center = calc_rebin(rebin, freq_center, iserr=False)
        # whire_noise
        white_noise = calc_rebin(rebin, white_noise, iserr=False)
        # v_noise
        v_noise = calc_rebin(rebin, v_noise, iserr=False)
        # integral_velocity
        integral_velocity = calc_rebin(rebin, integral_velocity, iserr=False)
        # y
        y = calc_rebin(rebin, y, iserr=False)
        # y_err
        y_err = calc_rebin(rebin, y_err, iserr=True)
        pass;

    if doPlot :
        # plot
        fig, axs = plt.subplots(2,2);
        fig.set_size_inches(12,9);
        plt.subplots_adjust(wspace=0.3, hspace=0.5, left=0.10, right=0.95,bottom=0.15, top=0.95)
 
        # plot spectrum
        ax = axs[0][0];
        ax.plot(freq_center, y, color='tab:blue'              , label='signal+white noise', marker='o', markersize=1., linestyle='-', linewidth=1.);
        ax.plot(freq_center, noise_floor+white_noise, color='tab:red', label='white noise', marker='o', markersize=1., linestyle='-', linewidth=1.);
        ax.set_xscale('linear');
        ax.set_yscale('linear');
        ax.set_xlabel('Frequency [GHz]', fontsize=15);
        ax.set_ylabel(f'Power [W / {freq_binwidth*1.e+6:.0f} kHz]', fontsize=15);
        ax.grid(True);
        ax.legend();
        ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
        ax.tick_params(axis='y', labelsize=12);
        ax.minorticks_on();
        #ax.set_xlim([freq[0],freq[-1]]);
        #ax.set_ylim([min(y)/2., max(y)*1.1]);
        ax.xaxis.get_major_formatter().set_useOffset(False);
        ax.yaxis.get_major_formatter().set_useOffset(False);
        ax.locator_params(axis='x', nbins=5);
        ax.locator_params(axis='y', nbins=5);
        xlim  = ax.get_xlim();
        xwidth= xlim[1]-xlim[0];
        x0    = xlim[0]+xwidth*0.05;
        ylim  = ax.get_ylim();
        ywidth= ylim[1]-ylim[0];
        y0    = ylim[1]-ywidth*0.05;
        yscale=0.08
        ax.text(x0, y0-ywidth*yscale*1 , fontsize=10, s='$A_{eff}$'+f' = {A_eff*1.e+4:.1f} [cm^2]');
        ax.text(x0, y0-ywidth*yscale*2 , fontsize=10, s='Gain = {:+.0f} [dB]'.format(np.log10(gain)*10.));
        ax.text(x0, y0-ywidth*yscale*3 , fontsize=10, s='$T_{noise}$'+f' = {T_noise:.1f} [K]');
        ax.text(x0, y0-ywidth*yscale*4 , fontsize=10, s=r'$\Delta_{\nu}$'+f' = {freq_binwidth/freq_conv*1.e-3:.1f} [kHz]');
        ax.text(x0, y0-ywidth*yscale*5 , fontsize=10, s='$P_{DP}$ after *gain'+f' = {P_DP:.1e} [W]');
        ax.text(x0, y0-ywidth*yscale*6.5, fontsize=8 , s='White noise ($\Delta N$) after *gain\n'+f'= {noise:.1e} [W]');
        ax.text(x0, y0-ywidth*yscale*8.0, fontsize=8 , s=f'Noise floor (offset) after *gain\n= {noise_floor:.1e} [W] ');
 
        # plot white noise histogram
        ax = axs[0][1];
        abs_max = max( abs(min(white_noise)), max(white_noise));
        ax.hist(white_noise, bins=20, range=[-abs_max, abs_max], histtype='stepfilled',
                align='mid', orientation='vertical', log=False, linewidth=0.5, linestyle='-', edgecolor='k',
                color='tab:blue', alpha=1.0, label='White noise fluctuation', stacked=False);
        ax.set_xscale('linear');
        ax.set_yscale('linear');
        ax.set_xlabel(f'Noise Power [W / {freq_binwidth*1.e+6:.0f} kHz]', fontsize=15);
        ax.set_ylabel(f'# of points', fontsize=15);
        ax.grid(True);
        ax.legend();
        ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
        ax.tick_params(axis='y', labelsize=12);
        ax.minorticks_on();
        ax.xaxis.get_major_formatter().set_useOffset(False);
        ax.yaxis.get_major_formatter().set_useOffset(False);
        ax.locator_params(axis='x', nbins=5);
        ax.locator_params(axis='y', nbins=5);
        xlim  = ax.get_xlim();
        xwidth= xlim[1]-xlim[0];
        x0    = xlim[0]+xwidth*0.05;
        ylim  = ax.get_ylim();
        ywidth= ylim[1]-ylim[0];
        y0    = ylim[1]-ywidth*0.05;
        yscale=0.08
        ax.text(x0, y0-ywidth*yscale*1 , fontsize=10, s='std'+' = {:.2e}'.format(np.std(white_noise))+r' [W]');
 
        # plot voltage noise histogram
        ax = axs[1][1];
        abs_max = max( abs(min(v_noise)), max(v_noise));
        ax.hist(v_noise, bins=20, range=[-abs_max, abs_max], histtype='stepfilled',
                align='mid', orientation='vertical', log=False, linewidth=0.5, linestyle='-', edgecolor='k',
                color='tab:blue', alpha=1.0, label='Voltage noise fluctuation', stacked=False);
        ax.set_xscale('linear');
        ax.set_yscale('linear');
        ax.set_xlabel(r'Voltage Noise [$\sqrt{W}$ / '+f'{freq_binwidth*1.e+6:.0f} kHz]', fontsize=15);
        ax.set_ylabel(f'# of points', fontsize=15);
        ax.grid(True);
        ax.legend();
        ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
        ax.tick_params(axis='y', labelsize=12);
        ax.minorticks_on();
        ax.xaxis.get_major_formatter().set_useOffset(False);
        ax.yaxis.get_major_formatter().set_useOffset(False);
        ax.locator_params(axis='x', nbins=5);
        ax.locator_params(axis='y', nbins=5);
        xlim  = ax.get_xlim();
        xwidth= xlim[1]-xlim[0];
        x0    = xlim[0]+xwidth*0.05;
        ylim  = ax.get_ylim();
        ywidth= ylim[1]-ylim[0];
        y0    = ylim[1]-ywidth*0.05;
        yscale=0.08
        ax.text(x0, y0-ywidth*yscale*1 , fontsize=10, s='std'+' = {:.2e}'.format(np.std(v_noise))+r' [$\sqrt{\mathrm{W}}$]');
        ax.text(x0, y0-ywidth*yscale*2 , fontsize=10, s=r'std$^2$'+' = {:.2e}'.format(np.std(v_noise)**2.)+r' [W]');
  
        # plot peak velocity integral histogram
        ax = axs[1][0];
        abs_max = max( abs(min(integral_velocity)), max(integral_velocity));
        ax.plot(freq_center, integral_velocity, linewidth=2, linestyle='-', color='tab:blue', label='Velocity PDF integral by freq. width');
        ax.set_xscale('linear');
        ax.set_yscale('linear');
        ax.set_xlabel(f'Frequency [GHz]', fontsize=15);
        ax.set_ylabel(f'Velocity PDF integral\nby freq. width (fraction)', fontsize=15);
        ax.grid(True);
        ax.legend();
        ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
        ax.tick_params(axis='y', labelsize=12);
        ax.minorticks_on();
        ax.xaxis.get_major_formatter().set_useOffset(False);
        ax.yaxis.get_major_formatter().set_useOffset(False);
        ax.locator_params(axis='x', nbins=5);
        ax.locator_params(axis='y', nbins=5);
  
        if verbosity>0: print(f'save: {outdir}/{outname}');
        fig.savefig(f'{outdir}/{outname}');
        pass;
 
    return freq_center/freq_conv, y, y_err;


def plot(x,y,y_err=None,fit_x=None,fit_y=None,label='Simulated spectrum',fit_label='Fitted spectrum',xlabel='x',ylabel='y',outdir='',outname='spectrum_fitted.pdf', texts=[]):
    # plot
    fig, axs = plt.subplots(1,1);
    fig.set_size_inches(8,6);
    plt.subplots_adjust(wspace=0.3, hspace=0.5, left=0.10, right=0.95,bottom=0.20, top=0.95)
 
    # plot spectrum
    ax = axs; # Only 1 axs
    ax.errorbar(x, y, xerr=None, yerr=y_err, color='tab:blue', label=label, marker='o', markersize=1., linestyle='-', linewidth=1.);
    if fit_x is not None and fit_y is not None:
        ax.plot(fit_x, fit_y, color='tab:red', label=fit_label, marker='', markersize=-0., linestyle='-', linewidth=1.);
        pass;
    ax.set_xscale('linear');
    ax.set_yscale('linear');
    ax.set_xlabel(xlabel, fontsize=15);
    ax.set_ylabel(ylabel, fontsize=15);
    ax.grid(True);
    ax.legend();
    ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
    ax.tick_params(axis='y', labelsize=12);
    ax.minorticks_on();
    #ax.set_xlim([freq[0],freq[-1]]);
    #ax.set_ylim([min(y)/2., max(y)*1.1]);
    ax.xaxis.get_major_formatter().set_useOffset(False);
    ax.yaxis.get_major_formatter().set_useOffset(False);
    ax.locator_params(axis='x', nbins=5);
    ax.locator_params(axis='y', nbins=5);
    xlim  = ax.get_xlim();
    xwidth= xlim[1]-xlim[0];
    x0    = xlim[0]+xwidth*0.05;
    ylim  = ax.get_ylim();
    ywidth= ylim[1]-ylim[0];
    y0    = ylim[1]-ywidth*0.05;
    yscale=0.05
    if texts is not None:
        for i, text in enumerate(texts) :
            ax.text(x0, y0-ywidth*yscale*(i+1) , fontsize=10, s=text);
            pass;
        pass;
 
    print(f'save: {outdir}/{outname}');
    fig.savefig(f'{outdir}/{outname}');
    return 0;
    


def fit(freq_binwidth = 5.e+3, # spectrum frequency bin width
        chi  = 2.e-10, # kumodes = 2.e-10
        T_noise = 10, # [K] System temperature
        gain = 1.e+6, # 1e+6 = +60dB
        time = 60., # [sec] 1min
        freq_0 = 20.e+9, # [Hz] peak position
        A_eff= (6.*1.e-2/2.)**2.*np.pi, # [m^2]
        rebin = 1,
        doPlot=True,outdir='',outname='spectrum_fitted.pdf',verbosity=2):
    freq_conv = 1.e-9; # Hz --> GHz

    def fitfunc(pars,x) :
        # pars[0]: peak position
        # pars[1]: P_DP (power of signal)
        # pars[2]: noise floor
        return peak_spectrum(x, freq_0=pars[0], P_DP=pars[1], a=0., b=pars[2], binwidth=freq_binwidth);

    init_pars  = [freq_0, 1.e-10, 1.e-10];
    limit_pars = [[10.e+9,30.e+9], [-1.,1.], [-1., 1.]];
    error_pars = [1.e-5,1.e-20,1.e-20];
    fix_pars   = [True,False,False];
    n_pars = len(init_pars);
    errdef = 1;

    x,y,y_err=create_spectrum(
        chi  = chi,  T_noise = T_noise, gain = gain, A_eff= A_eff,
        time = time, freq_0 = freq_0, freq_binwidth = freq_binwidth, rebin=rebin,
        doPlot=doPlot, verbosity=verbosity, outdir=outdir, outname='create_'+outname,
      );
    if verbosity>0: 
        print(f'x (size:{len(x)})={x}');
        print(f'y (size:{len(y)})={y}');
        print(f'y_err (size:{len(y_err)})={y_err}');
        pass;

    fitsquare = createfitsquare(x, y, func=fitfunc, nPar=n_pars, err=y_err);
    if verbosity>0:
        print('init_pars = ',init_pars);
        print('fitsquare return = ', fitsquare(*init_pars) );
        pass;
    result = minuitminosfit(fitsquare, init=init_pars, fix=fix_pars, limit=limit_pars, errordef=errdef, verbosity=verbosity);
    if verbosity>0: print('result =',result);
    pars = result[0][0];
    errs = result[0][3];
    freq0     = pars[0];
    freq0_err = errs[0];
    P_DP     = pars[1];
    P_DP_err = errs[1];
    noisefloor     = pars[2];
    noisefloor_err = errs[2];
    Ndof = len(x) - sum(fix_pars);

    if verbosity>0: printResult(result);

    if doPlot :
        # make plot with fitted result
        fit_x = np.linspace(x[0],x[-1],1000)
        fit_y = fitfunc(result[0][0], fit_x);
        #P_DP_digit = 10.**((int)(np.log10(abs(P_DP)))-1);
        #P_DP_text = f' = ({P_DP/P_DP_digit:.3f} '+r'$\pm$'+f' {P_DP_err/P_DP_digit:.3f})'+r'$\times 10^{'+'{:.0f}'.format(np.log10(P_DP_digit))+r'}$ [W]';
        P_DP_text = f' = {P_DP:.2e} '+r'$\pm$'+f' {P_DP_err:.2e} [W]';
        P_DP_limit = np.abs(P_DP) + P_DP_err;
        chi_limit  = calc_chi(P_DP_limit/gain, A_eff);
        #noisefloor_digit = 10.**((int)(np.log10(noisefloor))-1) if noisefloor>0. else 1.;
        #noisefloor_text = f' = ({noisefloor/noisefloor_digit:.3f} '+r'$\pm$'+f' {noisefloor_err/noisefloor_digit:.3f})'+r'$\times 10^{'+'{:.0f}'.format(np.log10(noisefloor_digit))+r'}$ [W]';
        noisefloor_text = f' = {noisefloor:.2e} '+r'$\pm$'+f' {noisefloor_err:.2e} [W]';
        plot(x*freq_conv,y,y_err,fit_x=fit_x*freq_conv,fit_y=fit_y,
                xlabel='Frequency [GHz]', 
                ylabel=f'Power [W / {freq_binwidth*1.e-3:.0f} kHz]',
                texts =[
                    'Fit result',
                    r'$f_{\mathrm{peak}}$'+f' = {freq0*1.e-9:.1f} '+ ' [GHz] (fixed)' if fix_pars[0] else( r'$\pm$'+f' {freq0_err*1.e-9:.1f} [GHz]') ,
                    r'$P_{\mathrm{DP}}$'+P_DP_text,
                    r'$P_{\mathrm{DP}}$ < '+f'{P_DP_limit:.1e}',
                    r'$\chi_{\mathrm{limit}}$ < '+f'{chi_limit:.1e}',
                    r'Noise floor'+noisefloor_text,
                    r'$\chi^2$'+f' = {result[0][1]:.1f}',
                    r'$N_{\mathrm{dof}}$'+f' = {Ndof}',
                    ],
                outdir=outdir, outname=outname);
        pass;

    return result[0];

def main(isTest=True,
        chi  = 2.e-10, # kumodes = 2.e-10
        T_noise = 10, # [K] System temperature
        gain = 1.e+6, # 1e+6 = +60dB
        time = 60., # [sec] 1min
        freq_0 = 20.e+9, # [Hz] peak position
        freq_binwidths = None, # [Hz] frequency bin width list
        A_phi= 6., # [cm] effective aperture phi
        nLoop=1000,
        rebin=1,
        outdir = 'figure',
        verbose = 0,
        outname = 'freq-bin-width_vs_fit-result.pdf',
        ):
    print(f'A_phi = {A_phi:.2f} cm');
    A_eff = (A_phi*1.e-2/2.)**2.*np.pi; # antenna

    if freq_binwidths is None :
        freq_binwidths= np.arange(1.0e+3,101.0e+3,10.e+3);
        pass;
    n_binwidth = len(freq_binwidths);
    if isTest: nLoop = 1;
    P_DP     = [];
    P_DP_err = [];
    for i, freq_binwidth in enumerate(freq_binwidths) :
        print(f'i={i}, freq_binwidth={freq_binwidth}');
        P_DP.append([]);
        P_DP_err.append([]);
        for n in range(nLoop):
            if n%50==0: print(f'n={n}');
            result = fit(freq_binwidth = freq_binwidth, # varied parameter
                    chi=chi, T_noise=T_noise, gain=gain, time=time, freq_0=freq_0, A_eff=A_eff, rebin=rebin, # fixed parameters
                    doPlot=True if n==0 else False,outdir=outdir, outname=f'spectrum_freqbinwidth{freq_binwidth*1.e-3:.1f}kHz_{outname}',verbosity=verbose+1 if n==0 else verbose);
            P_DP[i].append(result[0][1]);
            P_DP_err[i].append(result[3][1]);
            pass;
        pass;

    P_DP_mean = [ np.mean(array) for array in P_DP ];
    P_DP_err_mean = [ np.mean(array) for array in P_DP_err ];

    P_range = (np.min(P_DP), np.max(P_DP));

    # plot
    fig, axs = plt.subplots(2,2);
    fig.set_size_inches(12,8);
    plt.subplots_adjust(wspace=0.3, hspace=0.5, left=0.10, right=0.95,bottom=0.20, top=0.95)
    # plot P_DP_mean
    ax = axs[0][0];
    ax.errorbar(freq_binwidths*1.e-3, P_DP_mean, yerr=P_DP_err_mean, color='k', label='P_DP', capsize=5, alpha=0.3, marker='o', markersize=3., linestyle='', linewidth=1.);
    for i, freq_binwidth in enumerate(freq_binwidths) :
        ax.scatter([freq_binwidth*1.e-3]*nLoop, P_DP[i], color=colors[i], label=None, marker='o', s=1.);
        pass;
    ax.set_xscale('log');
    ax.set_yscale('linear');
    ax.set_xlabel('Frequency bin-width [kHz]', fontsize=15);
    ax.set_ylabel(r'Fitted $P_{\mathrm{DP}}$', fontsize=15);
    ax.grid(True);
    ax.legend();
    ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
    ax.tick_params(axis='y', labelsize=12);
    ax.minorticks_on();
    #ax.xaxis.get_major_formatter().set_useOffset(False);
    ax.yaxis.get_major_formatter().set_useOffset(False);
    #ax.locator_params(axis='x', nbins=5);
    ax.locator_params(axis='y', nbins=5);
    ax.set_ylim(P_range);
    # plot P_DP_err_mean
    ax = axs[0][1];
    ax.plot(freq_binwidths*1.e-3, P_DP_err_mean, color='tab:blue', label='P_DP', marker='o', markersize=2., linestyle='', linewidth=1.);
    ax.set_xscale('log');
    ax.set_yscale('linear');
    ax.set_xlabel('Frequency bin-width [kHz]', fontsize=15);
    ax.set_ylabel(r'Fitted error of $P_{\mathrm{DP}}$', fontsize=15);
    ax.grid(True);
    ax.legend();
    ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
    ax.tick_params(axis='y', labelsize=12);
    ax.minorticks_on();
    #ax.xaxis.get_major_formatter().set_useOffset(False);
    ax.yaxis.get_major_formatter().set_useOffset(False);
    #ax.locator_params(axis='x', nbins=5);
    ax.locator_params(axis='y', nbins=5);
    # plot P_DP histograms + Gaussian fit
    ax = axs[1][0];
    # Gaussian fit
    results = [];
    models  = [];

    for i, p in enumerate(P_DP):
        histo, bins = np.histogram(p, range=P_range, bins=20, density=False);
        bins_center = np.convolve(bins, [0.5,0.5], mode='valide');
        result, func = gauss_fit(bins_center, histo, y_err=None, stat_y_err=True);
        result=result[0];
        amp   = result[0][0];
        center= result[0][1];
        sigma = result[0][2];
        fit_x     = np.linspace(P_range[0], P_range[1], 1000);
        fit_gauss = func(fit_x, amp, center, sigma);
        models.append(func);
        results.append(result);
        ax.plot(fit_x, fit_gauss, color=colors[i], linestyle='--', linewidth=2, label='');
        pass;
    labels = [];
    for i, result in enumerate(results):
        sigma     = result[0][2]
        sigma_err = result[3][2];
        print(result);
        print(f'sigma = {sigma}+-{sigma_err}');
        label = r'$\Delta \nu$={:.2f} kHz'.format(freq_binwidths[i]*1e-3)+': gaussian $\sigma={:.1e}\pm{:.1e}$'.format(sigma, sigma_err);
        labels.append(label);
        pass; 
    ax.hist(P_DP, bins=20, range=P_range, stacked=False, align='mid', histtype='step', color=colors[0:len(P_DP)], label=labels);
    ax.set_xscale('linear');
    ax.set_yscale('linear');
    ax.set_xlabel(r'$P_{\mathrm{DP}}$ [W]', fontsize=15);
    ax.set_ylabel('# of fit', fontsize=15);
    ax.grid(True);
    ax.legend(mode = '',framealpha = 1,frameon = False,fontsize = 7,title='',borderaxespad=0.,labelspacing=1.2);
    ax.tick_params(axis='x', labelsize=12, labelrotation=60.);
    ax.tick_params(axis='y', labelsize=12);
    ax.minorticks_on();
    ax.xaxis.get_major_formatter().set_useOffset(False);
    ax.yaxis.get_major_formatter().set_useOffset(False);
    ax.locator_params(axis='x', nbins=5);
    ax.locator_params(axis='y', nbins=5);
 

    outfigname = 'aho.pdf' if isTest else f'{outdir}/{outname}';
    if verbose>0: print(f'save: {outfigname}');
    fig.savefig(outfigname);

    return 0;

if __name__=='__main__':
    isTest = False;
    nLoop=1000;
    chi  = 2.e-10; # kumodes = 2.e-10
    T_noise = 10; # [K] System temperature
    gain = 1.e+6; # 1e+6 = +60dB
    time = 60.; # [sec] 1min
    freq_0 = 20.e+9; # [Hz] peak position
    A_phi= 6.; # [cm] effective aperture phi
    rebin=1;
    outdir = 'figure';
    outname = 'freq-bin-width_vs_fit-result.pdf';
    verbose = 0;

    parser = argparse.ArgumentParser();
    parser.add_argument('--test', dest='test', type=int, default=(int)(isTest), help=f'Test or not (type=int, 0:False, 1:True | default: {isTest})');
    parser.add_argument('-n', '--nLoop', dest='nLoop', default=nLoop, type=int, help=f'Number of loop for each frequency bin-width (default: {nLoop} loops)');
    parser.add_argument('-c', '--chi', dest='chi', default=chi, type=float, help=f'chi (default: {chi:.2e} | kumodes limit = 2e-10 )');
    parser.add_argument('-T', '--T_noise', dest='T_noise', default=T_noise, type=float, help=f'Noise system temperature [K] (default: {T_noise:.1f} K)');
    parser.add_argument('-g', '--gain', dest='gain', default=gain, type=float, help=f'Gain by amplifiers (default: x{gain})');
    parser.add_argument('--gainDB', dest='gainDB', default=None, type=float, help=f'dB Gain by amplifiers (default: None)');
    parser.add_argument('-t', '--time', dest='time', default=time, type=float, help=f'Effective measurement time [sec] (default: {time} sec)');
    parser.add_argument('-A', '--A_phi', dest='A_phi', default=A_phi, type=float, help=f'Phi of effective circular aperture [cm] (A_eff=(A_phi/2)**2*pi [cm^2]) (default: {A_phi} cm)');
    parser.add_argument('-f', '--freq_0', dest='freq_0', default=freq_0, type=float, help=f'Peak frequency of signal [Hz] (default: {freq_0:.0e} Hz)');
    parser.add_argument('-b', '--binwidths', dest='freq_binwidths', default=None, type=str, help=f'Frequency bin width list [Hz] (default: None | ex.) "1e+3,2e+3,3e+3"');
    parser.add_argument('-r', '--rebin', dest='rebin', default=rebin, type=int, help=f'Number of rebinning for fitting (default: {rebin} bins)');
    parser.add_argument('--outdir', default=outdir, help=f'output directory (default: {outdir})');
    parser.add_argument('-o', '--outname', default=outname, help=f'output filename (default: {outname})');
    parser.add_argument('-v', '--verbose', default=verbose, type=int, help=f'Verbosity level (default: {verbose})');
    args = parser.parse_args();
    isTest = (args.test==1);
    nLoop  = args.nLoop;
    chi    = args.chi;
    T_noise= args.T_noise;
    gain   = args.gain;
    if args.gainDB is not None : gain = 10.**(args.gainDB*0.1);
    time   = args.time;
    A_phi  = args.A_phi
    freq_0 = args.freq_0;
    freq_binwidths = None if args.freq_binwidths is None else np.array([(float)(freq) for freq in args.freq_binwidths.split(',') ]);
    rebin = args.rebin
    outdir= args.outdir;
    outname= args.outname;
    verbose= args.verbose;

    print(f'isTest  = {isTest}');
    print(f'nLoop  = {nLoop}');
    print(f'chi     = {chi}');
    print(f'T_noise = {T_noise} [K]');
    print( 'gain    = {:+f} [dB]'.format(np.log10(gain)*10.));
    print(f'time    = {time} [sec]');
    print(f'A_phi   = {A_phi} [cm]');
    print(f'freq_0  = {freq_0*1.e-9:.0f} [GHz]');
    print(f'freq_binwidths  = {freq_binwidths} [Hz]');
    print(f'rebin   = {rebin} [bins]');
    print(f'outdir = {outdir}');
    print(f'outname = {outname}');
    print(f'verbose = {verbose}');

    main(isTest=isTest,
         chi  = chi,
         T_noise = T_noise,
         gain = gain,
         time = time, 
         freq_0 = freq_0, 
         freq_binwidths = freq_binwidths,
         A_phi= A_phi, 
         nLoop= nLoop,
         rebin= rebin,
         outdir = outdir,
         outname = outname,
         verbose = verbose,
         );
    pass;
