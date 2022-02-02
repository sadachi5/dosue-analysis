from iminuit import Minuit;
import numpy as np;
import time;
import matplotlib.pyplot as plt;

def drawcontour(
    minuit,
    fix,
    parIDs=[1,2],
    parLabels=None,
    numpoints=100,
    center=None,
    dxdy=None, #[dx,dy] (needs "center")
    outname='aho', 
    verbosity = 0,
    ) :

    if fix[parIDs[0]] or fix[parIDs[1]] :
        print('The contour parameters are fixed. Could not draw contours! --> Skip!!');
        return -1;
    
    contours = [];
    parname1 = 'x{}'.format(parIDs[0]);
    parname2 = 'x{}'.format(parIDs[1]);
    bound = 2.;
    if (not dxdy is None) and (not center is None) :
        bound = [[center[0]-dxdy[0],center[0]+dxdy[0]], 
                 [center[1]-dxdy[1],center[1]+dxdy[1]]];
        pass;
    x, y, fval = minuit.contour(parname1,parname2,size=numpoints,bound=bound);

    cont = plt.contour(x,y,fval);
    cont.clabel(fmt='%1.1f', fontsize=12);
    if not center is None :
        plt.scatter([center[0]],[center[1]],marker='*',label='center',color='k');
        plt.text(center[0],center[1],'(x,y)=({:.4},{:.3})'.format(center[0],center[1]));
        pass;
    xlabel = parname1 if parLabels is None else parLabels[0];
    ylabel = parname2 if parLabels is None else parLabels[1];
    plt.xlabel(xlabel,fontsize=16);
    plt.ylabel(ylabel,fontsize=16);
    plt.grid();
    plt.savefig('{}.png'.format(outname));
    plt.close();
    return 0;


def minuitminosfit(
    fitfunc,
    init,
    fix,
    limit,
    errordef = 1,
    verbosity=0
    ) :

    nPar = len(init);
    minuit = Minuit(
        fitfunc, 
        *init);
    if verbosity>0: print('minuit :\n{}'.format(minuit));
    names = minuit.parameters;
    minuit.errordef = Minuit.LEAST_SQUARES;
    for i in range(nPar) :
        minuit.fixed[i]=fix[i];
        minuit.limits[i]=limit[i];
        pass;
    
    minuit.migrad();
    minuit.hesse();
    minuit.minos();
    values = minuit.values;
    errs = minuit.errors;
    fmin = minuit.fval;
    nfcn = minuit.nfcn; 
    if verbosity>0: 
        print('values = {}'.format(values));
        print('errs = {}'.format(errs));
        print('fmin = {}'.format(fmin));
        print('nfcn = {}'.format(nfcn));
        print('covariance :\n{}'.format(minuit.covariance));
        print('is_valid   = {}'.format(minuit.valid   ));
        print('is_accurate= {}'.format(minuit.accurate));
        print('nfcn = {}'.format(nfcn));
        for i in range(nPar) :
            print('par{} = {}'.format(i, values[i]));
            pass;
        pass;
 
    result = [values, fmin, nfcn, errs];
 
    #del fitfunc, nPar, names, fmin, minuit;
    del fitfunc, nPar, names, fmin;
    return result, minuit;


def createfitfunc(func, x, y, err=None) : 

    def fitfunc(pars) : 
        return func(pars,x,y,err);
 
    return fitfunc;


def createfitsquare(x, y, func, nPar, err=None, errx=None, n_fix_pars=0, verbosity=0) :

    # initialize error
    #print(err);
    if err is None and len(err)==1 :
        print('err = None --> err  = [1,1,...]');
        err = [ 1 for i in range(len(x)) ];
        pass;
 
    def fitsquare(*pars) :
        y_fit = func(pars, x);
        diff = y - y_fit;
        errsquare = err*err if errx==None else err*err + errx*errx ;
        square = np.sum( np.power(diff,2.) / errsquare );
        return square;
        ## !NOTE!: You should consider which calculation is correct!
        #nof  = len(x) - len(pars)-n_fix_pars-1;
        #return np.sum( np.power( diff / err, 2.) )/nof;
 
    return fitsquare;


def truncateX(x, data, xlim=[None, None], err=None) :
    xmin = min(x);
    xmax = max(x);
    if xlim[0]==None : xlim[0] = xmin;
    if xlim[1]==None : xlim[1] = xmax;
    x_truncate    = [];
    data_truncate = [];
    err_truncate  = [];
    for i in range(len(x)):
        if x[i] < xlim[0] or x[i] > xlim[1] : continue;
        x_truncate   .append(x   [i]);
        data_truncate.append(data[i]);
        if len(err)>i : err_truncate .append(err [i]);
        pass;
    return np.array(x_truncate), np.array(data_truncate), np.array(err_truncate);


def printResult(result, parlabels=None, verbosity=0) :
    #result = [values, fmin, nfcn, errs];
    result = result[0];
    print('fmin (minimized value) = {}'.format(result[1]));
    print('nfcn (# of calls)      = {}'.format(result[2]));
    errs    = result[3];
    for n in range(len(result[0])) :
        if parlabels is None : parlabel = 'par[{}]'.format(n);
        else                 : parlabel = parlabels[n];
        print('{:30s} = {:.2e} +- {:.2e}'.format(parlabel, result[0][n],errs[n]));
        pass;
    return 0;

def gaussian(x,amp,center,sigma):
    return (amp/(np.sqrt(2.*np.pi)*sigma)) * np.exp(-(1.0*x-center)**2 / (2*sigma**2));

def gauss_fit(x,y,y_err=None,stat_y_err=False) :
    def fitfunc(pars, x):
        return gaussian(x, pars[0], pars[1], pars[2]);
    nPars = 3;
    if y_err is None : y_err = np.full(len(y), 1.);
    if stat_y_err : y_err = np.where( y>0., np.sqrt(y), 1.e+100 );
    print(f'stat_y_err={stat_y_err} / y_err = {y_err}');
    fitsquare = createfitsquare(x, y, func=fitfunc, nPar=nPars, err=y_err);

    print(f'x ({len(x)})={x}');
    print(f'y ({len(y)})={y}');
    center0= np.average(x,weights=y);
    height0= np.max(y);
    Ntotal = np.sum(y)
    sigma0 = np.sqrt( np.sum( np.power(x - center0,2.)*y ) / Ntotal );
    amp0   = height0 * np.sqrt(2.*np.pi)*sigma0;
    print(f'init center = {center0:.2e}');
    print(f'init sigma  = {sigma0:.2e}');
    print(f'init height = {height0:.2e}');
    print(f'init amp    = {amp0:.2e}');
    init_pars = [amp0, center0, sigma0];
    fix_pars  = [False,False,False];
    limit_pars= [[0.,amp0*100.], [center0-sigma0*100.,center0+sigma0*100], [0., sigma0*100.]];
    errdef    = 1;
    result    = minuitminosfit(fitsquare, init=init_pars, fix=fix_pars, limit=limit_pars, errordef=errdef, verbosity=0);
    return result, gaussian;


if __name__=='__main__' :
    t = np.linspace(0.,200.,2000);
    theta = t/100.;
    err = 0.01;
    errdef = 1;
    y = np.sin(2.*np.pi*theta) + 0.2 + np.random.randn()*err;
    #y = np.sin(2.*np.pi*theta) + 0.2;
    y_err = np.full(len(t), err);

    def fitfunc(pars,x,y) :
        return np.multiply(np.sin(np.multiply(x, 2.*np.pi*pars[1]) + pars[2]),  pars[0])  + pars[3];

    init_pars  = [1., 0.01, 0., 0.];
    limit_pars = [[-10.,10.], [-10.,10.], [-np.pi, np.pi], [-10.,10.]];
    #fix_pars   = [False,True,True,True];
    fix_pars   = [False,False,False,False];

    t_truncate, y_truncate, err_truncate = truncateX(t, y, [50., 150.], y_err);

    fitsquare = createfitsquare(t_truncate, y_truncate, fitfunc, len(init_pars), err_truncate);
    result    = minuitminosfit(fitsquare, init=init_pars, fix=fix_pars, limit=limit_pars, errordef=errdef, verbosity=2);
    errs      = result[3];

    printResult(result);

    t_fitrange = np.linspace(50., 150., 1000);
    fitted_y   = fitfunc(result[0], t_fitrange, None);

    plt.errorbar(t, y, y_err);
    plt.plot(t_fitrange, fitted_y, color='r');
    plt.savefig('aho.png');

    pass;
