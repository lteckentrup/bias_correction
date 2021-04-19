import numpy as np
import itertools
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d

def qdm(var, obs, pred, cor, threshold, nquantiles):

    eps = np.finfo(float).eps

    if (np.isnan(obs).all()==True and np.isnan(pred).any()==True and
        np.isnan(cor).any()==True):
        yout = itertools(np.nan, len(cor))

    else:
        ### For ratio data, treat exact zeros as left censored values less than
        ### threshold. Original R script adds jitter; ignore
        if var == 'prec':
            obs[np.where((obs<threshold) & (obs!=np.nan))] = \
            np.random.uniform(low=eps, high=threshold,
                              size=np.nansum(obs<threshold))
            pred[np.where((pred<threshold) & (pred!=np.nan))] = \
            np.random.uniform(low=eps, high=threshold,
                              size=np.nansum(pred<threshold))
            cor[np.where((cor<threshold) & (cor!=np.nan))] = \
            np.random.uniform(low=eps, high=threshold,
                              size=np.nansum(cor<threshold))

        # Calculate empirical quantiles using Weibull plotting position
        n = max(len(obs), len(pred), len(cor))
        if nquantiles == None or nquantiles == 0:
            nquant = n
        else:
            nquant = nquantiles

        tau = (np.linspace(1./(n+1), n/(n+1), num=nquant))

        ## precision of mquantile??
        quantobs = mquantiles(obs, prob = tau, alphap=0, betap=0)
        quantpred = mquantiles(pred, prob = tau, alphap=0, betap=0)
        quantcor = mquantiles(cor, prob = tau, alphap=0, betap=0)

        # Apply QDM bias correction original code had rule = 2??
        p2o = interp1d(quantcor, tau, kind='linear', bounds_error=False)
        taucor = p2o(cor)

        if var == 'prec':
            p2o_pred = interp1d(tau, quantpred, kind='linear', bounds_error=False)
            ### code might break if data that needs to be correct is zero so add eps
            delta = cor/(p2o_pred(taucor)+np.finfo(float).eps)
            p2o_yout = interp1d(tau, quantobs)
            yout = p2o_yout(taucor)*delta
        else:
            p2o_pred = interp1d(tau, quantpred, kind='linear', bounds_error=False)
            delta = cor - p2o_pred(taucor)
            p2o_yout = interp1d(tau, quantobs)
            yout = p2o_yout(taucor) + delta

        # For precip data, set values less than threshold to zero
        if var == 'prec':
            yout[yout < threshold] = 0

    return(yout)
