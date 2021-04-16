import numpy as np
import itertools
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d

def qdm(var, obs, pred, cor, threshold, nquantiles):
    if (np.isnan(obs).all()==True) and (np.isnan(pred).any()==True) and (np.isnan(cor).any()==True):
        yout = itertools(np.nan, len(cor))

    else:
        # Apply a small amount of jitter to accomodate ties due to limited measurement precision -> ignore for now
        #o <- jitter(o, jitter.factor)
        #p <- jitter(p, jitter.factor)
        #s <- jitter(s, jitter.factor)

        # For ratio data, treat exact zeros as left censored values less than pr.threshold
        if var == 'prec':
            obs[np.where((obs<2) & (obs!=np.nan))] = np.random.uniform(low=0,
                                                                       high=threshold,
                                                                       size=np.nansum(obs<threshold))
            pred[np.where((pred<2) & (pred!=np.nan))] = np.random.uniform(low=0,
                                                                          high=threshold,
                                                                          size=np.nansum(pred<threshold))
            cor[np.where((cor<2) & (cor!=np.nan))] = np.random.uniform(low=0,
                                                                       high=threshold,
                                                                       size=np.nansum(cor<threshold))

        # Calculate empirical quantiles using Weibull plotting position
        n = max(len(obs), len(pred), len(cor))
        if nquantiles == 0:
            nquantiles = n

        tau = (np.linspace(1./(n+1), n/(n+1), num=12))
        ## precision of mquantile??
        quantobs = mquantiles(obs, prob = tau, alphap=0, betap=0)
        quantpred = mquantiles(pred, prob = tau, alphap=0, betap=0)
        quantcor = mquantiles(cor, prob = tau, alphap=0, betap=0)

        # Apply QDM bias correction original code had rule = 2??
        p2o = interp1d(quantcor, tau, kind='linear', bounds_error=False)
        taucor = p2o(cor)

        if var == 'prec':
            p2o_pred = interp1d(tau, quantpred, kind='linear', bounds_error=False)
            delta = cor/p2o_pred(taucor)
            p2o_yout = interp1d(tau, quantobs)
            yout = p2o_yout(taucor)
        else:
            p2o_pred = interp1d(tau, quantpred, kind='linear', bounds_error=False)
            delta = cor - p2o_pred(taucor)
            p2o_yout = interp1d(tau, quantobs)
            yout = p2o_yout(taucor) + delta

        # For precip data, set values less than threshold to zero
        if var == 'prec':
            yout[yout < threshold] = 0

    return(yout)
