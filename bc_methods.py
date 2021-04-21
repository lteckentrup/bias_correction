import numpy as np
import itertools
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression
import math
from ppt_adjustement import adjustPrecipFreq
from statsmodels.distributions.empirical_distribution import ECDF

def delta(obs, pred, cor):
    c = obs + (np.nanmean(cor) - np.nanmean(pred))
    return(c)

def scaling_add(obs, pred, cor):
    c = cor - np.nanmean(pred) + np.nanmean(obs)
    return(c)

def scaling_multi(obs, pred, cor):
    c = (cor/np.nanmean(pred)) * np.nanmean(obs)
    return(c)

def mva(obs, pred, cor):
    c = (cor - np.nanmean(pred)) + np.nanstd(obs)/np.nanstd(pred) + np.nanmean(obs)
    return(c)

def dpm(var, obs, pred, cor, threshold, nquantiles, detrend):

    eps = np.finfo(float).eps
    if (np.isnan(obs).all()==True and np.isnan(pred).any()==True and
        np.isnan(cor).any()==True):
        yout = itertools(np.nan, len(cor))

    else:
        ### Original method adds jitter; ignore
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

        obs_mean = np.mean(obs)
        pred_mean = np.nanmean(pred)

        if var == 'prec':
            cor = cor/pred_mean*obs_mean
        else:
            cor = cor-pred_mean+obs_mean

        if detrend == True:
            lr = LinearRegression()
            lr.fit(np.arange(1,len(cor)+1).reshape(-1,1), cor)
            cor_mean = lr.predict(np.arange(1,len(cor)+1).reshape(-1,1))
        else:
            cor_mean = list(itertools.repeat(obs_mean, len(cor)))
            cor_mean = np.array(cor_mean)

        if nquantiles == None or nquantiles == 0:
            nquant = max(len(obs), len(pred))
        else:
            nquant = nquantiles

        tau = np.arange(1,nquant+1)/(nquant+1)
        tau = np.insert(tau, 0, 0)
        tau = np.append(tau, 1)

        if obs.any() > math.sqrt(eps):
            x = mquantiles(pred/pred_mean, prob=tau, alphap=1, betap=1)
            y = mquantiles(obs/obs_mean, prob=tau, alphap=1, betap=1)

            ### in R rule=2:1, looks like result is identical to this?
            p2o= interp1d(x, y, kind='linear', bounds_error=False)

            yout=p2o(cor/cor_mean)
            extrap=np.isnan(yout)
            yout[extrap]=np.nanmax(obs/obs_mean)*((cor[np.isnan(yout)]/
                                                   cor_mean[np.isnan(yout)]))/ \
                                                  (np.nanmax(pred/pred_mean))
            yout = yout*cor_mean
        elif obs.any() < math.sqrt(eps):
            x = mquantiles(pred/pred_mean, prob=tau, alphap=1, betap=1)
            y = mquantiles(obs/obs_mean, prob=tau, alphap=1, betap=1)
            p2o= interp1d(x, y, kind='linear', bounds_error=False)
            yout=p2o(cor/cor_mean)

            cond1 = np.isnan(yout)
            cond2 = cor/cor_mean < np.nanmin(pred/pred_mean)
            cond3 = cor/cor_mean > np.nanmax(pred/pred_mean)

            extrap_lower = cond1 & cond2
            extrap_upper = cond1 & cond3

            yout[extrap_lower] = np.nanmin(obs/obs_mean)*((cor[extrap_lower]/
                                                           cor_mean[extrap_lower])/
                                                          np.nanmin(pred/pred_mean))
            yout[extrap_upper] = np.nanmax(obs/obs_mean)*((cor[extrap_upper]/
                                                           cor_mean[extrap_upper])/
                                                          np.nanmax(pred/pred_mean))
            yout=yout*cor_mean

        else:
            x = mquantiles(pred-pred_mean, prob=tau, alphap=1, betap=1)
            y = mquantiles(obs-obs_mean, prob=tau, alphap=1, betap=1)
            p2o= interp1d(x, y, kind='linear', bounds_error=False)
            yout=p2o(cor-cor_mean)

            cond1 = np.isnan(yout)
            cond2 = cor-cor_mean < np.nanmin(pred-pred_mean)
            cond3 = cor-cor_mean > np.nanmax(pred-pred_mean)

            extrap_lower = cond1 & cond2
            extrap_upper = cond1 & cond3

            yout[extrap_lower] = np.nanmin(obs-obs_mean)+((cor[extrap_lower]-
                                                           cor_mean[extrap_lower])-
                                                          np.nanmin(pred-pred_mean))
            yout[extrap_upper] = np.nanmax(obs-obs_mean)+((cor[extrap_upper]-
                                                           cor_mean[extrap_upper])-
                                                          np.nanmax(pred-pred_mean))

            yout = yout+cor_mean

        if var == 'prec':
            yout[np.where(yout < math.sqrt(eps))] = 0

    return(yout)

def eqm(var, obs, pred, cor, threshold, nquantiles, extrapolate):
    if var == 'prec':
        if np.isnan(obs).any() == False:
            p, nPo, nPp, Pth = adjustPrecipFreq(obs, pred, threshold)
        else:
            nP = None

        cor_map = list(itertools.repeat(np.nan,len(cor)))

        if (np.isnan(obs).any()==False) and (np.isnan(pred).any()==False):
            if (len(np.where(p>Pth)[0])>0):
                noRain = np.where(cor<=Pth)
                rain = np.where(cor>Pth)
                drizzle = np.where((cor>Pth) & (cor<=np.min(p[np.where(p>Pth)])))

                if len(rain) > 0:
                    ecdf = ECDF(cor[rain])
                    eFrc =  ecdf.x
                    eFrc  = eFrc[~np.isinf(eFrc)]

                    if nquantiles == None or nquantiles == 0:
                        nquant = len(p)
                    else:
                        nquant = nquantiles

                    nbins = nquant
                    binmid = np.arange((1./nbins), 1., 1./nbins)

                    ### alphap=1, betap=1: p(k) = (k-1)/(n-1):
                    ### p(k) = mode[F(x[k])]. (R type 7, R default)
                    qo = mquantiles(obs[np.where(obs>threshold)], prob=binmid,
                                    alphap=1, betap=1)
                    qp = mquantiles(p[np.where(p>Pth)], prob=binmid, alphap=1,
                                    betap=1)
                    p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
                    cor_map = cor.copy()
                    cor_map[rain]=p2o(cor[rain])

                    if extrapolate is None:
                        cor_map_rain = cor_map[rain]
                        cor_map_rain[np.where(cor[rain]>np.nanmax(qp))] = qo[len(qo)-1]
                        cor_map_rain[np.where(cor[rain]<np.nanmin(qp))] = qo[0]

                    elif extrapolate == 'constant':
                        cor_map_rain = cor_map[rain]
                        cor_map_rain[np.where(cor[rain]>np.nanmax(qp))] = \
                        cor[rain][np.where(cor[rain]>np.nanmax(qp))]+ \
                        (qo[len(qo)-1]-qp[len(qo)-1])
                        cor_map_rain[np.where(cor[rain]<np.nanmin(qp))] = \
                        cor[rain][np.where(cor[rain]<np.nanmin(qp))]+(qo[0]-qp[0])

                    cor_map[rain] = cor_map_rain
                else:
                    cor_map[rain] = list(itertools.repeat(0,len(cor)))

                if len(drizzle) > 0:
                    ### print(eFrc[1]) ### equals cor_map[drizzle]; original R:
                    ### probs = eFrc(s[drizzle])
                    cor_map[drizzle] = mquantiles(cor[cor>np.nanmin(p[np.where(p>Pth)])],
                                                  prob=0.1428571, alphap=0,
                                                  betap=1)
                    cor_map = np.array(cor_map)
                    cor_map[noRain] = 0
                else:
                    cor_map = cor

    else:
        if np.isnan(obs).all() == True:
            cor_map = list(itertools.repeat(np.nan,len(cor)))
        elif np.isnan(pred).all() == True:
            cor_map = list(itertools.repeat(np.nan,len(cor)))
        elif (np.isnan(obs).any() == False) and (np.isnan(pred).any() == False):
            if nquantiles == None or nquantiles == 0:
                nquant = len(pred)
            else:
                nquant = nquantiles

            nbins = nquant
            binmid = np.arange((1./nbins), 1., 1./nbins)
            qo = mquantiles(obs, prob=binmid, alphap=1, betap=1)
            qp = mquantiles(pred, prob=binmid, alphap=1, betap=1)
            p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
            cor_map=p2o(cor)
            if extrapolate == 'constant':
                cor_map[np.where(cor>np.nanmax(qp))] = \
                cor[np.where(cor>np.nanmax(qp))]+(qo[len(qo)-1]-qp[len(qo)-1]) ### fix
                cor_map[np.where(cor<np.nanmin(qp))] = \
                cor[np.where(cor<np.nanmin(qp))]+(qo[0]-qp[0])
            else:
                cor_map[np.where(cor>np.nanmax(qp))] = qo[len(qo)-1]
                cor_map[np.where(cor<np.nanmin(qp))] = qo[0]

    return(cor_map)

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
            yout[yout<threshold] = 0

    return(yout)
