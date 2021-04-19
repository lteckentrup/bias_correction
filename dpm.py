
import numpy as np
import itertools
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
from sklearn.datasets import make_regression
from sklearn.linear_model import LinearRegression

def qdm(var, obs, pred, cor, threshold, nquantiles, detrend):
    if (np.isnan(obs).all()==True) and (np.isnan(pred).any()==True) and (np.isnan(cor).any()==True):
        yout = itertools(np.nan, len(cor))

    else:
        ### Original method adds jitter; ignore
        if var == 'prec':
            obs[np.where((obs<2) & (obs!=np.nan))] = np.random.uniform(low=np.finfo(float).eps,
                                                                       high=threshold,
                                                                       size=np.nansum(obs<threshold))
            pred[np.where((pred<2) & (pred!=np.nan))] = np.random.uniform(low=np.finfo(float).eps,
                                                                          high=threshold,
                                                                          size=np.nansum(pred<threshold))
            cor[np.where((cor<2) & (cor!=np.nan))] = np.random.uniform(low=np.finfo(float).eps,
                                                                       high=threshold,
                                                                       size=np.nansum(cor<threshold))

        obs_mean = np.mean(obs)
        pred_mean = np.nanmean(pred)

        if var == 'prec':
            cor_copy = cor.copy()
            cor_copy = cor_copy/pred_mean*obs_mean
        else:
            cor_copy = cor.copy()
            cor_copy = cor_copy-pred_mean+obs_mean

        if detrend == True:
            lr = LinearRegression()
            lr.fit(np.arange(1,len(cor)+1).reshape(-1,1), cor)
            obs_mean = lr.predict(np.arange(1,len(cor)+1).reshape(-1,1))
        else:
            cor_mean=obs_mean

        if nquantiles == 0:
            nquantiles = max(len(obs), len(pred))

        tau = np.arange(1,nquantiles+1)/(nquantiles+1)
        tau = np.insert(tau, 0, 0)
        tau = np.append(tau, 1)

        if obs.any() > np.finfo(float).eps:
            x = mquantiles(pred/pred_mean, prob=tau, alphap=1, betap=1)
            y = mquantiles(obs/obs_mean, prob=tau, alphap=1, betap=1)

            ### in R rule=2:1, looks like result is identical to this?
            p2o= interp1d(x, y, kind='linear', bounds_error=False)

            yout=p2o(cor/cor_mean)
            extrap=np.isnan(yout)
            yout[extrap]=np.nanmax(obs/obs_mean)*((cor[np.isnan(yout)]/cor_mean))/(np.nanmax(pred/pred_mean))
            yout = yout*cor_mean
        elif obs.any() < np.finfo(float).eps:
            x = mquantiles(pred/pred_mean, prob=tau, alphap=1, betap=1)
            y = mquantiles(obs/obs_mean, prob=tau, alphap=1, betap=1)
            p2o= interp1d(x, y, kind='linear', bounds_error=False)
            yout=p2o(cor/cor_mean)]

            cond1 = np.isnan(yout)
            cond2 = cor/cor_mean < np.nanmin(pred/pred_mean)
            cond3 = cor/cor_mean > np.nanmax(pred/pred_mean)

            extrap_lower = cond1 & cond2
            extrap_upper = cond1 & cond3

            yout[extrap_lower] = np.nanmin(obs/obs_mean)*((cor[extrap_lower]/cor_mean)/np.nanmin(pred/pred_mean))
            yout[extrap_upper] = np.nanmax(obs/obs_mean)*((cor[extrap_upper]/cor_mean)/np.nanmax(pred/pred_mean))
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

            yout[extrap_lower] = np.nanmin(obs-obs_mean)+((cor[extrap_lower]-cor_mean)-np.nanmin(pred-pred_mean))
            yout[extrap_upper] = np.nanmax(obs-obs_mean)+((cor[extrap_upper]-cor_mean)-np.nanmax(pred-pred_mean))

            yout = yout+cor_mean

        if var == 'prec':
            yout[np.where(yout < math.sqrt(np.finfo(float).eps))] = 0
            
    return(yout)
