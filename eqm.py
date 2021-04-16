import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
from ppt_adjustement import adjustPrecipFreq
import itertools
from statsmodels.distributions.empirical_distribution import ECDF

def eqm(var, obs, pred, s, threshold, nquantiles):
    if var == 'prec':
        if np.isnan(obs).any() == False:
            p, nPo, nPp, Pth = adjustPrecipFreq(obs, pred, threshold)
        else:
            nP = None

        smap = list(itertools.repeat(np.nan,len(s)))

        if (np.isnan(obs).any() == False) and (np.isnan(pred).any() == False):
            if (len(np.where(p > Pth)[0]) > 0):
                noRain = np.where(s <= Pth)
                rain = np.where(s > Pth)
                drizzle = np.where((s>Pth) & (s<= np.min(p[np.where(p > Pth)])))

                if (len(rain) > 0):
                    ecdf = ECDF(s[rain])
                    eFrc =  ecdf.x
                    eFrc  = eFrc[~np.isinf(eFrc)]
                    nquantiles = None
                    if nquantiles == None:
                        nquantiles = len(p)
                    else:
                        pass

                    nbins = nquantiles
                    binmid = np.arange((1./nbins), 1., 1./nbins)
                    qo = mquantiles(obs[np.where(obs>2)], prob=binmid, alphap=1, betap=1)
                    qp = mquantiles(p[np.where(p>Pth)], prob=binmid, alphap=1, betap=1)
                    p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
                    smap = s.copy()
                    smap[rain]=p2o(s[rain])

                    extrapolate = 'constant'
                    if extrapolate is None:
                        smap_rain = smap[rain]
                        smap_rain[np.where(s[rain]>np.nanmax(qp))] = qo[len(qo)-1]
                        smap_rain[np.where(s[rain]<np.nanmin(qp))] = qo[0]

                    elif extrapolate == 'constant':
                        smap_rain = smap[rain]
                        smap_rain[np.where(s[rain]>np.nanmax(qp))] = s[rain][np.where(s[rain] > np.nanmax(qp))] + (qo[len(qo)-1] - qp[len(qo)-1])
                        smap_rain[np.where(s[rain]<np.nanmin(qp))] = s[rain][np.where(s[rain] < np.nanmin(qp))] + (qo[0] - qp[0])

                    smap[rain] = smap_rain
                else:
                    smap[rain] = list(itertools.repeat(0,len(s)))


                if len(drizzle) > 0:
                    # print(eFrc[1]) ### equals smap[drizzle]; original R: probs = eFrc(s[drizzle])
                    smap[drizzle] = mquantiles(s[s > np.nanmin(p[np.where(p > Pth)])],
                                               prob=0.1428571, alphap=0, betap=1)
                    smap = np.array(smap)
                    smap[noRain] = 0

                else:
                    smap = s
    return(smap)

    else:
        if (all(is.na(o))):
            smap = list(itertools.repeat(np.nan,len(s)))
        elif (all(is.na(p))):
            smap = list(itertools.repeat(np.nan,len(s)))
        elif (any(!is.na(p)) & any(!is.na(o))):
            if (is.null(n.quantiles)) n.quantiles <- length(p):
                nbins = nquantiles
                binmid = np.arange((1./nbins), 1., 1./nbins)
                qo = mquantiles(obs[np.where(obs>2)], prob=binmid, alphap=1, betap=1)
                qp = mquantiles(p[np.where(p>Pth)], prob=binmid, alphap=1, betap=1)
                p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
                smap=p2o(s)
                if extrapolation == "constant":
                    smap[np.where(s>np.nanmax(qp))] = s[np.where(s>np.nanmax(qp))]+(qo[len(qo)-1]-qp[len(qo)-1])
                    smap[np.where(s<np.nanmin(qp))] = s[np.where(s < np.nanmin(qp))] + (qo[0] - qp[0])
                else:
                    smap[np.where(s>np.nanmax(qp))] = qo[len(qo)-1]
                    smap[np.where(s<np.nanmin(qp))] = qo[0]

    return(smap)
