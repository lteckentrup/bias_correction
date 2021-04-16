import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
from ppt_adjustement import adjustPrecipFreq
import itertools
from statsmodels.distributions.empirical_distribution import ECDF

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

                if (len(rain) > 0):
                    ecdf = ECDF(cor[rain])
                    eFrc =  ecdf.x
                    eFrc  = eFrc[~np.isinf(eFrc)]

                    if nquantiles == None:
                        nquantiles = len(p)
                    else:
                        pass

                    nbins = nquantiles
                    binmid = np.arange((1./nbins), 1., 1./nbins)

                    ### alphap=1, betap=1: p(k) = (k-1)/(n-1):
                    ### p(k) = mode[F(x[k])]. (R type 7, R default)
                    qo = mquantiles(obs[np.where(obs>2)], prob=binmid, alphap=1, betap=1)
                    qp = mquantiles(p[np.where(p>Pth)], prob=binmid, alphap=1, betap=1)
                    p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
                    cor_map = cor.copy()
                    cor_map[rain]=p2o(cor[rain])

                    if extrapolate is None:
                        cor_map_rain = cor_map[rain]
                        cor_map_rain[np.where(cor[rain]>np.nanmax(qp))] = qo[len(qo)-1]
                        cor_map_rain[np.where(cor[rain]<np.nanmin(qp))] = qo[0]

                    elif extrapolate == 'constant':
                        cor_map_rain = cor_map[rain]
                        print(cor_map_rain)
                        cor_map_rain[np.where(cor[rain]>np.nanmax(qp))] = cor[rain][np.where(cor[rain]>np.nanmax(qp))]+(qo[len(qo)-1]-qp[len(qo)-1])
                        cor_map_rain[np.where(cor[rain]<np.nanmin(qp))] = cor[rain][np.where(cor[rain]<np.nanmin(qp))]+(qo[0]-qp[0])

                    cor_map[rain] = cor_map_rain
                else:
                    cor_map[rain] = list(itertools.repeat(0,len(cor)))

                if len(drizzle) > 0:
                    # print(eFrc[1]) ### equals cor_map[drizzle]; original R: probs = eFrc(s[drizzle])
                    cor_map[drizzle] = mquantiles(cor[cor>np.nanmin(p[np.where(p>Pth)])],
                                               prob=0.1428571, alphap=0, betap=1)
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
            if nquantiles == None:
                nquantiles = len(pred)
            else:
                pass
            nbins = nquantiles
            binmid = np.arange((1./nbins), 1., 1./nbins)
            qo = mquantiles(obs, prob=binmid, alphap=1, betap=1)
            qp = mquantiles(pred, prob=binmid, alphap=1, betap=1)
            p2o = interp1d(qp, qo, kind='linear', bounds_error=False)
            cor_map=p2o(cor)
            if extrapolate == 'constant':
                cor_map[np.where(cor>np.nanmax(qp))] = cor[np.where(cor>np.nanmax(qp))]+(qo[len(qo)-1]-qp[len(qo)-1])
                cor_map[np.where(cor<np.nanmin(qp))] = cor[np.where(cor<np.nanmin(qp))]+(qo[0]-qp[0])
            else:
                cor_map[np.where(cor>np.nanmax(qp))] = qo[len(qo)-1]
                cor_map[np.where(cor<np.nanmin(qp))] = qo[0]

    return(cor_map)
