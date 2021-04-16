import numpy as np
from scipy.stats.mstats import mquantiles
from scipy.interpolate import interp1d
from ppt_adjustement import adjustPrecipFreq
import itertools
from statsmodels.distributions.empirical_distribution import ECDF

obs = np.array([448.2, 172.0881, 118.9816, 5.797349, 2, 0.7, 0.7, 0.1, 0.7,
                14, 41.78181, 94.99255])
pred = np.array([80.2817878635324, 268.832150047484, 270.459391305064,
                 113.710020941381, 28.1730174899666, 2.60631106038364,
                 2.81899801580404, 1.4109973622908, 1.23972919673748,
                 8.74956907080104, 29.8733589632633,   276.032919084638])
s = np.array([80.2817878635324, 268.832150047484, 270.459391305064,
              113.710020941381, 28.1730174899666, 2.60631106038364,
              2.81899801580404, 1.4109973622908, 1.23972919673748,
              8.74956907080104, 29.8733589632633,   276.032919084638])

if np.isnan(obs).any() == False:
    p, nPo, nPp, Pth = adjustPrecipFreq(obs, pred, 2)
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
                smap[rain][np.where(s[rain]>np.nanmax(qp))] = qo[len(qo)-1]
                smap[rain][np.where(s[rain]<np.nanmin(qp))] = qo[0]

            elif extrapolate == 'constant':
                smap[rain][np.where(s[rain]>np.nanmax(qp))] = s[rain][np.where(s[rain] > np.nanmax(qp))] + (qo[len(qo)-1] - qp[len(qo)-1])
                smap[rain][np.where(s[rain]<np.nanmin(qp))] = s[rain][np.where(s[rain] < np.nanmin(qp))] + (qo[0] - qp[0])

        else:
            smap[rain] = list(itertools.repeat(0,len(s)))


        if len(drizzle) > 0:
            # print(eFrc[1]) ### equals smap[drizzle]
            smap[drizzle] = mquantiles(s[s > np.nanmin(p[np.where(p > Pth)])], prob=0.1428571, alphap=0, betap=1)
            smap = np.array(smap)
            smap[noRain] = 0
            # print(smap)

        else:
            smap = s
