import numpy as np
import math
import itertools
from scipy.stats import gamma

def adjustPrecipFreq(obs, pred, threshold):
    o = obs[np.isfinite(obs)]
    p = pred[np.isfinite(pred)]

    # Number of dry days in 'o'
    nPo = len(o[o<=threshold])

    # Number of dry days that must be in 'p' to equal precip frequency in 'o'
    nPp = math.ceil(len(p)*nPo/len(o))

    # Index and values of ordered 'p'
    ix = np.argsort(p)
    Ps = sorted(p)
    
    # Code fails if threshold is zero 
    if nPo == 0:
        Pth = 0
    else:
        Pth = np.nanmax(Ps[nPo-1:nPo+1]) # in case nPp == length(Ps)

    print(Pth)
    # Themeßl (Themessl) modification (simulating rain for model dry days)
    inddrzl = np.where(np.array(Ps[nPp:len(Ps)]) < threshold)[0]
    if len(inddrzl) > 0:
        Os = sorted(o)
        indO = math.ceil(len(Os) * (nPp + max(inddrzl))/len(Ps))
        auxOs = Os[nPo:indO+1]
        if len(np.unique(np.array(auxOs))) > 6:
            # simulate precip for 'p' with a gamma adjusted in 'o' for values between
            shape, loc, scale = gamma.fit(auxOs, floc=0)
            Ps[nPp:(nPp+max(inddrzl))+1] = np.random.gamma(shape, scale, 8)
        else:
            Ps[nPp:(nPp+max(inddrzl))+1] = list(itertools.repeat(np.nanmean(auxOs),
                                                                 len(Ps[nPp:(nPp+max(inddrzl))+1])))

        # order 'Ps' after simulation
        Ps = sorted(Ps)

    # Make 0-s
    if nPo > 0:
        ind = min(nPp, len(p))
        Ps[:ind] = list(itertools.repeat(0, ind))

    p[ix] = Ps

    pred[np.isfinite(pred)] = p

    return(p, nPo, nPp, Pth)
