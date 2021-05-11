'''
Functions for compressing high vertical resolution turbulence profiles to a smaller 
number of layers for use in Monte Carlo simulation.

Three methods here:

    Equivalent Layers - simplest, conserves isoplanatic angle (Fusco 1999, DOI: 10.1117/12.363606)

    Optimal Grouping - optimal for tomographic AO (Saxenhuber 2017, DOI: 10.1364/AO.56.002621)

    Generalised Conservation of Turbulence Moments (GCTM) - conserves arbitrary 
        turbulence moments (Saxenhuber 2017, DOI: 10.1364/AO.56.002621)

Optimal Grouping in particular relies on numba to run fast.
'''
import numpy
from scipy.optimize import minimize
from numba import njit
from numpy import float64, int64, ndarray
from typing import List, Tuple

def equivalent_layers(h: ndarray, p: ndarray, L: int) -> Tuple[ndarray, ndarray]:
    '''
    Equivalent layers method of profile compression (Fusco 1999).

    Splits the profile into L "slabs", then sets the height of each slab as the 
    effective height ((integral cn2(h) * h^{5/3} dh) / integral cn2(h) dh)^(3/5)
    and the cn2 as the sum of cn2 in that slab.

    Parameters
        h (numpy.ndarray): heights of input profile layers
        p (numpy.ndarray): cn2dh values of input profile layers
        L (int): number of layers to compress down to
    
    Returns
        h_L (numpy.ndarray): compressed profile heights
        cn2_L (numpy.ndarray): compressed profile cn2dh per layer 

    '''
    h_el = numpy.zeros(L)
    cn2_el = numpy.zeros(L)

    hstep = (h.max()-h.min())/L
    alt_bins = numpy.arange(h.min(), h.max(), hstep)
    ix = numpy.digitize(h, alt_bins)
    for i in range(L):
        ix_tmp = ix==i+1
        cn2_el[i] = p[ix_tmp].sum()
        h_el[i] = ((p[ix_tmp] * h[ix_tmp]**(5/3)).sum() / p[ix_tmp].sum())**(3/5)
    return h_el, cn2_el

def optimal_grouping(R: int, L: int, h: ndarray, p: ndarray) -> Tuple[ndarray, ndarray]:
    '''
    Python implementation of algorithm 2 from Saxenhuber et al (2017). Performs the 
    "optimal grouping" algorithm , which finds the grouping that minimises the 
    cost function given in Eq. 7 of that paper. Requires numba.

    Parameters
        R (int): number of random starting groupings (recommended 10?)
        L (int): number of layers to compress down to
        h (numpy.ndarray): input profile heights
        p (numpy.ndarray): input profile cn2dh per layer

    Returns
        h_L (numpy.ndarray): compressed profile heights
        cn2_L (numpy.ndarray): compressed profile cn2dh per layer
    '''
    N = len(p)

    # set initial best grouping to be (approx) equal splits 
    gamma_best = numpy.linspace(0,N,L+1,dtype=int)[1:-1]
    gamma_best, G_best = _optGroupingMinimization(gamma_best, h, p)
    for r in range(R):
        gamma_new, G_new = _optGroupingMinimization(_random_grouping(N,L), h, p)
        if G_new < G_best:
            gamma_best = gamma_new
            G_best = G_new

    cn2 = []
    hmin_best = _G(_convert_splits_to_groups(gamma_best,N), h, p, return_hmin=True)[1]
    for groups in _convert_splits_to_groups(gamma_best,N):
        cn2.append(p[groups].sum())

    return numpy.array(hmin_best), numpy.array(cn2)

def GCTM(h: ndarray, p: ndarray, L: int, h_scaling: float=10000., cn2_scaling: float=100e-15) -> Tuple[ndarray, ndarray]:
    '''
    Generalised Conservation of Turbulence Moments compression method, from 
    Saxenhuber et al 2017. This compresses the profile whilst conserving 2L-1 
    moments of the profile, where a moment N is ((integral cn2(h) h^N)/integral cn2(h))^{1/N}.

    Parameters
        h (numpy.ndarray): heights of input profile layers
        p (numpy.ndarray): cn2dh values of input profile layers
        L (int): number of layers to compress down to
        h_scaling[optional] (float): scaling of heights to avoid overflow
        cn2_scaling[optional] (float): scaling of cn2 values to avoid overflow

    Returns
        h_L (numpy.ndarray): compressed profile heights
        cn2_L (numpy.ndarray): compressed profile cn2 
    '''
    mom0 = _moments(h/h_scaling, p/cn2_scaling, L)
    bounds = [(0,None) for i in range(2*L)]
    
    # calculate first guess using equivalent layers
    guess_h, guess_cn2 = equivalent_layers(h,p,L) 
    x0 = numpy.hstack([guess_h/h_scaling, guess_cn2/cn2_scaling]) # scale so that we don't have overflow issues
    res = minimize(_moments_minfunc, x0, args=(L, mom0), bounds = bounds)['x']
    return res[:L]*h_scaling, res[L:] * cn2_scaling

def _random_grouping(N: int, L: int) -> ndarray:
    '''
    Creates a random grouping of N elements into L groups
    '''
    options = numpy.arange(0,N-2)
    splits = numpy.sort(numpy.random.choice(options, size=L-1, replace=False))
    return splits 

def _vicinity(grouping: ndarray, N: int) -> List[ndarray]:
    '''
    Calculate the vicinity of a grouping, see Saxenhuber (2017)
    '''
    grouping_borders = grouping.copy()
    grouping_borders = numpy.insert(grouping_borders,0,-1)
    grouping_borders = numpy.append(grouping_borders,N-1)

    # calculate all possible splittings for each grouping
    pre_merge = []
    for i in range(len(grouping_borders)-1):
        for j in range(grouping_borders[i]+1, grouping_borders[i+1]):
            grouping_tmp = grouping.copy()
            grouping_tmp = numpy.insert(grouping_tmp,i,j)
            pre_merge.append(grouping_tmp)
    
    # of these groupings, merge every two neighbouring groupings in each
    merged = []
    for i in range(len(pre_merge)):
        for j in range(len(pre_merge[i])):
            merge_tmp = pre_merge[i].copy()
            merged.append(numpy.delete(merge_tmp,j))

    return merged 

def _convert_splits_to_groups(splits: ndarray, N: int) -> List[ndarray]:
    '''
    Convert indices of splits into explicit groupings
    '''
    out = []
    for i in range(len(splits)):
        if i == 0:
            out.append(numpy.arange(0,splits[i]+1))
        if i == len(splits)-1:
            out.append(numpy.arange(splits[i]+1, N))
            break
        out.append(numpy.arange(splits[i]+1, splits[i+1]+1))
    return out

def _G(grouping: List[ndarray], h: ndarray, p: ndarray, return_hmin: bool=False) -> Tuple[float64, List[int64]]:
    '''
    Optimal Grouping cost function (slow version)
    '''
    cost_function_value = 0
    hmin = []
    for group in grouping:
        ht_array = numpy.tile(group,(len(group),1)).T
        cost_function = (p[group] * numpy.abs(h[ht_array] - h[group])).sum(1)
        cost_function_value += cost_function.min()
        if return_hmin:
            hmin.append(h[group[numpy.argmin(cost_function)]])
    if return_hmin:
        return cost_function_value, hmin
    return cost_function_value

@njit
def _Gjit(splits, h, p, return_hmin=False):
    '''
    Optimal Grouping cost function (fast numba version)
    '''
    N = len(p)
    # this is just convert_splits_to_groups, numba doesnt like it outside in its own function
    grouping = []
    for i in range(len(splits)):
        if i == 0:
            grouping.append(numpy.arange(0,splits[i]+1))
        if i == len(splits)-1:
            grouping.append(numpy.arange(splits[i]+1, N))
            break
        grouping.append(numpy.arange(splits[i]+1, splits[i+1]+1))
        
    cost_function_value = 0.
    hmin = []
    for group in grouping:
        ptmp = p[group]
        htmp = h[group]
        cost_function = numpy.zeros(len(group)) 
        for i,g in enumerate(group):
            cost_function[i] = ((ptmp * numpy.abs(htmp-h[g])).sum())
        cost_function_value += numpy.min(cost_function)
        if return_hmin:
            hmin.append(h[group[numpy.argmin(cost_function)]])
    return cost_function_value

def _optGroupingMinimization(start_grouping: ndarray,h: ndarray,p: ndarray,maxiter: int=200) -> Tuple[ndarray, float64]:
    '''
    Minimisation of G (slow version)
    '''
    gamma_new = start_grouping
    N = len(p)
    for i in range(maxiter):
        gamma_old = gamma_new
        V = numpy.array(_vicinity(gamma_old, N))
        Gvals = numpy.zeros(len(V))
        for j,v in enumerate(V):
            Gvals[j] = _Gjit(v,h,p)
        Gopt_ix = numpy.argmin(Gvals)
        G_new = Gvals[Gopt_ix]
        gamma_new = V[Gopt_ix]
        if (gamma_new == gamma_old).all():
            break
    return gamma_new, G_new

def _moments(h: ndarray, p: ndarray, L: int) -> ndarray:
    '''
    Generalised turbulence moments function for GCTM
    '''
    mom = numpy.array([p * h**(i) for i in range(2*L-1)])
    return mom.sum(1)

def _moments_minfunc(args: ndarray, L: int, mom0: ndarray) -> float64:
    '''
    Minimisation function for GCTM. Minimises the sum of square difference between 
    2L-1 moments (see Saxenhuber et al 2017).
    '''
    h = args[:L]
    cn2 = args[L:]
    return ((_moments(h, cn2, L) - mom0)**2).sum()
