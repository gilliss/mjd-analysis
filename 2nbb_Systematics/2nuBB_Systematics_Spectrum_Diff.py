"""
Class that interfaces with the mjdsim results and ROOT objects.
This class should be able to convert root into pandas or numpy objects and return them
"""

import sys
from ROOT import TFile, TCanvas, TH1D, gPad, gApplication
import numpy as np
import matplotlib.pyplot as plt

def GetBinnedData(inFile, hName, nBinsX):
    f = TFile(inFile, 'READ')
    h = f.Get(hName)

    if (h.GetNbinsX() != nBinsX) or (h.GetNcells() != nBinsX + 2):
        print('Error', '  GetBinnedData: Unexpected binning', inFile, hName)

    hArray = np.zeros(nBinsX + 2, dtype = float)
    if h.GetEntries() > 0:
        # print('Debug', '  Working with:', h.GetName(), h.GetTitle(), h.GetNbinsX(), h.GetEntries(), h.Integral())
        hArray = np.frombuffer(h.GetArray(), dtype = 'float', count = nBinsX + 2, offset = 0) # getting array of data from PyDoubleBuffer object
        hArray = hArray[1:-1] # trimming off underflow and overflow bins
        hArray = np.copy(hArray)
        # print('Debug', '  np.sum(hArray) matches h.Integral():', np.sum(hArray) == h.Integral())
    else:
        print('Debug', '  0 entries. Returning hArray = zeros.')

    f.Close()
    del f
    del h

    return hArray

def GetSimData(inFile):
    return np.load(inFile)

def bin_ndarray(ndarray, new_shape, operation='sum'): # https://gist.github.com/derricw/95eab740e1b08b78c03f
    """
    Bins an ndarray in all axes based on the target shape, by summing or
        averaging.
    Number of output dimensions must match number of input dimensions.
    Example
    -------
    >>> m = np.arange(0,100,1).reshape((10,10))
    >>> n = bin_ndarray(m, new_shape=(5,5), operation='sum')
    >>> print(n)
    [[ 22  30  38  46  54]
     [102 110 118 126 134]
     [182 190 198 206 214]
     [262 270 278 286 294]
     [342 350 358 366 374]]
    """
    if not operation.lower() in ['sum', 'mean', 'average', 'avg']:
        raise ValueError("Operation {} not supported.".format(operation))
    if ndarray.ndim != len(new_shape):
        raise ValueError("Shape mismatch: {} -> {}".format(ndarray.shape,
                                                           new_shape))
    compression_pairs = [(d, c//d) for d, c in zip(new_shape,
                                                   ndarray.shape)]
    flattened = [l for p in compression_pairs for l in p]
    ndarray = ndarray.reshape(flattened)
    for i in range(len(new_shape)):
        if operation.lower() == "sum":
            ndarray = ndarray.sum(-1*(i+1))
        elif operation.lower() in ["mean", "average", "avg"]:
            ndarray = ndarray.mean(-1*(i+1))
    return ndarray

def FindBin(val, binCenters, binWidth):
    """
    Returns the bin index corresponding to a value along the binned axis
    """
    for i in range(len(binCenters)):
        if val >= (binCenters[i] - 0.5*binWidth) and val < (binCenters[i] + 0.5*binWidth): # included left edge, excluded right edge
            return i

from scipy.stats import chi2

def _count(a, axis=None):
    """
    Count the number of non-masked elements of an array.
    This function behaves like np.ma.count(), but is much faster
    for ndarrays.
    Ref: https://github.com/scipy/scipy/blob/v0.14.0/scipy/stats/stats.py#L3467
    """
    if hasattr(a, 'count'):
        num = a.count(axis=axis)
        if isinstance(num, np.ndarray) and num.ndim == 0:
            # In some cases, the `count` method returns a scalar array (e.g.
            # np.array(3)), but we want a plain integer.
            num = int(num)
    else:
        if axis is None:
            num = a.size
        else:
            num = a.shape[axis]
    return num

def chstwo(bins1, bins2, ddof=0, axis=0):
    """
    Chi-square test for difference between two data sets. Return the statistic and the p-value.
    Uses _count() to drop from the chi-square sum any entries for which both values are 0; the degrees of freedom are decremented for each dropped case.

    Comments on relation to NRC's chstwo() and SciPy's chi2.sf():
        -bins1 :: f_obs
        -bins2 :: f_exp
        -ddof :: adjustment to dof ... related to knstrn ... dof = num_obs - 1 - ddof ... see NRC discussion of arguments to chstwo algorithm.
            ... if data sets are of equal integral (perhaps normalized) then knstrn = 1, ddof = 0, dof = num_obs - 1 - 0
            ... if data sets are not of equal integral then knstrn = 0 ... ddof = -1, dof = num_obs - 1 - (-1) = num_obs
            ... could essentially rewrite as dof = num_obs - 1 - ddof --> dof = num_obs - knstrn - ddof, where ddof becomes any adjustment to dof beyond that of knstrn
        -Evaluating the prob from the chi2 distribution:
            ... NRC defines gammq as 1 - P, where P is probability that the observed chi2 for a correct model should be less than a value chi2.
                gammq(0.5*df, 0.5*chi2)
            ... For scipy, we'd have gammq = 1 - scipy.stats.chi2.cdf(x, df, loc=0, scale=1)
                or gammq = scipy.stats.chi2.sf(x, df, loc=0, scale=1) ... sf is survival function (also defined as 1 - cdf, but sf is sometimes more accurate)
    """
    # check inputs
    if len(bins1) != len(bins2):
        return 'Error: chstwo: len(bins1) != len(bins2)'
    # where bins1[i]=bins2[i]=0, mask entry i
    bins1, bins2 = np.ma.masked_where(condition = [(bins1 == bins2) & (bins1==0), (bins2 == bins1) & (bins2==0)], a = [bins1, bins2])
    # Do the test
    terms = (bins1 - bins2)**2 / (bins1 + bins2) # Terms with division by zero have been masked out. Terms evaluating to zero are kept.
    stat = terms.sum(axis=axis)
    num_obs = _count(terms, axis=axis) # returns number of non-masked terms in stat
    ddof = np.asarray(ddof)
    p = chi2.sf(stat, num_obs - 1 - ddof)
    # print('chi2.sf(stat = %f, dof = %d, ddof = %d)' % (stat, num_obs - 1 - ddof, ddof))
    return stat, p

###########
# Args
###########
data_config = str(sys.argv[1]) # ['0',...,'53','6']
sim_config = data_config[0]
bothMs = True
if int(data_config) < 5:
    bothMs = False

###########
# Useful vars and info
###########

# Dict for open data
Open_Enr_Exp_Dict = { # 'val': M1 + M2 (kg*dy), 'unc': sigma
'0': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'1': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'2': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'3': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'4': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'51': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'52': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'53': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'6': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)}
}

# Dict for blind data
Blind_Enr_Exp_Dict = { # 'val': M1 + M2 (kg*dy), 'unc': sigma
'1': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'2': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'53': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
'6': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)}
}

Open_DS_Enr_Exposure = Open_Enr_Exp_Dict[data_config]['val']
Open_DS_Enr_Exposure_Unc = Open_Enr_Exp_Dict[data_config]['unc']
Blind_DS_Enr_Exposure = Blind_Enr_Exp_Dict[data_config]['val']
Blind_DS_Enr_Exposure_Unc = Blind_Enr_Exp_Dict[data_config]['unc']
roiLo, roiHi = 1000, 1400

xmin_data = 0 # low edge (included)
xmax_data = 10000 # upper edge (excluded in ROOT, included in np and plt)
xBinWid_data = 1
nBinsX_data = int((xmax_data-xmin_data)/xBinWid_data)

xLimLo, xLimHi = 0, 3000

###########
# 20 keV Bins
###########

# histogram parameters
xmin = 0 # low edge (included)
xmax = 10000 # upper edge (excluded in ROOT, included in np and plt)
xBinWid = 20
nBinsX = int((xmax-xmin)/xBinWid)

xBinCenters = np.arange(xmin + 0.5*xBinWid, xmax + 0.5*xBinWid, xBinWid) # to be used as list of bin centers
roiLo_i, roiHi_i = FindBin(roiLo, xBinCenters, xBinWid), FindBin(roiHi, xBinCenters, xBinWid)

# get the open data
inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
inFile += '/open/' + ('DS%s_M1_CutScheme3.root' % data_config)
hName = 'h_Enr_CutScheme3'
dataOpen = GetBinnedData(inFile, hName, nBinsX_data)
if bothMs == True:
    inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
    inFile += '/' + 'open' + '/' + ('DS%s_M2_CutScheme3.root' % data_config)
    hName = 'h_Enr_CutScheme3'
    dataOpen += GetBinnedData(inFile, hName, nBinsX_data)
dataOpen = bin_ndarray(dataOpen, (nBinsX,), operation='sum')
dataOpenUnc = np.sqrt(dataOpen)
nCountsOpen = np.sum(dataOpen[roiLo_i:roiHi_i])

# get the blind data
inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
inFile += '/blind/' + ('DS%s_M1_CutScheme3.root' % data_config)
hName = 'h_Enr_CutScheme3'
dataBlind = GetBinnedData(inFile, hName, nBinsX_data)
if bothMs == True:
    inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
    inFile += '/' + 'blind' + '/' + ('DS%s_M2_CutScheme3.root' % data_config)
    hName = 'h_Enr_CutScheme3'
    dataBlind += GetBinnedData(inFile, hName, nBinsX_data)
dataBlind = bin_ndarray(dataBlind, (nBinsX,), operation='sum')
dataBlindUnc = np.sqrt(dataBlind)
nCountsBlind = np.sum(dataBlind[roiLo_i:roiHi_i])

# prep figure
fig = plt.figure() # clear current figure to prevent any previous figures or axes from persisting
plt.xlabel('keV')
plt.ylabel('c/kg/dy/keV (%d keV bins)' % xBinWid)
plt.xlim(xLimLo, xLimHi)

# scale data
dataOpen = dataOpen * (1/Open_DS_Enr_Exposure) / xBinWid
dataBlind = dataBlind * (1/Blind_DS_Enr_Exposure) / xBinWid

# plot data
plt.step(xBinCenters, dataOpen, where = 'mid', color = 'tab:blue', label = 'open')
plt.step(xBinCenters, dataBlind, where = 'mid', color = 'tab:orange', label = 'blind')

# plot text
plt.text(0.3, 0.70, 'Counts Open [%d,%d) = %.0f' % (roiLo, roiHi, nCountsOpen), transform=plt.gca().transAxes, fontsize=10)
plt.text(0.3, 0.65, 'Counts Blind [%d,%d) = %.0f' % (roiLo, roiHi, nCountsBlind), transform=plt.gca().transAxes, fontsize=10)

# plot legend
plt.legend()

# save figure
fig.savefig('DS%s_SpectrumDiff_%dkeVBins.pdf' % (data_config, xBinWid))

###########
# Show plots
###########
plt.show()
