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
dataOpenUncTot = np.sqrt( ( (1/Open_DS_Enr_Exposure) * np.sqrt(dataOpen) )**2 + ( (-dataOpen/(Open_DS_Enr_Exposure**2)) * Open_DS_Enr_Exposure_Unc )**2 ) / xBinWid
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
dataBlindUncTot = np.sqrt( ( (1/Blind_DS_Enr_Exposure) * np.sqrt(dataBlind) )**2 + ( (-dataBlind/(Blind_DS_Enr_Exposure**2)) * Blind_DS_Enr_Exposure_Unc )**2 ) / xBinWid
nCountsBlind = np.sum(dataBlind[roiLo_i:roiHi_i])

# scale data
dataOpen = dataOpen * (1/Open_DS_Enr_Exposure) / xBinWid
dataBlind = dataBlind * (1/Blind_DS_Enr_Exposure) / xBinWid

# prep figure
fig = plt.figure('fig_array3d_sigma')
import matplotlib.gridspec as gridspec
gs = gridspec.GridSpec(6, 1) # row, col
ax = plt.subplot(gs[0:3, 0])
ax_resid = plt.subplot(gs[3, 0])
ax_OB = plt.subplot(gs[4, 0])
ax_BO = plt.subplot(gs[5, 0])
plt.subplots_adjust(hspace = 0., top = 0.98, bottom = 0.1)

ax.set_ylabel('c/kg/dy/keV (%d keV bins)' % xBinWid)
ax.set_xlim(xLimLo, xLimHi)
ax.set_ylim(0., 0.035)
ax.set_xticks([])

# plot uncerts
lowerBlindUncTot = dataBlind-dataBlindUncTot
np.place(lowerBlindUncTot, lowerBlindUncTot<0, 0) # clip < 0
lowerOpenUncTot = dataOpen-dataOpenUncTot
np.place(lowerOpenUncTot, lowerOpenUncTot<0, 0) # clip < 0
ax.fill_between(xBinCenters, dataBlind+dataBlindUncTot, lowerBlindUncTot, step = 'mid', color = 'tab:orange', alpha = 0.2)
ax.fill_between(xBinCenters, dataOpen+dataOpenUncTot, lowerOpenUncTot, step = 'mid', color = 'tab:blue', alpha = 0.2)

# plot data
ax.step(xBinCenters, dataOpen, where = 'mid', color = 'tab:blue', label = 'open')
ax.step(xBinCenters, dataBlind, where = 'mid', color = 'tab:orange', label = 'blind')

# plot text
ax.text(0.3, 0.70, 'Counts Open [%d,%d) = %.0f' % (roiLo, roiHi, nCountsOpen), transform=plt.gca().transAxes, fontsize=10)
ax.text(0.3, 0.65, 'Counts Blind [%d,%d) = %.0f' % (roiLo, roiHi, nCountsBlind), transform=plt.gca().transAxes, fontsize=10)

# plot residual
dataResid = dataOpen - dataBlind
dataResidUncTot = np.sqrt( dataOpenUncTot**2 + dataBlindUncTot**2 )
bin100keV = FindBin(100,xBinCenters,xBinWid)
bin2039keV = FindBin(2039,xBinCenters,xBinWid)
ax_resid.fill_between(xBinCenters[bin100keV:], (dataResid+dataResidUncTot)[bin100keV:], (dataResid-dataResidUncTot)[bin100keV:], step = 'mid', color = 'lightgrey')
ax_resid.scatter(xBinCenters[bin100keV:], dataResid[bin100keV:], s = 1, color = 'k')
ax_resid.set_ylabel('O - B')
ax_resid.set_xlabel('keV')
ax_resid.set_xlim(xLimLo, xLimHi)
ax_resid.set_xticks([])
# plot the distribution of the residuals
hist, bin_edges = np.histogram(a = dataResid[bin100keV:bin2039keV], bins = 'auto', density = True)
hist_binWid = np.full(len(hist), bin_edges[1] - bin_edges[0])
hist_left = np.zeros(len(hist))
hist_left = np.full(len(hist), xLimHi)
hist  = -1 * hist
ax_resid.barh(bin_edges[:-1], hist, align = 'edge', height = hist_binWid, left = hist_left, color = 'r', alpha = 0.3, zorder = 0)
# ax_resid.hist(dataResid[bin100keV:bin2039keV], bins=30, orientation="horizontal", color = 'r', alpha = 0.2)

# plot O-B/sigO
dataSigO = np.abs(dataOpen - dataBlind)/dataOpenUncTot
ax_OB.scatter(xBinCenters[bin100keV:], dataSigO[bin100keV:], s = 1, color = 'k')
ax_OB.set_xlim(xLimLo, xLimHi)
ax_OB.set_ylim(0, 10)
ax_OB.set_ylabel('|$\Delta$|/$\sigma_{O}$')
ax_OB.set_xticks([])

# plot B-O/sigB
dataSigB = np.abs(dataOpen - dataBlind)/dataBlindUncTot
ax_BO.scatter(xBinCenters[bin100keV:], dataSigB[bin100keV:], s = 1, color = 'k')
ax_BO.set_xlim(xLimLo, xLimHi)
ax_BO.set_ylim(0, 10)
ax_BO.set_ylabel('|$\Delta$|/$\sigma_{B}$')
ax_BO.set_xlabel('keV')

# plot legend
ax.legend()

# save figure
# fig.savefig('DS%s_SpectrumDiff_%dkeVBins.pdf' % (data_config, xBinWid))

###########
# Show plots
###########
plt.show()
