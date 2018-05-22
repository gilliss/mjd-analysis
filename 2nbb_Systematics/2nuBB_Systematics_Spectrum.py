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



###########
# Args
###########
data_config = str(sys.argv[1]) # ['0',...,'53','6']
sim_config = data_config[0]
dataType = 'open'
bothMs = True
if int(data_config) < 5:
    bothMs = False

###########
# Useful vars and info
###########
# Dict for open data
Enr_Exp_Dict = { # 'val': M1 + M2 (kg*dy), 'unc': sigma
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
# Enr_Exp_Dict = { # 'val': M1 + M2 (kg*dy), 'unc': sigma
# '1': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
# '2': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
# '53': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)},
# '6': {'val': 0. + 0., 'unc': np.sqrt(0.**2 + 0.**2)}
# }

DS_Enr_Exposure = Enr_Exp_Dict[data_config]['val']
DS_Enr_Exposure_Unc = Enr_Exp_Dict[data_config]['unc']
dyPerYr = 365.25
roiLo, roiHi = 1000, 1400

xmin_data = 0 # low edge (included)
xmax_data = 10000 # upper edge (excluded in ROOT, included in np and plt)
xBinWid_data = 1
nBinsX_data = int((xmax_data-xmin_data)/xBinWid_data)

xLimLo, xLimHi = 0, 3000

eff_2nbb = 0.9

###########
# 1 keV Bins
###########

# histogram parameters
xmin = 0 # low edge (included)
xmax = 10000 # upper edge (excluded in ROOT, included in np and plt)
xBinWid = 1
nBinsX = int((xmax-xmin)/xBinWid)

xBinCenters = np.arange(xmin + 0.5*xBinWid, xmax + 0.5*xBinWid, xBinWid) # to be used as list of bin centers
roiLo_i, roiHi_i = FindBin(roiLo, xBinCenters, xBinWid), FindBin(roiHi, xBinCenters, xBinWid)

# get the data and sim
inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
inFile += '/' + dataType + '/' + ('DS%s_M1_CutScheme3.root' % data_config)
hName = 'h_Enr_CutScheme3'
data = GetBinnedData(inFile, hName, nBinsX_data)
if bothMs == True:
    inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
    inFile += '/' + dataType + '/' + ('DS%s_M2_CutScheme3.root' % data_config)
    hName = 'h_Enr_CutScheme3'
    data += GetBinnedData(inFile, hName, nBinsX_data)
dataUnc = np.sqrt(data)
nCounts = np.sum(data[roiLo_i:roiHi_i])

inFile = '/global/homes/g/gilliss/BuildSpectra_Output'
inFile += '/' + ('DS%s/cut2/l4' % sim_config)
inFile += '/' + ('2v_2_DS%s.npy' % sim_config)
sim = GetSimData(inFile)

# prep figure
plt.figure() # clear current figure to prevent any previous figures or axes from persisting
plt.xlabel('keV')
plt.ylabel('c/kg/dy/keV (%d keV bins)' % xBinWid)
plt.xlim(xLimLo, xLimHi)

# plot data
data = data * (1/DS_Enr_Exposure) / xBinWid
dataUnc = dataUnc * (1/DS_Enr_Exposure) / xBinWid
plt.step(xBinCenters, data, where = 'mid', color = 'k', label = 'data')

# plot sim
sim = sim / dyPerYr # (cnt/kg/yr) / (dy/yr) = (cnt/kg/dy)
sim = eff_2nbb * sim
plt.step(xBinCenters, sim, where = 'mid', color = 'darkorange', label = '2v sim')

# plot text
plt.text(0.3, 0.70, 'Counts Data [%d,%d) = %.0f' % (roiLo, roiHi, nCounts), transform=plt.gca().transAxes, fontsize=10)
plt.text(0.3, 0.65, 'Rates [%d,%d) (c/kg/dy):' % (roiLo, roiHi), transform=plt.gca().transAxes, fontsize=10)
integral = np.sum(data[roiLo_i:roiHi_i]) * xBinWid
sigma_N = np.sqrt(nCounts)
sigma_MT = DS_Enr_Exposure_Unc
sigma_tot = np.sqrt( ( (1/DS_Enr_Exposure) * sigma_N )**2 + ( (-nCounts/(DS_Enr_Exposure**2)) * sigma_MT )**2 )
plt.text(0.3, 0.60, '    Data = %.2f +/- %.4f ($\sigma_{N}$:%.2f, $\sigma_{MT}$:%.2f)' % (integral, sigma_tot, sigma_N, sigma_MT), transform=plt.gca().transAxes, fontsize=10)
integral = np.sum(sim[roiLo_i:roiHi_i]) * xBinWid
plt.text(0.3, 0.55, '    Sim = %.2f' % integral, transform=plt.gca().transAxes, fontsize=10)

# plot legend
plt.legend()

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

# get the data and sim
inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
inFile += '/' + dataType + '/' + ('DS%s_M1_CutScheme3.root' % data_config)
hName = 'h_Enr_CutScheme3'
data = GetBinnedData(inFile, hName, nBinsX_data)
if bothMs == True:
    inFile = '/global/homes/g/gilliss/2nuBB_Systematics/Update_May2018/output_data'
    inFile += '/' + dataType + '/' + ('DS%s_M2_CutScheme3.root' % data_config)
    hName = 'h_Enr_CutScheme3'
    data += GetBinnedData(inFile, hName, nBinsX_data)
data = bin_ndarray(data, (nBinsX,), operation='sum')
dataUnc = np.sqrt(data)
dataUncTot = np.sqrt( ( (1/DS_Enr_Exposure) * np.sqrt(data) )**2 + ( (-data/(DS_Enr_Exposure**2)) * DS_Enr_Exposure_Unc )**2 ) / xBinWid
nCounts = np.sum(data[roiLo_i:roiHi_i])

inFile = '/global/homes/g/gilliss/BuildSpectra_Output'
inFile += '/' + ('DS%s/cut2/l4' % sim_config)
inFile += '/' + ('2v_2_DS%s.npy' % sim_config)
sim = GetSimData(inFile)
sim = bin_ndarray(sim, (nBinsX,), operation='sum')

# prep figure
fig20 = plt.figure()
plt.xlabel('keV')
plt.ylabel('c/kg/dy/keV (%d keV bins)' % xBinWid)
plt.xlim(xLimLo, xLimHi)

# plot data
data = data * (1/DS_Enr_Exposure) / xBinWid
dataUnc = dataUnc * (1/DS_Enr_Exposure) / xBinWid
plt.step(xBinCenters, data, where = 'mid', color = 'k', label = 'data')

# plot uncert
lowerUnc = data-dataUnc
np.place(lowerUnc, lowerUnc<0, 0) # clip < 0
lowerUncTot = data-dataUncTot
np.place(lowerUncTot, lowerUncTot<0, 0) # clip < 0
plt.fill_between(xBinCenters, data+dataUncTot, lowerUncTot, step = 'mid', color = 'lightblue', label = 'tot. unc.')
# plt.fill_between(xBinCenters, data+dataUnc, lowerUnc, step = 'mid', color = 'deepskyblue', label = 'stat. unc.')
# plt.fill_between(xBinCenters, data+dataUnc, lowerUnc, step = 'mid', color = 'lightblue', label = 'stat. unc.')

# plot sim
sim = sim / dyPerYr / xBinWid # (cnt/kg/yr) / (dy/yr) = (cnt/kg/dy)
sim = eff_2nbb * sim
plt.step(xBinCenters, sim, where = 'mid', color = 'darkorange', label = '2v sim')

# plot text
plt.text(0.3, 0.70, 'Counts Data [%d,%d) = %.0f' % (roiLo, roiHi, nCounts), transform=plt.gca().transAxes, fontsize=10)
plt.text(0.3, 0.65, 'Rates [%d,%d) (c/kg/dy):' % (roiLo, roiHi), transform=plt.gca().transAxes, fontsize=10)
integral = np.sum(data[roiLo_i:roiHi_i]) * xBinWid
sigma_N = np.sqrt(nCounts)
sigma_MT = DS_Enr_Exposure_Unc
sigma_tot = np.sqrt( ( (1/DS_Enr_Exposure) * sigma_N )**2 + ( (-nCounts/(DS_Enr_Exposure**2)) * sigma_MT )**2 )
plt.text(0.3, 0.60, '    Data = %.2f +/- %.4f ($\sigma_{N}$:%.2f, $\sigma_{MT}$:%.2f)' % (integral, sigma_tot, sigma_N, sigma_MT), transform=plt.gca().transAxes, fontsize=10)
integral = np.sum(sim[roiLo_i:roiHi_i]) * xBinWid
plt.text(0.3, 0.55, '    Sim = %.2f' % integral, transform=plt.gca().transAxes, fontsize=10)

# plot legend
plt.legend()

# save figure
fig20.savefig('DS%s_Spectrum_%dkeVBins.pdf' % (data_config, xBinWid))

###########
# Show plots
###########
plt.show()
