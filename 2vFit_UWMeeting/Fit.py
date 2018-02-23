# This script reads in a text file of a clean MJD spectrum as well as a text file of a 2vBB pdf.
# Once text files are read in, they are turned into histograms and plotted.
# Next, a spectral fit is performed. Starting off simply. ..

# ISSUE is likely that I need to do the fit in terms of counts, not in terms of very small normalized rates of counts. For starters, the k parameter of the poisson distribution needs to be an integer ... For instance poisson.pmf(10.25,10) evaluates to 0 b/c 10.25 is not an integer.
# Imagine X+Y as histograms in units of counts, and lowercase letters as normalization factors:                          c(X+Y) = aX + bY ... c = (aX+bY)/(X+Y)

# Currently getting 1.87910814374e+21 compared to 1.926e21. Including deadlayer efficiency would bring up my guess... Also, folding in BG sources will bring up my guess (i.e. right now, unaccounted-for BG counts are giving the 2vBB region more counts and driving down the halflife).
# peak_x  32700
# 1-alpha  0.68
# upperBound  33350
# lowerBound  32050
# T12  1.87910814374e+21
# T12+  3.81098375486e+19
# T12-  3.66242966547e+19
#### 1.88 (+/- 0.04) x 10^21 y where +/- is smallest 68% interval 

###############################
# IMPORTS
###############################

from scipy.stats import rv_continuous, rv_discrete, poisson, norm, uniform
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

###############################
# PREP PLOTS
###############################

xlo = 0.
xhi = 3000.
wid = 10.
nbins = (xhi-xlo)/wid
xscale_lo = 0
xscale_hi = 3000
yscale_lo = 0.5
yscale_hi = 3e2

fit_lo = 500 # energy fit range limits
fit_hi = 2039
lo_S = 28000 # scale factor search limits
hi_S = 38000
scaleStep = 10 # step size for scale factor search

fig1, (ax1, ax2) = plt.subplots(2, 1)
fig1.subplots_adjust(top=None,bottom=None,hspace=0.5)

fig2, (ax3, ax4) = plt.subplots(2, 1)
fig2.subplots_adjust(top=None,bottom=None,hspace=0.5)

fig3, (ax5, ax6) = plt.subplots(2, 1)
fig3.subplots_adjust(top=None,bottom=None,hspace=0.5)

fig4, (ax7) = plt.subplots(1, 1)
fig4.subplots_adjust(top=None,bottom=None,hspace=0.5)

###############################
# DATA
###############################

# READ IN FILE
dataFilePath = "/Users/tomgilliss/Desktop/UNC/ENAP/Simulation/2vBB_Fit/2vBB_Fit/HistDataForFit/DS34_AllCuts_HistContents.txt"
data_x, data_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

# PLOT DATA
ax1.step(data_x, data_y, color='k') #ax1.bar(data_x, data_y, width=wid, color='k', fill=False)
ax1.set_title('MJD Data, DS3+DS4')
ax1.set_xlabel('Energy (keV)')
ax1.set_ylabel('Cnt/10keV')
ax1.set_yscale('log',nonposy='clip')
ax1.set_ylim(yscale_lo,yscale_hi)
ax1.set_xlim(xscale_lo,xscale_hi)

ax7.step(data_x, data_y, color='k') #ax1.bar(data_x, data_y, width=wid, color='k', fill=False)
ax7.set_title('MJD Data, DS3+DS4')
ax7.set_xlabel('Energy (keV)')
ax7.set_ylabel('Cnt/10keV')
ax7.set_yscale('log',nonposy='clip')
ax7.set_ylim(yscale_lo,yscale_hi)
ax7.set_xlim(xscale_lo,xscale_hi)

###############################
# 2vBB PDF
###############################

# READ IN FILE
dataFilePath = "/Users/tomgilliss/Desktop/UNC/ENAP/Simulation/2vBB_Fit/2vBB_Fit/HistDataForFit/76Ge2vBBSums_HistContents.txt"
pdf_x, pdf_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

# PLOT DATA
ax2.step(pdf_x, pdf_y, color='r')
ax2.set_title('2vBB PDF')
ax2.set_xlabel('Energy (keV)')
ax2.set_ylabel('Cnt/10keV')
ax2.set_yscale('log',nonposy='clip')
ax2.set_xlim(xscale_lo,xscale_hi)

###############################
###############################
###############################
###############################
# Fit
###############################
###############################
###############################
###############################

###############################
# DEFINITIONS FOR FIT
###############################

# Flat prior distribution for scaling factor S (note, this might not be the most uninformative distr for this model).

def PriorS(s):
    return uniform.pdf(s, loc=lo_S, scale=hi_S);
def vals_prior_S(n):
    return uniform.rvs(loc=lo_S, scale=hi_S, size=n)

# Likelihood function. Input is data and output is product of likelihoods for each data point. Poisson likelihood, with 2vBB model as expected value, and data as input.
def Likelihood(S,lo_E,hi_E): # this version of the function uses meshgrid magic
    product = 1;
    for i in range(0,int(nbins)):
        if(data_x[i]<hi_E and data_x[i]>lo_E and pdf_x[i]<hi_E and pdf_x[i]>lo_E):
            temp = poisson.pmf(data_y[i],S*pdf_y[i])
            product = product*temp
            if(temp == 0):
                print (data_x[i],data_y[i],S*pdf_y[i],temp)
    return product

###############################
# DO FIT
###############################

# Grid search over S from lo_S to hi_S in steps of 10
scaleVals = np.arange(lo_S,hi_S,scaleStep)
results = np.zeros(scaleVals.size)
for i in range(0,scaleVals.size):
    results[i] = Likelihood(scaleVals[i],fit_lo,fit_hi)

# PLOT RESULTS
#ax3.step(scaleVals, results, color='r')
#ax3.set_title('Likelihood')
#ax3.set_xlabel('Scale Factor')
#ax3.set_ylabel('Likelihood')

###############################
# EXTRACT THE BEST FIT
###############################

S_bestfit = scaleVals[np.argmax(results)]
print "S_bestfit %d" % (S_bestfit) # 32700 ... compare to 31923.8310602 GERDA data and some efficiencies

# Normalize the likelihood distribution
integral_Likelihood = 0.
for i in range(0,results.size):
    integral_Likelihood = integral_Likelihood + results[i]*scaleStep
results_norm = results*(1/integral_Likelihood)

# Check that normalization worked
print "non-normalized likelihood integral ", integral_Likelihood
integral_Likelihood = 0.
for i in range(0,results_norm.size):
    integral_Likelihood = integral_Likelihood + results_norm[i]*scaleStep
print "normalized likelihood integral ", integral_Likelihood

# PLOT NORMALIZED LIKELIHOOD DISTRIBUTION
ax3.step(scaleVals, results_norm, color='r')
ax3.set_title('Likelihood')
ax3.set_xlabel('Scale Factor')
ax3.set_ylabel('Likelihood')
ax3.set_ylim(bottom=0)

###############################
# TRANSLATE TO T1/2
###############################

# INTEGRAL OF COUNTS FROM 2vBB PDF [J. Kotila and F. Iachello, Phys. Rev. C 85, 034316 (19 March 2012)]
pdfCnts = 0.1

# GENERAL
N_A = 6.022140857e23 # (1/mol)
m_A = 75.6e-3 # (kg/mol) molar mass of 76Ge

# FROM GERDA
T12 = 1.926e21 # (yr); GERDA: Eur. Phys. J. C, 75:416, 2015

# FROM MJD
f_enr = 0.88 # (enrichment fraction) MJD
f_av = 1 # (active volume fraction) ****this is currently absorbed in mAct_EnrDets_DS*
e_av = 0.8 # https://docs.google.com/spreadsheets/d/1wk8nY7PH4ndDnZ6AJ_myhlDYI0_f7WBPH6ILLOjAbNI/edit#gid=1829051574**** a guess; not sure about this

runtime_DS3 = 2585130/(60.*60.*24.*365.25) # (yr) **** should I use the mean Julian year?
mAct_EnrDets_DS3 = 12.631 # (kg) active enriched mass from isGood Enr dets; includes f_av
mActEnr_EnrDets_DS3 = f_enr*mAct_EnrDets_DS3;

runtime_DS4 = 2047670/(60.*60.*24.*365.25) # (yr) **** should I use the mean Julian year?
mAct_EnrDets_DS4 = 5.471 # (kg) active enriched mass from isGood Enr dets; includes f_av
mActEnr_EnrDets_DS4 = f_enr*mAct_EnrDets_DS4;

# CALCULATE
SumObjects = e_av*runtime_DS3*mActEnr_EnrDets_DS3 + e_av*runtime_DS4*mActEnr_EnrDets_DS4
N_obs = S_bestfit*pdfCnts
T12 = np.log(2)*(N_A/m_A)*(1/N_obs)*SumObjects
print "T12 ", T12

###############################
# GET UNCERTAINTY FROM LIKELIHOOD
###############################

# Determine minimal interval of width 1-alpha
stepSize = scaleStep
distr_x = scaleVals
distr_y = results_norm
peak_x_i = np.argmax(distr_y)
peak_x = distr_x[peak_x_i]
cumulativeP = 0.
lowerBound = peak_x
upperBound = peak_x

alpha = 1 - 0.68 # 68% central interval
print "1-alpha ", 1-alpha

# FOR DATA THAT COVERS EVENLY SPACED INTEGERS
if(distr_y[peak_x_i] > (1-alpha)):
    lowerBound = peak_x
    upperBound = peak_x
else:
    intervWidth = 0
    while cumulativeP < (1-alpha):
        cumulativeP = 0
        if(peak_x < intervWidth - np.floor(intervWidth/2.)):
            lowerBound = 0
            upperBound = intervWidth
        else:
            lowerBound = peak_x - intervWidth*stepSize
            upperBound = peak_x + intervWidth*stepSize
        for k in range(lowerBound,upperBound+1,stepSize):
            cumulativeP = cumulativeP + stepSize*distr_y[np.where(distr_x==k)][0]
        intervWidth = intervWidth + 1
#print cumulativeP
uB_raw = upperBound-stepSize
lB_raw = lowerBound+stepSize
#print "upperBound ", upperBound-stepSize
#print "lowerBound ", lowerBound+stepSize

# UPPER/LOWER BOUNDS based on numbers calc'd from Interval.py with S_bestfit = 32700
uB_T12 = np.log(2)*(N_A/m_A)*(1/(lB_raw*pdfCnts))*SumObjects # swap +/- bounds b/c of inverse relationship of S and T12
lB_T12 = np.log(2)*(N_A/m_A)*(1/(uB_raw*pdfCnts))*SumObjects

print "T12+ ", uB_T12 - T12
print "T12- ", T12 - lB_T12

# SHADE THE CONFIDENCE INTERVAL
fill_confint_x = []
fill_confint_y = []
for i in range(0,distr_x.size):
    if(distr_x[i] >= lB_raw and distr_x[i] <= uB_raw):
        fill_confint_x.append(distr_x[i])
        fill_confint_y.append(distr_y[i])
ax3.fill_between(fill_confint_x,fill_confint_y,0,facecolor='k', alpha=0.1)


###############################
# PLOT BEST FIT TO DATA
###############################

fill_fitrange_x = []
fill_fitrange_y = []
for i in range(0,data_x.size):
    if(data_x[i] >= fit_lo and data_x[i] <= fit_hi):
        fill_fitrange_x.append(data_x[i])
        fill_fitrange_y.append(yscale_hi) #fill_fitrange_y.append(data_y[i])
#ax4.fill_between(fill_fitrange_x,fill_fitrange_y,0,facecolor='k', alpha=0.1)
ax4.step(data_x, data_y, color='k')
ax4.step(pdf_x, S_bestfit*pdf_y, color='r')
ax4.set_title('Best Fit')
ax4.set_xlabel('Energy (keV)')
ax4.set_ylabel('Cnt/10keV')
ax4.set_yscale('log',nonposy='clip')
ax4.set_ylim(yscale_lo,yscale_hi)
ax4.set_xlim(xscale_lo,xscale_hi)

ax7.step(data_x, data_y, color='k')
ax7.step(pdf_x, S_bestfit*pdf_y, color='r')
ax7.set_title('Best Fit Scaling to DS3+DS4')
ax7.set_xlabel('Energy (keV)')
ax7.set_ylabel('Cnt/10keV')
ax7.set_yscale('log',nonposy='clip')
ax7.set_ylim(yscale_lo,yscale_hi)
ax7.set_xlim(xscale_lo,xscale_hi)

###############################
# GET RESIDUALS OR DATA/MODEL RATIO
###############################

ratio = data_y/(S_bestfit*pdf_y)
ax5.plot(data_x, ratio, color='k')
ax5.set_xlim(xscale_lo,xscale_hi)
ax5.set_ylim(0,10)
ax5.set_xlabel('Energy (keV)')
ax5.set_ylabel('Ratio Data/Model')

residual = data_y - (S_bestfit*pdf_y)
ax6.plot(data_x, residual, color='k')
ax6.set_xlim(xscale_lo,xscale_hi)
ax6.set_ylim(-50,100)
ax6.set_xlabel('Energy (keV)')
ax6.set_ylabel('Residual Data-Model')

###############################
# GET CONFIDENCE INTERVALS ON DATA, RESIDUALS, OR RATIO
###############################

# For each conf interv, get array of y-values for upper and lower bounds. Then turn those arrays into arrays of ratios relative to the corresponding data y-value
interv_lo = np.zeros(data_x.size)
interv_hi = np.zeros(data_x.size)
interv_68_lo = np.zeros(data_x.size)
interv_68_hi = np.zeros(data_x.size)
interv_95_lo = np.zeros(data_x.size)
interv_95_hi = np.zeros(data_x.size)
interv_995_lo = np.zeros(data_x.size)
interv_995_hi = np.zeros(data_x.size)

# Determine minimal interval of width 1-alpha
expected_y = S_bestfit*pdf_y

# FOR .68 CONFIDENCE
alpha = 1 - 0.68 # 68% central interval
cumulativeP = 0.
for i in range(0,data_x.size):
    mu = expected_y[i]
    dmu = np.floor(mu)
    cumulativeP = 0.
    if(poisson.pmf(dmu,mu) > (1-alpha)): # FOR POISSON DATA AND MODEL
        lowerBound = dmu
        upperBound = dmu
    else:
        intervWidth = 0
        while cumulativeP < (1-alpha):
            cumulativeP = 0
            if(dmu < intervWidth - np.floor(intervWidth/2.)):
                lowerBound = 0
                upperBound = intervWidth
            else:
                lowerBound = int(dmu - np.floor(intervWidth/2.))
                upperBound = int(dmu + np.ceil(intervWidth/2.))
            for k in range(int(lowerBound),int(upperBound+1)):
                cumulativeP = cumulativeP + poisson.pmf(k,mu)
            intervWidth = intervWidth + 1
    interv_68_hi[i] = upperBound
    interv_68_lo[i] = lowerBound
ax7.fill_between(data_x,interv_68_hi,interv_68_lo,facecolor='g', alpha=0.2)


#ax4.plot(data_x, interv_68_hi, color='b')
#ax4.plot(data_x, interv_68_lo, color='g')

# CONFIDENCE INTERVALS FOR RATIOS AND RESIDUALS
#ratio_68_hi = data_y/interv_68_hi
#ratio_68_lo = data_y/interv_68_lo
#ax5.fill_between(data_x,ratio_68_hi,ratio_68_lo,facecolor='k', alpha=0.1)
#ax5.plot(data_x, ratio_68_hi, color='b')
#ax5.plot(data_x, ratio_68_lo, color='g')



# FOR .95 CONFIDENCE
alpha = 1 - 0.95 # 68% central interval
cumulativeP = 0.
for i in range(0,data_x.size):
    mu = expected_y[i]
    dmu = np.floor(mu)
    cumulativeP = 0.
    if(poisson.pmf(dmu,mu) > (1-alpha)):
        lowerBound = dmu
        upperBound = dmu
    else:
        intervWidth = 0
        while cumulativeP < (1-alpha):
            cumulativeP = 0
            if(dmu < intervWidth - np.floor(intervWidth/2.)):
                lowerBound = 0
                upperBound = intervWidth
            else:
                lowerBound = int(dmu - np.floor(intervWidth/2.))
                upperBound = int(dmu + np.ceil(intervWidth/2.))
            for k in range(int(lowerBound),int(upperBound+1)):
                cumulativeP = cumulativeP + poisson.pmf(k,mu)
            intervWidth = intervWidth + 1
    interv_95_hi[i] = upperBound
    interv_95_lo[i] = lowerBound
ax7.fill_between(data_x,interv_95_hi,interv_68_hi,facecolor='y', alpha=0.2)
ax7.fill_between(data_x,interv_68_lo,interv_95_lo,facecolor='y', alpha=0.2)

# FOR .995 CONFIDENCE
alpha = 1 - 0.995 # 68% central interval
cumulativeP = 0.
for i in range(0,data_x.size):
    mu = expected_y[i]
    dmu = np.floor(mu)
    cumulativeP = 0.
    if(poisson.pmf(dmu,mu) > (1-alpha)):
        lowerBound = dmu
        upperBound = dmu
    else:
        intervWidth = 0
        while cumulativeP < (1-alpha):
            cumulativeP = 0
            if(dmu < intervWidth - np.floor(intervWidth/2.)):
                lowerBound = 0
                upperBound = intervWidth
            else:
                lowerBound = int(dmu - np.floor(intervWidth/2.))
                upperBound = int(dmu + np.ceil(intervWidth/2.))
            for k in range(int(lowerBound),int(upperBound+1)):
                cumulativeP = cumulativeP + poisson.pmf(k,mu)
            intervWidth = intervWidth + 1
    interv_995_hi[i] = upperBound
    interv_995_lo[i] = lowerBound
ax7.fill_between(data_x,interv_995_hi,interv_95_hi,facecolor='r', alpha=0.2)
ax7.fill_between(data_x,interv_95_lo,interv_995_lo,facecolor='r', alpha=0.2)

###############################
# SHOW PLOTS
###############################

# SHOW
plt.show()

# Resources
# https://docs.scipy.org/doc/numpy/reference/generated/numpy.loadtxt.html
# https://matplotlib.org/api/pyplot_api.html#matplotlib.pyplot.step
#