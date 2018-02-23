"""

L ~ Poisson(data|scale*sim)
scale ~ Flat()
sim = 2nbb theoretical spectrum, normalized
data = DS52 w/ cuts, in terms of counts

Beware having "zero" terms in the likelihood product if the fit range crosses into that territory ...

"""

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, rcParams
import corner
#%matplotlib inline # is this a jupyter command?

###############
# Get data and sim
###############

# READ IN DATA FILE
dataFilePath = "./text_30keVBins/DS52cutScheme22.txt"
data_x, data_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

# SET FIT RANGE
fit_lo = 300#500 # energy fit range limits
fit_hi = 2630#1735 #1935 # would choose 2039, but not enough stats
fit_index_lo = 0
fit_index_hi = data_x.size - 1

for i in range(0,data_x.size):
    if data_x[i] < fit_lo + 15 and data_x[i] >= fit_lo - 15:
        fit_index_lo = i
    if data_x[i] <= fit_hi + 15 and data_x[i] > fit_hi - 15:
        fit_index_hi = i

# READ IN SIM FILES
dataFilePath = "./text_30keVBins/Sim2vDS52cutScheme1Det2.txt"
pdf_2v_x, pdf_2v_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

dataFilePath = "./text_30keVBins/SimThDS52cutScheme1Det2.txt"
pdf_Th_x, pdf_Th_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

dataFilePath = "./text_30keVBins/SimUDS52cutScheme1Det2.txt"
pdf_U_x, pdf_U_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

dataFilePath = "./text_30keVBins/SimKDS52cutScheme1Det2.txt"
pdf_K_x, pdf_K_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

dataFilePath = "./text_30keVBins/SimCoDS52cutScheme1Det2.txt"
pdf_Co_x, pdf_Co_y = np.loadtxt(dataFilePath, delimiter=" ", usecols=(0,1), unpack=True)

# SANITY CHECKS
print("fit index range: ", fit_index_lo, fit_index_hi)
print("fit energy range: ", data_x[fit_index_lo], data_x[fit_index_hi])
print("checking some data and sim shapes:")
print("shape data_y ", data_y.shape)
print("shape pdf_2v_y ", pdf_2v_y.shape)
print("shape pdf_Th_y ", pdf_Th_y.shape)


###############
# Specifiy model
###############

#### Can implement fit range via indexing in the L definition, i.e. sim_y[indexFor500:indexFor2039+1]
#fit_lo = 500 # energy fit range limits
#fit_hi = 2039
#
#scale_lo = 20000.
#scale_hi = 40000.

from pymc3 import Model, Normal, Poisson, Uniform

# Create instance of Model class
basic_model = Model()

# Specify model components
with basic_model:
    # Priors for unknown model parameters (Stochastic random vars)
    s = Uniform('s', lower=0, upper=200000, shape=5) # scaling of bg component
    # Expected value of outcome (Deterministic var)
    sim_2v_y = s[0]*pdf_2v_y # scaled bg component
    sim_Th_y = s[1]*pdf_Th_y
    sim_U_y = s[2]*pdf_U_y
    sim_K_y = s[3]*pdf_K_y
    sim_Co_y = s[4]*pdf_Co_y
    model_y = sim_2v_y + sim_Th_y + sim_U_y + sim_K_y + sim_Co_y
    # Likelihood (sampling distribution) of observations (Observed Stochastic var)
    L = Poisson('L', mu=model_y[fit_index_lo:fit_index_hi], observed=data_y[fit_index_lo:fit_index_hi])

###############
# Fit model and get posterior estimates for parameters
###############

# Import a sampler
from pymc3 import Metropolis, HamiltonianMC, sample

# Setup sampler within the context of the model
with basic_model:
    # Set some starting value guesses
    start = {'s': [0.,0.,0.,0.,0.]}
    # Instantiate sampler
    step = HamiltonianMC([s])
    # Draw 2000 posterior samples
    trace = sample(draws=30000, step=step, start=start)
    #trace = sample(draws=30000, step=step)

from pymc3 import summary, traceplot, autocorrplot

summary(trace)
traceplot(trace)
autocorrplot(trace)

burnin = 20000
print("trace shape:", trace['s'].shape)
print("trace shape post-burnin: ", trace['s'][burnin:,:].shape)
s_trace_postburn = trace['s'][burnin:,:]
print("trace shape post-burnin copy: ", s_trace_postburn.shape)
s_mean = np.mean(s_trace_postburn, axis=0)
print("s_mean ",s_mean)

# CORNER PLOTS
import corner

s_trace_postburn_vstack = np.vstack([s_trace_postburn[:,0], s_trace_postburn[:,1], s_trace_postburn[:,2], s_trace_postburn[:,3], s_trace_postburn[:,4]])
print("s_trace_postburn_vstack shape: ", s_trace_postburn_vstack.shape)

figure = corner.corner(s_trace_postburn,
                       labels=[r"$2v$", r"$Th$", r"$U$", r"$K$", r"$Co$"],
                       quantiles=[0.16, 0.5, 0.84],
                       show_titles=True,
                       title_kwargs={"fontsize": 12})
#figure.savefig("corner.png", dpi=300)

###############################
# PLOT TRACES
###############################

fig_traces, (ax_tr0, ax_tr1, ax_tr2, ax_tr3, ax_tr4) = plt.subplots(5,1)
fig_traces.subplots_adjust(top=None,bottom=None,hspace=.45)
ax_tr_x = np.linspace(0, 999, 1000)
print("ax_tr_x shape: ", ax_tr_x.shape)
ax_tr0.scatter(ax_tr_x, trace['s'][0:1000,0], s=10, facecolors='none', edgecolors='k')
ax_tr1.scatter(ax_tr_x, trace['s'][0:1000,1], s=10, facecolors='none', edgecolors='k')
ax_tr2.scatter(ax_tr_x, trace['s'][0:1000,2], s=10, facecolors='none', edgecolors='k')
ax_tr3.scatter(ax_tr_x, trace['s'][0:1000,3], s=10, facecolors='none', edgecolors='k')
ax_tr4.scatter(ax_tr_x, trace['s'][0:1000,4], s=10, facecolors='none', edgecolors='k')
print("s0 trace: ", trace['s'][0:20,0])
print("s1 trace: ", trace['s'][0:20,1])
print("s2 trace: ", trace['s'][0:20,2])
print("s3 trace: ", trace['s'][0:20,3])
print("s4 trace: ", trace['s'][0:20,4])


###############
# Plots
###############

xlo = 0.
xhi = 3000.
wid = 10.
nbins = (xhi-xlo)/wid
xscale_lo = 0
xscale_hi = 3000
yscale_lo = 0.5
yscale_hi = 5e3

#fig, (ax3, ax2) = plt.subplots(2, 1)
#fig.subplots_adjust(top=None,bottom=None,hspace=.41)

#fig = plt.figure(figsize=(6, 8))
#gs = gridspec.GridSpec(2, 1, width_ratios=[3, 1])
#ax3 = plt.subplot(gs[0])
#ax2 = plt.subplot(gs[1])

fig, (ax3, ax2) = plt.subplots(2,1, gridspec_kw = {'height_ratios':[2, 1]})
fig.subplots_adjust(top=None,bottom=None,hspace=.45)

# PLOT DATA
#ax1.step(data_x, data_y, color='k') #ax1.bar(data_x, data_y, width=wid, color='k', fill=False)
#ax1.set_title('MJD Data, DS5b')
#ax1.set_xlabel('Energy (keV)')
#ax1.set_ylabel('Cnt/30keV')
#ax1.set_yscale('log',nonposy='clip')
##ax1.set_ylim(yscale_lo,yscale_hi)
#ax1.set_xlim(xscale_lo,xscale_hi)

# PLOT SIM
ax2.step(pdf_2v_x, pdf_2v_y, color='b')
ax2.step(pdf_Th_x, pdf_Th_y, color='tab:orange')
ax2.step(pdf_U_x, pdf_U_y, color='g')
ax2.step(pdf_K_x, pdf_K_y, color='r')
ax2.step(pdf_Co_x, pdf_Co_y, color='m')
ax2.set_title('Sim PDFs')
ax2.set_xlabel('Energy (keV)')
ax2.set_ylabel('Cnt/30keV')
ax2.set_yscale('log',nonposy='clip')
#ax2.set_ylim(yscale_lo,yscale_hi)
ax2.set_xlim(xscale_lo,xscale_hi)

###############################
# PLOT BEST FIT TO DATA
###############################

fill_fitrange_x = []
fill_fitrange_y = []
fill_fitrange_yy = []
for i in range(0,data_x.size):
    if(data_x[i] >= fit_lo and data_x[i] <= fit_hi):
        fill_fitrange_x.append(data_x[i])
        fill_fitrange_y.append(yscale_hi)
        fill_fitrange_yy.append(data_y[i])
        #fill_fitrange_y.append(np.amax(data_y))#fill_fitrange_y.append(data_y[i])
ax3.fill_between(fill_fitrange_x,fill_fitrange_y,0,facecolor='k', alpha=0.08) # alpha=.1
#ax3.fill_between(fill_fitrange_x,fill_fitrange_yy,fill_fitrange_y,facecolor='k', alpha=0.08)
#ax3.fill_between(data_x,data_y,0,facecolor='k', alpha=0.08)
#ax3.step(data_x, data_y, color='k')
ax3.scatter(data_x-15, data_y, s=10, facecolors='none', edgecolors='k') #c='k', marker='o')
ax3.step(pdf_2v_x, s_mean[0]*pdf_2v_y, color='b', linewidth=1.0)
ax3.step(pdf_Th_x, s_mean[1]*pdf_Th_y, color='tab:orange', linewidth=1.0)
ax3.step(pdf_U_x, s_mean[2]*pdf_U_y, color='g', linewidth=1.0)
ax3.step(pdf_K_x, s_mean[3]*pdf_K_y, color='r', linewidth=1.0)
ax3.step(pdf_Co_x, s_mean[4]*pdf_Co_y, color='m', linewidth=1.0)
all_y = s_mean[0]*pdf_2v_y + s_mean[1]*pdf_Th_y + s_mean[2]*pdf_U_y + s_mean[3]*pdf_K_y + s_mean[4]*pdf_Co_y
ax3.step(pdf_Th_x, all_y, color='c')
ax3.set_title('DS5b, Fit')
#ax3.set_xlabel('Energy (keV)')
ax3.set_ylabel('Cnt/30keV')
ax3.set_yscale('log',nonposy='clip')
ax3.set_ylim(yscale_lo,yscale_hi)
ax3.set_xlim(xscale_lo,xscale_hi)

#ratio = data_y/all_y
#ax5.plot(data_x, ratio, color='k')
#ax5.set_xlim(xscale_lo,xscale_hi)
#ax5.set_ylim(0,10)
#ax5.set_xlabel('Energy (keV)')
#ax5.set_ylabel('Ratio Data/Model')

###############################
# PRINT T1/2
###############################

pdfCnts = 1.0

# GENERAL
N_A = 6.022140857e23 # (1/mol)
m_A = 75.6e-3 # (kg/mol) molar mass of 76Ge

# FROM GERDA
#T12 = 1.926e21 # (yr); GERDA: Eur. Phys. J. C, 75:416, 2015

# FROM MJD
e_av = 0.8
exposureDS5bEnr = (492.149+182.187)/365.25 # kg-yr
exposureDS5bNat = (138.285+197.6)/365.25 # kg-yr

# CALCULATE
tot_Exp_Eff = (0.88)*exposureDS5bEnr + (0.078)*exposureDS5bNat
N_obs = s_mean[0]*pdfCnts
T12 = np.log(2)*(N_A/m_A)*(1/N_obs)*tot_Exp_Eff
print ("T12 ", T12)

###############
# Show figures
###############
plt.show()
