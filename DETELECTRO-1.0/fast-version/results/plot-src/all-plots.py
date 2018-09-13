# Determinant QMC data visualization

import numpy as np
import matplotlib.pyplot as plt
import os
import warnings
cwd = os.getcwd()
import seaborn as sns
sns.set()
sns.set_palette("Blues_r")
sns.set_style("white")
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

## Are you testing or making plots?

SAVEPLOTS = True

NSITES = 2
dt = 0.015625
beta = 2.
L = 128
t = 1.
U = 4.
mu = 0.
totalMCSweeps = 5000
freq = 4
intsize = 32
geom = 1
ny = 0

mainDir = ('../data/' + str(NSITES) + \
           'sites_L=' + str(L) + \
           '_beta=' + str(beta) + \
           '_dt_' + str(dt) + '_t_' + \
           str(t) + '_U_'+ str(U) + '_mu_' + str(mu))

simDir = (mainDir + '/data-to-reproduce/' + \
          'totalMCSweeps_' + str(totalMCSweeps) + \
          '_freq_' + str(freq) + '_intsize_' + str(intsize) + \
          '_geom_' + str(geom) + '_ny_' + str(ny) )

plotDir = (mainDir + '/plots/' + \
           'totalMCSweeps_' + str(totalMCSweeps) + \
           '_freq_' + str(freq) + '_intsize_' + str(intsize) + \
           '_geom_' + str(geom) + '_ny_' + str(ny) )

# Load weights to plot

weights = np.loadtxt(simDir + '/Log-weights.csv')

nLatSweeps = weights.size

latSweeps = np.arange(nLatSweeps) + 1  #measured in lattice sweeps

# Metropolis Sampling convergence

plt.figure(1)
plt.scatter(latSweeps / L, weights, s = 0.3) #show time in space-time sweeps
plt.xlabel(r"Space-time sweep")
plt.ylabel(r'$\log \frac{| P(\mathbf{h}) | }{ | P(\mathbf{h_0}) | } $')
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/Log-weights.png', dpi = 600)

# Sign problem

sweeps = np.arange(totalMCSweeps) + 1
signs = np.loadtxt(simDir + '/Local-av-sign.csv')
avSign = np.mean(signs) * np.ones(len(sweeps))

fig = plt.figure(2)
ax = fig.add_subplot(111)
plt.xlabel(r'Space-time sweep')
plt.ylabel(r'$\left\langle sign [P(\mathbf{h})] \right\rangle$ ')
ax.scatter(sweeps, signs, s = 1, color = "#34495e", label = 'Sign of statistical weight')
ax.plot(sweeps, avSign, linewidth = 1, color = "#e74c3c", label = 'Average sign')
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/Local-av-sign.png', dpi = 600,
                bbox_extra_artists=(lgd,), bbox_inches='tight')

# Measurements

## Electron density

electronDensity = np.loadtxt(simDir + '/electronDensity.csv')
plt.figure(3)
plt.scatter(sweeps, electronDensity, s = 1)
plt.xlabel(r'Space-time sweep')
plt.ylabel(r'$\rho$')
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/electronDensity.png', dpi = 600)

## Double occupancy

doubleOc = np.loadtxt(simDir + '/doubleOccupancy.csv')
plt.figure(4)
plt.scatter(sweeps, doubleOc, s = 2)
plt.xlabel(r'Space-time sweep')
plt.ylabel(r'$\left\langle \hat{n}_\uparrow \hat{n}_\downarrow \right\rangle$')
if SAVEPLOTS == 1:
    plt.savefig(plotDir + '/doubleOccupancy.png', dpi = 600)

## Auto-correlation time for double occupancy $\left\langle n_\uparrow n_\downarrow \right \rangle$

W = 100
m1 = np.mean(doubleOc[W:])
tMax = 500
chi_nUp_nDown = np.zeros(tMax)

for t in range(tMax):
    chi_nUp_nDown[t] = np.sum( ( doubleOc[W:tMax + W] - m1 ) \
                              * ( doubleOc[t + W:t + tMax + W]\
                                 - np.mean(doubleOc[t + W:] ) ) )

plt.figure(5)
plt.scatter(np.arange(tMax) + W, chi_nUp_nDown, s = 2)
plt.xlabel(r'Space-time sweep')
plt.ylabel(r'$\chi_{_{\left\langle \hat{n}_\uparrow \hat{n}_\downarrow \right\rangle}}$')
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/doubleOcAutoCorr.png', dpi = 600)

## Average double occupancy and error (variance)

A = 4 # auto-correlation time in sweeps
meanDoubleOc = 0
meanDoubleOcSq = 0
for idx, D in enumerate(doubleOc[W:]):
    if idx % A == 0:
        meanDoubleOc += D
        meanDoubleOcSq += D**2

meanDoubleOc *= A / ( totalMCSweeps - W )
meanDoubleOcSq *= A / ( totalMCSweeps - W )
varDoubleOc = np.sqrt( (meanDoubleOcSq - meanDoubleOc**2) / (NSITES - 1) )
print('<n_up n_dw>: ', meanDoubleOc)
print('U <n_up n_dw>: ', U * meanDoubleOc )
print('variance: ', varDoubleOc)

# Magnetic structure factor $S(\mathbf q) = \frac{1}{N} \sum_{i, j} e^{i \mathbf q \cdot  (\mathbf i - \mathbf j)} \left\langle \mathbf S_{\mathbf i} \cdot \mathbf S_{\mathbf j} \right\rangle $

## Auto-correlation time in the measurement of the correlation function

magCorrMeas = np.loadtxt(simDir + '/EqTimeSzCorrelations.csv')

tMax = 100
chiMag = np.zeros(tMax)
corrZeroZero = np.zeros(totalMCSweeps)
for m in range(totalMCSweeps - W):
    corrZeroZero[m] = magCorrMeas[(m + W)*NSITES:(m + W + 1)*NSITES][0, 0]

meanCorr = np.mean(corrZeroZero[W:])

for t in range(tMax):
    chiMag[t] = np.sum( ( corrZeroZero[W:tMax + W] - meanCorr )\
                       * ( corrZeroZero[t + W:t + tMax +W]\
                          - np.mean(corrZeroZero[t + W:]) ) )

plt.figure(6)
plt.scatter(np.arange(tMax) + W, chiMag, s = 2)
plt.xlabel(r'Space-time sweep')
plt.ylabel(r'$\chi_{_{<\mathbf{S}_0^{\quad 2}>}}$')

## Average spin-spin correlation function

magCorr = np.zeros((NSITES, NSITES))
for m in range(totalMCSweeps - W):
    if m % A == 0:
        magCorr += ( magCorrMeas[(m+W)*NSITES:(m + 1+W)*NSITES] \
                    - magCorr ) / ( ( m + 1) / A )

plt.figure(7)
plt.scatter( NSITES - np.arange(NSITES), magCorr[int(NSITES/2), :], s = 2, marker = 'o')
plt.plot( NSITES - np.arange(NSITES), magCorr[int(NSITES/2), :], linewidth = 0.5)
plt.xlabel(r'$j$')
plt.ylabel(r'$\left\langle S_i S_j \right\rangle$, $i = $' + str(int(NSITES/2)))
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/magCorr.png', dpi = 600)

## Fourier transform to obtain the structure factor

n_qs = 500
qMax = 2*np.pi
qs = np.arange(0, qMax + qMax/n_qs, qMax/n_qs)
S = np.zeros(n_qs+1)
dist = 0

for idx, q in enumerate(qs):
    for x in range(NSITES):
        for y in range(NSITES):
            if ( x - y ) > NSITES/2:
                dist = NSITES - abs(x - y)
            else:
                dist = abs(x - y)
            S[idx] += np.cos( dist * q ) * magCorr[x, y]

S /=  NSITES # factor of 16 comes from 1/2 spins and overcounting

plt.figure(8)
plt.scatter(qs/2/np.pi, S, s = 1.2, color = 'salmon')
plt.plot(qs/2/np.pi, S, linewidth = 0.8)
plt.xlabel(r'$q$')
plt.ylabel(r'$S(q)$')
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/S(q).png', dpi = 600)

### It has a maximum at $q = \pi$ signaling quasi-AF order

print("S has a maximum at q = pi signaling quasi-AF order", qs[np.argmax(S)])

## Compare with the results for 60 sites obtained by Yi et al.

s_compare = np.loadtxt("../data/Yi1997/s_compare.csv")
fig = plt.figure(9)
ax = fig.add_subplot(111)
plt.xlabel(r'$q$')
plt.ylabel(r'$S(q)$')
ax.scatter(s_compare[:, 0], s_compare[:, 1], s = 4, color = 'salmon', marker = 'x', label = 'H.Yi et al. (1997)')
ax.plot(qs/2/np.pi, S, linewidth = 0.5, label = 'This work')
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/s_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 600)

# Magnetic susceptibility $\chi(\mathbf q) = \frac{1}{N} \sum_{i, j}
#e^{i \mathbf q \cdot  (\mathbf i - \mathbf j)} \int_0^\beta \left\langle \mathbf
#S_{\mathbf i}(\tau) \cdot \mathbf S_{\mathbf j}(0) \right\rangle d\tau $

UneqMagCorrMeas = np.loadtxt(simDir + '/UneqTimeSzCorrelations.csv', skiprows = 1)

UneqMagCorr = np.zeros((NSITES, NSITES))

# Average unequal time spin-spin correlation function

for m in range(totalMCSweeps - W):
    if m % 100 == 0:
        UneqMagCorr += ( UneqMagCorrMeas[(m + W)*NSITES:(m + W + 1)*NSITES]\
                        - UneqMagCorr ) / ( ( m + 1) / 100 )

# Compute structure factor

n_qs = 500
qMax = 2*np.pi
threshold = 0
qs = np.arange(threshold * qMax/n_qs, qMax - (threshold - 1) * qMax/n_qs , qMax/n_qs)
MagSus = np.zeros(n_qs - 2 * threshold + 1)

for idx, q in enumerate(qs):
    for x in range(int(np.sqrt(NSITES))):
        for y in range(int(np.sqrt(NSITES))):
            if ( x - y ) > NSITES/2:
                dist = NSITES - (x - y)
            else:
                dist = x - y
            MagSus[idx] += np.cos( dist * q) * UneqMagCorr[x, y]

MagSus /= (NSITES) / dt # factor of 4 comes from 1/2 spins

plt.figure(10)
plt.scatter(qs / 2 / np.pi, MagSus, s = 1.2, color = 'salmon')
plt.plot(qs / 2 / np.pi, MagSus, linewidth = 0.8)

plt.xlabel(r'q')
plt.ylabel(r'$\chi(q)$')
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/chi(q).png', bbox_inches='tight', dpi = 600)

## Compare with the results for 60 sites obtained by Yi et al.

chi_compare = np.loadtxt("../data/Yi1997/chi_compare.csv")
fig = plt.figure(11)
ax = fig.add_subplot(111)
plt.xlabel(r'$q$')
plt.ylabel(r'$\chi(q)$')
ax.scatter(chi_compare[:, 0], chi_compare[:, 1], s = 4, color = 'salmon', marker = 'x', label = 'H.Yi et al. (1997)')
ax.plot(qs/2/np.pi, MagSus, linewidth = 0.5, label = 'This work')
lgd = ax.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/chi_compare.png', bbox_extra_artists=(lgd,), bbox_inches='tight', dpi = 600)

# Again, it has a maximum at $q = \pi$ signaling quasi-AF order

threshold =30
print("Chi has a maximum at q = pi signaling quasi-AF order", qs[np.argmax(MagSus[threshold:n_qs-threshold]) + threshold])

## Auto-correlation time in the measurement of the unequal-time correlation function

autoCorr = np.zeros(tMax)
Uneq = np.zeros(totalMCSweeps-W)
for m in range(totalMCSweeps - W):
    Uneq[m] = UneqMagCorrMeas[(m + W)*NSITES:(m + W + 1)*NSITES][0, 1]
av = np.mean(Uneq[W:])
for t in range(tMax):
    autoCorr[t] = np.sum( ( Uneq[W:tMax+W] - av ) * ( Uneq[t+W:t + tMax+W] - np.mean(Uneq[t+W:]) ) )
plt.figure(12)
plt.scatter(np.arange(tMax), autoCorr, s = 1.5)
plt.xlabel(r'Space-time sweep')
plt.ylabel(r'$\chi_{_{<\mathbf{S}_0 (\tau) \cdot \mathbf{S}_1 (0) >}}$')
if SAVEPLOTS == True:
    plt.savefig(plotDir + '/chiAutoCorr.png', dpi = 600)

