import numpy as np
import matplotlib.pyplot as plt
import os
cwd = os.getcwd()
import seaborn as sns
sns.set()
sns.set_style("whitegrid")
sns.set_palette("GnBu_d")
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']

SAVEDIR = "../plots/single-site-hubbard"
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# $\rho(\mu)$
### Electron density $\rho$ for varying chemical potential
### $\mu$ and temperature $T = \beta^{-1}$, but fixed $U = 4$.
### As the temperature decreases, a Mott plateau sets in.
### The Mott insulating gap already seen here is an important feature of the Hubbard model.

# Set the parameters of the Hubbard model and temperature

U = 4
mu = np.arange(-6, 6, 0.001)
rho = np.zeros(np.size(mu))

t = np.array([0.25, 0.5, 0.75, 1., 1.5, 2.])

fig1 = plt.figure(1)
ax1 = fig1.add_subplot(111)
plt.xlabel(r'$\mu$')
plt.ylabel(r'$\rho$')

for T in t:
    for r in range(len(rho)):
        rho[r] = 2 * ( np.exp( (U/4 + mu[r]) / T )
                    + np.exp( (2*mu[r] - U/4) / T ) ) / ( np.exp(- U / 4 / T)
                      + 2 * np.exp( (U/4 + mu[r]) / T ) + np.exp( (2*mu[r] - U/4) / T ) )
    ax1.plot(mu, rho, label = 'T = ' + str(T), linewidth = 0.5)

    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #       ncol=2, mode="expand", borderaxespad=0.)
    lgd = ax1.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(SAVEDIR + "/rhoVsMu.png", dpi = 600, bbox_extra_artists=(lgd,), bbox_inches='tight')

# $\left\langle m^2 \right\rangle (U)$

t = np.array([0.25, 0.5, 0.75, 1., 1.5, 2.])

fig2 = plt.figure(2)
ax2 = fig2.add_subplot(111)

mu = 0
U = np.arange(0, 12, 0.001)
mSq = np.zeros(np.size(U))

plt.ylabel(r'$\left\langle m^2 \right\rangle$')
plt.xlabel(r'$U$')

for T in t:
    for u in range(len(U)):
        mSq[u] = 2 * ( np.exp( (U[u]/4 + mu) / T )\
        + np.exp( (2*mu - U[u]/4) / T ) ) \
        / ( np.exp(- U[u] / 4 / T )+ 2 * np.exp( (U[u] /4 + mu) / T ) + np.exp( (2*mu - U[u]/4) / T ) ) \
        - 2 * np.exp( (2*mu - U[u]/4) / T ) \
        / ( np.exp(- U[u] / 4 / T )+ 2 * np.exp( (U[u] /4 + mu) / T ) + np.exp( (2*mu - U[u]/4) / T ) )
    ax2.plot(U, mSq, label = 'T = ' + str(T), linewidth = 0.5)

    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #       ncol=2, mode="expand", borderaxespad=0.)
    lgd = ax2.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)

plt.savefig(SAVEDIR + "/mSqVsU.png", dpi = 600, bbox_extra_artists=(lgd,), bbox_inches='tight')

# $\left\langle m^2 \right\rangle (T)$ for varying $\mu$

U = 4
mu_s = np.arange(0, 7, 1)
T = np.arange(0.1, 6, 0.001)
mSq = np.zeros(np.size(T))

fig3 = plt.figure(3)
ax3 = fig3.add_subplot(111)

plt.ylabel(r'$\left\langle m^2 \right\rangle$')
plt.xlabel(r'$T$')

for mu in mu_s:
    for t in range(len(T)):
        mSq[t] = 2 * ( np.exp( (U/4 + mu) / T[t] )\
        + np.exp( (2*mu - U/4) / T[t] ) ) \
        / ( np.exp(- U / 4 / T[t] )+ 2 * np.exp( (U /4 + mu) / T[t] ) + np.exp( (2*mu - U/4) / T[t] ) ) \
        - 2 * np.exp( (2*mu - U/4) / T[t] ) \
        / ( np.exp(- U / 4 / T[t] )+ 2 * np.exp( (U /4 + mu) / T[t] ) + np.exp( (2*mu - U/4) / T[t] ) )
    if mu != 0:
        ax3.plot(T, mSq, label = r'$\mu = \pm$' + str(mu), linewidth = 0.5)
    else:
        ax3.plot(T, mSq, label = r'$\mu$ = ' + str(mu), linewidth = 0.5)

    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #       ncol=2, mode="expand", borderaxespad=0.)
    lgd = ax3.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(SAVEDIR + "/mSqVsT.png", dpi = 600, bbox_extra_artists=(lgd,), bbox_inches='tight')

# $\left\langle m^2 \right\rangle (T)$ for varying $U$

u = np.arange(1, 14, 2)
mu = 0
T = np.arange(0.005, 6, 0.005)
mSq = np.zeros(np.size(T))

fig4 = plt.figure(4)
ax4 = fig4.add_subplot(111)

plt.ylabel(r'$\left\langle m^2 \right\rangle$')
plt.xlabel(r'$T$')

for U in u:
    for t in range(len(T)):
        mSq[t] = 2 * ( np.exp( (U/4 + mu) / T[t] )\
        + np.exp( (2*mu - U/4) / T[t] ) ) \
        / ( np.exp(- U / 4 / T[t] )+ 2 * np.exp( (U /4 + mu) / T[t] ) + np.exp( (2*mu - U/4) / T[t] ) ) \
        - 2 * np.exp( (2*mu - U/4) / T[t] ) \
        / ( np.exp(- U / 4 / T[t] )+ 2 * np.exp( (U /4 + mu) / T[t] ) + np.exp( (2*mu - U/4) / T[t] ) )
    ax4.plot(T, mSq, label = r'$U$ = ' + str(U), linewidth = 0.5)

    #plt.legend(bbox_to_anchor=(0., 1.02, 1., .102), loc=3,
    #       ncol=2, mode="expand", borderaxespad=0.)
    lgd = ax4.legend(bbox_to_anchor=(1.05, 1), loc=2, borderaxespad=0.)
plt.savefig(SAVEDIR + "/mSqVsT_and_U.png", dpi = 600, bbox_extra_artists=(lgd,), bbox_inches='tight')
