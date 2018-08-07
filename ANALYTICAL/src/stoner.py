import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import os
cwd = os.getcwd()
sns.set()
sns.set_style("white")
from matplotlib import rcParams
rcParams['font.family'] = 'sans-serif'
rcParams['font.sans-serif'] = ['Arial']
sns.set_palette("GnBu_d")

SAVEDIR = "../plots/stoner"
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)

# Density of states for a 1D tight binding model

def DOS(E):
    return 1 / ( np.pi * np.sqrt( 4 - E**2 ) )

dE = 0.001
Emax = 1.9
Es = np.arange(-Emax, Emax+dE, dE)
Efermi = 1
deltaE = 0.1

fig1, ax1 = plt.subplots()

ax1.plot(Es, DOS(Es), linewidth = 2)
ax1.fill_between(Es, DOS(Emax), DOS(Es)\
                 , where = Es < Efermi, facecolor='aliceblue')
ax1.fill_between(Es[int( (Emax+ Efermi - deltaE) / dE):]\
                 , DOS(Emax), DOS(Es[int( (Emax+ Efermi - deltaE) / dE):]),\
                 where = Es[int( (Emax+ Efermi - deltaE) / dE):] < Efermi + deltaE, facecolor='lightyellow')
ax1.plot(np.ones((Es.size))*Efermi,\
         np.linspace(DOS(Efermi), DOS(Emax), Es.size), color = 'seagreen', linewidth = 2.5)
ax1.plot(np.ones((Es.size))*(Efermi + deltaE),\
         np.linspace(DOS(Efermi + deltaE), DOS(Emax), Es.size), color = 'orange', linewidth = 2.5)
ax1.plot(np.ones((Es.size))*(Efermi - deltaE),\
         np.linspace(DOS(Efermi - deltaE), DOS(Emax), Es.size), color = 'red', linewidth = 2.5)
plt.xlabel(r'$E$')
plt.ylabel(r'$N_\downarrow(E)$')
plt.savefig(SAVEDIR + "/DOSminus.png", dpi = 600)

fig2, ax2 = plt.subplots()

ax2.plot(Es, DOS(Es), linewidth = 2)
ax2.fill_between(Es, DOS(Emax), DOS(Es)\
                 , where = Es < Efermi - deltaE, facecolor='aliceblue')
ax2.plot(np.ones((Es.size))*Efermi,\
         np.linspace(DOS(Efermi), DOS(Emax), Es.size), color = 'seagreen', linewidth = 2.5)
ax2.plot(np.ones((Es.size))*(Efermi + deltaE),\
         np.linspace(DOS(Efermi + deltaE), DOS(Emax), Es.size), color = 'orange', linewidth = 2.5)
ax2.plot(np.ones((Es.size))*(Efermi - deltaE),\
         np.linspace(DOS(Efermi - deltaE), DOS(Emax), Es.size), color = 'red', linewidth = 2.5)
plt.xlabel(r'$E$')
plt.ylabel(r'$N_\uparrow(E)$')
plt.savefig(SAVEDIR + "/DOSplus.png", dpi = 600)


