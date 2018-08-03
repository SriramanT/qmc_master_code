import numpy as np
import copy
import os
import warnings
cwd = os.getcwd()

# Retrieve simulation parameters

simulation = np.genfromtxt('simulationParameters.txt')

simulationParameters = simulation[:, 1]

NSITES = int(simulationParameters[0])
dt = simulationParameters[1]
beta = simulationParameters[2]
L = int(simulationParameters[3])
t = simulationParameters[4]
U = simulationParameters[5]
mu = simulationParameters[6]
totalMCSweeps = int(simulationParameters[7])
freq = int(simulationParameters[8])
intsize = int(simulationParameters[9])

# Load configuration weights

weights = np.loadtxt('LOGweights.txt', skiprows = 1)

observables = np.loadtxt('measurementsScalars.txt', skiprows = 1)

electronDensity = observables[:, 0]
doubleOc = observables[:, 1]

magCorrMeas = np.loadtxt('SpinSpinCorrelations.txt', skiprows = 1)

#UneqMagCorrMeas = np.loadtxt('UneqTimeSpinSpinCorrelations.txt', skiprows = 1)

directory1 = (str(NSITES) + \
              'sites_L=' + str(L) + \
              '_beta=' + str(beta) + \
              '_dt_' + str(dt) + '_t_' + \
              str(t) + '_U_'+ str(U) + '_mu_' + str(mu))

directory2 = (directory1 + '/data-to-reproduce')


if not os.path.exists(directory1):
  os.makedirs(directory1)

if not os.path.exists(directory2):
    os.makedirs(directory2)

    np.savetxt(directory2 + '/weights_' + \
               'totalMCSweeps_' + str(totalMCSweeps) + \
               '_freq_' + str(freq) + '_intsize_' + str(intsize) + '.txt', (weights))
    np.savetxt(directory2 + '/simulationParameters_' + \
              'totalMCSweeps_' + str(totalMCSweeps) + \
              '_freq_' + str(freq) + '_intsize_' + str(intsize) + '.txt', (simulationParameters))

    np.savetxt(directory2 + '/electronDensity' + str(totalMCSweeps) + \
               '_freq_' + str(freq) + '_intsize_' + str(intsize) + '.txt', (electronDensity))

    np.savetxt(directory2 + '/doubleOc' + str(totalMCSweeps) + \
               '_freq_' + str(freq) + '_intsize_' + str(intsize) + '.txt', (doubleOc))

    np.savetxt(directory2 + '/magCorr' + str(totalMCSweeps) + \
               '_freq_' + str(freq) + '_intsize_' + str(intsize) + '.txt', (magCorrMeas))

#    np.savetxt(directory2 + '/UneqMagCorr' + str(totalMCSweeps) + \
#               '_freq_' + str(freq) + '_intsize_' + str(intsize) + '.txt', (UneqMagCorrMeas))
