import numpy as np
import copy
import os
import warnings
cwd = os.getcwd()

# Retrieve simulation parameters

simulation = np.genfromtxt('../temp-data/simulationParameters.csv', delimiter=',')

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
geom = int(simulationParameters[10])
ny = int(simulationParameters[11])

# Load configuration weights

weights = np.loadtxt('../temp-data/Log-weights.csv', skiprows = 1)

scalarMeasurements = np.loadtxt('../temp-data/MeasurementsScalars.csv', usecols = (1,), delimiter=',')

directory1 = ('data/' + str(NSITES) + \
              'sites_L=' + str(L) + \
              '_beta=' + str(beta) + \
              '_dt_' + str(dt) + '_t_' + \
              str(t) + '_U_'+ str(U) + '_mu_' + str(mu))

directory2 = (directory1 + '/data-to-reproduce/' + \
              'totalMCSweeps_' + str(totalMCSweeps) + \
              '_freq_' + str(freq) + '_intsize_' + str(intsize) + \
              '_geom_' + str(geom) + '_ny_' + str(ny) )

directory3 = (directory1 + '/plots/' + \
              'totalMCSweeps_' + str(totalMCSweeps) + \
              '_freq_' + str(freq) + '_intsize_' + str(intsize) + \
              '_geom_' + str(geom) + '_ny_' + str(ny) )


if not os.path.exists(directory1):
    os.makedirs(directory1)

if not os.path.exists(directory2):
    os.makedirs(directory2)

if not os.path.exists(directory3):
    os.makedirs(directory3)

np.savetxt(directory2 + '/Log-weights.csv', (weights), header = "Configuration log weight")

paramNames  = np.array(['Number of sites,', 'dtau,', 'beta,', 'L,',\
't,', 'U,', 'mu,', 'Number of MC Sweeps,', 'Frequency of recomputing G,',\
'Number of multiplied Bs after stabilization,', 'Geometry,', 'Ny,'])
namedParams = np.zeros(paramNames.size, dtype=[('pnames', 'U50'), ('params', float)])
namedParams['pnames'] = paramNames
namedParams['params'] = simulationParameters

np.savetxt(directory2 + '/simulationParameters.csv', namedParams, fmt="%50s% 10.7f" )

measNames  = np.array(['Electron density <n>,', 'd<n>,', 'Average sign <s>,'])
namedMeas = np.zeros(measNames.size, dtype=[('mnames', 'U50'), ('meas', float)])
namedMeas['mnames'] = measNames
namedMeas['meas'] = scalarMeasurements

np.savetxt(directory2 + '/MeasurementsScalars.csv', namedMeas, fmt="%50s%10.7f" )
