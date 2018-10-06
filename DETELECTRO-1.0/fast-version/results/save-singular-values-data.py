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
Dup = np.loadtxt('../temp-data/Dup.csv', delimiter = ',')
Ddown = np.loadtxt('../temp-data/Ddown.csv', delimiter = ',')
kappas = np.loadtxt('../temp-data/conditionNumbers.csv', delimiter = ',')

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

np.savetxt(directory2 + '/Dup.csv', (Dup), \
header = "Dup")

np.savetxt(directory2 + '/Ddown.csv', (Ddown), \
header = "Ddown")

np.savetxt(directory2 + '/kappaUpDown.csv', (kappas) )
