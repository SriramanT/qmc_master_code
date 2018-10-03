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

scalarMeasurements = np.loadtxt('../temp-data/MeasurementsScalars.csv', usecols = 1, delimiter=',')

magCorrZZMeas = np.loadtxt('../temp-data/EqTimeSzCorrelations.csv', skiprows = 1, delimiter = ',')

magCorrZZMeasError = np.loadtxt('../temp-data/EqTimeSzCorrelationsError.csv', skiprows = 1, delimiter = ',')

#magCorrXXMeas = np.loadtxt('../temp-data/EqTimeSxCorrelations.csv', skiprows = 1, delimiter = ',')
#
#magCorrXXMeasError = np.loadtxt('../temp-data/EqTimeSxCorrelationsError.csv', skiprows = 1, delimiter = ',')

try:
    UneqMagCorrMeasZZ = np.loadtxt('../temp-data/UneqTimeSzCorrelations.csv', skiprows = 1, delimiter = ',')
    UneqMagCorrMeasZZError = np.loadtxt('../temp-data/UneqTimeSzCorrelationsError.csv', skiprows = 1, delimiter = ',')
#    UneqMagCorrMeasXX = np.loadtxt('../temp-data/UneqTimeSxCorrelations.csv', skiprows = 1, delimiter = ',')
#    UneqMagCorrMeasXXError = np.loadtxt('../temp-data/UneqTimeSxCorrelationsError.csv', skiprows = 1, delimiter = ',')
except IOError:
    print("\nExisting data contains only equal time measurements")

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

measNames  = np.array(['Electron density <n>,', 'd<n>,', 'Double occupancy <n+ n->,',\
'd<n+ n->,', 'ZZ AF Structure Factor,', 'Sq.Magnetization <m^2>,', 'd<m^2>,', 'Hkin,', 'dHkin,',\
'Hint,', 'E,', 'Average sign <s>,'])
namedMeas = np.zeros(measNames.size, dtype=[('mnames', 'U50'), ('meas', float)])
namedMeas['mnames'] = measNames
namedMeas['meas'] = scalarMeasurements

np.savetxt(directory2 + '/MeasurementsScalars.csv', namedMeas, fmt="%50s%10.7f" )

np.savetxt(directory2 + '/EqTimeSzCorrelations.csv', (magCorrZZMeas), \
header = "<Sz_i Sz_j >")

np.savetxt(directory2 + '/EqTimeSzCorrelationsError.csv', (magCorrZZMeasError), \
header = "d<Sz_i Sz_j >")

#np.savetxt(directory2 + '/EqTimeSxCorrelations.csv', (magCorrXXMeas), \
#header = "<Sx_i Sx_j >")
#
#np.savetxt(directory2 + '/EqTimeSxCorrelationsError.csv', (magCorrXXMeasError), \
#header = "d<Sx_i Sx_j >")

try:
    np.savetxt(directory2 + '/UneqTimeSzCorrelations.csv', (UneqMagCorrMeasZZ), \
    header = "int_0^beta dt <Sz_i (t) Sz_j (0) >")
    np.savetxt(directory2 + '/UneqTimeSzCorrelationsError.csv', (UneqMagCorrMeasZZError), \
    header = "d int_0^beta dt <Sz_i (t) Sz_j (0) >")
#    np.savetxt(directory2 + '/UneqTimeSxCorrelations.csv', (UneqMagCorrMeasXX), \
#    header = "int_0^beta dt <Sx_i (t) Sx_j (0) >")
#    np.savetxt(directory2 + '/UneqTimeSxCorrelationsError.csv', (UneqMagCorrMeasXXError), \
#    header = "d int_0^beta dt <Sx_i (t) Sx_j (0) >")
except NameError:
    print("\nIf you want unequal time measurements as well recompile the code \
and run the simulation again using\n\nmake clean\n\n\
make nsites=<Number of sites> dt_inv=<Trotter error> beta=<Inverse temperature> green_afresh_freq=<Frequency of Recomputing G from scratch> eq_or_uneq=src/mainUneqTime.cpp object\
=src/mainUneqTime.o\n\nand \n\n./simulation <t> <U> <mu> <geom> <Ny> <Total Number of Sweeps (Space-Time)> <Number of Warm-up Sweeps (Space-Time)>  <Number of Auto-correlation Sweeps (Space-Time)>\n")
