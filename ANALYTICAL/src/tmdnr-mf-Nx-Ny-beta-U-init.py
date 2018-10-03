
# coding: utf-8

# In[19]:


import numpy as np
import numpy.linalg as la
import os
import copy
import warnings
import scipy.optimize as opt


# # Save data to local files or from hard drive?

# In[20]:


localData = True


# # Define the model

# In[28]:


nOrb = 3

nHole = 0

tmd = 'MoS_2'


if tmd == 'MoS_2' :

    abs_t0 = 0.184

    e1 = 1.046 / abs_t0
    e2 = 2.104 / abs_t0
    t0 = - 1
    t1 = 0.401 / abs_t0
    t2 = 0.507 / abs_t0
    t11 = 0.218 / abs_t0
    t12 = 0.338 / abs_t0
    t22 = 0.057 / abs_t0

E0 = np.array([[e1, 0, 0],
               [0, e2, 0],
               [0, 0, e2]])

E1 = np.array([[t0, t1, t2],
               [-t1, t11, t12],
               [t2, -t12, t22]])

E4 = np.array([[t0, -t1, t2],
               [t1, t11, -t12],
               [t2, t12, t22]])

E2 = np.array([[t0, 0.5 * t1 - np.sqrt(3) / 2 * t2, - np.sqrt(3) / 2 * t1 - 0.5 * t2],
               [-0.5 * t1 - np.sqrt(3) / 2 * t2, 0.25 * ( t11 + 3 * t22 ), np.sqrt(3) / 4 * ( t22 - t11 ) - t12],
               [np.sqrt(3) / 2 * t1 - 0.5 * t2, np.sqrt(3) / 4 * ( t22 - t11 ) + t12, ( 3 * t11 + t22) / 4 ]])

E5 = np.array([[t0, - 0.5 * t1 - np.sqrt(3) / 2 * t2, np.sqrt(3) / 2 * t1 - 0.5 * t2],
               [0.5 * t1 - np.sqrt(3) / 2 * t2, 0.25 * ( t11 + 3 * t22 ), np.sqrt(3) / 4 * ( t22 - t11 ) + t12],
               [-np.sqrt(3) / 2 * t1 - 0.5 * t2, np.sqrt(3) / 4 * ( t22 - t11 ) - t12, ( 3 * t11 + t22) / 4 ]])

E3 = np.array([[t0, - 0.5 * t1 + np.sqrt(3) / 2 * t2, -np.sqrt(3) / 2 * t1 - 0.5 * t2],
               [0.5 * t1 + np.sqrt(3) / 2 * t2, 0.25 * ( t11 + 3 * t22 ), -np.sqrt(3) / 4 * ( t22 - t11 ) + t12],
               [np.sqrt(3) / 2 * t1 - 0.5 * t2, -np.sqrt(3) / 4 * ( t22 - t11 ) - t12, ( 3 * t11 + t22) / 4 ]])

E6 = np.array([[t0, 0.5 * t1 + np.sqrt(3) / 2 * t2, np.sqrt(3) / 2 * t1 - 0.5 * t2],
               [-0.5 * t1 + np.sqrt(3) / 2 * t2, 0.25 * ( t11 + 3 * t22 ), -np.sqrt(3) / 4 * ( t22 - t11 ) - t12],
               [-np.sqrt(3) / 2 * t1 - 0.5 * t2, -np.sqrt(3) / 4 * ( t22 - t11 ) + t12, ( 3 * t11 + t22) / 4 ]])

hoppings = np.array([E0, E1, E2, E3, E4, E5, E6])

inftyCutOff = 50 # above it, beta is practically infinity

def fermi(e, mu, beta):
    '''
        For zero temperature, set beta = inftyCutOff
        '''
    if beta == inftyCutOff:
        return (e < mu).astype(int)
    else:
        return 1 / ( 1 + np.exp( beta * ( e - mu ) ) )

def Hribbon(k, Ny):
    Hrib = np.zeros((3 * Ny, 3 * Ny), dtype=np.complex64)
    
    h1 = np.array([[e1 + 2 * t0 * np.cos(k),
                2.j * np.sin(k) * t1,
                2 * t2 * np.cos(k)],
               
               [-2.j * np.sin(k) * t1,
                e2 + 2 * t11 * np.cos(k),
                2.j * np.sin(k) * t12],
               
               [2 * t2 * np.cos(k),
                -2.j * np.sin(k) * t12,
                e2 + 2 * t22 * np.cos(k)]
               
               ], dtype=np.complex64)
        
    h2 = np.array([
              
          [ 2 * t0 * np.cos(k/2) ,
           1.j * np.sin(k/2) * ( t1 - np.sqrt(3) * t2 ) ,
           -1. * np.cos(k/2) * ( np.sqrt(3) * t1 + t2 )] ,
          
          [ -1.j * np.sin(k/2) * ( t1 + np.sqrt(3) * t2 ),
           0.5 * np.cos(k/2) * ( t11 + 3 * t22 ),
           1.j * np.sin(k/2) * ( np.sqrt(3) / 2 * ( t22 - t11 ) - 2 * t12 ) ],
          
          [ np.cos(k/2) * ( np.sqrt(3) * t1 - t2 ),
           1.j * np.sin(k/2) * ( np.sqrt(3)/2 * ( t22 - t11 ) + 2 * t12 ),
           0.5 * np.cos(k/2) * ( 3 * t11 + t22 ) ]
          
          ], dtype=np.complex64)
   
    for y in range(1, Ny-1):
       Hrib[3*y:3*(y+1), 3*y:3*(y+1)] = h1
       Hrib[3*(y-1):3*y, 3*y:3*(y+1)] = (h2.conj()).T
       Hrib[3*(y+1):3*(y+2), 3*y:3*(y+1)] = h2

    Hrib[0:3, 0:3] = h1
    Hrib[3*(Ny-1):3*(Ny), 3*(Ny-1):3*(Ny)] = h1
    Hrib[3*(Ny-2):3*(Ny-1), 3*(Ny-1):3*Ny] = (h2.conj()).T
    Hrib[3:6, 0:3] = h2

    return Hrib

def solve_self_consistent(Nx, Ny, invTemp, U, initCond, beta0):
    '''
        Nx : Number of ks;
        Ny : Number of atoms along the ribbon's transverse direction
        intTemp : Inverse Temperature Beta
        U : On-site interaction
        initCond : 1: Ferromagnetic, 2: AF, 3: Paramagnetic
        beta0 = 1.1 # must be > 1, otherwise beta decreases
        
        Returns nUp, nDown, energies, itSwitch, lastIt, eUp, eDown, wUp, wDown
    '''
    
    
    N = nOrb * Nx * Ny # number of sites (orbital + spatial)

    betaTarget = invTemp # inftyCutOff in fermi means infty , i.e. T = 0
    beta = beta0 # beta starts as beta0

    itMax = 100
    it = 0
    lbda = 0.5 / (1.1 * itMax) # the factor multiplied by itMax impedes P ( I ) < \delta
    itSwitch = 0
    dampFreq = 1

    np.random.seed(1)

    if initCond == 1:
        
        nUp = np.ones(3*Ny) - 0.01 * np.random.rand(3*Ny)
        nDown = np.zeros(3*Ny) + 0.01 * np.random.rand(3*Ny)
    
    if initCond == 2:
        
        nUp = np.zeros(3*Ny)
        nDown = np.zeros(3*Ny)

        spinFlipper = 0.5

        for i in range(3*Ny):
            if i % (3) == 0:
                spinFlipper *= -1
            nUp[i] = 0.5 + spinFlipper *.1 * np.random.rand()
            nDown[i] = 0.5 - spinFlipper *.1 * np.random.rand()

    if initCond == 3:
        nUp = np.ones(3*Ny) - 0.01 * np.random.rand(3*Ny)

        nDown = np.ones(3*Ny) + 0.01 * np.random.rand(3*Ny)

    # Initialize energies
    energies = np.zeros(itMax)

    # Tolerance

    delta = 1e-6

    deltaUp = delta + 1
    deltaDown = delta + 1

    eUp = np.zeros((Nx, 3*Ny))
    wUp = np.zeros((Nx, 3*Ny, 3*Ny), dtype=np.complex64)
    eDown = np.zeros((Nx, 3*Ny))
    wDown = np.zeros((Nx, 3*Ny, 3*Ny), dtype=np.complex64)

    ks = np.linspace(-np.pi,np.pi, num=Nx, endpoint=False)

    while (it < itMax and deltaUp > delta and deltaDown > delta):
        
       # Annealing
       
        if (beta < inftyCutOff             and beta < betaTarget) : # > infty: zero temperature case
            beta = beta0 ** it
            if beta > betaTarget:
                itSwitch = it
                beta = betaTarget
        else:
            beta = betaTarget

        print('beta: ', beta)

        for kCount, k in enumerate(ks):

            C = - nUp * nDown
    
            K = Hribbon(k, Ny)
    
            Hup = K + U * np.eye(3*Ny) * ( nDown + C / 2 )
            Hdown = K + U * np.eye(3*Ny) * ( nUp + C / 2 )
        
            eUp[kCount, :], wUp[kCount, :, :] = la.eigh(Hup)
            eDown[kCount, :], wDown[kCount, :, :] = la.eigh(Hdown)
    
        nUpOld = nUp.copy()
        nDownOld = nDown.copy()

        def rootToChem(chemPot):
            return ( np.sum(fermi(eUp, chemPot, beta ))                     +  np.sum(fermi(eDown, chemPot, beta )) ) / Nx / Ny - (2 - nHole)

        mu = opt.bisect(rootToChem, -50, 50)
    
        for i in range(3*Ny):
            nUp[i] = 0
            nDown[i] = 0
            for n in range(3*Ny):
                for q in range(Nx):
                    nUp[i] += abs(wUp[q, i, n])**2 * fermi(eUp[q, n].real, mu , beta)
                    nDown[i] += abs(wDown[q, i, n])**2 * fermi(eDown[q, n].real, mu, beta)
            nUp[i] /= Nx
            nDown[i] /= Nx

        # Damping
        if it % dampFreq == 0:
            nUp = ( 1 / 2 + lbda * it ) * nUp            + ( 1 / 2 - lbda * it) * nUpOld
            nDown = ( 1 / 2 + lbda * it ) * nDown            + ( 1 / 2 - lbda * it) * nDownOld
    
        deltaUp = np.dot(nUp - nUpOld, nUp - nUpOld) / np.dot(nUpOld, nUpOld)
        deltaDown = np.dot(nDown - nDownOld, nDown - nDownOld) / np.dot(nDownOld, nDownOld)

        # To check convergence
        print('delta nUp: ', deltaUp)
        print('delta nDown: ', deltaDown)
        # Check if chemical potential is imposing
        # the right number of particles
        print('<n>: ', (nUp.sum() + nDown.sum() ) / ( 3 * Ny ) )
        
        if beta < inftyCutOff:
            energies[it] = U * np.dot(nUp, nDown) + mu * (nUp + nDown).sum()             - 1 / beta / Nx * ( ( np.log( 1 + np.exp( - beta * ( eUp - mu ) ) ) ).sum() +                           (np.log( 1 + np.exp( - beta * ( eDown - mu ) ) ) ).sum() )
        else:
            energies[it] = U * np.dot(nUp, nDown) + mu * (nUp + nDown).sum()             
        it += 1
            
    lastIt = it

    return nUp, nDown, energies, itSwitch, lastIt, eUp, eDown, mu,np.absolute(wUp.flatten('C'))**2, np.absolute(wDown.flatten('C'))**2

def savedata(SAVESUBDIR, data):
    np.savetxt(SAVESUBDIR + "nUp.txt", data[0])
    np.savetxt(SAVESUBDIR + "nDown.txt", data[1])
    np.savetxt(SAVESUBDIR + "energies.txt", data[2])
    np.savetxt(SAVESUBDIR + "itSwitch_lastIt_mu.txt", (data[3], data[4], data[7]))
    np.savetxt(SAVESUBDIR + "eUp.txt", (data[5]))
    np.savetxt(SAVESUBDIR + "eDown.txt", (data[6]))
    np.savetxt(SAVESUBDIR + "wUp.txt", (data[8]))
    np.savetxt(SAVESUBDIR + "wDown.txt", (data[9]))
    np.savetxt(SAVESUBDIR + "modelParams.txt",              (abs_t0, e1, e2, t0, t1, t2, t11, t12, t22))


# In[29]:


cwd = os.getcwd()

if localData == False:
    
    ## SAVE FILES TO HARD DRIVE (THEY ARE TOO BIG!). SET PATH BELOW.
    SAVEDIR = "../../../../../../../Volumes/ADATA HD710/MScData/MeanFieldTMDnanoribbon/"
else:
    ## SAVE FILES IN plots FOLDER
    SAVEDIR = "../plots/MeanFieldTMDnanoribbon/"


# ## Sweep U parameter for T = 0 (zero temperature)

# In[ ]:


Nx = 512
Ny = 16
beta0 = 1.15
beta = 8
initCond = 1

#Us = [1, 3, 5, 7, 10, 15, 20]

#Us = np.arange(18, 20)

#Us = np.arange(14.8, 16.2, 0.2)

Us = [10, 12, 13, 13.1, 13.2, 13.3, 13.4, 13.5,\
      14, 14.8, 15, 15.2, 15.22, 15.24, 15.26, 15.28, 15.3,15.32,15.34,15.38, 15.4, 15.42, 15.44, 15.46, 15.48, 15.5, 15.6, 15.8, 16, 17, 18, 19]

#Us = [20]

for U in Us:
    SAVESUBDIR = SAVEDIR + "/Nx=" + str(Nx) +        "_Ny=" + str(Ny) + "_U=" + str(U) + "_beta=" + str(beta) + "/"
    if not os.path.exists(SAVESUBDIR):
        os.makedirs(SAVESUBDIR)

    data = solve_self_consistent(Nx, Ny, beta, U, initCond, beta0)
    savedata(SAVESUBDIR, data)

