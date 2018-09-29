
# coding: utf-8

# In[3]:


import numpy as np
import numpy.linalg as la
import os
import copy
import warnings
import scipy.optimize as opt
import sys


##      DEFINE THE MODEL    ##


nOrb = 3

nHole = 0

tmd = 'MoS_2'


# In[78]:


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

inftyCutOff = 100 # above it, beta is practically infinity

def fermi(e, mu, beta):
    '''
        For zero temperature, set beta = inftyCutOff
        '''
    if beta == inftyCutOff:
        return (e < mu).astype(int)
    else:
        return 1 / ( 1 + np.exp( beta * ( e - mu ) ) )

def iTriang(x, y, Nx, Ny):
    return Nx * y + x

def triangularNano(Nx, Ny, nOrb, hoppings):
    T = np.zeros((nOrb*Nx*Ny, nOrb*Nx*Ny))
    for x in range(Nx):
        for y in range(Ny):
            # Diagonal term
            T[ iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny) + 1) * nOrb,              iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny) + 1) * nOrb ]            = hoppings[0]
            
            # E1
            T[ iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny) + 1)*nOrb,              iTriang( (x + 1) % Nx , y, Nx, Ny)*nOrb:(iTriang( (x + 1) % Nx , y, Nx, Ny) + 1)*nOrb ]              = hoppings[1]
            
            # E4
            T[ iTriang( (x + 1) % Nx , y, Nx, Ny)*nOrb:(iTriang( (x + 1) % Nx , y, Nx, Ny) + 1)*nOrb              , iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb ] = hoppings[4]
            
            if y == 0:
                T[ iTriang(x, 0, Nx, Ny)*nOrb:(iTriang(x, 0, Nx, Ny)+1)*nOrb,                  iTriang( x, 1, Nx, Ny)*nOrb:(iTriang( x, 1, Nx, Ny)+1)*nOrb ]                  = hoppings[6]
                
                T[ iTriang(x, 1, Nx, Ny)*nOrb:(iTriang(x, 1, Nx, Ny)+1)*nOrb,                  iTriang( x, 0, Nx, Ny)*nOrb:(iTriang( x, 0, Nx, Ny)+1)*nOrb ]                  = hoppings[3]
                
                if x == 0:
                    T[iTriang(x, 0, Nx, Ny)*nOrb:(iTriang(x, 0, Nx, Ny) + 1)*nOrb,                      iTriang( Nx - 1, 1, Nx, Ny)*nOrb:(iTriang( Nx - 1, 1, Nx, Ny) + 1)*nOrb]                      = hoppings[5]
                    T[iTriang( Nx - 1, 1, Nx, Ny)*nOrb:(iTriang( Nx - 1, 1, Nx, Ny)+1)*nOrb,                      iTriang(x, 0, Nx, Ny)*nOrb:(iTriang(x, 0, Nx, Ny)+1)*nOrb] = hoppings[2]
                else:
                    T[iTriang(x, 0, Nx, Ny)*nOrb:(iTriang(x, 0, Nx, Ny)+1)*nOrb,                      iTriang( x - 1, 1, Nx, Ny)*nOrb:(iTriang( x - 1, 1, Nx, Ny)+1)*nOrb] = hoppings[5]
                    T[iTriang(x - 1, 1, Nx, Ny)*nOrb:(iTriang(x - 1, 1, Nx, Ny)+1)*nOrb,                      iTriang(x, 0, Nx, Ny)*nOrb:(iTriang(x, 0, Nx, Ny)+1)*nOrb] = hoppings[2]
        else:
            if y == Ny - 1:
                T[iTriang(x, Ny - 1 , Nx, Ny)*nOrb:(iTriang(x, Ny - 1 , Nx, Ny) + 1)*nOrb,                  iTriang( (x + 1) % Nx , Ny - 2, Nx, Ny)*nOrb:(iTriang( (x + 1) % Nx , Ny - 2, Nx, Ny) + 1)*nOrb]= hoppings[2]
                T[iTriang( (x + 1) % Nx , Ny - 2, Nx, Ny)*nOrb:(iTriang( (x + 1) % Nx , Ny - 2, Nx, Ny)+1)*nOrb,                  iTriang(x, Ny - 1, Nx, Ny)*nOrb:(iTriang(x, Ny - 1, Nx, Ny)+1)*nOrb] = hoppings[5]
                T[iTriang(x, Ny - 1, Nx, Ny)*nOrb:(iTriang(x, Ny - 1, Nx, Ny)+1)*nOrb,                  iTriang( x, Ny - 2, Nx, Ny)*nOrb:(iTriang( x, Ny - 2, Nx, Ny)+1)*nOrb] = hoppings[3]
                T[iTriang(x, Ny - 2, Nx, Ny)*nOrb:(iTriang(x, Ny - 2, Nx, Ny)+1)*nOrb,                  iTriang( x, Ny - 1, Nx, Ny)*nOrb:(iTriang( x, Ny - 1, Nx, Ny)+1)*nOrb] = hoppings[6]
                
            else:
                T[iTriang(x, y , Nx, Ny)*nOrb:(iTriang(x, y , Nx, Ny)+1)*nOrb,                  iTriang( (x + 1) % Nx , y - 1, Nx, Ny)*nOrb:(iTriang( (x + 1) % Nx , y - 1, Nx, Ny)+1)*nOrb] = hoppings[2]
                T[iTriang( (x + 1) % Nx , y - 1, Nx, Ny)*nOrb:(iTriang( (x + 1) % Nx , y - 1, Nx, Ny)+1)*nOrb,                  iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb] = hoppings[5]
                T[iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb,                  iTriang( x, y - 1, Nx, Ny)*nOrb:(iTriang( x, y - 1, Nx, Ny)+1)*nOrb] = hoppings[3]
                T[iTriang(x, y - 1, Nx, Ny)*nOrb:(iTriang(x, y - 1, Nx, Ny)+1)*nOrb,                  iTriang( x, y, Nx, Ny)*nOrb:(iTriang( x, y, Nx, Ny)+1)*nOrb] = hoppings[6]
                if x == 0:
                    T[iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb,                      iTriang( Nx - 1, y + 1, Nx, Ny)*nOrb:(iTriang( Nx - 1, y + 1, Nx, Ny)+1)*nOrb] = hoppings[5]
                    T[iTriang( Nx - 1, y + 1, Nx, Ny)*nOrb:(iTriang( Nx - 1, y + 1, Nx, Ny)+1)*nOrb,                      iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb] = hoppings[2]
                else:
                    T[iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb,                      iTriang(x - 1, y + 1, Nx, Ny)*nOrb:(iTriang(x - 1, y + 1, Nx, Ny)+1)*nOrb] = hoppings[5]
                    T[iTriang(x - 1, y + 1, Nx, Ny)*nOrb:(iTriang(x - 1, y + 1, Nx, Ny)+1)*nOrb,                      iTriang(x, y, Nx, Ny)*nOrb:(iTriang(x, y, Nx, Ny)+1)*nOrb] = hoppings[2]
    return T

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


def solve_self_consistent(Nx, Ny, invTemp, U, initCond):
    '''
        Nx : Number of ks;
        Ny : Number of atoms along the ribbon's transverse direction
        intTemp : Inverse Temperature Beta
        U : On-site interaction
        initCond : 1: Ferromagnetic, 2: AF, 3: Paramagnetic
        
        Returns nUp, nDown, energies, itSwitch, lastIt, eUp, eDown, wUp, wDown
    '''
    
    K = triangularNano(Nx, Ny, nOrb, hoppings)
    
    N = nOrb * Nx * Ny # number of sites (orbital + spatial)

    beta0 = 1.5 # must be > 1, otherwise beta decreases
    betaTarget = invTemp # inftyCutOff in fermi means infty , i.e. T = 0
    beta = beta0 # beta starts as beta0

    itMax = 100
    it = 0
    lbda = 0.5 / (1.1 * itMax) # the factor multiplied by itMax impedes P ( I ) < \delta
    itSwitch = 0
    dampFreq = 1

    np.random.seed(1)

    if initCond == 1:
        
        nUp = np.ones(3*Ny) - 0.1 * np.random.rand(3*Ny)
        nDown = np.zeros(3*Ny) + 0.1 * np.random.rand(3*Ny)
    
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
        nUp = np.ones(3*Ny) - 0.1 * np.random.rand(3*Ny)

        nDown = np.ones(3*Ny) + 0.1 * np.random.rand(3*Ny)

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
       
        if (beta < inftyCutOff            and beta < betaTarget) : # > infty: zero temperature case
           beta = beta0 ** it
           if beta > betaTarget:
               itSwitch = it
               print(itSwitch)
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
                
        energies[it] = U * np.dot(nUp, nDown) / 3 / Ny + mu * (nUp + nDown).sum()         - 1 / beta * ( ( np.log( 1 + np.exp( - beta * ( eUp - mu ) ) ) ).sum() +                       (np.log( 1 + np.exp( - beta * ( eDown - mu ) ) ) ).sum() )
            
        it += 1
            
    lastIt = it
    print(lastIt)

    return nUp, nDown, energies, itSwitch, lastIt, eUp, eDown, wUp, wDown, mu

def savedata(SAVESUBDIR, data):
    np.savetxt(SAVESUBDIR + "nUp.txt", data[0])
    np.savetxt(SAVESUBDIR + "nDown.txt", data[1])
    np.savetxt(SAVESUBDIR + "energies.txt", data[2])
    np.savetxt(SAVESUBDIR + "itSwitch_lastIt_mu.txt", (data[3], data[4], data[9]))
    np.savetxt(SAVESUBDIR + "eUp.txt", (data[5]))
    np.savetxt(SAVESUBDIR + "eDown.txt", (data[6]))
    np.savetxt(SAVESUBDIR + "wUp.txt", (data[7]).flatten('C'))
    np.savetxt(SAVESUBDIR + "wDown.txt", (data[8]).flatten('C'))
    np.savetxt(SAVESUBDIR + "modelParams.txt",              (abs_t0, e1, e2, t0, t1, t2, t11, t12, t22))


# In[79]:


cwd = os.getcwd()
SAVEDIR = "../plots/MeanFieldTMDnanoribbon/"
if not os.path.exists(SAVEDIR):
    os.makedirs(SAVEDIR)


# ## Small ribbon

# In[93]:


Nx = int(sys.argv[1])
Ny = int(sys.argv[2])
beta = int(sys.argv[3])
U = int(sys.argv[4])
initCond = int(sys.argv[5])


SAVESUBDIR = SAVEDIR + "/Nx=" + str(Nx) +\
    "_Ny=" + str(Ny) + "_U=" + str(U) + "_beta=" + str(beta) + "/"
if not os.path.exists(SAVESUBDIR):
    os.makedirs(SAVESUBDIR)
    
data = solve_self_consistent(Nx, Ny, beta, U, initCond)
savedata(SAVESUBDIR, data)


