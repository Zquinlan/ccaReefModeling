import numpy as np
import pandas as pd
from multiprocessing import Pool
from SimulateFunctions import run_model

#Initial parameters which will be used in the model
dt = 1 # time step (year)
NUMSTEPS = 1000
NUMYEARS = int(NUMSTEPS/dt)

## No Alkalinity parameters
#Growth rates
r = 0.2 #coral growth rate
rA = 0.4 #CCA growth rate
gamma = 0.7 #algal growth rate

#overgrowth rates
a = 0.3 #algal overgrowth onto corals
algalCCA = 0.3 #algal overgrowth onto CCA
coralCCA = 0.1 #coral overgrowth onto CCA

#death, grazing
mu = 0.02 #coral background mortality
g0 = 0.2
g1 = 1 #maximum grazing rate

#settlement 
f = 0.1 #fecundity
phi = 1
kappa = -0.4 #Value used to scale settlement to increased settlement present on free space rather than cca tissue. Approximated from R-R Williams et al. (2010)

parameters_dict = {'dt':dt, 
                   'NUMYEARS': NUMYEARS, 
                    'r': r, 
                    'a': a,
                    'mu': mu,
                    'gamma': gamma,
                    'g0': g0,
                    'g1': g1,
                    'rA': rA,
                    'algalCCA': algalCCA,
                    'coralCCA': coralCCA,
                    'f': f,
                    'kappa': kappa,
                    'phi': phi
                    }

#Defining kappa and mortality which will be adjusted throughout the model
kappaSteps = 0.05
kappa0 = -1
kappaMax = 1

muSteps= 0.005
mu0 = 0.01
muMax = 0.2

kappaPerms = int((kappaMax-kappa0)/kappaSteps)
muPerms = int((muMax-mu0)/muSteps)

class simulatePhi:
    def __init__(self, parameterDictionary, xName, xMin, xMax, xPerms, yName, yMin, yMax, yPerms, totalPermutations):
        super().__init__()
        """
        parameterDictionary is all of the parameters which the model takes as input
        xName is the name of the first parameter you want to change as a string
        xMin is the minimmum value of the first parameter you want to change
        xMax is the maximium value of the first parameter you want to change
        xPerms is the number of steps to take to get from xMin to xMax
        yName is the name of the first parameter you want to change as a string
        yMin is the minimmum value of the first parameter you want to change
        yMax is the maximium value of the first parameter you want to change
        yPerms is the number of steps to take to get from yMin to yMax
        
        totalPermutations is the number of variations of C0, A0, and M0 which the model will test for each x and y value
        """
        simulations = pd.DataFrame({'xParam' : [],
                                         'yParam' : [],
                                         'c0' : [],
                                         'm0' : [],
                                         'a0' : [],
                                         'cEnd' : [],
                                         'mEnd' : [],
                                         'aEnd' : []})
        for g in np.linspace(-1, 1, 20):
            for i in np.linspace(xMin, xMax, xPerms):
                for t in np.linspace(yMin, yMax, yPerms):
                    parameterDictionary['phi'] = g
                    parameterDictionary[xName] = i
                    parameterDictionary[yName] = t
                    for c in [0.01, 0.05,0.5,0.75]:
                        for m in [0.01, 0.05,0.5,0.75]:
                            for a in [0.01, 0.05,0.5,0.75]:
                                if (c + m + a <= 1) : 
                                    C_array, M_array, A_array = run_model(c,m,a, parameterDictionary)
                                    dfTemp = pd.DataFrame({'xParam' : [i],
                                                        'yParam' : [t],
                                                        'phi' : parameterDictionary['phi'],
                                                        'g0' : [g],
                                                        'c0' : [c],
                                                        'm0' : [m],
                                                        'a0' : [a],
                                                        'cEnd' : [C_array[-1]],
                                                        'mEnd' : [M_array[-1]],
                                                        'aEnd' : [A_array[-1]]})

                                    simulations = pd.concat([simulations, dfTemp], ignore_index = True)
        self.simulations = simulations


def simulatePhi_wrapper(item):
    return simulatePhi(parameterDictionary=parameters_dict, 
                        xName='mu', xMin=mu0, xMax=muMax, xPerms=20,
                        yName='kappa', yMin=kappa0, yMax=kappaMax, yPerms=20,
                        totalPermutations=1000).simulations


workers = 300 # Setting the number of processors. 16 are available on my machine. Running with 14 processors took x mins

with Pool(workers) as pool:
    # phiSims = pool.map(simulatePhi_wrapper, parameters_dict.values())
    grazingSims = pool.map(simulatePhi_wrapper, parameters_dict.values())

# phiDF = pd.concat(phiSims, ignore_index=True)
grazingDF = pd.concat(grazingSims, ignore_index=True)


# phiDF.to_csv('phiData.csv')
grazingDF.to_csv('grazingData.csv')