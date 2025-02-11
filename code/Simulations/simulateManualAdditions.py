import numpy as np
import pandas as pd
from multiprocessing import Pool
from simulateEach import dualSimulate

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
muMax = 0.08

kappaPerms = int((kappaMax-kappa0)/kappaSteps)
muPerms = int((muMax-mu0)/muSteps)

workers = 20 # Setting the number of processors. 16 are available on my machine. Running with 14 processors took x mins

models = {}

for m in ['macroMinus', 'yearlyAdditions', 'yearlyCCAMacro' 'yearlyMacroCoral', 'yearlyCoral']:
    def simulationWrapper(item):
        return dualSimulate(parameterDictionary=parameters_dict, 
                            xName='mu', xMin=mu0, xMax=muMax, xPerms=muPerms,
                            yName='kappa', yMin=kappa0, yMax=kappaMax, yPerms=kappaPerms, model = m).simulations 
    
    with Pool(workers) as pool:
    # phiSims = pool.map(simulatePhi_wrapper, parameters_dict.values())
        sims = pool.map(simulationWrapper, parameters_dict.values())

    df = pd.concat(sims, ignore_index=True)

    models[m] = df
    modelName = str(m) + '.csv'
    models[m].to_csv(modelName)
