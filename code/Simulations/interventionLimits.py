import numpy as np
import pandas as pd
from multiprocessing import Pool
from simulateEach import interventionLimits

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
mu = 0.04 #coral background mortality
g0 = 0.2
g1 = 1 #maximum grazing rate

#settlement 
f = 0.1 #fecundity
phi = 1
kappa = 0.75 #Value used to scale settlement to increased settlement present on free space rather than cca tissue. Approximated from R-R Williams et al. (2010)

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


workers = 20 # Setting the number of processors. 16 are available on my machine. Running with 14 processors took x mins

models = {}

# for m in ['baseModel', 'macroMinus', 'yearlyAdditions', 'yearlyCCAMacro', 'yearlyMacroCoral', 'yearlyCoral']:
# for m in ['macroMinus', 'yearlyAdditions']:
for m in ['yearlyCCACoral', 'yearlyMacroCoral', 'yearlyCoral']:
    def simulationWrapper(item):
        return interventionLimits(parameterDictionary=parameters_dict, model = m).simulations 
    
    with Pool(workers) as pool:
    # phiSims = pool.map(simulatePhi_wrapper, parameters_dict.values())
        sims = pool.map(simulationWrapper, parameters_dict.values())

    df = pd.concat(sims, ignore_index=True)

    models[m] = df
    modelName = str(m) + '.csv'
    models[m].to_csv(modelName)
