import numpy as np
import pandas as pd


#Define the initial equations, runge-kutta methods and run them together with run_model()
#Defining the differential equations for Coral, macroalgae and CCA
def dNdt(C,M,A,P):
    
    dt = P['dt']
    NUMYEARS = P['NUMYEARS']
    r = P['r']
    a = P['a']
    mu = P['mu']
    gamma = P['gamma']
    g0 = P['g0']
    g1 = P['g1']
    rA = P['rA']
    algalCCA = P['algalCCA']
    coralCCA = P['coralCCA']
    f=P['f']
    kappa=P['kappa']
    phi=P['phi']
    
    # Calculating base larvae looking to settle
    Fr = (1-C-M-A)
    beta = f*C*A #Base amount of swimming larvae cued to settle
    pA = (1+kappa)/((1 + kappa)+(1-kappa))
    pF = (1-kappa)/((1 + kappa)+(1-kappa))
    gM = (g0 + (g1 - g0)*(C**0.5))*M
    # gM = g0*M
    
    #! Calculate the derivative    
    dC = r*C*Fr - a*C*M - mu*C + coralCCA*C*A + beta*phi*(pA*A + pF*Fr) #Coral equation with CCA recruitment and overgrowth incorporated
    dM = gamma*M*Fr + a*C*M - gM + algalCCA*A*M #Macroalgae equation with CCA overgrowth incorporated
    dA = rA*A*Fr - coralCCA*C*A - algalCCA*A*M - beta*phi*pA*A # AA equation 

    # dC = r*C*Fr - a*C*M - mu*C #Coral equation with CCA recruitment and overgrowth incorporated
    # dM = gamma*M*Fr + a*C*M - gM #Macroalgae equation with CCA overgrowth incorporated

    return dC, dM, dA
    # return dC, dM
#Defining 4th-order Runge-Kutta equations
def RK4(C, M, A, P): #4th-order Runge-Kutta
    # Initial values
    C_init = C
    M_init = M
    A_init = A

    # Step 1
    k1_C, k1_M, k1_A = dNdt(C, M, A, P)

    # Step 2
    C1 = C + 0.5 * k1_C
    M1 = M + 0.5 * k1_M
    A1 = A + 0.5 * k1_A
    k2_C, k2_M, k2_A = dNdt(C1, M1, A1, P)

    # Step 3
    C2 = C + 0.5 * k2_C
    M2 = M + 0.5 * k2_M
    A2 = A + 0.5 * k2_A
    k3_C, k3_M, k3_A = dNdt(C2, M2, A2, P)

    # Step 4
    C3 = C + k3_C
    M3 = M + k3_M
    A3 = A + k3_A
    k4_C, k4_M, k4_A = dNdt(C3, M3, A3, P)

    # Calculate weighted average
    dC = (k1_C + 2 * k2_C + 2 * k3_C + k4_C) / 6
    dM = (k1_M + 2 * k2_M + 2 * k3_M + k4_M) / 6
    dA = (k1_A + 2 * k2_A + 2 * k3_A + k4_A) / 6

    # Update values
    C = C_init + dC
    M = M_init + dM
    A = A_init + dA

    return C, M, A


def run_model(INIT_C, INIT_M, INIT_A, addition, P):
   
    NUMYEARS = P['NUMYEARS']
    
    C = np.zeros((NUMYEARS+1))
    M = np.zeros((NUMYEARS+1))
    A = np.zeros((NUMYEARS+1))
    
    C[0] = INIT_C
    M[0] = INIT_M
    A[0] = INIT_A
    
    
    for year in np.arange(0,NUMYEARS):
            C[year+1],M[year+1],A[year+1] = RK4(C[year],M[year],A[year],P)

            if (M[year +1] > C[year+1]) & (M[year+1] >=0.5):
                A[year+1] = A[year+1] + (1- C[year+1]- M[year+1]- A[year+1])*addition
    
    return C, M, A

#dual Simulate takes two model
class dualSimulate:
    def __init__(self, parameterDictionary, xName, xMin, xMax, xPerms, yName, yMin, yMax, yPerms):
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
                                         'additions': [],
                                         'cEnd' : [],
                                         'mEnd' : [],
                                         'aEnd' : []})
        
        for i in np.linspace(xMin, xMax, xPerms):
            for t in np.linspace(yMin, yMax, yPerms):
                parameterDictionary[xName] = i
                parameterDictionary[yName] = t
                for x in np.linspace(0.01, 1, 10):
                    for c in [0.05,0.45,0.75]:
                        for m in [0.05,0.45,0.75]:
                            for a in [0.05,0.45,0.75]:
                                if (c + m + a <= 1) : 
                                    C_array, M_array, A_array = run_model(c,m,a,x, parameterDictionary)
                                    dfTemp = pd.DataFrame({'xParam' : [i],
                                                        'yParam' : [t],
                                                        'c0' : [c],
                                                        'm0' : [m],
                                                        'a0' : [a],
                                                        'additions' : [x],
                                                        'cEnd' : [C_array[-1]],
                                                        'mEnd' : [M_array[-1]],
                                                        'aEnd' : [A_array[-1]]})

                                    simulations = pd.concat([simulations, dfTemp], ignore_index = True)
                                                


        self.simulations = simulations