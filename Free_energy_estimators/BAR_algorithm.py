# BAR_algortihm : Bennett Acceptance Ratio algorithm
#####################################################
# BAR original algorithm applied to the Crooks Theorem (non-equilibrium simulation) 
# The algorithm is solved self-sonsistently through interations
# Units are F=[kJ/mol] and T=[K]
# The W values are given as numpy array as input to the function 

import numpy as np
toll = 1e-6  # Tollerance for the algorithm 0.000001
kB = 0.00831446261815324 # [kJ/(mol*K)]
T = 298
beta = 1./(kB*T)

def bennett(v0,v1,T,toll=toll):
    '''
    The function applies iteratively the BAR algorithm solving it self-consistently
    The function returns the value of the difference of Free Energy in KbT units
    Args:
    v0: forward work vector (np array)
    v1: backward work vector (np array)
    T: T in K
    toll:tollerance for the algorithm 
    '''
    kB = 0.00831446261815324 # [kJ/(mol*K)]
    beta = 1./(kB*T)
    n0 = len(v0) 
    n1 = len(v1)        
    it = 0 #counter for the wile loop 

    df_old = 0
    den = 1
    num = 0
    eta = 0.1 #eta for the gradient descendent minimization
    
    C = df_old+np.log(n1/n0)
    ratio = num/den 

    while(abs(ratio-1)>toll):
    #denominator forward
        arg = beta*v0-C
        den = np.sum(1.0/(1.0+np.exp(arg)))
    #numerator backward
        arg = -beta*v1+C
        num = np.sum(1.0/(1.0+np.exp(arg)))
      
        ratio = num/den     
        df = eta*np.log(ratio)+df_old      
        df_old = df            
        C = df+np.log(n1/n0)
        it = it+1
    return df
