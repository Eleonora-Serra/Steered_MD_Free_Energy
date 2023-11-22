import pandas as pd
import scipy.integrate as spi
import numpy as np
import scipy as sc

#Compute the Free Energy Difference as ratio of the patition functions of the bound and unbound states
#The partition functions are computed as integral over the PMF

#CONSTANTS

T = 298 #in k
kb = sc.constants.Boltzmann # 1.380649e-23 J
Na= sc.constants.Avogadro # 6.02214076e+23 mol-1
kbT_kJ_mol = kb*0.001*Na*T
kbT_kcal_mol = kbT_kJ_mol/4.184
beta = (1/(kbT_kcal_mol))
#beta= (1/0.5922)

def binding_free_energy(file_FreeEnergy,
                                  discrim_frame,
                                  invert,
                                  ):
    data=pd.read_csv(file_FreeEnergy, sep=',', header=0, names=['S(x)', 'FreeEnergy'])
    x = data['S(x)'].tolist()
    FreeEnergy_list = data['FreeEnergy'].tolist()
    y = np.exp(-beta*np.array(FreeEnergy_list))
    if invert == False:
        bound_region = spi.simpson(y[:discrim_frame],x[:discrim_frame])
        unbound_region = spi.simpson(y[discrim_frame:],x[discrim_frame:])
    if invert == True:
        unbound_region = spi.simpson(y[:discrim_frame],x[:discrim_frame])
        bound_region = spi.simpson(y[discrim_frame:],x[discrim_frame:])
    Free_Energy_Binding = -1/beta*np.log(bound_region/unbound_region) 
    return Free_Energy_Binding

