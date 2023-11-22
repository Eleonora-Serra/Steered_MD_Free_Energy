import scipy as sc
import numpy as np

T = 298 # K
kb = 1.380649*10**(-23) # J K-1
Na= 6.02214076*10**(+23) # mol-1
kbT_kJ_mol = kb*0.001*Na*T # J*K-1*mol-1*K*10^-3 = kJ mol-1
cal_to_J = 4.184 # kJ / kcal
kbT_kcal_mol = kbT_kJ_mol/cal_to_J # kJ mol-1 / (kJ / kcal) = kcal mol-1
beta = (1/(kbT_kcal_mol)) # mol kcal-1

def volume_correction(V_unbound):
    V_std = 1661 #A^3
    V_unb_V_std = V_unbound/V_std
    return -1/beta*np.log(V_unb_V_std)  # 1/mol kcal-1 = kcal mol-1


def absolute_binding_free_energy(V_unbound, binding_free_energy):
    return volume_correction(V_unbound) + binding_free_energy 

### returns absolute binding free energy in kcal mol-1###s