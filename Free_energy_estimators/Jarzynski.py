# Jarzynski_estimator : Jarzynski estimator applied to a list of W values
# cumulant_expansion: approximation of the J. exponential average
#####################################################
# Units are F=[kJ/mol] and T=[K]
# The W values are given as numpy array as input to the function 


import scipy as sc
import math
import statistics
import numpy as np

from Free_energy_estimators.Utils import Convert_kJ_kcal

#Jarzynski identity 
def Jarzyinski_function(list_work:list,
                        invert=False):
    '''
    Function to apply the J. equality to a list of Work values. ΔF=−β^−1ln⟨exp−βW⟩
    If inverte = False (default), the sign of the computed Free Energy is inverted and we obrain the unbinding Free Energy.
    If invert = True the function gives back the binding Free Energy.
    By default, inverted is = False.
    It returns the value of the binding Free Energy and print it.
    '''
    # Define the constants:
    T = 298 #in k
    kb = sc.constants.Boltzmann # 1.380649e-23 J
    Na = sc.constants.Avogadro # 6.02214076e+23 mol-1
    kbT_kJ_mol = kb*0.001*Na*T
    beta = (1/(kbT_kJ_mol))
    lista = []
    #strb = 'unbinding'
    strb = 'original sign'
    for w in list_work:
        lista.append(math.exp(-w*beta))
    
    mean=statistics.mean(lista)
    Free_Energy = -np.log(mean)/beta
        
    if invert == True:
        Free_Energy = -Free_Energy
        #strb = 'binding'
        strb ='inverted sign'
            
    #print('Jarzynski Free Energy [kJ/mol] from %s runs:'%strb, Free_Energy)
    #print('Jarzynski Free Energy [kJ/mol] with %s:'%strb, Free_Energy)

    Free_Energy_kcal = Convert_kJ_kcal(Free_Energy)
    #print('Jarzynski Free Energy [kcal/mol] from %s runs:'%strb, Free_Energy_kcal)
    #print('Jarzynski Free Energy [kcal/mol]  with %s:'%strb, Free_Energy_kcal)
    
        
    return Free_Energy_kcal

# Cumulant expansion of the exponential average truncated at the second order    
def cumulant_expansion(list_work, invert=False):
    '''
    Computing the exponential average of the Jarzynski equation as cumulant expansion truncated at the second order.
    Change sign for the binding with the inverte parameter.
    '''
    # Define the constants:
    T = 298 #in k
    kb = sc.constants.Boltzmann # 1.380649e-23 J
    Na= sc.constants.Avogadro # 6.02214076e+23 mol-1
    kbT_kJ_mol = kb*0.001*Na*T
    beta = (1/(kbT_kJ_mol))
    list_work = np.array(list_work)
    #strb = 'unbinding'
    strb = 'original sign'

    df = np.mean(list_work)-(beta/2)*np.var(list_work)
    
    if invert == True:
        df = -df
        #strb ='binding'
        strb='inverted sign'
        
    #print('Jarzynski Cumulant expansion: Free Energy [kJ/mol] from %s runs:'%strb, df)
    #print('Jarzynski Cumulant expansion: Free Energy [kJ/mol] with %s:'%strb, df)

    Free_Energy_kcal = Convert_kJ_kcal(df)
    #print('Jarzynski Cumulant expansion: Free Energy [kcal/mol] from %s runs:'%strb, Free_Energy_kcal)
    #print('Jarzynski Cumulant expansion: Free Energy [kcal/mol] with %s:'%strb, Free_Energy_kcal)
        
    return Free_Energy_kcal

