
import pandas as pd
import numpy as np

from Free_energy_estimators.Utils import extract_list_work, Convert_kJ_kcal, Convert_list
from Free_energy_estimators.BAR_algorithm import bennett

# INITIAL STEP AND FINAL STEP ARE DEFINED IN CRESCENT ORDER
# function to calculate the free energy difference using BAR associated to each value of S
def BAR_for_PMF(initial_s,
                final_s,
                file_bind,
                file_unbind,
                file_prefix):
    
    T = 298
    toll = toll = 1e-6  # 0.000001
    kB = 0.00831446261815324 # [kJ/(mol*K)]
    beta = 1./(kB*T)

    BAR_values_list = []

    for i in range(initial_s,final_s+1):
     
        bind=f"{file_bind}{file_prefix}{i}"
        unbind=f"{file_unbind}{file_prefix}{i}"
                
        list_work_unbind = extract_list_work(unbind,'list_work_unbind')
        list_work_bind = extract_list_work(bind,'list_work_bind')
        list_work_bind = Convert_list(list_work_bind)
        
        # converting list to array
        v0 = np.array(list_work_unbind)  
        v1 = np.array(list_work_bind)
        
        
        df = bennett(v0,v1,T,toll)
        df = (1/beta)*df
        df = Convert_kJ_kcal(df)
      
        BAR_values_list.append(df)
    return BAR_values_list

