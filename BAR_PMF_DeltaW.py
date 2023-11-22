import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from BAR_PMF_S_cumulativeW import BAR_for_PMF
from Free_energy_estimators.Utils import Convert_list

#directory with all colvar files, name of files, name of output


import pandas as pd
import numpy as np

from Free_energy_estimators.Utils import extract_list_work, Convert_kJ_kcal, Convert_list
from Free_energy_estimators.BAR_algorithm import bennett

#INITIAL STEP AND FINAL STEP ARE DEFINED IN CRESCENT ORDER

def BAR_for_PMF_DeltaW(final_step:int,
                      file_bind,
                    file_unbind,
                      file_prefix,
                    initial_step,
                        ):
    
    T = 298
    toll = toll = 1e-6  # 0.000001
    kB = 0.00831446261815324 # [kJ/(mol*K)]
    beta = 1./(kB*T)

    BAR_values_list = []

    for i in range(initial_step,final_step):
     
        bind=f"{file_bind}/{file_prefix}_{i+1}-{i}"
   
        unbind=f"{file_unbind}/{file_prefix}_{i+1}-{i}"
                
        list_work_unbind = extract_list_work(unbind,'list_work_unbind')
        list_work_bind = extract_list_work(bind,'list_work_bind')
        list_work_bind = Convert_list(list_work_bind)
        
        # converting list to array
        v0 = np.array(list_work_unbind)  
        v1 = np.array(list_work_bind)
        
        #print('UNBINDING',v0)
        #print('BINDING',v1)
        
        df = bennett(v0,v1,T,toll)
        #print(df)
        df = (1/beta)*df

        df = Convert_kJ_kcal(df)
        print('Free Energy [kcal/mol]:', df)
        
      
        BAR_values_list.append(df)
    return BAR_values_list





