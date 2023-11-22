import pandas as pd
import plumed
import numpy as np
import math
import os

from Free_energy_estimators.Utils import extract_list_work
from Free_energy_estimators.Jarzynski import Jarzyinski_function


# function to calculate the free energy difference using Jarzynski associated to each value of S
def Jarzynsky_for_PMF(initial_s,
                      final_s,
                      output_files_path,
                      output_prefix,
                      invert
                        ):
    
    Jarzynski_values_list = []
    
    #define the range to read work_s_ files
    if (invert == False):
        start = initial_s
        stop = final_s+1
        step = +1
    else:
        start = initial_s
        stop = final_s-1
        step = -1

    #for each S, read work_s_ file and calculate free energy of binding using Jarzynski
    for i in np.arange(start,stop,step):
        file=f"{output_files_path}/{output_prefix}{i}"
        list_work= extract_list_work(file,'list_work')    
        if invert == False:    
            Free_Energy_binding = Jarzyinski_function(list_work=list_work, invert = False)
            Jarzynski_values_list.append(Free_Energy_binding)
        if invert == True:
            Free_Energy_binding = Jarzyinski_function(list_work=list_work, invert = True)
            Jarzynski_values_list.append(Free_Energy_binding)
    return Jarzynski_values_list
