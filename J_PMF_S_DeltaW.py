import pandas as pd
import plumed
import numpy as np
import os

#directory with all colvar files, name of files, name of output
from Free_energy_estimators.Utils import extract_list_work
from Free_energy_estimators.Jarzynski import Jarzyinski_function


#step used is 0.5 for s and 1 for time

#loop over s


def deltaW_definition(path_dir,
                      input_prefix,
                      output_file_path,
                      output_name,
                      initial_step,
                      final_step,
                      invert,
                      ):
                      
    #create directory for the output files
    os.mkdir(output_file_path)

    if invert == False:
        for s in range(initial_step,final_step,1):
    #read 2 work columns (each with work associated to colvar from 1 to 50) for s and s+0.5
            work_file_s = f"{path_dir}/{input_prefix}{s}"
            work_file_s_plus_half = f"{path_dir}/{input_prefix}{s+1}"
            df_s = pd.read_csv(work_file_s, sep=',', header=None, names=['colvar_num_file', 'work_s'])
            df_s_plus_half = pd.read_csv(work_file_s_plus_half, sep=',', header=None, names=['colvar_num_file', 'work_s_plus_half'])
    #work_diff = None
        
    #calculate the difference betweeen 2 columns and create a new dataframe with it
            work_diff = round(df_s_plus_half['work_s_plus_half'] - df_s['work_s'], 4)
            delta_work = pd.DataFrame([df_s_plus_half['colvar_num_file'], work_diff]).transpose()
    #print(delta_work)
            delta_work.to_csv(f'{output_file_path}/{output_name}_{s+1}-{s}', header=None, index=None, sep=',')

    if invert == True:
        for s in range(initial_step,final_step,-1):
            work_file_s = f"{path_dir}/{input_prefix}{s}"
            work_file_s_minus_one = f"{path_dir}/{input_prefix}{s-1}"
            df_s = pd.read_csv(work_file_s, sep=',', header=None, names=['colvar_num_file', 'work_s'])
            df_s_minus_one = pd.read_csv(work_file_s_minus_one, sep=',', header=None, names=['colvar_num_file', 'work_s_minus_one'])
    #work_diff = None
        
    #calculate the difference betweeen 2 columns and create a new dataframe with it
            work_diff = round(df_s['work_s'] - df_s_minus_one['work_s_minus_one'],4)
            delta_work = pd.DataFrame([df_s_minus_one['colvar_num_file'], work_diff]).transpose()
    #print(delta_work)
            delta_work.to_csv(f'{output_file_path}/{output_name}_{s}-{s-1}', header=None, index=None, sep=',')



def Jarzynsky_for_PMF_deltaW(final_step:int,
                      wdir,
                      file_prefix,
                    initial_step,
                    invert,
                        ):
    
    Jarzynski_values_list = []
    if invert == False:
        for i in range(initial_step,final_step):
            print(i)
            file=f"{wdir}/{file_prefix}_{i+1}-{i}"
            print(file)
            list_work= extract_list_work(file,'list_work')
            Free_Energy_binding = Jarzyinski_function(list_work=list_work, invert= False)
            Jarzynski_values_list.append(Free_Energy_binding)
    if invert == True:
        for i in range(final_step,initial_step):
            file=f"{wdir}/{file_prefix}_{i+1}-{i}"
            print(file)

            list_work= extract_list_work(file,'list_work')
            Free_Energy_binding = Jarzyinski_function(list_work=list_work, invert= True)
            Jarzynski_values_list.append(Free_Energy_binding)

    print(Jarzynski_values_list)
    return Jarzynski_values_list

