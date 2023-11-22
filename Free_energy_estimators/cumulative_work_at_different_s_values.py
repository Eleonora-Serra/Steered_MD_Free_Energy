import pandas as pd
import plumed
import numpy as np
import os



#Function to read the step to be used from the plumed.dat file
def read_step_datfile(file:str,
             itial_line_toskip=8,
             final_line_tiskip=3):
    output_list=[]
    with open(file, 'r') as filin:
        lines = filin.readlines()[itial_line_toskip:-final_line_tiskip]
        for line in lines:
        #print(line)
            line_str= line.split(' ')
            last= line_str[-1]
            last2 = last.split('=')
            last2=last2[-1].replace("\n", "")
            output_list.append(last2)
        output_list = output_list[1:]
        output_list.insert(0,'0')
        output_list=[float(x) for x in output_list]
    return output_list


#Function to read the Colvar files and create a file with replicas work values for each value of s reading the time (Cumulative W at each step)
###STEP = 1###
#N.B. work_at_different_s_values reads the time in the colvar files at which different values of s were defined in the plumed.dat
#and creates one file for each time where s was defined with 2 columns: colvar_num_file and cumulative work at that time
def work_at_different_s_values(replicas,
                                input_files_path,
                                colvar_file,
                                output_files_path,
                                output_prefix,
                                initial_s,
                                final_s,
                                plumed_dat_file,
                                invert
                               ):
    
    ''' N.B.: if invert = True, then S decreases with the time (default for binding)
    the final work associated with final time of simulation is subtracted from work values at each time
    as a result, is possible to compare binding and unbinding simulations
    (W=0 is always associated with the lower value of S)
    '''
    
    #create the directory for all the output files
    os.mkdir(output_files_path)
    
    #read steps of plumed file and transform then in time in ps with dt=0.002
    dt = 0.002
    s = initial_s
    
    plumed_steps=read_step_datfile(file= f'{plumed_dat_file}')   
    time_ps = [item * dt for item in plumed_steps]

    #write in a csv the file the plumed_steps, time_ps and the value of s
    df_time = pd.DataFrame({'steps':plumed_steps,'time':time_ps})

    if invert == False:
        df_time['S(x)'] = range(initial_s,final_s+1, 1)

    if invert == True:
        df_time['S(x)'] = range(initial_s,final_s-1, -1)

    df_time.to_csv(f'{output_files_path}/s_and_reference_time', header=True, sep=',', index = False)
    '''

    df_time = pd.read_csv(f'{output_files_path}/s_and_reference_time', header=0, sep=',', names=['Plumed_step', 'time', 'S(x)' ])
    time_ps = df_time['time']
    '''

    #read all colvar files
    dfs = [plumed.read_as_pandas(f"{input_files_path}/{colvar_file}{i}") for i in range(1,replicas+1)]

    #if invert = True (default for binding), the time_initial correspond to the lower value of S
    #and the corresponding value of work for each df, work_initial, is save in all_work_initial to be subtracted to every W values
    if invert == True:
        time_initial = df_time.iloc[-1]['time']
        
        final_work_initial = []
         #loop on all times of the list, produce one file for each time (i.e. s)
        for i,df in enumerate(dfs):
            if not df.loc[df['time'] >= time_initial, 'restraint.work'].empty:
                work_initial = df.loc[df['time'] >= time_initial, 'restraint.work'].iloc[0]
                df_work_initial = pd.DataFrame([f'{"./covlar_"}{i+1}', work_initial]).transpose()
                df_work_initial.columns = ['colvar_name', 'work_initial']
                final_work_initial.append(df_work_initial)
            else:
                print(f'ERROR: NO WORK VALUE IN THE RANGE REQUIRED {i}')

        all_works_initial = pd.concat(final_work_initial)
        all_works_initial.to_csv(f'{output_files_path}{output_prefix}{final_s}', header=None, index=None, sep=',')

    
    #create a file for each value of s corresponding to the time in time_ps
    for time in time_ps:

        final_work = []
        for i,df in enumerate(dfs):

            #find the work corresponding to specific time (or for the first greater time)
            if not df.loc[df['time'] >= time, 'restraint.work'].empty:
                work_in_range = df.loc[df['time'] >= time, 'restraint.work'].iloc[0]
                df_work = pd.DataFrame([f'{"./covlar_"}{i+1}', work_in_range]).transpose()
                df_work.columns = ['colvar_name', 'work']
                final_work.append(df_work)
            else:
                print(f'ERROR: no work value found for colvar: {i} and time: {time}')

    #concatenate all work values of each colvar and write them in a csv

        all_final_works = pd.concat(final_work)
        # if invert = True, then the difference from all_works_initial and all_final_works is calculated
        if invert == True:
            all_final_works['work'] =  all_works_initial['work_initial'] - all_final_works['work']
        all_final_works.to_csv(f'{output_files_path}{output_prefix}{s}', header=None, index=None, sep=',')

        #increment the name of file by 1 (or decrease it by 1 for invert = True)
        if invert == False:
            s += 1
        else:
            s -= 1


'''
#Function to create work files reading the values of s at different steps!
###STEP = CAN BE SELECTED!###
#This function creates one file for each step of s with 2 columns: colvar_num_file and work for the first appearance of that s value

def work_at_different_s_values_first_appearance (replicas,
                                input_files_path,
                                colvar_file,
                                output_files_path,
                                output_prefix,
                                initial_s,
                                final_s,
                                invert,
                                step=0.5,
                               ):

'''# N.B.: if invert = True, then S decreases with the time (default for binding)
    #the final work associated with final time of simulation is subtracted from work values at each time
    #as a result, is possible to compare binding and unbinding simulations
    #(W=0 is always associated with the lower value of S)
'''
    
    #create the directory for all the output files
    os.mkdir(output_files_path)
    
    # Generally, for unbinding invert = False, S increases (start =  s_min, stop = s_max+step, step = +step)
    #for binding invert = True, S decreases (start = s_max, stop = s_min-step, step = -step)
    if invert == False:
        start = initial_s
        stop = final_s+step
        step = +step
    else:
        start = initial_s
        stop = final_s-step
        step = -step



    #read all colvar files
    dfs = [plumed.read_as_pandas(f"{input_files_path}/{colvar_file}{i}") for i in range(1,replicas+1)]
    #directory with all colvar files, name of files, name of output

    #loop on all times of the list, produce one file for each step (i.e. s)
    for s in np.arange(start, stop, step):
        final_work = []
        #loop to read all colvars from different replicas
        for i, df in enumerate(dfs):
            #find s_real (read in p.sss) higher than s (theorical) for every s (theoretical) and create df with it
            df_range = df.loc[(df['p.sss'] >= s) & (df['p.sss'] < s+step)]
            
            #if the df is not empty, then create a df with the corresponding work
            if not df_range.empty:
                work_in_range = df_range['restraint.work'].iloc[0]
            #else consider the expanded condition
            else:
                df_range = df.loc[(df['p.sss'] >= s) & (df['p.sss'] < s+(step*2))]
                #if the df is not empty, then create a df with the corresponding work and print the warning that a different conditions has been considered for specific s and i
                if not df_range.empty:
                    print(f"warning! for s = {s} and i = {i+1} there is no s+{step}, s+{step*2} is considered")
                    work_in_range = df_range['restraint.work'].iloc[0]
                #else expand the conditions a second time
                else:
                    df_range = df.loc[(df['p.sss'] >= s) & (df['p.sss'] < s+(step*3))]
                    #if the df is not empty, then create a df with the corresponding work and print the warning that a different conditions has been considered for specific s and i
                    if not df_range.empty:
                        print(f"warning! for s = {s} and i = {i+1} there is no s+{step*2}, s+{step*3} is considered")
                        work_in_range = df_range['restraint.work'].iloc[0]
                    #else warning
                    else:
                        df_range = df.loc[(df['p.sss'] >= s) & (df['p.sss'] < s+(step*4))]
                        if not df_range.empty:
                            print(f"WARNING! for s = {s} and i = {i+1} there is no s+{step*3}, s+{step*4} is considered")
                            work_in_range = df_range['restraint.work'].iloc[0]
                        else:
                            df_range = df.loc[(df['p.sss'] >= s) & (df['p.sss'] < s+(step*5))]
                            if not df_range.empty:
                                print(f"ATTENTION! for s = {s} and i = {i+1} there is no s+{step*4}, s+{step*5} is considered")
                                work_in_range = df_range['restraint.work'].iloc[0]
                            else:
                                df_range = df.loc[(df['p.sss'] >= s) & (df['p.sss'] < s+(step*6))]
                                work_in_range = df_range['restraint.work'].iloc[0]
                                print(f"ATTENTION!!! for s = {s} and i = {i+1} there is no s+{step*5}, s+{step*6} is considered")
            #if work values are determined, create a df with them and append them to a list
            if work_in_range is not None:
                df_work = pd.DataFrame([f'{"./covlar_"}{i+1}', work_in_range]).transpose()
                final_work.append(df_work)
    #create a final df for each s value and save a file for each - concatenate all work values of each colvar and write them in a csv
    #! the name of the file is the s value
        all_final_works = pd.concat(final_work)
        all_final_works.to_csv(f'{output_files_path}{output_prefix}{s}', header=None, index=None, sep=',')

'''
