
import pandas as pd
import plumed
import numpy as np
import os

#######################################################
dt = 0.002
######################################################

#Function to read the Colvar files and create a file with replicas work values for each value of s reading the first appeareance for the S value
#N.B. work_at_different_s_values_first_app uses the work related to the first appearance of the S value
#and creates one file for each S with 2 columns: colvar_num_file and cumulative work at that time

def work_at_different_s_values_first_app(replicas,

                      input_files_path,

                      colvar_file,

                      output_files_path,

                      output_prefix,

                      initial_s,

                      final_s,

                      invert):
        


        os.mkdir(output_files_path)

        

        if invert == False:

                s_list = np.array(range(initial_s,final_s+1, 1))

        if invert == True:

                s_list = np.array(range(initial_s,final_s-1, -1))


        #read all colvar files

        dfs = [plumed.read_as_pandas(f"{input_files_path}/{colvar_file}{i}") for i in range(1,replicas+1)]


        #if invert = True (default for binding), the time_initial correspond to the lower value of S

        #and the corresponding value of work for each df, work_initial, is save in all_work_initial to be subtracted to every W values

        if invert == True:

                final_work_initial = []

                #loop on all S of the list, produce one file for each time (i.e. s)

                for i, df in enumerate(dfs):

                        if not df.loc[df['p.sss'] == final_s, 'restraint.work'].empty:

                                work_initial = df.loc[df['p.sss'] == final_s, 'restraint.work'].iloc[0]

                                df_work_initial = pd.DataFrame([f'{"covlar_"}{i+1}', work_initial]).transpose()

                                df_work_initial.columns = ['colvar_name', 'work_initial']

                                final_work_initial.append(df_work_initial)

                        else:

                                final_s_extended = final_s + 0.01

                                #print(f'{i+1}, {final_s_extended}')

                                while df.loc[(df['p.sss'] - final_s_extended) < 0.005, 'restraint.work'].empty:

                                        final_s_extended += 0.01

                                        #print(f'{i+1}, {final_s_extended}')

                                work_initial = df.loc[(df['p.sss'] - final_s_extended) < 0.005, 'restraint.work'].iloc[0]

                                df_work_initial = pd.DataFrame([f'{"covlar_"}{i+1}', work_initial]).transpose()

                                df_work_initial.columns = ['colvar_name', 'work_initial']

                                final_work_initial.append(df_work_initial)

                                print(f'initial S for colvar_{i+1} used for binding: {final_s_extended} instead of {final_s}')

                                print(f'please consider to use a value of initial_s close to {final_s_extended} in input')


                        all_works_initial = pd.concat(final_work_initial)

                all_works_initial.to_csv(f'{output_files_path}/{output_prefix}{final_s}', header=None, index=None, sep=',')

        #create a file for each value of s 

        for s in s_list:


                final_work = []

                for i, df in enumerate(dfs):

                       

                        #find the work corresponding to specific time (or for the first greater time)

                        

                        if invert==True and s==initial_s:

                                work_in_range = df.loc[0, 'restraint.work']

                                print(f'{s}, {work_in_range}')

                                df_work = pd.DataFrame([f'{"covlar_"}{i+1}', work_in_range]).transpose()

                                df_work.columns = ['colvar_name', 'work']

                                final_work.append(df_work)

                        

                        elif invert==False and s==final_s:

                                work_in_range = df.loc[df.index[-1], 'restraint.work']

                                print(f'{s}, {work_in_range}')

                                df_work = pd.DataFrame([f'{"covlar_"}{i+1}', work_in_range]).transpose()

                                df_work.columns = ['colvar_name', 'work']

                                final_work.append(df_work)

                        

                        else:

                                if not df.loc[df['p.sss'] == s, 'restraint.work'].empty:

                                        work_in_range = df.loc[df['p.sss'] == s, 'restraint.work'].iloc[0]

                                        print(f'{s}, {work_in_range}')

                                        df_work = pd.DataFrame([f'{"covlar_"}{i+1}', work_in_range]).transpose()

                                        df_work.columns = ['colvar_name', 'work']

                                        final_work.append(df_work)

                                

                                else:

                                        s_extended = s - 0.01

                                        while df.loc[np.abs(df['p.sss'] - s_extended) < 0.005, 'restraint.work'].empty:

                                                s_extended += 0.01

                                                #print(f'{i+1}, {final_s_extended}')

                                        work_in_range = df.loc[np.abs(df['p.sss'] - s_extended) < 0.005, 'restraint.work'].iloc[0]

                                        print(f'{s}, {s_extended}, {work_in_range}')

                                        df_work = pd.DataFrame([f'{"covlar_"}{i+1}', work_in_range]).transpose()

                                        df_work.columns = ['colvar_name', 'work']

                                        final_work.append(df_work)


                        #concatenate all work values of each colvar and write them in a csv

                        all_final_works = pd.concat(final_work)

                

                # if invert = True, then the difference from all_works_initial and all_final_works is calculated

                if invert == True:

                        all_final_works['work'] =  all_works_initial['work_initial'] - all_final_works['work']

                all_final_works.to_csv(f'{output_files_path}/{output_prefix}{s}', header=None, index=None, sep=',')

