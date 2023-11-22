#import plumed
import matplotlib.pyplot as plt
from scipy.stats import norm
import pandas as pd
from pathlib import Path
import numpy as np
import six
import seaborn as sns


#Function to read the Colvar file without plumed API
def read_colvar(fname, verbose=True):
    '''
    The funtion takes as argument the path of the colvar file to be red and it returns a dataframe
    containing the colvar data. 
    '''
    if verbose: print(fname)
    f = open(fname) if isinstance(fname, six.string_types) else fname

    with f:
        line = f.readline()
        assert line[:2] == '#!', 'Missing or incorrect header'
        columns = line.split()[2:]
        f.seek(0)
        df = pd.read_csv(f, delim_whitespace=True, comment="#", skiprows=1, header=None)

    df.columns = columns
    df.time = map(np.round, df.time)
    return df

# Function to plot the W profile (work accumulated during the simulation) as function of time 
def plt_work_vs_time(nrun:int,
                    wdir:Path,
                    file_prefix:str,
                    column_name = 'restraint.work'):
    '''
    Funtion to plot the W profile as function of t. 
    The name of the Colvar_file from different replica is composed by: file_prefix+nrun
    Arguments:
    nrun: total number of replicas considered. This argument is important to plot the W profile from multiple 
    replica together. The function plot from 1 to nrun. 
    wdir: working directory where all the Colvar_file from multiple replica are placed. 
    file_prefix: prefix of the Colvar_file name 
    column_name: name of the column of the work, by default is restraint.work
    '''
    for i in range(1,nrun+1):
        data=read_colvar(f"{wdir}/{file_prefix}{i}")
        data.rename(columns = {f'{column_name}':'work'}, inplace = True)
        plt.plot(data['time'].values ,data['work'].values)
        plt.title('work values as function of time')
        plt.xlabel("time[ps]")
        plt.ylabel("work [Kj/mol]")
    plt.show()

#Function to plot the work accumulated during the simulation as function of the Cv s(x)
def plt_work_vs_s(nrun:int,
                 wdir:Path,
                 file_prefix:str,
                column_name='restraint.work'):
    '''
    Funtion to plot the W profile as function of the CV s. 
    The name of the Colvar_file from different replica is composed by: file_prefix+nrun
    Arguments:
    nrun: total number of replicas considered. This argument is important to plot the W profile from multiple 
    replica together. The function plot from 1 to nrun. 
    wdir: working directory where all the Colvar_file from multiple replica are placed. 
    file_prefix: prefix of the Colvar_file name 
    column_name: name of the column of the work, by default is restraint.work
    '''
    for i in range(1,nrun+1):
        data=read_colvar(f"{wdir}/{file_prefix}{i}")
        data.rename(columns = {f'{column_name}':'work'}, inplace = True)
        data.rename(columns = {f'p.sss':'sss'}, inplace = True)
        plt.plot(data['sss'].values ,data['work'].values)
        plt.title('work values as function of the collective variable S')
        plt.xlabel("Collectiva variable:s")
        plt.ylabel("work [Kj/mol]")
    plt.show()

# Function to Plot the CV profile as function of time 
def plt_s_profile(nrun:int,
                 wdir:Path,
                file_prefix:str):
    '''
    Funtion to plot the CV profile as function of time. 
    The name of the Colvar_file from different replica is composed by: file_prefix+nrun
    Arguments:
    nrun: total number of replicas considered. This argument is important to plot the W profile from multiple 
    replica together. The function plot from 1 to nrun. 
    wdir: working directory where all the Colvar_file from multiple replica are placed. 
    file_prefix: prefix of the Colvar_file name 
    '''
    for i in range(1,nrun+1):
        data=read_colvar(f"{wdir}/{file_prefix}{i}")
        data.rename(columns = {f'p.sss':'sss'}, inplace = True)
        plt.plot(data['time'].values ,data['sss'].values)
        plt.title('Collective variable S profile')
        plt.xlabel("time [ps]")
        plt.ylabel("Collectiva variable:s")
    plt.show()
    
    #Ispiration Article:
#Free energy calculation of single molecular interaction using Jarzynski’s identity method: the case of HIV-1 protease inhibitor system


#Function to be used in plot_work_and_distribution function
def get_params(arr):
    mean = np.mean(arr)
    variance = np.var(arr)
    sigma = np.sqrt(variance)
    return mean, variance, sigma

#FUnction to plot the W profile of multiple replica and the W distribution in the same plot
def plot_work_and_distribution(nrun:int,
                               wdir:Path,
                               file_prefix:str,
                               list_final_work:list,
                               column_name='restraint.work'
                               ):
    
    '''
    Function to plot the W profile of multiple replica as function of time + W distribution in the same plot.
    The name of the Colvar_file from different replica is composed by: file_prefix+nrun
    Arguments:
    nrun: total number of replicas considered. This argument is important to plot the W profile from multiple 
    replica together. The function plot from 1 to nrun. 
    wdir: working directory where all the Colvar_file from multiple replica are placed. 
    file_prefix: prefix of the Colvar_file name 
    list_final_work: list of the final work to be used for the W distribution plot.
    column_name: name of the column of the work, by default is restraint.work
    '''

    fig, (a0,a1) = plt.subplots(1,2, sharey=True, gridspec_kw = {'width_ratios': [3,1]})
    #fig = plt.figure(figsize=(6, 6))
    #fig.add_subplot(1,2,1)
    
    for i in range(1,nrun+1):
        data=read_colvar(f"{wdir}/{file_prefix}{i}")
        data.rename(columns = {f'{column_name}':'work'}, inplace = True)
        a0.plot(data['time'].values ,data['work'].values)
        #a0.plot(data['time'].values ,data['work'].values, label=f"work_run{i}")

        a0.set_xlabel("time[ps]")
        a0.set_ylabel("work [kJ/mol]")
        #a0.legend()

    #fig.add_subplot(1,2,2)
    mean, variance, sigma = get_params(np.array(list_final_work))
    y = np.linspace(min(list_final_work), max(list_final_work), 300)

    a1.hist(list_final_work, density=True, histtype='stepfilled',
            orientation='horizontal', color='gray')
    a1.plot(norm.pdf(y, mean, sigma), y, color='red')
    
    a1.set_xlabel('Probability')

    plt.show()