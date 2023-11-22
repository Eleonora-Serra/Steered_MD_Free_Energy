# General functios

import scipy as sc
import numpy as np
import scipy.special as scp
from pathlib import Path

# Function to extract W values from a file text with 2 columns 
def extract_list_work (file:Path, 
                       output_list:str):
    '''
    The function takes as arguments the file text containing the W values and the name of the output file 
    as string. 
    It returns a list of work values 
    '''
    output_list = []
    with open(file, 'r') as filin:
            for line in filin:
                line_str= line.split(',')
                last= line_str[-1]
                last=float(last)
                output_list.append(last)
    return (output_list)

# Function to convert and invert the sign of a list
def Convert_list(lst:list):
    return [ -i for i in lst ]

def Convert_kJ_kcal(kJ_value:float):
     return kJ_value/4.184

# Quicly define the W distribution properties
def Work_distribution_properties (list_work:list,
                                    binding=True):
    '''
    Function to quicly define the W distribution properties. By default binding = True.
    It takes a list as input and convert it into a numpy array.
    '''
    list_work = np.array(list_work)
    mean =  np.mean(list_work)
    var = np.var(list_work)
    if binding==True:
         print("Total Work from binding simulations:\nMean:",mean, "Variance:", var)
    if binding ==False:
         print("Total Work from unbinding simulations:\nMean:",mean, "Variance:", var)
