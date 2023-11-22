#Application of the BAR algorithm 
#The W values are converted into np arrays 
#The binding W values are converted: opposite sign

import numpy as np
from pathlib import Path

from Utils import extract_list_work
from Utils import Convert_list
from Utils import Convert_kJ_kcal

from BAR_algorithm import bennett


ligand = 'CB8-G8'
path_root = '/home/eserra@iit.local/Crooks/HST-GST/%s'%ligand
simTime = '150ns'

file_binding = Path(f"{path_root}/binding/work_binding_{simTime}")
file_unbinding = Path(f"{path_root}/unbinding/work_unbinding_{simTime}")

list_work_bind = extract_list_work(file_binding, 'list_work_bind')
print(list_work_bind)
list_work_ubind = extract_list_work(file_unbinding, 'list_work_ubind')
print(list_work_ubind)

# CONVERT THE BINDING VALUES
converted_work = Convert_list(list_work_bind)
print('list_converted_work_binding:', converted_work)

# converting list to array
v0 = np.array(list_work_ubind) #foreword W: unbinding 
v1 = np.array(converted_work) #backward W: binding 


T = 298
toll = toll = 1e-6  # 0.000001
kB = 0.00831446261815324 # [kJ/(mol*K)]
beta = 1./(kB*T)

df = bennett(v0,v1,T,toll)
df = (1/beta)*df

df = Convert_kJ_kcal(df)
print('Free Energy [kcal/mol]:', df)
# It returns the value of the Free Energy in kcal/mol with the wrong sign