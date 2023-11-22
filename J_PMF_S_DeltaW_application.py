import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from J_PMF_S_DeltaW import deltaW_definition
from J_PMF_S_DeltaW import Jarzynsky_for_PMF_deltaW
#from Free_energy_estimators.Utils import Convert_list


strb =  'unbinding'

ligand = 'CB8-G8'
path_root = '/home/eserra@iit.local/Crooks/HST-GST/%s'%ligand
simTime = '100ns'


path_dir =path_root+'/%s/Colvar_%s_%s/work_different_s_time'%(strb,simTime,strb)
output_file_path = path_root+'/%s/Colvar_%s_%s/delta_work'%(strb,simTime,strb)
input_prefix = 'work_s_'
output_name = 'delta_work_s'

# For the Unbinding initial_step = 3, final_step = 45, crescent = True. FOr Binding initial_step = 45 final_step = 3.

initial_step=3
final_step=45


if (strb=='binding'):
    t = final_step
    final_step = initial_step
    initial_step = t


binding = False
if (strb=='binding'):
    binding=True

# WRITE THE FILES WITH THE DELTA WORK
deltaW_definition(path_dir=path_dir,
                      input_prefix=input_prefix,
                      output_file_path=output_file_path,
                      output_name=output_name,
                      initial_step=initial_step,
                      final_step=final_step,
                      crescent =not binding ) 

wdir = output_file_path
file_prefix= output_name
Jarzynski_values_list = Jarzynsky_for_PMF_deltaW(final_step=final_step,
                                            wdir=wdir,
                                            file_prefix=file_prefix,
                                            initial_step=initial_step,
                                            crescent = not binding)
print(Jarzynski_values_list)

# PLOT THE PMF 

file = path_root+'/%s/Colvar_%s_%s/s_and_reference_time'%(strb,simTime,strb)

data=pd.read_csv(file, sep=',', header=1, names=['Plumed_step', 'time_ns', 'S(x)' ])
s_list = data['S(x)'].tolist()
x = s_list

if strb=='binding':
    x = x[::-1]


Jarzynski_values_list_sum = []
running_sum = 0

for i in range(len(Jarzynski_values_list)):
    running_sum += Jarzynski_values_list[i]
    Jarzynski_values_list_sum.append(running_sum)
print(Jarzynski_values_list_sum)

y = Jarzynski_values_list_sum

print(x)

plt.plot(x, y, color='green', linestyle='dashed', linewidth = 2,
         marker='o', markerfacecolor='blue', markersize=5)
  
plt.xlabel('S(x)')
plt.ylabel('PMF [kcal/mol]')
plt.title('Jarzynski Free Energy profile along the Collective Variable S(x)')
plt.savefig(f'{wdir}/PMF_Jarzynski_DeltaWork{strb}_{simTime}.png') 

plt.show()

                      