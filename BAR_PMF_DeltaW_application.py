import pandas as pd
import numpy as np
import matplotlib.pyplot as plt


from BAR_PMF_S_cumulativeW import BAR_for_PMF

ligand = 'Trypsine_benzamidine'
path_root = '/home/eserra@iit.local/work/%s'%ligand
simTime = '100ns'


path_bind =path_root+'/SteeredMD/binding/Colvar_%s_binding/work_different_s_time/'%(simTime)
path_unbind =path_root+'/SteeredMD/unbinding/Colvar_%s_unbinding/work_different_s_time/'%(simTime)
file_prefix = 'work_s_'
initial_step=3
final_step=54




#FUNCTION TO COMPUTE THE BAR Free ENERGY IN THE DIFFERENT WINDOWS WITH THE CUMULATIVE WORK 

BAR_values_list = BAR_for_PMF(final_s=final_step,
                            file_bind=path_bind,
                            file_unbind=path_unbind,
                            file_prefix=file_prefix,
                            initial_s=3
                                            )
print(BAR_values_list)

# PLOT THE PMF --> PLOT SALVATO NELLA DIRECTORY DI UNBINDING
file = f'{path_root}/SteeredMD/unbinding/Colvar_{simTime}_unbinding/work_different_s_time/s_and_reference_time'
data=pd.read_csv(file, sep=',', header=1, names=['Plumed_step', 'time_ns','S(x)' ])
s_list = data['S(x)'].tolist()
x = s_list

BAR_values_list_sum = []
running_sum = 0

for i in range(len(BAR_values_list)):
    running_sum += BAR_values_list[i]
    BAR_values_list_sum.append(running_sum)
print(BAR_values_list_sum)

y = BAR_values_list_sum

plt.plot(x, y, color='green', linestyle='dashed', linewidth = 2,
         marker='o', markerfacecolor='blue', markersize=5)
  
plt.xlabel('S(x)')
plt.ylabel('PMF [kcal/mol]')
plt.title('Free Energy profile applying BAR along the Collective Variable S(x)')
plt.savefig(f'{path_unbind}/PMF_BAR_DeltaW.png') 

plt.show()