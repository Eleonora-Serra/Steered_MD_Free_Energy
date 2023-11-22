import pandas as pd
import matplotlib.pyplot as plt

from Free_energy_estimators.cumulative_work_at_different_s_values import work_at_different_s_values
from BAR_PMF_S_cumulativeW import BAR_for_PMF
from Ratio_partition_functions import binding_free_energy
from volume_correction import absolute_binding_free_energy

####################################################
###INPUTS - INFO ABOUT SIMULATIONS PERFORMED###

#V_unbound = 1336.03 #A^3
V_unbound = 1651.2
calculate_absolute_binding_free_energy = False

#ligand of simulation
ligand = 'Abl_Gleevec'
#ligand = 'Trypsine_benzamidine'
#ligand = 'CB8-G8'

#time of simulation
sim_time = 200

#number of simulation replicas performed
replicas = 30

#frame discriminating between the bound and unbound regions
discrim_frame = 22

##################################################
###INPUT AND OUTPUT DIRECTORIES AND PREFIX###
wdir = '/home/eserra@iit.local/work/%s/SteeredMD'%(ligand)
#wdir = '/home/eserra@iit.local/work/Crooks/HST-GST/%s/'%(ligand)


input_files_path_binding = wdir+'/binding/Colvar_%sns_binding_newPath/'%(sim_time)

input_files_path_unbinding = wdir+'/unbinding/Colvar_%sns_unbinding_newPath/'%(sim_time)
colvar_files_prefix = 'colvar_'

output_files_path_binding = wdir+'/binding/Colvar_%sns_binding_newPath/work_different_s_time/'%(sim_time)
output_files_path_unbinding = wdir+'/unbinding/Colvar_%sns_unbinding_newPath/work_different_s_time/'%(sim_time)


output_prefix = 'work_s_'



#plumed file used as input for simulations
plumed_dat_file_binding=input_files_path_binding+'abl_gleevec_binding_200ns.dat'
plumed_dat_file_unbinding=input_files_path_unbinding+'abl_gleevec_unbinding_200ns.dat'

#plumed_dat_file_binding=wdir+'Utils/G8_bound_%sns.dat'%(sim_time)
#plumed_dat_file_unbinding=wdir+'/G8_unbound_%sns.dat'%(sim_time)


##################################################
###VALUES OF PCV S DURING SIMULATION###
'''N.B.: UNBINDING: S INCREASES; BINDING: S DECREASES
we assume that in unbinding simulations the value of s increases,
while it decreases in binding simulations, thus initial and final values of s are inverted by default'''

#initial and final values of s (PCV)
initial_s = 1
final_s = 33

##################################################
# to make work_s_ files with cumulative W set make_cumulative_work_files = True
make_cumulative_work_files = False

##################################################

# create a work_s_ file with W values from 50 replicas for each s
if make_cumulative_work_files == True:
    #file for binding simulation
    work_at_different_s_values(replicas=replicas,
                            input_files_path=input_files_path_binding,
                            colvar_file=colvar_files_prefix,
                            output_files_path=output_files_path_binding,
                            output_prefix=output_prefix,
                            plumed_dat_file = plumed_dat_file_binding,
                            initial_s=final_s,
                            final_s=initial_s,
                            invert=True, 
                            )
    
    #file for unbinding simulation
    work_at_different_s_values(replicas=replicas,
                            input_files_path=input_files_path_unbinding,
                            colvar_file=colvar_files_prefix,
                            output_files_path=output_files_path_unbinding,
                            output_prefix=output_prefix,
                            plumed_dat_file = plumed_dat_file_unbinding,
                            initial_s=initial_s,
                            final_s=final_s,
                            invert=False, 
                            ) 
    

# function to calculate the free energy difference associated to each value of S using BAR on cumulative W values

BAR_values_list = BAR_for_PMF(initial_s=initial_s,
                              final_s=final_s,
                              file_bind=output_files_path_binding,
                              file_unbind=output_files_path_unbinding,
                              file_prefix=output_prefix,
                                            )
print(BAR_values_list)
print(len(BAR_values_list))

# # Plot the PMF obtained using BAR and create a csv with PMF (BAR binding free energy vs S)
file = f'{output_files_path_unbinding}/s_and_reference_time'
data=pd.read_csv(file, sep=',', header=0, names=['Plumed_step', 'time_ns','S(x)' ])
s_list = data['S(x)'].tolist()
x = [(i-min(s_list))/(max(s_list)-min(s_list)) for i in s_list]
y = BAR_values_list

file_BAR_FreeEnergy = output_files_path_unbinding +'BAR_FreeEnergy_%s.csv'%(sim_time)
BAR_df = pd.DataFrame({'S(x)':x, 'BAR_FreeEnergy':BAR_values_list})
BAR_df.to_csv(file_BAR_FreeEnergy, header = True, sep=',', index = None)

plt.plot(x, y, color='green', linestyle='dashed', linewidth = 2,
         marker='o', markerfacecolor='blue', markersize=5)
  
plt.xlabel('S(x)')
plt.ylabel('PMF [kcal/mol]')
plt.title('BAR Free Energy profile along the Collective Variable S(x)')
plt.savefig(f'{output_files_path_unbinding}/PMF_BAR_FreeEnergy_{sim_time}.png') 

# compute the BAR free energy difference as the ratio of the partition function of the two different states

BAR_Free_Energy_Binding = binding_free_energy(file_FreeEnergy=file_BAR_FreeEnergy,
                                                    discrim_frame=discrim_frame,
                                                    invert=False)
print(f'BAR Free Energy difference = {round(BAR_Free_Energy_Binding,3)} kcal/mol')

# compute the BAR standard free energy difference 

if calculate_absolute_binding_free_energy == True:
    BAR_Absolute_Free_Energy_Binding = absolute_binding_free_energy(V_unbound=V_unbound, binding_free_energy=BAR_Free_Energy_Binding)
    print(f'BAR Absolute Free Energy difference = {round(BAR_Absolute_Free_Energy_Binding,3)} kcal/mol')