import pandas as pd
import matplotlib.pyplot as plt

from Free_energy_estimators.cumulative_work_at_different_s_values import work_at_different_s_values
from J_PMF_S_cumulativeW import Jarzynsky_for_PMF
from Ratio_partition_functions import binding_free_energy
from volume_correction import absolute_binding_free_energy

###INPUTS - INFO ABOUT SIMULATIONS PERFORMED###

#V_unbound = 1336.03 #A^3
V_unbound = 1651.2
calculate_absolute_binding_free_energy = False

#ligand of simulation
ligand = 'Abl_Gleevec'
#ligand = 'CB8-G8'

#time of simulation
sim_time = 200

#type of simulation (binding or unbinding)
simulation = 'binding'

#number of simulation replicas performed
replicas = 30

#frame discriminating between the bound and unbound regions
discrim_frame = 23

##################################################
###INPUT AND OUTPUT DIRECTORIES AND PREFIX###

wdir = '/home/eserra@iit.local/work/%s/SteeredMD'%(ligand)
#wdir = '/home/eserra@iit.local/work/Crooks/HST-GST/%s/'%(ligand)


input_files_path = wdir+'/%s/Colvar_%sns_%s_newPath/'%(simulation,sim_time,simulation)

colvar_files_prefix = 'colvar_bind_correct_'

output_files_path = wdir+'/%s/Colvar_%sns_%s_newPath/work_different_s_time/'%(simulation,sim_time,simulation)

output_prefix = 'work_s_'


#plumed file used as input for simulations
plumed_dat_file=input_files_path+'abl_gleevec_binding_200ns.dat'

#plumed_dat_file=wdir+'Utils/G8_bound_%sns.dat'%(sim_time)

##################################################
###VALUES OF PCV S DURING SIMULATION###
'''N.B.: UNBINDING: S INCREASES; BINDING: S DECREASES
we assume that in unbinding simulations the value of s increases,
while it decreases in binding simulations, thus initial and final values of s are inverted by default'''

#initial and final values of s (PCV)

initial_s = 1
final_s = 33

#if invert = True, then values of final_s and initial_s are inverted
invert = False
#by default, invert = True for binding simulations - comment this line if this is not the case
if simulation == 'binding':
    invert = True

##################################################
# to make work_s_ files with cumulative W reading time set make_cumulative_work_files_reading_time = True
# to make work_s_ files with cumulative W reading s set make_cumulative_work_files_reading_s = True
make_cumulative_work_files_reading_time = True
'''make_cumulative_work_files_reading_s = True

#if make_cumulative_work_files_reading_s is used, the STEP must be defined!'''
step_reading_s = 0.5
    
##################################################

#inversion of final and initial values of s (default for binding simulations)
if (invert == True):
    t = final_s
    final_s = initial_s
    initial_s = t


# create a work_s_ file with W values from 50 replicas for each s

#if make_cumulative_work_files_reading_time is used, the STEP is always 1
if make_cumulative_work_files_reading_time == True:
    step = 1
    work_at_different_s_values(replicas=replicas,
                            input_files_path=input_files_path,
                            colvar_file=colvar_files_prefix,
                            output_files_path=output_files_path,
                            output_prefix=output_prefix,
                            plumed_dat_file = plumed_dat_file,
                            initial_s=initial_s,
                            final_s=final_s,
                            invert=invert, 
                            )
'''   
#if make_cumulative_work_files_reading_s is used, the STEP must be defined! and different output_files_path is used    
if make_cumulative_work_files_reading_s == True: 
    step = step_reading_s
    output_files_path = wdir+'/%s-%sns/work_different_s_first_appearance/'%(simulation,sim_time)   
    work_at_different_s_values_first_appearance (replicas=replicas,
                                input_files_path=input_files_path,
                                colvar_file=colvar_files_prefix,
                                output_files_path=output_files_path,
                                output_prefix=output_prefix,
                                initial_s=initial_s,
                                final_s=final_s,
                                invert=invert,
                                step=step,
                               )'''


# function to calculate the free energy difference associated to each value of S using Jarzynski on cumulative W values

Jarzynski_values_list = Jarzynsky_for_PMF(initial_s=initial_s,
                                          final_s=final_s,
                                          output_files_path=output_files_path,
                                          output_prefix=output_prefix,
                                          invert=invert)
print(Jarzynski_values_list)


# Plot the PMF obtained using Jazrynski and create a csv with PMF (Jarzynski binding free energy vs S)
file = f'{output_files_path}/s_and_reference_time'
data = pd.read_csv(file, sep=',', header=0, names=['Plumed_step', 'time_ns', 'S(x)' ])
s_list = data['S(x)'].tolist()
x = [(i-min(s_list))/(max(s_list)-min(s_list)) for i in s_list]
y = Jarzynski_values_list

file_Jarzynski_FreeEnergy = output_files_path +'Jarzynski_FreeEnergy_%s_%s.csv'%(simulation,sim_time)
Jarzynski_df = pd.DataFrame({'S(x)':x, 'Jarzynski_FreeEnergy':Jarzynski_values_list})
Jarzynski_df.to_csv(file_Jarzynski_FreeEnergy, header = True, sep=',', index = None)

plt.plot(x, y, color='green', linestyle='dashed', linewidth = 2,
         marker='o', markerfacecolor='blue', markersize=5)
  
plt.xlabel('S(x)')
plt.ylabel('PMF [kcal/mol]')
plt.title('Jarzynski Free Energy profile along the Collective Variable S(x)')
plt.savefig(f'{output_files_path}/PMF_Jarzynski_FreeEnergy_{simulation}_{sim_time}.png') 
plt.close()


# compute the Jarzynski free energy difference as the ratio of the partition function of the two different states

Jarzynski_Free_Energy_Binding = binding_free_energy(file_FreeEnergy=file_Jarzynski_FreeEnergy,
                                                    discrim_frame=discrim_frame,
                                                    invert=invert)
print(f'Jarzynski Free Energy difference = {round(Jarzynski_Free_Energy_Binding,3)} kcal/mol')

# compute the Jarzynski standard free energy difference 

if calculate_absolute_binding_free_energy == True:
    BAR_Absolute_Free_Energy_Binding = absolute_binding_free_energy(V_unbound=V_unbound, binding_free_energy=Jarzynski_Free_Energy_Binding)
    print(f'Jarzynski Absolute Free Energy difference = {round(BAR_Absolute_Free_Energy_Binding,3)} kcal/mol')


### to compute free energy for several frames ###

#free_energy_frames = []
#i_frame = 40
#f_frame = 120

'''for frame in range(i_frame, f_frame):
    Jarzynski_Free_Energy_Binding = binding_free_energy(file_FreeEnergy=file_Jarzynski_FreeEnergy,
                                                    discrim_frame=frame,
                                                    invert=invert)    
    free_energy_frames.append(Jarzynski_Free_Energy_Binding)

frames = list(range(i_frame, f_frame))
discrim_x = [(frame-initial_s)/(final_s-initial_s) for frame in frames]
file_Jarzynski_FreeEnergy_vs_frame = output_files_path +'/Jarzynski_FreeEnergy_%s_%s_vs_frame.csv'%(simulation,sim_time)
Jarzynski_df = pd.DataFrame({'S(x)*':discrim_x, 'x*': frames, 'Jarzynski_FreeEnergy':free_energy_frames})
Jarzynski_df.to_csv(file_Jarzynski_FreeEnergy_vs_frame, header = True, sep=',', index = None)
print(f'x:{discrim_x}, free_energy:{free_energy_frames}')
plt.plot(discrim_x,free_energy_frames, color='green', linewidth = 2, marker='o', markerfacecolor='blue', markersize=2)
plt.xlabel('discriminating frame')
plt.ylabel('Free Energy [kcal/mol]')
plt.savefig(f'{output_files_path}/Jarzynski_FreeEnergy_{simulation}_{sim_time}_vs_frame.png')'''