#Utilization of the Plot functions

from pathlib import Path
from Colvar_plot import plot_work_and_distribution
from Colvar_plot import plt_work_vs_s, plt_work_vs_time, plt_s_profile

from Free_energy_estimators.Utils import extract_list_work, Convert_list


#wdir = '/home/eserra@iit.local/work/Trypsine_benzamidine/SteeredMD/unbinding/Colvar_100ns_unbinding/'
wdir = '/home/eserra@iit.local/work/Abl_Gleevec/SteeredMD/unbinding/Colvar_200ns_unbinding_newPath/'

nrun=30
#file_prefix = 'colvar_'
file_prefix='colvar_unb_correct_'


# Plot the W profiles for multiple replica and plot the W distribution in the same plot
file_unb_100 = f'{wdir}Unique_colv_unbinding_200.txt'
list_work_unbind_100 = extract_list_work(file_unb_100, 'list_work_bind_100')
#plot_work_and_distribution(nrun=50, wdir=wdir, file_prefix=file_prefix, column_name='restraint.work', list_final_work=list_work_unbind_100)

#Plot the W profile as function of CV:s
plt_work_vs_s(nrun=nrun, wdir=wdir, file_prefix=file_prefix, column_name='restraint.work')

#Plot the W profile as function of time:s
plt_work_vs_time(nrun=nrun, wdir=wdir, file_prefix=file_prefix, column_name='restraint.work')

#Plot the CV profile as function of time 
plt_s_profile(nrun=nrun, wdir=wdir, file_prefix=file_prefix)

