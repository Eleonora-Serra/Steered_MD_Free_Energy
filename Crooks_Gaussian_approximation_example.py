from pathlib import Path
from Free_energy_estimators.Utils import extract_list_work
from Free_energy_estimators.Utils import Convert_list
from Free_energy_estimators.Utils import Work_distribution_properties

from Free_energy_estimators.Crooks_Gaussian_approximation import real_distribution_plot
from Free_energy_estimators.Crooks_Gaussian_approximation import gaussian_distribution_plot
from Free_energy_estimators.Crooks_Gaussian_approximation import gaussian_distributions
from Free_energy_estimators.Crooks_Gaussian_approximation import gaussian_intersection_plot_location

#directory with all colvar files, name of files, name of output

ligand = 'Abl_Gleevec'
path_root = '/home/eserra@iit.local/work/%s/SteeredMD'%(ligand)
simTime = '50ns'

file_binding = Path(f"{path_root}/binding/Colvar_50ns_binding/Unique_colv_binding_{simTime}.txt")
file_unbinding = Path(f"{path_root}/unbinding/Colvar_50ns_unbinding/Unique_colv_unbinding_{simTime}.txt")

list_work_bind = extract_list_work(file_binding, 'list_work_bind')
list_work_unbind = extract_list_work(file_unbinding, 'list_work_unbind')

#print(list_work_bind)
#print(list_work_ubind)


real_distribution_plot(list_work=list_work_bind, binding=True)
#Work_distribution_properties(list_work=list_work_ubind, binding=False)

gaussian_distribution_plot(list_work=list_work_bind, binding=True, save = False)

works_list = [list_work_bind, list_work_unbind]
gaussian_distributions(works_list=works_list,save = False)

gaussian_intersection_plot_location (works_list=works_list)

