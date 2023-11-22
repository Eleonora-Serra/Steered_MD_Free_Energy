from pathlib import Path
from Utils import extract_list_work
from Utils import Convert_list
from Utils import Work_distribution_properties

from Jarzynski import Jarzyinski_function
from Jarzynski import cumulant_expansion


#directory with all colvar files, name of files, name of output

ligand = 'CB8-G8'
path_root = '/home/eserra@iit.local/Crooks/HST-GST/%s'%ligand
simTime = '150ns'

#file_binding = Path(f"{path_root}/binding/work_binding_{simTime}")
file_binding = '/home/eserra@iit.local/Crooks/HST-GST/CB8-G8/binding/Colvar_150ns_binding/work_different_s_time/work_s_4'

#file_unbinding = Path(f"{path_root}/unbinding/work_unbinding_{simTime}")

list_work_bind = extract_list_work(file_binding, 'list_work_bind')
print(list_work_bind)
#list_work_ubind = extract_list_work(file_unbinding, 'list_work_ubind')
#print(list_work_ubind)

#QUICK INFORMATIONS ABOUT THE W DISTRIBUTIONS
Work_distribution_properties(list_work=list_work_bind, binding=True)
#Work_distribution_properties(list_work=list_work_ubind, binding=False)

# APLLY JARZYNSKI
Free_Energy_binding = Jarzyinski_function(list_work =list_work_bind, invert= False)
#Free_Energy_ubinding = Jarzyinski_function(list_work =list_work_ubind, invert= False)

# EXPONENTIAL AVERAGE COMPUTED AS CUMULANT EXPANSION
Free_Energy_binding_cumulant = cumulant_expansion(list_work =list_work_bind, invert= False)
#Free_Energy_unbinding_cumulant = cumulant_expansion(list_work =list_work_ubind, invert= False)