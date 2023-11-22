import numpy as np
import scipy.integrate as spi
import random
import scipy as sc

from Free_energy_estimators.Utils import extract_list_work, Convert_kJ_kcal, Convert_list
from Free_energy_estimators.Jarzynski import Jarzyinski_function
from Free_energy_estimators.BAR_algorithm import bennett


### constants ###
T = 298
kb = sc.constants.Boltzmann # 1.380649e-23 J
Na = sc.constants.Avogadro # 6.02214076e+23 mol-1
kbT_kJ_mol = kb*0.001*Na*T
beta = (1/(kbT_kJ_mol))
toll = 1e-6 

### bootrstap_randomize used in bootstrap_Jarzynski_binding_energy and bootstrap_BAR_binding_energy_test ###
### randomize two list of work values ###
def bootstrap_randomize(list_work):
    vector = []
    randomized_list = []
    for i in range(0,len(list_work)):
        vector.append(random.randint(0, len(list_work)-1))
        randomized_list=[list_work[x] for x in vector]
    return randomized_list

### calculate the bootstrap binding energy for Jarzynski estimator ### 
def bootstrap_Jarzynski_binding_energy(i, start, stop, step, n_bootstrap_sampling, wdir, name_file, discrim_frame, invert):
    s_list = range(start,stop,step)
    Free_Energy_Binding_list = []
    for n in range(0, n_bootstrap_sampling):
        PMF_list = []
        for s in range(start,stop,step):
            file = f'{wdir}/{name_file}{s}'
            work_list = extract_list_work(file, 'work_list')
            randomized_list = bootstrap_randomize(list_work=work_list)[:i+1]
            PMF_list.append(Jarzyinski_function(list_work=randomized_list, invert=invert))

        x = [(i-min(s_list))/(max(s_list)-min(s_list)) for i in s_list]
        y = np.exp(-beta*np.array(PMF_list))
        
        bound_region = spi.simpson(y[:discrim_frame],x[:discrim_frame])
        unbound_region = spi.simpson(y[discrim_frame:],x[discrim_frame:])
        Free_Energy_Binding = -1/beta*np.log(bound_region/unbound_region) 
        Free_Energy_Binding_list.append(Free_Energy_Binding)
    return np.mean(Free_Energy_Binding_list), np.std(Free_Energy_Binding_list)


### double_bootrstap_randomize used in bootstrap_BAR_binding_energy ###
### randomize two list of work values ###
def double_bootstrap_randomize(list_work1, list_work2):
    vector1 = []
    vector2 = []
    randomized_list_1 = []
    randomized_list_2 = []
    for i in range(0,len(list_work1)):
        vector1.append(random.randint(0, len(list_work1)-1))
        randomized_list_1=[list_work1[x] for x in vector1]
    for i in range(0,len(list_work2)):
        vector2.append(random.randint(0, len(list_work2)-1))
        randomized_list_2=[list_work2[x] for x in vector2]
    return randomized_list_1,randomized_list_2


### calculate the bootstrap binding energy for BAR estimator ### 
def bootstrap_BAR_binding_energy(i, start, stop, step, n_bootstrap_sampling, wdir, name_file_binding, name_file_unbinding, discrim_frame):
    s_list = range(start,stop,step)
    Free_Energy_Binding_list = []

    for n in range(0,n_bootstrap_sampling):
        PMF_list = []
        for s in range(start,stop,step):
            work_unbind_list = extract_list_work(f'{wdir}/{name_file_unbinding}{s}', 'work_bind_list')
            work_bind_list = Convert_list(extract_list_work(f'{wdir}/{name_file_binding}{s}', 'work_unbind_list'))
            v0 = np.array(work_unbind_list) #forward W: unbinding 
            v1 = np.array(work_bind_list) #backward W: binding 
            randomized_list_1, randomized_list_2 = double_bootstrap_randomize(list_work1=v0, list_work2=v1)
            randomized_list_1 = np.array(randomized_list_1)[:i+1]
            randomized_list_2 = np.array(randomized_list_2)[:i+1]
            PMF = (1/beta)*bennett(randomized_list_1,randomized_list_2,T)
            PMF_list.append(Convert_kJ_kcal(PMF))

        x = [(i-min(s_list))/(max(s_list)-min(s_list)) for i in s_list]
        y = np.exp(-beta*np.array(PMF_list))
        bound_region = spi.simpson(y[:discrim_frame],x[:discrim_frame])
        unbound_region = spi.simpson(y[discrim_frame:],x[discrim_frame:])
        Free_Energy_Binding = -1/beta*np.log(bound_region/unbound_region) 
        Free_Energy_Binding_list.append(Free_Energy_Binding)
    return np.mean(Free_Energy_Binding_list), np.std(Free_Energy_Binding_list)


### calculate the bootstrap binding energy for BAR estimator but using bootstrap_randomize ### 
def bootstrap_BAR_binding_energy_test(i, start, stop, step, n_bootstrap_sampling, wdir, name_file_binding, name_file_unbinding, discrim_frame):
    s_list = range(start,stop,step)
    Free_Energy_Binding_list = []

    for n in range(0,n_bootstrap_sampling):
        PMF_list = []
        for s in range(start,stop,step):
            file_bind = f'{wdir}/{name_file_binding}{s}'
            file_unbind = f'{wdir}/{name_file_unbinding}{s}'
            work_bind_list = extract_list_work(file_bind, 'work_bind_list')
            work_unbind_list = extract_list_work(file_unbind, 'work_unbind_list')
            converted_work_bind = Convert_list(work_bind_list)
            v0 = np.array(work_unbind_list) #forward W: unbinding 
            v1 = np.array(converted_work_bind) #backward W: binding 
            randomized_list_1 = np.array(bootstrap_randomize(list_work=v0))[:i+1]
            randomized_list_2 = np.array(bootstrap_randomize(list_work=v1))[:i+1]
            PMF = bennett(randomized_list_1,randomized_list_2,T,toll)
            PMF_binding = (1/beta)*PMF
            PMF_list.append(Convert_kJ_kcal(PMF_binding))

        x = [(i-min(s_list))/(max(s_list)-min(s_list)) for i in s_list]
        y = np.exp(-beta*np.array(PMF_list))
        
        bound_region = spi.simpson(y[:discrim_frame],x[:discrim_frame])
        unbound_region = spi.simpson(y[discrim_frame:],x[discrim_frame:])
        Free_Energy_Binding = -1/beta*np.log(bound_region/unbound_region) 
        Free_Energy_Binding_list.append(Free_Energy_Binding)
    return np.mean(Free_Energy_Binding_list), np.std(Free_Energy_Binding_list)