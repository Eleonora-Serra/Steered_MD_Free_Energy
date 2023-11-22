from bootstrap import bootstrap_Jarzynski_binding_energy
import matplotlib.pyplot as plt
import numpy as np

####################################################
###INPUTS - INFO ABOUT SIMULATIONS PERFORMED###
#time of simulation
sim_time = 150

#number of bootstrap sampling to be performed
n_bootstrap_sampling = 50

#number of simulation replicas performed
replicas = 50

#initial and final values of s (PCV)
initial_s = 3
final_s = 45

#frame discriminating between the bound and unbound region
discrim_frame=10
##################################################


#define work directories and name of files (usally not to be modified)
wdir = '/Users/alessiaghidini/Desktop/steered_md_on_predetermined_pcvs/steered-md'
name_file = 'cb8-g8/unbinding-%dns/work_different_s_time/work_s_'%sim_time


#define start, stop, and step values of s (PCV)
start = initial_s
stop = final_s+1
step = +1

invert = False

'''
if process == 'binding':
    
    start = initial_step
    stop = final_step-1
    step = -1
    invert = True
'''

### calculate bootstrap binding free energy with Jarzynski estimator ### 

mu, sigma = bootstrap_Jarzynski_binding_energy(i=replicas,
                                               start=start,
                                               stop=stop,
                                               step=step,
                                               n_bootstrap_sampling=n_bootstrap_sampling,
                                               wdir=wdir,
                                               name_file=name_file,
                                               discrim_frame=discrim_frame,
                                               invert=invert
                                               )
print(f'Bootstrap binding free energy = {round(mu,2)} +/-  {round(3*sigma,2)} kcal/mol')


'''
y_mean, y_std = zip(*[bootstrap_Jarzynski_binding_energy(i,
                                                        start=start,
                                                        stop=stop,
                                                        step=step,
                                                        n_bootstrap_sampling=n_bootstrap_sampling,
                                                        wdir=wdir,
                                                        name_file=name_file,
                                                        discrim_frame=discrim_frame,
                                                        invert=invert) for i in range(replicas+1)])

plt.xlabel('number of replicas')
plt.ylabel('Bootstrap binding Free Energy [kcal/mol]')
plt.errorbar(np.arange(replicas+1), y_mean, yerr=y_std, fmt='.k')
plt.savefig(f'bootstrap_n_replicas_{process}_{sim_time}.png') 

plt.show()
'''