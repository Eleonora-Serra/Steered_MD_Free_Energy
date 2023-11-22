import matplotlib.pyplot as plt
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import norm
import statistics
import seaborn as sns

# Import Function to extract W values from a file and make a list of them
# Import function to convert the list 

#Function to calculate mean and standard deviation of W values and plot the real distribution 

def real_distribution_plot (list_work:list, binding = True):
    '''
    Function to plot the "real distributions" as histrograms without the Gaussian assumption.
    It uses the default parameters selected by sns as the ones for the bin sized
    '''
    #fig = plt.figure()
    sns.displot(list_work, kde=True , legend=True)
    plt.xlabel('Work values (kJ/mol)')
    plt.ylabel('Count')
    if binding == True:
        plt.legend(labels = ["Binding"])
    if binding == False:
        plt.legend(labels = ["Unbinding"])
    plt.show()

def gaussian_distribution_plot(list_work:list, binding=False, save = False):
    '''
    Fuction that takes as input a list of Work values and plot the normal distribution.
    binding and save are by default False. 
    if binding= True we are considering binding W values and the Legend of the plot changes.
    if save = True the plot will be saved as imagine 
    '''
    mean = statistics.mean(list_work)
    sd = statistics.stdev(list_work)

    x_axis = np.arange(-1000, 1000, 0.01) + mean
    ics = x_axis
    ipsilon = norm.pdf(ics, mean, sd)
    plt.plot(ics,ipsilon)
    plt.xlabel('Work values (kJ/mol)')
    plt.ylabel('Frequency')
    if binding == True:
        plt.legend(labels = ['Binding'])
        if save == True:
            plt.savefig('Gaussian_distribution_Binding.png')
    if binding == False:
        plt.legend(labels=['Unbinding'])
        if save == True:
            plt.savefig('Gaussian_distribution_Unbinding.png')
    plt.show()



def gaussian_distributions (works_list, save = False):
    '''
    Function that plots the 2 W distributions together. The W values are supplied in a list composed by
    the 2 W lists. The first list supplied is the binding_list.
    Arg:works_list = List of the 2 lists of work, i.e. [list_work_bind, list_work_ubind]
        To have the correct legend the first list passed to the function is the binding one.
        save= False by default 
    '''
    lista =['Binding', 'Unbinding']
    for i,w in enumerate(works_list):
        x_axis = np.arange(-1000, 1000, 0.01) 
        ics = x_axis
        mean = statistics.mean(w)
        sd = statistics.stdev(w)
        ipsilon = norm.pdf(ics, mean, sd)
        plt.plot(ics, ipsilon, label=f'{lista[i]}' )
    plt.xlabel('Work values (kJ/mol)')
    plt.ylabel('Frequency')
    plt.legend()
    if save == True:
            plt.savefig('Gaussian_distributions.png')
    plt.show()


def solve(m1,m2,std1,std2):
    '''
    Find the intersectio between 2 gausian distributions: solve the system and return the solutions
    '''
    a = 1/(2*std1**2) - 1/(2*std2**2)
    b = m2/(std2**2) - m1/(std1**2)
    c = m1**2 /(2*std1**2) - m2**2 / (2*std2**2) - np.log(std2/std1)
    return np.roots([a,b,c])

def gaussian_intersection(works_list):
    '''
    It takes as input the list of the W value lists [list_work_bind, list_work_ubind].
    It returns the intersection points between the 2 gaussians: It returns values and not a plot. 
    One of the 2 solutions is the binding Free Energy with the opposite sign.
    '''
    m = []
    std = []
    for w in works_list:
        x_axis = np.arange(-1000, 1000, 0.01) 
        m.append(statistics.mean(w))
        std.append(statistics.stdev(w))
    
    return solve(m[0],m[1],std[0],std[1])

def gaussian_intersection_plot(works_list, save = False):
    '''
    It plots the intersection points between the 2 distributions
    Save = False by default
    '''
    m1 = statistics.mean(works_list[1])
    std1 = statistics.stdev(works_list[1])
    lista =['Binding', 'Unbinding']
    for i,w in enumerate(works_list):
        x_axis = np.arange(-1000, 1000, 0.01) 
        ics = x_axis
        mean = statistics.mean(w)
        sd = statistics.stdev(w)
        ipsilon = norm.pdf(ics, mean, sd)
        plt.plot(ics, ipsilon, label=f'{lista[i]}' )

    plt.plot(gaussian_intersection(works_list=works_list), norm.pdf(gaussian_intersection(works_list=works_list),m1,std1), 'o')
    plt.xlabel('Work values (kJ/mol)')
    plt.ylabel('Frequency')
    plt.legend()
    if save == True:
            plt.savefig('Gaussian_distributions.png')
    plt.show()

def get_intersection_locations(y1,y2,test=False,x=None): 
    """
    return indices of the intersection point/s.
    """
    idxs=np.argwhere(np.diff(np.sign(y1 - y2))).flatten()
    if test:
        x=range(len(y1)) if x is None else x
        plt.figure(figsize=[10,10])
        ax=plt.subplot()
        ax.plot(x,y1,color='r',label='Binding',alpha=0.5)
        ax.plot(x,y2,color='b',label='Unbinding',alpha=0.5)
        _=[ax.axvline(x[i],color='k') for i in idxs]
        _=[ax.text(x[i],ax.get_ylim()[1],f"{x[i]:1.1f}",ha='center',va='bottom') for i in idxs]
        ax.legend(bbox_to_anchor=[1,1])
        ax.set(xlabel='Work values (kJ/mol)',ylabel='density')
        plt.show()
    return idxs

def gaussian_intersection_plot_location (works_list):
    '''
    It plots the 2 distributions, the intersection points and the precise location of them on the plot.
    '''
    x = np.arange(-1000, 1000, 0.01) 
    m = []
    std = []
    for w in works_list: 
        m.append(statistics.mean(w))
        std.append(statistics.stdev(w))
    y1 = norm.pdf(x, m[0], std[0])
    y2 = norm.pdf(x, m[1], std[1])
    get_intersection_locations(y1=y1,y2=y2,x=x,test=True) # returns indice/s array([10173])
 

