# # IMP!: before running this file, create a new conda environment with the command: `conda env create --file environment.yaml --name emergence_patterns`
# # and switch to that environment with `conda activate emergence_patterns` or running this file.

import numpy as np
import os
from scipy.spatial import ConvexHull
from itertools import combinations
from math import factorial
from time import perf_counter
from multiprocessing import Pool
from numba import njit
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib import rcParams

def equilibrium2(payoff_vector, indices):
    vector = np.zeros(payoff_vector.shape)
    common_factor = 1/np.sum(1/payoff_vector[indices])
    vector[indices] = common_factor/payoff_vector[indices]
    return vector[1:]

@njit
def equilibrium(payoff_vector, indices):
    common_denominator = 0
    vector = np.zeros(payoff_vector.shape)
    for i in indices:
        term = -1
        for j in indices:
            if i==j:
                continue
            else:
                if term==-1:
                    term = payoff_vector[j]
                else:
                    term *= payoff_vector[j]
        vector[i] = term
        common_denominator += term
    for i in indices:
        vector[i] /= common_denominator
    return vector[1:]

def fractional_volume(M,w_range,index):
    frac_vols = np.zeros(w_range.size)
    for count,w in enumerate(w_range):
        A = 1 - w/2 + w/2*np.sin(4*np.pi*np.arange(M)/M)
        A = np.abs(A*(A>0) + 0.001*(A<=0))
        total_volume = 1/factorial(M-1)
        A_shifted = np.roll(A,-index)
        
        num_points = 2**(M-1)
        points = np.zeros([num_points,M-1])
        counter = 1
        
        indices = np.arange(1,M,1)
        for k in range(1,M,1):
            for comb in combinations(indices,k):
                points[counter,:] = equilibrium(A_shifted,(0,*comb))
                counter += 1
        
        if np.all(points==0.):
            frac_vols[count] = 0.
        else:
            hull = ConvexHull(points)
            frac_vols[count] = hull.volume/total_volume
    return frac_vols

if __name__=="__main__":

    M = 7  # number of total possible types
    w_numpts = 100 # number of points of w, the strength of selection
    w_min = 0 # minimum value of selection strength
    w_max = 1.5 # maximum value of selection strength; limit: 1
    alive_threshold = 0.01 # threshold in fractional volume to consider the type present in the system at max selection strength
    parallel_process = True
    
    processes = 9 
    
    w_range = np.linspace(w_min,w_max,w_numpts)
    volumes = np.zeros([w_numpts,M])
    
    indices_forward = []
    indices_backward = np.zeros([M],dtype=np.int64)
    A = np.sin(4*np.pi*np.arange(M)/M)
    for p in range(M):
        if np.any(np.isclose(A[p],A[:p])):
            indices_backward[p] = indices_forward.index(np.where(np.isclose(A[p],A[:p]))[0][0])
            continue
        else:
            indices_forward.append(p)
            indices_backward[p] = len(indices_forward) - 1
    
    num_calcs = len(indices_forward)
    volume_part = np.zeros([w_numpts,num_calcs])
    start = perf_counter()
    
    if parallel_process:
        parameters = zip([M]*num_calcs,[w_range]*num_calcs,indices_forward)
        with Pool(processes=processes) as p:
            results = p.starmap_async(fractional_volume,parameters)
            p.close()
            p.join()
        
        if hasattr(results, 'get'):
            results = results.get()
        
        for i,arr in enumerate(results):
            volume_part[:,i] = arr
    else:
        for count,p in enumerate(indices_forward):
            volume_part[:,count] = fractional_volume(M,w_range,p)
    
    volumes = volume_part[:,indices_backward]
    
    end = perf_counter()
    print(f"Time taken: {round(end-start,2)} s")
    
    # rcParams['figure.dpi'] = 600
    rcParams['font.family'] = 'Poppins'
    title_pars = {'fontsize':16, 'fontweight':'semibold'}
    label_pars = {'fontsize':14, 'fontweight':'book'}
    tick_pars = {'labelsize':13} 
    cb_label_pars = {'fontsize':13, 'fontweight':'normal'}
    cb_tick_pars = {'labelsize':12}
    
    c = np.sin(4*np.pi*np.arange(M)/M)
    norm = mpl.colors.Normalize(vmin=c.min(),vmax=c.max())
    cmap = mpl.cm.ScalarMappable(norm=norm,cmap=mpl.cm.rainbow)
    cmap.set_array([])
    fig,ax= plt.subplots()
    
    for i in range(M):
        ax.plot(w_range,volumes[:,i],linewidth=2,c=cmap.to_rgba(c[i]))
    plt.xlabel("Strength of Selection w",**label_pars)
    plt.ylabel("Relative Size of Basin of Attraction",**label_pars)
    plt.title(f"Number of types: {M}",**title_pars)
    plt.tick_params(**tick_pars)
    cbar = fig.colorbar(cmap)
    cbar.ax.set_ylabel('Environmental Fitness Benefit',**cb_label_pars)
    cbar.ax.tick_params(**cb_tick_pars)
    plt.tight_layout()

    try:
        os.mkdir("Figures")
    except:
        pass
    
    plt.savefig("Figures/figure_S1.pdf")
    plt.show()
    
    possible_types = np.sum(volumes[-1,:]>0.01)
    print(f"Number of types alive: {possible_types}/{M}")