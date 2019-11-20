#Benchmarking scripts for CW3
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

import matplotlib.pyplot as plt
import time
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from gillespie import *



def bench(y, params):
    if isinstance(params,list):
        out = [params[0]]
        out.append(list(y*np.asarray(params[1:])))
    else:
        out = params
    return out

scale = np.round(10**np.arange(0,3.3,0.2),0)

tspan = [0, 1000]
y0 = 0
repeats = 4
t_direct = []
t_first = []
for i, s in enumerate(scale):
    t_direct.append([])
    t_first.append([])
    print(s)
    for r in range(repeats):
        transition_rules = [100]
        params = np.random.rand(int(s))
        transition_rules += [-1]*(len(params)-1)
        tic = time.clock()
        gillespie(lambda y: bench(y, params), transition_rules, tspan, y0,method='direct')
        toc = time.clock()
        t_direct[i].append(toc-tic)
        tic = time.clock()
        gillespie(lambda y: bench(y, params), transition_rules, tspan, y0,method='first')
        toc = time.clock()
        t_first[i].append(toc-tic)
avg_first = [np.mean(test) for test in t_first]
avg_direct = [np.mean(test) for test in t_direct]
var_first = [np.var(test) for test in t_first]
var_direct = [np.var(test) for test in t_first]
    
#%% Plotting
    
fig, ax = plt.subplots()
ax.errorbar(scale,avg_first, yerr=var_first, label='First Reaction Method')
ax.errorbar(scale,avg_direct, yerr=var_direct, label='Direct Method')
ax.set_xscale('log')
ax.set_yscale('log')
ax.grid()
ax.set_xlabel('Model complexity / number of reactions')
ax.set_ylabel('Execution time / s')
ax.set_title('Computational scaling of Gillespie algorithms\n' + 'Error bars represent variance')
ax.legend()