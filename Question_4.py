import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

import matplotlib.pyplot as plt
import math
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

from gillespie import *

def propensities(y, k0, k1):
    return [k0, k1*y]


k0 = 0.2
k1 = 0.01
tspan = [0, 1000]
y0 = 0
transition_rules = [1, -1]


#%% Sample traces for cell division model

y0 = 0

fig2, ax2 = plt.subplots(figsize=(10,5))
cmap = plt.get_cmap("tab10")
division_time = 20*60
n_generations = 3
n_runs = 5


t2 = [0]*n_runs
y2 = [0]*n_runs
for i in range(n_runs):
    (t2[i],y2[i]) = cell_partition(lambda y: propensities(y, k0, k1), transition_rules, division_time, n_generations, y0, segregation='independent')
    if n_runs < 10:
        ax2.plot(t2[i],y2[i],zorder=-1, linewidth = 0.2,c=cmap(i))
        div_idx = np.where(t2[i] % division_time == 0)
        #Plot scatter of points after division
        ax2.scatter(t2[i][div_idx], y2[i][div_idx],marker='o',facecolors='none',edgecolors = cmap(i),zorder=2)
        #Plot scatter of points before division
        ax2.scatter(t2[i][np.subtract(div_idx,1)], y2[i][np.subtract(div_idx,1)],marker='o',facecolors='none',edgecolors=cmap(i),zorder=2)
    if (100*i/n_runs) % 1 == 0:
        print("\r{:.0f}%".format(100*(i+1)/n_runs),end='\r')
(t_avg, y_avg) = closest_mean(t2, y2, [0, division_time*n_generations], 1)
ax2.plot(t_avg, y_avg, zorder=1, c = 'tab:cyan')
ax2.set_xlabel("t/s")
ax2.set_ylabel("mRNA Count")
ax2.set_title("mRNA sample traces, average in cyan\n"+ 'Partitioning between sets of same-colour markers')

plt.show()


#%% Mean and Fano factor, with cell division

division_time = 20*60
n_generations = 1000
y0 = 10
(t,y) = cell_partition(lambda y: propensities(y, k0, k1), transition_rules, division_time, n_generations, y0, segregation='independent')

print("Mean = {:.2f}, variance = {:.2f}, Fano factor = {:.2f}".format(time_normalised_mean(t,y), time_normalised_var(t,y), time_normalised_var(t,y)/time_normalised_mean(t,y)))

generations = range(1,n_generations+1)
fano = np.zeros(n_generations)
mean = np.zeros(n_generations)
variance = np.zeros(n_generations)

partition_vals = []
partition_mean = np.zeros(n_generations)
partition_variance = np.zeros(n_generations)
partition_fano = np.zeros(n_generations)
for i in range(n_generations):
    idx = int((i+1)*(len(y)-1)/(n_generations+1))
    mean[i] = time_normalised_mean(t[0:idx],y[0:idx])
    variance[i] = time_normalised_var(t[0:idx],y[0:idx])
    fano[i] = variance[i]/mean[i]
    
    partition_idx = t == division_time*i
    partition_vals.append(y[partition_idx])
    partition_mean[i] = np.mean(partition_vals)
    partition_variance[i] = np.var(partition_vals)
    partition_fano[i] = partition_variance[i]/partition_mean[i]

fig,[[ax1,ax2],[ax3,ax4],[ax5,ax6]]=plt.subplots(nrows=3,ncols=2,figsize=(7,5.5))
fig.suptitle(r'Statistic convergence as number of cell generations increases' + '\nLog and linear scales')
ax1.plot([1,n_generations],[mean[-1],mean[-1]],c='r')
ax1.plot(generations,mean)
ax1.set_xscale('log')
ax1.set_xlim([1,n_generations])
ax1.set_ylabel(r'Mean')

ax2.plot([1,n_generations],[mean[-1],mean[-1]],c='r')
ax2.plot(generations,mean)
ax2.set_xscale('linear')
ax2.set_xlim([1,n_generations])

ax3.plot([1,n_generations],[variance[-1],variance[-1]],c='r')
ax3.plot(generations,variance)
ax3.set_xscale('log')
ax3.set_xlim([1,n_generations])
ax3.set_ylabel('Variance')

ax4.plot([1,n_generations],[variance[-1],variance[-1]],c='r')
ax4.plot(generations,variance)
ax4.set_xscale('linear')
ax4.set_xlim([1,n_generations])

ax5.plot([1,n_generations],[fano[-1],fano[-1]],c='r')
ax5.plot(generations,fano)
ax5.set_xscale('log')
ax5.set_xlim([1,n_generations])
ax5.set_xlabel('Generations')
ax5.set_ylabel('Fano factor')

ax6.plot([1,n_generations],[fano[-1],fano[-1]],c='r')
ax6.plot(generations,fano)
ax6.set_xscale('linear')
ax6.set_xlim([1,n_generations])
ax6.set_xlabel('Generations')

fig2,[[ax7,ax8],[ax9,ax10],[ax11,ax12]]=plt.subplots(nrows=3,ncols=2,figsize=(7,5.5))
fig2.suptitle(r'Partition statistic convergence as number of cell generations increases' + '\nLog and linear scales')
ax7.plot([1,n_generations],[10,10],c='r')
ax7.plot(generations,partition_mean)
ax7.set_xscale('log')
ax7.set_xlim([1,n_generations])
ax7.set_ylabel(r'Partition mean')

ax8.plot([1,n_generations],[10,10],c='r')
ax8.plot(generations,partition_mean)
ax8.set_xscale('linear')
ax8.set_xlim([1,n_generations])

ax9.plot([1,n_generations],[10,10],c='r')
ax9.plot(generations,partition_variance)
ax9.set_xscale('log')
ax9.set_xlim([1,n_generations])
ax9.set_ylabel('Partition variance')

ax10.plot([1,n_generations],[10,10],c='r')
ax10.plot(generations,partition_variance)
ax10.set_xscale('linear')
ax10.set_xlim([1,n_generations])

ax11.plot([1,n_generations],[1,1],c='r')
ax11.plot(generations,partition_fano)
ax11.set_xscale('log')
ax11.set_xlim([1,n_generations])
ax11.set_xlabel('Generations')
ax11.set_ylabel('Partition Fano factor')

ax12.plot([1,n_generations],[1,1],c='r')
ax12.plot(generations,partition_fano)
ax12.set_xscale('linear')
ax12.set_xlim([1,n_generations])
ax12.set_xlabel('Generations')