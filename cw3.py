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

#%% 5 simulations
n_sim=5
tspan=[0,1000]
y0=0
t = [0]*n_sim
y = [0]*n_sim
fig, ax = plt.subplots()

for i in range(n_sim):
    (t[i],y[i]) = gillespie(lambda y: propensities(y, k0, k1), transition_rules, tspan, y0,method='direct')
    print("Mean = {:.2f}, variance = {:.2f}".format(np.mean(y[i]), np.var(y[i])))

    ax.plot(t[i],y[i])
    
ax.set_xlabel('t/s')
ax.set_ylabel('mRNA number')
ax.set_title('Sample traces of mRNA production')
plt.show()


#%% Fano factor convergence, n_sims -> inf

tspan = [0, 1000]
y0 = 20
n_tot = 30000
y = np.zeros(n_tot)
for i in range(n_tot):
    (thold,yhold) = gillespie(lambda y: propensities(y, k0, k1), transition_rules, tspan, y0,method='direct')
    y[i] = yhold[-1]
    if (100*i/n_tot) % 1 == 0:
        print("\r{:.0f}%".format(100*(i+1)/n_tot),end='\r')
print("Mean = {:.2f}, variance = {:.2f}, Fano factor = {:.2f}".format(np.mean(y), np.var(y), np.var(y)/np.mean(y)))


#%% Distribution plots

fig, ax = plt.subplots()
(ints, freq) = integer_histogram(y)
ax.bar(ints,freq)
ax.set_xlabel('mRNA number')
ax.set_ylabel('Relative abundance')
ax.set_title('mRNA distribution without cell cycle')
#Poisson overlay
P=20
poisson = np.zeros(len(ints))
for i, integer in enumerate(ints):
    poisson[i] = np.exp(-P)*P**integer/math.factorial(integer)
ax.plot(np.array(ints),poisson,c='r')
ax.set_title('mRNA distribution without cell cycle\nOverlaid with Poisson(20)')


#%% Statistic convergence plots

n_points = n_tot
n_sims = np.zeros(n_points)
fano = np.zeros(n_points)
mean = np.zeros(n_points)
variance = np.zeros(n_points)
simspan = [1,n_tot]
for i in range(n_points):
    idx = int((i+1)*(len(y)-1)/(n_points+1))
    n_sims[i] = idx
    fano[i] = np.var(y[0:idx])/np.mean(y[0:idx])
    mean[i] = np.mean(y[0:idx])
    variance[i] = np.var(y[0:idx])
fig,[[ax1,ax2],[ax3,ax4],[ax5,ax6]]=plt.subplots(nrows=3,ncols=2,figsize=(7,5.5))
fig.suptitle(r'Statistic convergence as simulation count increases' + '\nLog and linear scales')
ax1.plot(simspan,[20,20],c='r')
ax1.plot(n_sims,mean)
ax1.set_xscale('log')
ax1.set_xlim(simspan)
ax1.set_ylabel(r'Mean')

ax2.plot(simspan,[20,20],c='r')
ax2.plot(n_sims,mean)
ax2.set_xscale('linear')
ax2.set_xlim(simspan)

ax3.plot(simspan,[20,20],c='r')
ax3.plot(n_sims,variance)
ax3.set_xscale('log')
ax3.set_xlim(simspan)
ax3.set_ylabel('Variance')

ax4.plot(simspan,[20,20],c='r')
ax4.plot(n_sims,variance)
ax4.set_xscale('linear')
ax4.set_xlim(simspan)

ax5.plot(simspan,[1,1],c='r')
ax5.plot(n_sims,fano)
ax5.set_xscale('log')
ax5.set_xlim(simspan)
ax5.set_xlabel('Simulation count')
ax5.set_ylabel('Fano factor')

ax6.plot(simspan,[1,1],c='r')
ax6.plot(n_sims,fano)
ax6.set_xscale('linear')
ax6.set_xlim(simspan)
ax6.set_xlabel('Simulation count')
    

#%% Fano factor convergence, simulation time -> infinity, using time-normalised statistics

tspan = [0, 1]
y0 = 20
(t,y) = gillespie(lambda y: propensities(y, k0, k1), transition_rules, tspan, y0,method='direct')
print("Mean = {:.2f}, variance = {:.2f}, Fano factor = {:.2f}".format(time_normalised_mean(t,y), time_normalised_var(t,y), time_normalised_var(t,y)/time_normalised_mean(t,y)))


n_points = 1000
t_points = np.zeros(n_points)
fano = np.zeros(n_points)
mean = np.zeros(n_points)
variance = np.zeros(n_points)
for i in range(n_points):
    idx = int((i+1)*(len(y)-1)/(n_points+1))
    t_points[i] = t[idx]
    mean[i] = time_normalised_mean(t[0:idx],y[0:idx])
    variance[i] = time_normalised_var(t[0:idx],y[0:idx])
    fano[i] = variance[i]/mean[i]

fig,[[ax1,ax2],[ax3,ax4],[ax5,ax6]]=plt.subplots(nrows=3,ncols=2,figsize=(7,5.5))
fig.suptitle(r'Statistic convergence as simulation length increases' + '\nLog and linear scales')
ax1.plot(tspan,[20,20],c='r')
ax1.plot(t_points,mean)
ax1.set_xscale('log')
ax1.set_xlim(tspan)
ax1.set_ylabel(r'Mean')

ax2.plot(tspan,[20,20],c='r')
ax2.plot(t_points,mean)
ax2.set_xscale('linear')
ax2.set_xlim(tspan)

ax3.plot(tspan,[20,20],c='r')
ax3.plot(t_points,variance)
ax3.set_xscale('log')
ax3.set_xlim(tspan)
ax3.set_ylabel('Variance')

ax4.plot(tspan,[20,20],c='r')
ax4.plot(t_points,variance)
ax4.set_xscale('linear')
ax4.set_xlim(tspan)

ax5.plot(tspan,[1,1],c='r')
ax5.plot(t_points,fano)
ax5.set_xscale('log')
ax5.set_xlim(tspan)
ax5.set_xlabel('Simulation time')
ax5.set_ylabel('Fano factor')

ax6.plot(tspan,[1,1],c='r')
ax6.plot(t_points,fano)
ax6.set_xscale('linear')
ax6.set_xlim(tspan)
ax6.set_xlabel('Simulation time')

#%% Incorrect convergence illustration 

n_points = 1000
t_points = np.zeros(n_points)
fano = np.zeros(n_points)
mean = np.zeros(n_points)
variance = np.zeros(n_points)
for i in range(n_points):
    idx = int((i+1)*(len(y)-1)/(n_points+1))
    t_points[i] = t[idx]
    mean[i] = np.mean(y[0:idx])
    variance[i] = np.var(y[0:idx])
    fano[i] = variance[i]/mean[i]

fig,[[ax1,ax2],[ax3,ax4],[ax5,ax6]]=plt.subplots(nrows=3,ncols=2,figsize=(7,5.5))
fig.suptitle(r'Incorrect statistic convergence as simulation length increases' + '\nLog and linear scales')
ax1.plot(tspan,[20,20],c='r')
ax1.plot(t_points,mean)
ax1.set_xscale('log')
ax1.set_xlim(tspan)
ax1.set_ylabel(r'Mean')

ax2.plot(tspan,[20,20],c='r')
ax2.plot(t_points,mean)
ax2.set_xscale('linear')
ax2.set_xlim(tspan)

ax3.plot(tspan,[20,20],c='r')
ax3.plot(t_points,variance)
ax3.set_xscale('log')
ax3.set_xlim(tspan)
ax3.set_ylabel('Variance')

ax4.plot(tspan,[20,20],c='r')
ax4.plot(t_points,variance)
ax4.set_xscale('linear')
ax4.set_xlim(tspan)

ax5.plot(tspan,[1,1],c='r')
ax5.plot(t_points,fano)
ax5.set_xscale('log')
ax5.set_xlim(tspan)
ax5.set_xlabel('Simulation time')
ax5.set_ylabel('Fano factor')

ax6.plot(tspan,[1,1],c='r')
ax6.plot(t_points,fano)
ax6.set_xscale('linear')
ax6.set_xlim(tspan)
ax6.set_xlabel('Simulation time')

#%% Cell division with binomial partitioning (binned mean)

y0 = 0

fig2, ax2 = plt.subplots()
cmap = plt.get_cmap("tab10")
division_time = 20*60
n_generations = 6
n_runs = 1000


t2 = [0]*n_runs
y2 = [0]*n_runs
for i in range(n_runs):
    (t2[i],y2[i]) = cell_partition(lambda y: propensities(y, k0, k1), transition_rules, division_time, n_generations, y0, segregation='independent')
    if n_runs < 10:
        ax2.plot(t2[i],y2[i],zorder=-1, linewidth = 0.2,c='gold')
        div_idx = np.where(t2[i] % division_time == 0)
        #Plot scatter of points after division
        ax2.scatter(t2[i][div_idx], y2[i][div_idx],marker='o',facecolors='none',edgecolors = cmap(i),zorder=2)
        #Plot scatter of points before division
        ax2.scatter(t2[i][np.subtract(div_idx,1)], y2[i][np.subtract(div_idx,1)],marker='o',facecolors='none',edgecolors=cmap(i),zorder=2)
    if (100*i/n_runs) % 1 == 0:
        print("\r{:.0f}%".format(100*(i+1)/n_runs),end='\r')
(t_avg, y_avg) = binned_mean(t2, y2, [0, division_time*n_generations], 1)
ax2.plot(t_avg, y_avg, zorder=1, c = 'tab:cyan')
ax2.set_xlabel("t/s")
ax2.set_ylabel("Average mRNA Count")
ax2.set_title("Average mRNA with cell cycle, 1000 runs, 1s bin width")

plt.show()

#%% Same as above, but with closest mean

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

#%% Approximate Bayesian Computation

k0 = 0.2
k1 = 0.01

y0 = 0
division_time = 20*60
n_generations = 10
tspan = [0, n_generations*division_time]
#Data holders
params = [0]*2
similarity = [0]*2
#Parameters
param_ranges = np.asarray([[0, 0.2], [0, 0.02]])
d_stats = [10,10]
n_samples = 2000
epsilon = 0.5

#data_stats is given as [mean, var]
(params[0], similarity[0]) = approximate_bayesian(lambda params: gillespie(lambda y: propensities(y, params[0], params[1]), transition_rules, tspan, y0), param_ranges, d_stats, n_samples,epsilon)
(params[1], similarity[1]) = approximate_bayesian(lambda params: cell_partition(lambda y: propensities(y, params[0],params[1]), transition_rules, division_time, n_generations, y0, segregation='independent'), param_ranges, d_stats, n_samples,epsilon)

#%% Similarity histograms
fig,ax = plt.subplots()
ax.hist(similarity[0],bins=20,alpha=0.7,edgecolor=None,facecolor='b',label='Original',density=True)
ax.hist(similarity[1],bins=20,alpha=0.7,edgecolor=None,facecolor='r',label='Cell Division',density=True)
ax.set_xlabel('Similarity Distance')
ax.set_ylabel('Relative Frequency')
ax.set_title('Overlay of similarity measure histograms\n'+ str(n_samples) + ' simulations')
ax.legend()
plt.show()

#%% Posterior histograms

fig3, [[ax3,ax4],[ax5,ax6]] =plt.subplots(2,2)
fig3.suptitle('Posterior distributions for original (top) and cell cycle (bottom) models\n' + r'{} simulations, $\epsilon={}$'.format(n_samples, epsilon))

ax3.hist(params[0][:,0],bins=20,density=True,facecolor='b')
ax3.set_ylabel("Relative Frequency")

ax4.hist(params[0][:,1],bins=20,density=True,facecolor='b')

ax5.hist(params[1][:,0],bins=20,density=True,facecolor='b')
ax5.set_xlabel("$k_0$")
ax5.set_ylabel("Relative Frequency")

ax6.hist(params[1][:,1],bins=20,density=True,facecolor='b')
ax6.set_xlabel("$k_1$")

plt.show()

#%% Multimodel Approximate Bayesian Computation

k0 = 0.2
k1 = 0.01
y0 = 0
division_time = 20*60
n_generations = 6
tspan = [0, n_generations*division_time]
#Data holders
params = [0]*2
similarity = [0]*2
#Parameters
param_ranges = np.asarray([[0, 0.2], [0, 0.02]])
d_stats = [10,10]
n_samples = 1000
epsilon = 0.5

model_set = [lambda params: gillespie(lambda y: propensities(y, params[0], params[1]), transition_rules, tspan, y0),\
             lambda params: cell_partition(lambda y: propensities(y, params[0],params[1]), transition_rules, division_time, n_generations, y0, segregation='independent')]

#data_stats is given as [mean, var]
(multi_params, multi_similarity, multi_model) = approximate_bayesian(model_set, param_ranges, d_stats, n_samples,epsilon)
#%%

fig4, (ax1,ax2,ax3) = plt.subplots(nrows=1,ncols=3)
fig4.set_size_inches(7,4)
fig4.suptitle("Posterior distributions and relative model selection frequency\nfor multiple-model ABC")
ax1.hist(multi_params[multi_model==0,0],bins=20,facecolor='b',alpha=0.7,density=True)
ax1.hist(multi_params[multi_model==1,0],bins=20,facecolor='r',alpha=0.7,density=True)
ax1.set_xlabel("$k_0$")
ax1.set_ylabel("Relative Frequency")

ax2.hist(multi_params[multi_model==0,1],bins=20,facecolor='b',alpha=0.7,density=True)
ax2.hist(multi_params[multi_model==1,1],bins=20,facecolor='r',alpha=0.7,density=True)
ax2.set_xlabel("$k_1$")

ax3.bar([1,2],[sum(multi_model==0)/len(multi_model),sum(multi_model==1)/len(multi_model)],color=['b','r'],alpha=0.7)
ax3.set_xticks([1,2])
ax3.set_xticklabels(['Original','Cell\nDivision'])
plt.show()
