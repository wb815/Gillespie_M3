import numpy as np
from collections import deque


def gillespie(fun, transition_rules, tspan, y0, method='first'):
    """
    Function for performing Gillespie simulations for arbitrary systems, with user-defined
    propensities and transition rules. Transition rules can correspond to reaction 
    stoichiometries, or reaction associations, in the case where species are used as templates.
    Parameters:
        fun : callable
            Propensity functions for the system.
            The calling signature is fun(y). Here y has shape (n,).
            Fun must return array_like with shape (p,).
            
        transition-rules : array_like, shape (p,n)
        Rules for changes in y upon each transition. Each row corresponds to the updates to all y for a specific
        transition. Ensure this formatting, especially if number of propensities = number of species
    
        t_span : 2-tuple of floats
        Interval of simulation (t0, tf). The solver starts with t=t0 and integrates until it reaches t=tf.
        
        y0 : array_like, shape (n,)
        Initial state.
        
        method : string: 'first' or 'direct'
        The algorithm to use for the simulation. They are statistically identical and possess the same 
        asymptotic complexity, but 'direct' has proven more efficient for small systems (<25 reactions),
        and 'first' for larger systems. The Next Reaction Method (Gibson and Bruck, 2000), is even more 
        efficient, but not implemented.
        
    
    Returns:
        t : ndarray, shape (n_points,)
        Time points.
        
        y : ndarray, shape (n_points,n)
        Values of the solution at t.
    """
    #Ignore divide by 0 errors (inevitable for first reaction method (ie sampling with propensity of 0)
    np.seterr(divide='ignore')
    #Initialise deques (used for speedy append)
    t = deque([tspan[0]])
    y = deque([y0])
    #Not really needed due to break clause, but using just because I can
    break_loop = False
    while not break_loop:
        #O(R)
        propensities = fun(y[-1])
        if method=='first':
            #O(R)
            # 1/alpha because of how Python defines exponential distributions
            taus = np.random.exponential(np.divide(1,propensities))
            #Get index of next reaction
            #O(R)
            mu = np.argmin(taus)
            #Calculate next  values
            t_next = t[-1] + taus[mu]
            y_next = np.add(y[-1], transition_rules[mu])
        elif method=='direct':
            #Sample the reaction "selector" from a distribution ~U([0,1))
            mu = np.random.rand()
            #O(R)
            S = sum(propensities)
            #Recursively subtract the relative propensities to identify where the selector landed
            #O(R)
            for i, p in enumerate(propensities):
                if mu < p/S:
                    #When selector has landed, set it to the index mu of the selected reaction
                    mu = i
                    break
                else:
                    mu -= p/S
            #Select the reaction time from an exponential distribution parametrised by the sum of the propensities
            # 1/S because of how Python defines exponential distributions
            tau = np.random.exponential(1/S)
            t_next = t[-1] + tau
            y_next = np.add(y[-1], transition_rules[mu])
        else:
            raise ValueError('Invalid method requested.')
        #Check for events (including timeout)
        if t_next > tspan[1]:
            break_loop = True
            #Break to avoid overrunning tspan
            break
        #Update time
        t.append(t_next)
        #Update abundance of reaction products
        y.append(y_next)
        
    t = np.asarray(t)
    y = np.asarray(y)
    
    return (t,y)#GillespieResult(t=t,y=y)

def cell_partition(fun, transition_rules, division_time, n_generations, y0, method='direct', segregation='independent'):
    """
    Function for performing Gillespie simulations for arbitrary cellular systems, including cell division, 
    with user-defined propensities and transition rules. Transition rules can correspond to reaction 
    stoichiometries, or reaction associations, in the case where species are used as templates.
    Currently, only independent segregation of species upon cell division is implemented.
    Parameters:
        fun : callable
            Propensity functions for the system.
            The calling signature is fun(y). Here y has shape (n,).
            Fun must return array_like with shape (p,).
            
        transition-rules : array_like, shape (p,n)
        Rules for changes in y upon each transition. Each row corresponds to the updates to all y for a specific
        transition. Ensure this formatting, especially if number of propensities = number of species
    
        t_span : 2-tuple of floats
        Interval of simulation (t0, tf). The solver starts with t=t0 and integrates until it reaches t=tf.
        
        y0 : array_like, shape (n,)
        Initial state.
        
        method : string: 'first' or 'direct'
        The algorithm to use for the simulation. They are statistically identical and possess the same 
        asymptotic complexity, but 'direct' has proven more efficient for small systems (<25 reactions),
        and 'first' for larger systems. The Next Reaction Method (Gibson and Bruck, 2000), is even more 
        efficient, but not implemented.
        
        segregation : string: 'independent'
        Placeholder parameter: currently, only independent partitioning is implemented, but this parameter
        is included to allow for expansion to ordered or disordered partitioning in the future.
        
    
    Returns:
        t : ndarray, shape (n_points,)
        Time points.
        
        y : ndarray, shape (n_points,n)
        Values of the solution at t.
    """
    for n in range(n_generations):
        tspan = [n*division_time, (n+1)*division_time]
        if n == 0:
            (t,y) = gillespie(fun, transition_rules, tspan, y0, method)
        else:
            if segregation == 'independent':
                y0 = np.random.binomial(y[-1],0.5)
            (th,yh) = gillespie(fun, transition_rules, tspan, y0, method)
            t = np.append(t, th, axis=0)
            y = np.append(y, yh, axis=0)
    return (t,y)

def integer_histogram(y):
    """
    Function returning a histogram for integer data.
    Parameters:
        y : array_like, shape (l,n)
        Input data
    
    Returns:
        ints : ndarray, shape (i,)
        Array of integer values
        
        freq : ndarray, shape (i,n)
        Relative frequency of ints
    """
    ints = range(int(np.min(y)),int(np.max(y)))
    freq = np.zeros(len(ints))
    y = list(y)
    for i,num in enumerate(ints):
        freq[i] = y.count(num)/len(y)
    return (ints, freq)
        
    

def binned_mean(t, y, tspan, resolution, include_start=True):
    """
    Function to calculate the binned mean of a number of time traces. 
    Parameters:
        t : list of array_like, length (s)
        List of input time traces
        
        y : list of array_like, length (s)
        List of input data
        
        tspan : array_like, shape (2,)
        Time span of the data to be averaged
        
        include_start : bool
        Whether or not to include the mean of the initial datapoint as the 
        start of the returned mean
    
    Returns:
        t_out : ndarray, shape (i,)
        Array of time points, one at the centre of each bin
        
        y_out : ndarray, shape (i,n)
        Array of mean y values within each bin, across all input traces
    """
    bin_centres = np.arange(tspan[0]+resolution/2, tspan[1], resolution)
    bin_divs =  np.arange(tspan[0], tspan[1]+resolution/2, resolution)
    if np.ndim(y[0]) == 1:
        bin_vals = np.zeros(len(bin_centres))
        lin = True
    else:
        bin_vals = np.zeros([len(bin_centres),np.shape(y[0])[1]])
        lin = False
    for i in range(len(bin_centres)):
         binned_idx = np.logical_and(t[0]>= bin_divs[i], t[0]<bin_divs[i+1])
         y_pool = y[0][binned_idx]
         for j in range(1,len(t)):
            binned_idx = np.logical_and(t[j]>= bin_divs[i], t[j]<bin_divs[i+1])
            if lin:
                y_pool = np.append(y_pool, y[j][binned_idx])
            else:
                y_pool = np.append(y_pool, y[j][binned_idx], axis=0)
         bin_vals[i] = np.mean(y_pool)
    if include_start:
        t_out = np.append(tspan[0], bin_centres)
        y0 = np.asarray(y[0][0])
        for i in range(1,len(y)):
            y0 = np.append(y0, np.asarray(y[i][0]))
        y0 = y0.mean(axis=0)
        y_out = np.append(y0, bin_vals)
    else:
        t_out = bin_centres
        y_out = bin_vals
    return (t_out,y_out)

def closest_mean(t,y,tspan,resolution):
    """
    Function to calculate the binned mean of a number of time traces. 
    Parameters:
        t : list of array_like, length (s)
        List of input time traces
        
        y : list of array_like, length (s)
        List of input data
        
        tspan : array_like, shape (2,)
        Time span of the data to be averaged
        
        resolution : float
        Sampling time period, will return one mean value per every "resolution"
        time units
    
    Returns:
        t_out : ndarray, shape (i,)
        Array of time points, one every "resolution" time units
        
        y_out : ndarray, shape (i,n)
        Array of mean y values across all input traces, with each value to be averaged
        taken as the closest in time to the respective time point in t_out
    """
    t_out = np.arange(tspan[0],tspan[1]+resolution,resolution)
    if np.ndim(y[0]) == 1:
        y_out = np.zeros(len(t_out))
        lin = True
    else:
        y_out = np.zeros([len(t_out),np.shape(y[0])[1]])
        lin = False
    for i,t_curr in enumerate(t_out):
        closest_idx = np.abs(t[0] - t_curr).argmin()
        y_pool = y[0][closest_idx]

        for j in range(1,len(t)):
            closest_idx = np.abs(t[j] - t_curr).argmin()

            if lin:
                y_pool = np.append(y_pool, y[j][closest_idx])
            else:
                y_pool = np.append(y_pool, y[j][closest_idx],axis=0)
        y_out[i] = np.mean(y_pool)
    return (t_out, y_out)

def time_normalised_mean(t,y):
    """
    Function to calculate the time-normalised mean of time series data
    Parameters:
        t : array_like, shape (n_points,)
        Time points of the input time series
        
        y : array_like, shape (n_points,n)
        Values of the input time series
    
    Returns:
        mu : 
        if n == 1:
            float
        else:
            array_like, shape (n,)
        Mean of y values
    """
    T = t[-1]-t[0]
    y = np.asarray(y)
    mu = np.sum(np.multiply(y[:-1],np.divide(t[1:]-t[:-1],T)),axis=0)
    return mu

def time_normalised_var(t,y):
    """
    Function to calculate the time-normalised variance of time series data
    Parameters:
        t : array_like, shape (n_points,)
        Time points of the input time series
        
        y : array_like, shape (n_points,n)
        Values of the input time series
    
    Returns:
        var : 
        if n == 1:
            float
        else:
            array_like, shape (n,)
        Variance of y values
    """
    T = t[-1]-t[0]
    y = np.asarray(y)
    n = len(t)-1
    mu = time_normalised_mean(t,y)
    var = np.divide(np.sum(np.multiply(np.power(np.subtract(y[:-1],mu),2),t[1:]-t[:-1]),axis=0),T-1)
    return var

def approximate_bayesian(model_set, param_ranges, d_stats, n_samples, epsilon):
    """
    Function to perform approximate Bayesian computation on candidate models. 
    Currently hard-coded to only accept experimental data in the form of mean and variance,
    and to utilise the similarity distance function in the coursework.
    Parameters:
        model_set : callable, or list of callables
        Model or set of models upon which to perform ABC. If list, multimodel ABC
        is performed
        
        param_ranges : array_like, shape (k,2)
        Prior ranges for parameters 
        
        d_stats : array_like, shape (2,)
        Mean (index 0) and variance (index 1) of experimental data. Hard-coded to only 
        accept these summary statistics at this time.
        
        n_samples : int
        Number of acceptable parameter sets to return
        
        epsilon : float
        Similarity distance threshold
    
    Returns:
        params : ndarray, shape (n_samples, k)
        Accepted parameter sets
        
        similarity : ndarray, shape (n_samples,)
        Similarity distance values for accepted parameter sets
        
        models : ndarray, shape (n_samples,)
        Model index for accepted parameter sets. Only returned if multimodel ABC
        is used
    """
    if isinstance(model_set,list):
        # Treat model set as priors
        params = np.ndarray([n_samples,param_ranges.shape[0]])
        similarity = np.ndarray(n_samples)
        models = np.ndarray(n_samples)
        for i in range(n_samples):
            satisfactory = False
            test = 0
            while satisfactory == False:
                candidate_params = np.random.uniform(param_ranges[:,0], param_ranges[:,1])
                candidate_model = np.random.randint(len(model_set))
                (t,y) = model_set[candidate_model](candidate_params)
                #Hard-coded similarity measures, maybe open up to general if we can be bothered
                candidate_similarity = (time_normalised_mean(t,y) - d_stats[0])**2 + (time_normalised_var(t,y) - d_stats[1])**2
                test += 1
                if candidate_similarity <= epsilon:
                    satisfactory = True
                    params[i,:] = candidate_params
                    models[i] = candidate_model
                    similarity[i] = candidate_similarity
                if test > 1000:
                    if i == 0:
                        print("No satisfactory condition found. Try increasing epsilon")
                    else:
                        print("Few satisfactory conditions found. Try increasing epsilon")
                    params = params[0:i,:]
                    models = models[0:i]
                    similarity = similarity[0:i]
                    return (params, similarity)
            if (100*i/n_samples) % 1 == 0:
                print("\r{:.0f}%".format(100*(i+1)/n_samples),end='\r')
        print("\r100%")
        return (params, similarity, models)
    else:
        # Single-model ABC
        params = np.ndarray([n_samples,param_ranges.shape[0]])
        similarity = np.ndarray(n_samples)
        model = model_set
        for i in range(n_samples):
            satisfactory = False
            test = 0
            while satisfactory == False:
                candidate_params = np.random.uniform(param_ranges[:,0], param_ranges[:,1])
                (t,y) = model(candidate_params)
                #Hard-coded similarity measures, maybe open up to general if we can be bothered
                candidate_similarity = (time_normalised_mean(t,y) - d_stats[0])**2 + (time_normalised_var(t,y) - d_stats[1])**2
                test += 1
                if candidate_similarity <= epsilon:
                    satisfactory = True
                    params[i,:] = candidate_params
                    similarity[i] = candidate_similarity
                if test > 1000:
                    if i == 0:
                        print("No satisfactory condition found. Try increasing epsilon")
                    else:
                        print("Few satisfactory conditions found. Try increasing epsilon")
                    params = params[0:i,:]
                    similarity = similarity[0:i]
                    return (params, similarity)
            if (100*i/n_samples) % 1 == 0:
                print("\r{:.0f}%".format(100*(i+1)/n_samples),end='\r')
        print("\r100%")
        return (params, similarity)