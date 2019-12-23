# Gillespie_M3

Python3 repository for the MRes in Systems and Synthetic Biology Coursework 3 (Theoretical Systems Biology) at Imperial College London.

CONTENTS
----------

benchmarking.py
----------
  Basic computational scaling testing for First Reaction and Direct implementations of the Gillspie algorithm


cw3.py
----------
  Working file for this coursework, containing all the code required to execute all functions and draw all graphs in this coursework (with a few exceptions)


gillespie.py
----------
  Gillespie simulation module, see below for contents


probability.py
----------
  Python file containing code required to calculate the probability distribution, mean and variance of the distribution derived in this coursework, which is described by the probability mass function
  
<img src="https://latex.codecogs.com/svg.latex?\Large&space;f(x;\lambda,p)=\sum\limits_{r=x}^{\infty}\frac{e^{-\lambda}\lambda^r}{x!(r-x)!}p^x(1-p)^{r-x}" title="\Large f(x;\lambda,p)=\sum\limits_{r=x}^{\infty}\frac{e^{-\lambda}\lambda^r}{x!(r-x)!}p^x(1-p)^{r-x}" />

Question_4.py
----------
  A subsection of cw3.py containing code required to answer Question 4 from this coursework. Copied in a seperate file for ease of marking
  
--------------------------------------------------

gillespie.py
----------

gillespie
----------

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
        
----------

cell_partition
----------

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
        
----------

integer_histogram
----------

    Function returning a histogram for integer data.
    
    Parameters:
        y : array_like, shape (l,n)
        Input data
    
    Returns:
        ints : ndarray, shape (i,)
        Array of integer values
        
        freq : ndarray, shape (i,n)
        Relative frequency of ints
        
----------

binned_mean
----------

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
        
----------

closest_mean
----------

    Function to calculate the mean of a number of time traces based on the nearest value to 
    desired sample points. 
    
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
        
----------

time_normalised_mean
----------


The formula for time-normalised mean is given by

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\mu(Y)=\frac{1}{\sum\limits_{j=1}^{N}(dt_j)}\sum\limits_{i=1}^{N}(y_idt_i)" title="\Large \mu(Y)=\frac{1}{\sum\limits_{j=1}^{N}(dt_j)}\sum\limits_{i=1}^{N}(y_idt_i)" />

where Y = {y<sub>1</sub>,..,y<sub>N</sub>} and dt<sub>i</sub> is the time spent in state y<sub>i</sub>.

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
        
----------

time_normalised_var
----------

The formula for time-normalised variance is given by

<img src="https://latex.codecogs.com/svg.latex?\Large&space;\sigma^2(Y)=\frac{N}{(N-1)\sum\limits_{j=1}^{N}(dt_j)}\sum\limits_{i=1}^{N}((y_i-\mu)^2dt_i)" title="\Large \sigma^2(Y)=\frac{N}{(N-1)\sum\limits_{j=1}^{N}(dt_j)}\sum\limits_{i=1}^{N}((y_i-\mu)^2dt_i)" />

where Y = {y<sub>1</sub>,..,y<sub>N</sub>} and dt<sub>i</sub> is the time spent in state y<sub>i</sub>.

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
        
----------

approximate_bayesian
----------

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
