import math
import matplotlib
matplotlib.rcParams['mathtext.fontset'] = 'stix'
matplotlib.rcParams['font.family'] = 'STIXGeneral'

def pmf(x, l, p, threshold = 1e-8):
    """
    Probability mass function for the derived distribution
    Parameters:
        x : int
        Value for which to calculate the pmf
        
        l : int
        Lambda (Poisson parameter)
        
        p : int
        Probability (Binomial parameter)
    
    Returns:
        f : float
        Value of PMF
    """
    x = int(x)
    l = int(l)
    r=x
    n = math.e**(-l)*l**r*p**x*(1-p)**(r-x)
    d = math.factorial(x)*math.factorial(r-x)
    f = n/d
    while True:
        r += 1
        n = math.e**(-l)*l**r*p**x*(1-p)**(r-x)
        d = math.factorial(x)*math.factorial(r-x)
        f += n/d
        if n < threshold*d:
            return f

def pmf_mean(fun, threshold=1e-6):
    """
    Function to calculate the expected value of an arbitrary distribution
    Parameters: 
        fun : callable
        Probability mass function for which to calculate the expected value
        
        threshold : float
        Cutoff threshold at which to stop iterating
        
    Returns:
        mu : float
        Expected value of PMF
    """
    x = 1
    contrib = fun(x)
    mu = contrib
    x += 1
    while contrib > threshold:
        contrib = x*fun(x)
        mu += contrib
        x += 1
    return mu

def pmf_var(fun, threshold=1e-6):
    """
    Function to calculate the variance of an arbitrary distribution
    Parameters: 
        fun : callable
        Probability mass function for which to calculate the variance
        
        threshold : float
        Cutoff threshold at which to stop iterating
        
    Returns:
        var : float
        Variance of PMF
    """
    mu = pmf_mean(fun)
    x = 1
    contrib = fun(x)
    var = contrib
    x += 1
    while contrib > threshold:
        contrib = x**2*fun(x)
        var += contrib
        x += 1
    var -= mu**2
    return var
    

l = 20
p = 0.5

ints = np.arange(round(1.5*l))
vals = np.zeros(len(ints))

for i, integer in enumerate(ints):
    vals[i] = pmf(integer, l, p)

fig, ax = plt.subplots()
ax.bar(ints,vals)
ax.set_xlabel("x")
ax.set_ylabel("P(X=x)")
ax.set_title("Probability distribution for the derived distribution\n" + "$\lambda = 20$, $p = 0.5$")

mu = pmf_mean(lambda x: pmf(x, l, p))

print("Mean = {}".format(mu))

var = pmf_var(lambda x: pmf(x, l, p))

print("Variance = {}".format(var))
