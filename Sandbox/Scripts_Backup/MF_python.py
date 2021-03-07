#!/usr/bin/env python3
#%%

""" The practical of regex section """
# Author: Jingkai Sun (ks3020@ic.ac.uk)

from lmfit import Minimizer, Parameters, report_fit
import numpy as np
import matplotlib.pylab as plt

t = np.arange(0, 24, 2)
N = np.array([32500, 33000, 38000, 105000, 445000, 1430000, 3020000, 4720000, 5670000, 5870000, 5930000, 5940000])

np.random.seed(1234) # Set random seed for reproducibility

N_rand = N * (1 + np.random.normal(scale = 0.1, size = len(N)))

plt.plot(t, N_rand, 'r+', markersize = 15, markeredgewidth = 2, label = 'Data')
plt.xlabel('t', fontsize = 20)
plt.ylabel(r'$N$', fontsize = 20)
plt.ticklabel_format(style="scientific", scilimits=[0, 3])

# Create object for storing parameters
params_linear = Parameters()

# Add parameters and initial values to it
params_linear.add('a', value = 1)
params_linear.add('b', value = 1)
params_linear.add('c', value = 1)
params_linear.add('d', value = 1)

# Write down the objective function that we want to minimize
def residuals_linear(params, t, data):
    """
    Calculate cubic growth and substract data
    """
    
    # Get an ordered dictionary of parameter values
    v = params.valuesdict()
    
    # Cubic model
    model = v['a'] * t ** 3 + v['b'] * t ** 2 + v['c'] * t + v['d']
    return model - data # Return residuals

# Create a Minimizer object
minner = Minimizer(residuals_linear, params_linear, fcn_args=(t, np.log(N_rand)))

# Perform the minimization
fit_linear_NLLS = minner.minimize()
report_fit(fit_linear_NLLS)

# Using OLS
fit_linear_OLS = np.polyfit(t, np.log(N_rand), 3)
print(fit_linear_OLS)

## Comparing the NLLS and OLS fits
fit_linear_NLLS.params
par_dict = fit_linear_NLLS.params.valuesdict().values()

## Transform these into an array
par = np.array(list(par_dict))

## Check the differences in the parameter values obtined with lmfit and polyfit
print(fit_linear_OLS - par)

# Calculating the residuals
## Construct the fitted polynomial equation
my_poly = np.poly1d(fit_linear_OLS)
## Compute predicted values
ypred = my_poly(t)
## Calculating residuals
residuals = ypred - np.log(N_rand)

residuals_NLLS = residuals_linear(fit_linear_NLLS.params, t, np.log(N_rand))
# %%
