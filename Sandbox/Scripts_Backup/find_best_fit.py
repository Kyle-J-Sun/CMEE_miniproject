#!/usr/bin/env python3

""" The function of finding the best fit """
# Script Name: find_best_fit.py
# Author: Jingkai Sun (ks3020@ic.ac.uk)

from numpy import exp
import numpy as np
from model_comparison import get_AIC, get_BIC
from lmfit import Minimizer, Parameters, report_fit

def add_params(params, values, min_ = None, max_ = None):
    ''' Function for fitting the non-linear model '''
    # Create object and add parameters and initial values to it.
    params_model = Parameters()
    for i in range(len(params)):
        if min_ != None and max_ != None:
            params_model.add(params[i], value=values[i],
                             vary=True, min = min_[i], max = max_[i])
        else:
            params_model.add(params[i], value=values[i],
                             vary=True, min=None, max=None)
    return params_model

def find_best_fit(x, y, params, test_number, residuals_model, verbose=False, method='uniform', 
                  value_range=None, mean_sd_norm=None, min_=None, max_=None, log_scale = False,
                  log_scale_residuals_model = None):
    values = []
    AIC_values = []
    BIC_values = []
    num = test_number
    k = 0
    while k <= num:
        k += 1
        if method == 'uniform':
            params_fit = add_params(
                params=params,
                values=[round(np.random.uniform(elem[0], elem[1]), ndigits=3)
                        for elem in value_range],
                min_= min_,
                max_= max_
            )
            
        elif method == 'normal':
            params_fit = add_params(
                params = params,
                values = [round(np.random.normal(loc=elem[0], scale=elem[1]), ndigits=3)
                          for elem in mean_sd_norm],
                min_=min_,
                max_=max_
            )

        try:
            minner = Minimizer(residuals_model, params_fit, fcn_args=(x, y))
            # Perform the minimisation
            fit_model = minner.minimize()
        except:
            # print("error")
            continue
        
        if log_scale == True:
            AIC_value = get_AIC(log_scale_residuals_model(fit_model.params, x, exp(y)),
                                y, fit_model.params)

            BIC_value = get_BIC(log_scale_residuals_model(fit_model.params, x, exp(y)),
                                y, fit_model.params)
        else:
            AIC_value = get_AIC(fit_model.residual,
                                y, fit_model.params)

            BIC_value = get_BIC(fit_model.residual,
                                y, fit_model.params)
            
        values.append(params_fit)
        AIC_values.append(AIC_value)
        BIC_values.append(BIC_value)
        
    minner = Minimizer(
        residuals_model, values[AIC_values.index(min(AIC_values))], fcn_args=(x, y))
    # Perform the minimisation
    fit_model = minner.minimize()
    
    min_AIC = min(AIC_values)
    min_BIC = BIC_values[AIC_values.index(min(AIC_values))]
    print("Minimum AIC: ", min_AIC)
    print("Minimum BIC: ", min_BIC, '\n')

    print("Initial values used:")
    for key, value in dict(values[AIC_values.index(min(AIC_values))].valuesdict()).items():
        print("%s: %f" % (key, value))
    print('\n')
    
    print("Converged Info:")
    fit_model.params.pretty_print(columns=['value', 'stderr', 'min', 'max'])
    print('\n')

    if verbose == True:
        print(report_fit(fit_model))
    return fit_model, min_AIC, min_BIC

    
