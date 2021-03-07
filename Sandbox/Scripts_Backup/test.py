import sys
import scipy as sc
from lmfit import Minimizer, Parameters, report_fit
from numpy import exp, log, pi
import numpy as np
import pandas as pd
from multiprocessing import Process, Pool
import matplotlib.pylab as plt
from model_residuals import *
from model_comparison import get_AIC, get_BIC
from model_visualisation import fitting_plot_NLLS, fitting_plot_OLS
from find_best_fit import find_best_fit


def main(argv):
    ''' The main function '''
    
    k = 8.617e-5


    d = pd.read_csv("../Data/ThermRespData.csv")
    d.columns.values
    # Drop all the values that equals to 0
    d = d[d['OriginalTraitValue'] > 0]
    d['LoggedOriginalTraitValue'] = log(d['OriginalTraitValue'])
    # d = d[d['LoggedOriginalTraitValue'] > -20]
    sub = d[d['ID'] == 2]
    pool = Pool(processes = 2)

    
    results = pool.apply_async(
        func=find_best_fit, 
        kwds={
            'x': sub['ConTemp'],
            'y': sub['OriginalTraitValue'],
            'params_add': ('B0', 'T0', 'Tm'),
            'method': 'normal',
            'mean_sd_norm': ([0.01, 2.1], [-10, 10], [36.85, 15.55]),
            'test_number': 12,
            'residuals_model': residuals_briere
        }
    )

    results2 = Process(
        target=find_best_fit,
        kwargs={
            'x': sub['ConTemp'],
            'y': sub['OriginalTraitValue'],
            'params_add': ('B0', 'E', 'E_h', 'T_h'),
            'method': 'uniform',
            'value_range': ([-3, 3], [0.01, 2], [0.03, 3], [26.85, 43.15]),
            'test_number': 20,
            'log_scale': True,
            'log_scale_residuals_model': residuals_school_original,
            'residuals_model': residuals_school_log,
            'min_': (None, 0, 0, None),
            'max_': (None, None, None, None)
        }
    )
    
    pool.close()
    pool.join()
    print(results)
    print(results2)
    return 0


if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)
