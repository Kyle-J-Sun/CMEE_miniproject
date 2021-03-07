#!/usr/bin/env python3

""" To fit the model for Thermal Performence Curves """
# Script Name: TPC_main.py
# Author: Jingkai Sun (ks3020@ic.ac.uk)

# Import packages
import sys
import scipy as sc
from lmfit import Minimizer, Parameters, report_fit
from numpy import exp, log, pi
import numpy as np
import pandas as pd
from multiprocessing import Process
import matplotlib.pylab as plt
from model_residuals import *
from model_comparison import get_AIC, get_BIC
from model_visualisation import fitting_plot_OLS
from find_best_fit import find_best_fit

def linear_fitting(sub):
    ''' Fitting three linear model using OLS method '''
    # fit_lin = np.polyfit(sub['ConTemp'], sub['OriginalTraitValue'], 1)
    # lin_residual = residual(
    #     fit_lin, sub['ConTemp'], sub['OriginalTraitValue'])

    fit_quad = np.polyfit(
        sub['ConTemp'], sub['OriginalTraitValue'], 2)
    quad_residual = residual(
        fit_quad, sub['ConTemp'], sub['OriginalTraitValue'])

    fit_poly3 = np.polyfit(
        sub['ConTemp'], sub['OriginalTraitValue'], 3)
    poly3_residual = residual(
        fit_poly3, sub['ConTemp'], sub['OriginalTraitValue'])

    # AIC_linear = get_AIC(
    #     lin_residual, sub['OriginalTraitValue'], fit_lin)
    # BIC_linear = get_BIC(
    #     lin_residual, sub['OriginalTraitValue'], fit_lin)

    AIC_quad = get_AIC(
        quad_residual, sub['OriginalTraitValue'], fit_quad)
    BIC_quad = get_BIC(
        quad_residual, sub['OriginalTraitValue'], fit_quad)

    AIC_cubic = get_AIC(
        poly3_residual, sub['OriginalTraitValue'], fit_poly3)
    BIC_cubic = get_BIC(
        poly3_residual, sub['OriginalTraitValue'], fit_poly3)

    fitting_plot_OLS(
        x=sub['ConTemp'],
        fit_model=(fit_quad, fit_poly3),
        name=("Quadratic", "Polynomial"),
        color=("red", "Orange")
    )

    # Print AIC,BIC for linear models
    print("Fitting Infomation of Linear Models:", '\n')
    # print('-'*80)
    ################################################################################
    # Linear
    # print("Linear AIC: ", AIC_linear)
    # print("Linear BIC: ", BIC_linear)

    ## Quadratic

    print("Quadratic AIC: ", AIC_quad)
    print("Quadratic BIC: ", BIC_quad)

    ## Polynomial

    print("Poly3 AIC: ", AIC_cubic)
    print("Poly3 BIC: ", BIC_cubic)
    print('\n')
    #################################################################################
    return AIC_quad, BIC_quad, AIC_cubic, BIC_cubic

def briere_fitting(sub, briere_test_number):
    ''' Fitting briere model with original trait value'''
    try:
        # Fitting Non-Linear model
        print("Fitting Information of Briere Model:")
        # np.random.seed(123)
        fit_briere, AIC_briere, BIC_briere = find_best_fit(
            sub['ConTemp'],
            sub['OriginalTraitValue'],
            params=('B0', 'T0', 'Tm'),
            method='normal',
            mean_sd_norm=([0.01, 10], [10, 10], [36.85, 30]),
            test_number=briere_test_number,
            residuals_model=residuals_briere
        )

        plt.rcParams['figure.figsize'] = [20, 15]
        result_model = sub['OriginalTraitValue'] + fit_briere.residual
        plt.plot(sub['ConTemp'], result_model, '.', markerfacecolor="blue",
            markeredgecolor="blue", markersize=15, label="Briere")
        # Get a smooth curve by plugging a time vector to the fitted non-linear model
        t_vec = np.linspace(
            min(sub['ConTemp']), max(sub['ConTemp']), 1000)
        N_vec = np.ones(len(t_vec))
        residual_smooth = residuals_briere(
            fit_briere.params, t_vec, N_vec)
        plt.plot(t_vec, residual_smooth + N_vec,
                    "blue", linestyle='--', linewidth=2)
    except:
        # Asign 999999 to School's AIC and BIC standing for fitting failure
        AIC_briere = 999999
        BIC_briere = 999999
    return AIC_briere, BIC_briere

def school_fitting(sub, school_test_number):
    '''Fitting logged schoolfield model with logged data'''
    try:
        print("Fitting Information of Logged SchoolField Model:")
        fit_school_log, AIC_school, BIC_school = find_best_fit(
            sub['ConTemp'],
            sub['LoggedOriginalTraitValue'],
            params=('B0', 'E', 'E_h', 'T_h'),
            method='uniform',
            # mean_sd_norm=([0.5, 2.1], [1.2, 0.7],
            #                 [45.21, 44.98], [36.85, 10.5]),
            value_range=([-20, 20], [0.01, 10],
                            [0.03, 30], [28.85, 45]),
            test_number=school_test_number,
            log_scale=True,
            log_scale_residuals_model=residuals_school_original,
            residuals_model=residuals_school_log,
            min_=(None, 0, 0, None),
            max_=(None, None, None, None)
        )

        result_model = sub['OriginalTraitValue'] + residuals_school_original(
                fit_school_log.params, sub['ConTemp'], sub['OriginalTraitValue'])
        plt.plot(sub['ConTemp'], result_model, '.', markerfacecolor="magenta",
                    markeredgecolor="magenta", markersize=15, label="Schoolfield")
        # Get a smooth curve by plugging a time vector to the fitted non-linear model
        t_vec = np.linspace(
            min(sub['ConTemp']), max(sub['ConTemp']), 1000)
        N_vec = np.ones(len(t_vec))
        residual_smooth = residuals_school_original(
            fit_school_log.params, t_vec, N_vec)
        plt.plot(t_vec, residual_smooth + N_vec,
                    "magenta", linestyle='--', linewidth=2)
    except:
        AIC_school = 999999  # Asign 999999 to School's AIC and BIC standing for fitting failure
        BIC_school = 999999
    return AIC_school, BIC_school


def save_results(AIC_quad, AIC_cubic, AIC_briere, AIC_school,
                 BIC_quad, BIC_cubic, BIC_briere, BIC_school, i, f, s):
    ''' Save aic and bic results as local csv files'''
    model_AICs = {
                "Quadratic": AIC_quad,
                "Cubic": AIC_cubic,
                "Briere": AIC_briere,
                "Schoolfield": AIC_school
                }

    model_BICs = {
                "Quadratic": BIC_quad,
                "Cubic": BIC_cubic,
                "Briere": BIC_briere,
                "Schoolfield": BIC_school
                }
            
    # Save AIC BIC infos
    f.write("%d,%3f,%3f,%3f,%3f,%s" % (i,AIC_quad, AIC_cubic,
                                        AIC_briere, AIC_school, min(model_AICs, key=model_AICs.get)) + '\n')
    s.write("%d,%3f,%3f,%3f,%3f,%s" % (i,BIC_quad, BIC_cubic,
                                        BIC_briere, BIC_school, min(model_BICs, key=model_BICs.get)) + '\n')
    
    return

def data_cleaning():
    d = pd.read_csv("../Data/ThermRespData.csv")
    # d.columns.values
    # Drop all the values that equals to 0
    d = d[d['OriginalTraitValue'] > 0]
    d['LoggedOriginalTraitValue'] = log(d['OriginalTraitValue'])
    
    # Write header
    f = open("../Results/fitting_info_aic.csv", 'w')
    f.write("ID,Quadratic,Cubic,Briere,Schoolfield,Best_fit" + '\n')

    # Write header
    s = open("../Results/fitting_info_bic.csv", 'w')
    s.write("ID,Quadratic,Cubic,Briere,Schoolfield,Best_fit" + '\n')
    return d, f, s

def main(argv):
    ''' The main function '''
    briere_test_number = int(argv[1])
    school_test_number = int(argv[2])

    global k
    k = 8.617e-5
    d, f, s = data_cleaning()
    
    for i in range(1, 15):
    # for i in list(d['ID'].unique()):
        sub = d[d['ID'] == i]
        print("Fitting the ID %d ..." % i)
        print("-" * 80)
        if len(sub) > 4:
            # Model Fitting
            
            AIC_quad, BIC_quad, AIC_cubic, BIC_cubic = linear_fitting(sub)
            
            AIC_briere, BIC_briere = briere_fitting(sub, briere_test_number)
                
            AIC_school, BIC_school = school_fitting(sub, school_test_number)
                
            # Plot data points
            plt.plot(sub['ConTemp'], sub['OriginalTraitValue'], 'r+', markersize=15,
                        markeredgewidth=2, label='Data')
            
            # Plot legend
            plt.legend(fontsize=20)
            plt.xlabel('Temperature(°C)', fontsize=20)
            plt.ylabel('%s\n%s' % (sub['OrignalTraitName'].unique()[0], sub['OriginalTraitUnit'].unique()[0]),
                        fontsize=20)
            plt.ticklabel_format(style='scientific', scilimits=[0, 3])
            plt.savefig("../Results/TPC_fitting%d.pdf" % i)
            plt.close()
            
            print("-"*80, '\n')
            
            save_results(AIC_quad, AIC_cubic, AIC_briere, AIC_school,
                         BIC_quad, BIC_cubic, BIC_briere, BIC_school, i, f, s)

        else:
            print("Fitting Failed! ID %d has only %d observation(s)" % (i, len(sub)))
     
    f.close()
    s.close()
    
    return 0

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

