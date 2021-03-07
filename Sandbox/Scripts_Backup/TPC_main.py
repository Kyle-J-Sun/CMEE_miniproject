#!/usr/bin/env python3

""" To fit the model for Thermal Performence Curves """
# Script Name: TPC_main.py
# Author: Jingkai Sun (ks3020@ic.ac.uk)

#TODO:
# 1. To drop all negative and zero values. (finished)
# 2. Define get_starting_values() function to get start value (finished)
# 3. To fit School model with apprioprate starting values (finished)
# 4. add the bounding (add_many) (finished)
# 5. beautify the printing resutls (finished)
# 6. Midify sampling implementation for finding more appropriate parameters (finished)
# 7. Analysis of the model for finding biological interpretation
# 8. Writing LaTeX
# 9. !!!Modified AIC and BIC of logged fieldmodel!!! (finished)

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
from model_visualisation import fitting_plot_NLLS, fitting_plot_OLS
from find_best_fit import find_best_fit
    

def main(argv):
    ''' The main function '''
    briere_test_number = int(argv[1])
    school_test_number = int(argv[2])
    
    global k
    k = 8.617e-5
    d = pd.read_csv("../Data/ThermRespData.csv")
    d.columns.values
    # Drop all the values that equals to 0
    d = d[d['OriginalTraitValue'] > 0]
    d['LoggedOriginalTraitValue'] = log(d['OriginalTraitValue'])
    
    # Write header
    f = open("../Results/fitting_info_aic.csv", 'w')
    f.write("ID,Linear,Quadratic,Cubic,Briere,Schoolfield,Best_fit" + '\n')

    # Write header
    s = open("../Results/fitting_info_bic.csv", 'w')
    s.write("ID,Linear,Quadratic,Cubic,Briere,Schoolfield,Best_fit" + '\n')

    
    for i in list(d['ID'].unique()):
        sub = d[d['ID'] == i]
        print("Fitting the ID %d ..." % i)
        print("-" * 80)
        if len(sub) > 4:
            # Model Fitting
            ''' Fitting three linear model using OLS method '''
            fit_lin = np.polyfit(sub['ConTemp'], sub['OriginalTraitValue'], 1)
            lin_residual = residual(
                fit_lin, sub['ConTemp'], sub['OriginalTraitValue'])

            fit_quad = np.polyfit(
                sub['ConTemp'], sub['OriginalTraitValue'], 2)
            quad_residual = residual(
                fit_quad, sub['ConTemp'], sub['OriginalTraitValue'])
            
            fit_poly3 = np.polyfit(
                sub['ConTemp'], sub['OriginalTraitValue'], 3)
            poly3_residual = residual(
                fit_poly3, sub['ConTemp'], sub['OriginalTraitValue'])
            
            fitting_plot_OLS(
                x=sub['ConTemp'],
                fit_model=(fit_lin, fit_quad, fit_poly3),
                name=("Linear", "Quadratic", "Polynomial"),
                color=("blue", "red", "Orange")
            )
            
            AIC_linear = get_AIC(lin_residual, sub['OriginalTraitValue'], fit_lin)
            BIC_linear = get_BIC(lin_residual, sub['OriginalTraitValue'], fit_lin)
            
            AIC_quad = get_AIC(quad_residual, sub['OriginalTraitValue'], fit_quad)
            BIC_quad = get_BIC(quad_residual, sub['OriginalTraitValue'], fit_quad)
            
            AIC_cubic = get_AIC(poly3_residual, sub['OriginalTraitValue'], fit_poly3)
            BIC_cubic = get_BIC(poly3_residual, sub['OriginalTraitValue'], fit_poly3)
            
            plt.plot(sub['ConTemp'], sub['OriginalTraitValue'], 'r+', markersize=15,
                        markeredgewidth=2, label='Data')
            # Plot legend
            plt.legend(fontsize=20)
            plt.xlabel('Temperature(°C)', fontsize=20)
            plt.ylabel('%s\n%s' % (sub['OrignalTraitName'].unique()[0], sub['OriginalTraitUnit'].unique()[0]), fontsize=20)
            plt.ticklabel_format(style='scientific', scilimits=[0, 3])
            plt.savefig(
                "../Results/Linear/TPC_fit_polymodels%d.pdf" % i)
            plt.close()
            
            try:
                # Fitting Non-Linear model
                ''' Fitting briere model with original trait value'''
                print("Fitting information of Briere Model:")
                # np.random.seed(123)
                fit_briere, AIC_briere, BIC_briere = find_best_fit(
                    sub['ConTemp'],
                    sub['OriginalTraitValue'],
                    params_add=('B0', 'T0', 'Tm'),
                    method = 'normal',
                    mean_sd_norm=([0.01, 2.1], [-10, 10], [36.85, 15.55]),
                    test_number=briere_test_number,
                    residuals_model=residuals_briere
                )
                
                fitting_plot_NLLS(
                    x=sub['ConTemp'],
                    y=sub['OriginalTraitValue'],
                    fit_model=(fit_briere,),
                    residuals_model=(residuals_briere,),
                    name=("Briere",),
                    color=("grey",)
                )
                
                # Plot data points
                plt.plot(sub['ConTemp'], sub['OriginalTraitValue'], 'r+', markersize=15,
                        markeredgewidth=2, label='Data')
                # Plot legend
                plt.legend(fontsize=20)
                plt.xlabel('Temperature(°C)', fontsize=20)
                plt.ylabel('%s\n%s' % (sub['OrignalTraitName'].unique()[0], sub['OriginalTraitUnit'].unique()[0]),
                        fontsize=20)
                plt.ticklabel_format(style='scientific', scilimits=[0, 3])
                plt.savefig("../Results/Briere/TPC_fit_briere%d.pdf" % i)
                plt.close()
            except:
                AIC_briere = 999999 # Asign 999999 to Briere's AIC and BIC standing for fitting failure
                BIC_briere = 999999
            try:
                '''Fitting logged schoolfield model with logged data'''
                # np.random.seed(111)
                print("Fitting information of Logged SchoolField Model:")
                fit_school_log, AIC_school, BIC_school = find_best_fit(
                    sub['ConTemp'],
                    sub['LoggedOriginalTraitValue'],
                    params_add=('B0', 'E', 'E_h', 'T_h'),
                    method='uniform',
                    # mean_sd_norm=([0.5, 2.1], [1.2, 0.7],
                    #                 [45.21, 44.98], [36.85, 10.5]),
                    value_range=([-5, 5], [0.01, 2], [0.03, 3], [26.85, 43.15]),
                    test_number=school_test_number,
                    log_scale=True,
                    log_scale_residuals_model=residuals_school_original,
                    residuals_model=residuals_school_log,
                    min_=(None, 0, 0, None),
                    max_=(None, None, None, None)
                )
                
                fitting_plot_NLLS(
                    x=sub['ConTemp'],
                    y=sub['OriginalTraitValue'],
                    fit_model=(fit_school_log,),
                    residuals_model=(residuals_school_original,),
                    name=("SchoolField",),
                    color=("magenta",),
                    log_scale=True
                )
                
                # Plot data points
                plt.plot(sub['ConTemp'], sub['OriginalTraitValue'], 'r+', markersize=15,
                         markeredgewidth=2, label='Data')
                # Plot legend
                plt.legend(fontsize=20)
                plt.xlabel('Temperature(°C)', fontsize=20)
                plt.ylabel('%s\n%s' % (sub['OrignalTraitName'].unique()[0], sub['OriginalTraitUnit'].unique()[0]),
                           fontsize=20)
                plt.ticklabel_format(style='scientific', scilimits=[0, 3])
                plt.savefig("../Results/School/TPC_fit_school%d.pdf" % i)
                plt.close()
                
                fitting_plot_NLLS(
                    x=sub['ConTemp'],
                    y=sub['LoggedOriginalTraitValue'],
                    fit_model=(fit_school_log,),
                    residuals_model=(residuals_school_log,),
                    name=("Log-transformed SchoolField",),
                    color=("magenta",)
                )

                plt.plot(sub['ConTemp'], sub['LoggedOriginalTraitValue'], 'r+', markersize=15,
                         markeredgewidth=2, label='Data')
                # Plot legend
                plt.legend(fontsize=20)
                plt.xlabel('Temperature(C)', fontsize=20)
                plt.ylabel('%s\n%s(Log-Transformed)' % (sub['OrignalTraitName'].unique()[0], 
                                                        sub['OriginalTraitUnit'].unique()[0]), fontsize=20)
                plt.ticklabel_format(style='scientific', scilimits=[0, 3])
                plt.savefig(
                    "../Results/Log_School/Thermal_fit_log_School%d.pdf" % i)
                plt.close()
            except:
                AIC_school = 999999  # Asign 999999 to School's AIC and BIC standing for fitting failure
                BIC_school = 999999      
            # Print AIC,BIC for linear models
            print('\n')
            print("Fitting Infomation of Linear Models:", '\n')
            # print('-'*80)
            ################################################################################
            # Linear

            print("Linear AIC: ", AIC_linear)
            print("Linear BIC: ", BIC_linear)
            
            ## Quadratic

            print("Quadratic AIC: ", AIC_quad)
            print("Quadratic BIC: ", BIC_quad)
            
            ## Polynomial

            print("Poly3 AIC: ", AIC_cubic)
            print("Poly3 BIC: ", BIC_cubic)
        
            #################################################################################
            print("-"*80, '\n\n')
            
            model_AICs = {"Linear": AIC_linear,
                            "Quadratic": AIC_quad,
                            "Cubic": AIC_cubic,
                            "Briere": AIC_briere,
                            "Schoolfield": AIC_school}
            
            model_BICs = {"Linear": BIC_linear,
                            "Quadratic": BIC_quad,
                            "Cubic": BIC_cubic,
                            "Briere": BIC_briere,
                            "Schoolfield": BIC_school}
            
            # Save AIC BIC infos
            f.write("%d,%3f,%3f,%3f,%3f,%3f,%s" % (i,AIC_linear, AIC_quad, AIC_cubic,
                                                AIC_briere, AIC_school, min(model_AICs, key=model_AICs.get)) + '\n')
            s.write("%d,%3f,%3f,%3f,%3f,%3f,%s" % (i,BIC_linear, BIC_quad, BIC_cubic,
                                                BIC_briere, BIC_school, min(model_BICs, key=model_BICs.get)) + '\n')
        else:
            print("Fitting Failed! ID %d has only %d observation(s)" % (i, len(sub)))
     
    f.close()
    s.close()
    
    return 0

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

           # # Plot data points
    # plt.plot(sub['ConTemp'], sub['OriginalTraitValue'], 'r+', markersize=15,
    #             markeredgewidth=2, label='Data')
    # # Plot legend
    # plt.legend(fontsize=20)
    # plt.xlabel('Temperature(C)', fontsize=20)
    # plt.ylabel('%s\n%s' % (sub['OrignalTraitName'].unique()[0], sub['OriginalTraitUnit'].unique()[0]),
    #             fontsize=20)
    # plt.ticklabel_format(style='scientific', scilimits=[0, 3])
    # plt.savefig("../Results/Main_plots/TPC_fit_python%d.pdf" % i)
    # plt.close()

    # SchoolField FITTING IN other paper:
    # B0(h-1) = 1.42 (-0.056 to 2.91)
    # Th(K) = 314.7 (314 to 315.3)
    # E(kJ) = -5.43 (-59.8 to 49.0)
    # Eh(kJ) = 687.9 (402.1 to 973.7)
    # El(kJ) = -141.1 (-182.2 to -100.1)
    # Tl(K) = 297.7 (286 to 309.3)
