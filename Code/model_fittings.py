#!/usr/bin/env python3

""" Functions of fitting models """
# Script Name: model_fittings.py
# Author: Jingkai Sun (ks3020@ic.ac.uk)

import numpy as np
import pandas as pd
from lmfit import minimize, Parameters, report_fit
from numpy import exp, log, poly1d, pi
import matplotlib.pylab as plt

# Define all residual functions
def resi(fit, t, data):
    ''' Function to compute resiudals of fitted model using OLS method'''
    # Construct the fitted polynomial equation
    my_fit = poly1d(fit)
    #Compute predicted values
    ypred = my_fit(t)
    # Calculating residuals
    return ypred - data


def resi_bri(params, t, data):
    ''' create a briere model and subtract data '''
    # Get an ordered dictionary of parameter values
    v = params.valuesdict()
    # Briere Model
    model = v['B0'] * t * (t - v['T0']) * (abs(v['Tm'] - t) ** (1/2)) * \
        ((t < v['Tm']).astype(float)) * ((t > v['T0']).astype(float))
    return model - data


def resi_sch_log(params, T, log_data):
    ''' create a SchoolField model and subtract data '''
    k = 8.617e-5
    # Get an ordered dictionary of parameter values
    v = params.valuesdict()
    B0 = v['B0']
    E = v['E']
    Th = v['Th']
    Eh = v['Eh']
    # SchoolField Model
    # model = log(B0) - E * (T + 1/(283.15 * k)) - log(1 + exp(Eh * (1/(k*Th) - T)))
    model = log(B0) - E/k * (1/T - 1/283.15) - \
        log(1 + exp(Eh/k * (1/Th - 1/T)))
    return model - log_data


def resi_sch_orig(params, T, data):
    ''' create a SchoolField model and subtract data '''
    k = 8.617e-5
    # Get an ordered dictionary of parameter values
    v = params.valuesdict()
    B0 = v['B0']
    E = v['E']
    Th = v['Th']
    Eh = v['Eh']
    # SchoolField Model
    model = log(B0) - E/k * (1/T - 1/283.15) - \
        log(1 + exp(Eh/k * (1/Th - 1/T)))
    return exp(model) - data


def get_AIC(residuals, data, params):
    ''' Functions to get AIC value'''
    n = len(data)
    RSS = sum(residuals ** 2)
    p_model = len(params)
    return n + 2 + n * log((2 * pi) / n) + n * log(RSS) + 2 * p_model


def get_BIC(residuals, data, params):
    ''' Functions to get BIC value'''
    n = len(data)
    RSS = sum(residuals ** 2)
    p_model = len(params)
    return n + 2 + n * log((2 * pi) / n) + n * log(RSS) + p_model * log(n)


def linear_models(id, data):
    ''' Fitting three linear model using OLS method '''
    T = np.asarray(data["ConTemp"][data['ID'] == id])  # Temperature in kelvin
    B = np.asarray(data["OriginalTraitValue"][data['ID'] == id])  # Unlogged Trait value

    # Get fitting info of quadratic model
    fit_quad = np.polyfit(T, B, 2)
    rss2,tss2,rsq2,aic2,bic2,aicc2 = get_fitting_info_linear(T, B, fit_quad)

    results_all_quad = {
        "id": [id],
        "n" : len(B),
        "quad.a": fit_quad[0],
        "quad.b1": fit_quad[1],
        "quad.b2": fit_quad[2],
        "quad.Rsq": rsq2,
        "quad.aic": aic2,
        "quad.bic": bic2,
        "quad.aicc": aicc2
    }

    # Get fitting info of cubic model
    fit_cubic = np.polyfit(T, B, 3)
    rss3,tss3,rsq3,aic3,bic3,aicc3 = get_fitting_info_linear(T,B,fit_cubic)

    results_all_cubic = {
        "id": [id],
        "n": len(B),
        "cubic.a": fit_cubic[0],
        "cubic.b1": fit_cubic[1],
        "cubic.b2": fit_cubic[2],
        "cubic.b3": fit_cubic[3],
        "cubic.Rsq": rsq3,
        "cubic.aic": aic3,
        "cubic.bic": bic3,
        "cubic.aicc": aicc3
    }

    return pd.DataFrame(results_all_quad), pd.DataFrame(results_all_cubic), fit_quad, fit_cubic


def Briere_model(id, data, min_round, max_round):
    """ Briere model fitting  """
    #Variables
    sub = data.loc[data["ID"] == id]
    T = np.asarray(data['ConTemp'][data['ID'] == id])  # Temperature
    B = np.asarray(data['OriginalTraitValue'][data['ID'] == id])  # Trait Value
    # B0 = sub.loc[sub["ConTemp"] == np.min(sub["ConTemp"]), "OriginalTraitValue"]
    B0_min = np.min(B)
    B0_max = np.max(B)
    
    # if len(B0) >= 2:
    #     B0 = np.mean(B0)

    # Starting values
    B0 = 0.01
    T0 = min(T)
    Tm = max(T)

    # Initiate an empty dictionary for results
    results_all = {"id": [id],
                   "n": len(B),
                   "bri.B0": [B0],
                   "bri.T0": [T0],
                   "bri.chisqr": [np.NAN],
                   "bri.Rsq": [np.NAN],
                   "bri.aic": [np.NAN],
                   "bri.bic": [np.NAN],
                   "bri.aicc": [np.NAN]
                   }

    best_p = 0
    round = 0
    while round < max_round:
        round += 1
        # If TPC converged within mean rounds, or max round is exceeded, break the loop
        if results_all['bri.aicc'] != [np.NaN] and round > min_round and results_all["bri.Rsq"] >= -1:
            break

        #Starting values defined for rounds
        #For the first round use setted starting values
        if round == 1:
            p = Parameters()
            p.add_many(("B0", B0), ("T0", T0), ("Tm", Tm))

        #For other rounds, assign B0 with random number distributed between 0 and B0 * 2
        else:
            p = Parameters()
            p.add("B0", value=np.random.normal(B0, 0.10 * B0))
            p.add("T0", value=np.random.normal(T0, 20))
            p.add("Tm", value=np.random.normal(Tm, 40))

        try:
            output = minimize(resi_bri, p, args=(T, B))
            rss,tss,rsq,aic,bic,aicc = get_fitting_info(T, B, output, resi_bri)

            # Replace the fitting results if aicc is smaller,
            # or if the previous fit didn't converge
            if aicc < results_all["bri.aicc"] or results_all["bri.aicc"] == [np.NaN] and rsq >= -1:
                best_p = output.params
                results_all = {
                    "id": [id],
                    "n": len(B),
                    "bri.B0": [output.params["B0"].value],
                    "bri.T0": [output.params["T0"].value],
                    "bri.Tm": [output.params["Tm"].value],
                    "bri.chisqr": [output.chisqr],
                    "bri.Rsq": rsq,
                    "bri.aic": aic,
                    "bri.bic": bic,
                    "bri.aicc": aicc
                }
        except ValueError:
            pass
        continue

    return pd.DataFrame(results_all), best_p

# Schoolfield model simplified, appropriate for without low termperature data


def Schoolfield_model(id, data, min_round, max_round):
    """ fitting simplified Schoolfield model with some rounds """
    sub = data.loc[data["ID"] == id]

    T = np.asarray(sub["ConTemp_Kelvin"])  # Temperature in kelvin
    B = np.asarray(sub["OriginalTraitValue"])  # Unlogged Trait value
    B_log = np.asarray(sub["LoggedOriginalTraitValue"])  # Logged Trait value
    
    # Th_min = np.min(T)
    # Th_max = np.max(T)
    B0_min = np.min(B)
    B0_max = np.max(B)
    B0_mean = np.mean(np.array(B0_min, B0_max))

    T_max_C = sub["ConTemp_Kelvin"][sub["LoggedOriginalTraitValue"] == np.max(sub["LoggedOriginalTraitValue"])]

    if len(T_max_C) >= 2:
        T_max_C = np.mean(T_max_C)

    # Starting values(matters on converge)
    B0 = sub.loc[sub["ConTemp_Kelvin"] == np.min(sub["ConTemp_Kelvin"]), "OriginalTraitValue"]
    
    if len(B0) >= 2:
        B0 = np.min(B0)
        
    B0 = float(B0)
    E, Eh = Arrhenius_model(id, data)
    Th = float(T_max_C)

    # Add starting values with bounding
    p = Parameters()
    p.add("B0", value=B0)
    p.add("E", value=E, min=0.000001, max=40)
    p.add("Th", value=Th, min=273.15 - 10, max=273.15 + 150)
    p.add("Eh", value=Eh, min=0.000001, max=100)

    # Initiate an empty dictionary for results
    results_all = {
        "id": [id],
        "n": len(B),
        "sch.B0": [B0],
        "sch.E": [E],
        "sch.Th": [Th],
        "sch.Eh": [Eh],
        "sch.chisqr": [np.NAN],
        "sch.Rsq": [np.NAN],
        "sch.aic": [np.NAN],
        "sch.bic": [np.NAN],
        "sch.aicc": [np.NAN],
        "sch.Rsq_log": [np.NAN],
        "sch.aic_log": [np.NAN],
        "sch.bic_log": [np.NAN],
        "sch.aicc_log": [np.NAN]
    }

    best_p = 0
    
    round = 0
    while round < max_round:

        round += 1
        if results_all["sch.aicc"] != [np.NaN] and round > min_round and results_all['sch.Rsq'] >= -1:
            break
        #Starting values defined for rounds
        #For the first round use set starting values
        if round == 1:
            p = Parameters()
            p.add("B0", value=B0)
            p.add("E", value=E, min=0.000001, max = 40)
            p.add("Th", value=Th, min=273.15 - 10, max=273.15 + 150)
            p.add("Eh", value=Eh, min=0.000001, max = 100)
        else:
            p = Parameters()
            p.add("B0", value=np.random.normal(B0, B0 * 0.15))
            p.add("E", value=np.random.normal(E, E * 0.15), min=0.000001, max=40)
            p.add("Th", value=Th, min=273.15 - 10, max=273.15 + 150)
            p.add("Eh", value=np.random.normal(Eh, Eh * 0.15), min=0.000001, max=100)

        #Results for model fitting on logged values
        try:
            output = minimize(resi_sch_log, p, args=(T, B_log))
            ss_log,tss_log,rsq_log,aic_log,bic_log,aicc_log = get_fitting_info(T, B_log, output, resi_sch_log)

            # Calculate unlogged fitting
            rss,tss,rsq,aic,bic,aicc = get_fitting_info(T, B, output, resi_sch_orig)

            if aicc < results_all["sch.aicc"] or results_all["sch.aicc"] == [np.NaN] and rsq >= -1:
                best_p = output.params
                results_all = {
                    "id": [id],
                    "n": [len(B)],
                    "sch.B0": [output.params["B0"].value],
                    "sch.E": [output.params["E"].value],
                    "sch.Th": [output.params["Th"].value],
                    "sch.Eh": [output.params["Eh"].value],
                    "sch.chisqr": [output.chisqr],
                    "sch.Rsq": rsq,
                    "sch.aic": aic,
                    "sch.bic": bic,
                    "sch.aicc": aicc,
                    "sch.Rsq_log": rsq_log,
                    "sch.aic_log": bic_log,
                    "sch.bic_log": aic_log,
                    "sch.aicc_log": aicc_log
                }
        except ValueError:
            pass
        continue

    return (pd.DataFrame(results_all), best_p)
    

def Arrhenius_model(id, data):
    ''' Functions to get Ea and Eh parameters'''
    # k = 8.617e-5
    sub = data[data["ID"] == id]

    B_log = np.asarray(sub["LoggedOriginalTraitValue"])
    Trans_T = np.asarray(sub["Transferred_ConTemp"])

    fit_arr = np.polyfit(Trans_T, B_log, 1)

    if np.isnan(fit_arr[1]):
        return 10, 50
    else:
        return abs(fit_arr[1]), abs(fit_arr)[1] * 10


def get_fitting_info(x, y, output, resid):
    ''' Getting all fitting infos'''
    resi = resid(output.params, x, y)
    rss = sum(resi ** 2)
    tss = sum((y - np.mean(y)) ** 2)
    rsq = 1 - (rss/tss)
    aic = get_AIC(resi, y, output.params)
    bic = get_BIC(resi, y, output.params)
    K = len(output.params)
    aicc = aic + 2 * K * ((K + 1)/(len(x) - K - 1))
    return rss, tss, rsq, aic, bic, aicc


def get_fitting_info_linear(x, y, output):
    resid = resi(output, x, y)
    rss = sum(resid ** 2)
    tss = sum((y - np.mean(y)) ** 2)
    rsq = 1 - (rss/tss)
    aic = get_AIC(resid, y, output)
    bic = get_BIC(resid, y, output)
    K = len(output)
    aicc = aic + 2 * K * ((K + 1)/(len(x) - K - 1))
    return rss, tss, rsq, aic, bic, aicc


def fitting_plot_OLS(x, fit_model=None, residuals_model=None, name=None, color=None):
    ''' Function to plot OLS fitted model '''
    for elem in range(len(fit_model)):
        plt.rcParams['figure.figsize'] = [20, 15]
        result_linear = np.poly1d(fit_model[elem])(x)
        plt.plot(x, result_linear, '.', markerfacecolor=color[elem],
                 markeredgecolor=color[elem], markersize=15, label=name[elem])
        # Get a smooth curve by plugging a time vector to the fitted linear model
        t_vec = np.linspace(min(x), max(x), 1000)
        ypred_smooth = np.poly1d(fit_model[elem])(t_vec)
        plt.plot(t_vec, ypred_smooth, color[elem], linestyle='--', linewidth=2)


def models_visualisation(id, d, params):
    """ To visualise 4 models """
    sub = d[d["ID"] == id]
    T = np.asarray(sub["ConTemp"])  # Temperature in kelvin
    T_kelvin = np.asarray(sub["ConTemp_Kelvin"])
    # Unlogged Trait value
    B = np.asarray(sub["OriginalTraitValue"])
    B_name = sub["OrignalTraitName"].iloc[0]
    B_unit = sub["OriginalTraitUnit"].iloc[0]

    fitting_plot_OLS(
        x=T,
        fit_model=(params[0], params[1]),
        name=("Quadratic", "Cubic"),
        color=("red", "Orange")
    )
    if params[2] != 0:
        # Briere model
        result_model = B + resi_bri(params[2], T, B)
        plt.plot(T, result_model, '.', markerfacecolor="blue",
                 markeredgecolor="blue", markersize=15, label="Briere")
        # Get a smooth curve by plugging a time vector to the fitted non-linear model
        t_vec = np.linspace(min(T), max(T), 1000)
        N_vec = np.ones(len(t_vec))
        residual_smooth = resi_bri(params[2], t_vec, N_vec)
        plt.plot(t_vec, residual_smooth + N_vec,
                 "blue", linestyle='--', linewidth=2)

    if params[3] != 0 and params[3] != None:
        # Simplified Schoolfield
        result_model = B + resi_sch_orig(params[3], T_kelvin, B)
        plt.plot(T, result_model, '.', markerfacecolor="magenta",
                 markeredgecolor="magenta", markersize=15, label="Schoolfield")
        # Get a smooth curve by plugging a time vector to the fitted non-linear model
        t_vec = np.linspace(min(T_kelvin), max(T_kelvin), 1000)
        N_vec = np.ones(len(t_vec))
        residual_smooth = resi_sch_orig(params[3], t_vec, N_vec)
        plt.plot(t_vec - 273.15, residual_smooth + N_vec,
                 "magenta", linestyle='--', linewidth=2)

    # Plot data points
    plt.plot(T, B, 'r+', markersize=15,
             markeredgewidth=2, label='Data')
    # Plot legend
    plt.legend(fontsize=20)
    plt.xlabel('Temperature(°C)', fontsize=20)
    plt.ylabel('%s\n%s' % (B_name, B_unit), fontsize=20)
    plt.ticklabel_format(style='scientific', scilimits=[0, 3])
