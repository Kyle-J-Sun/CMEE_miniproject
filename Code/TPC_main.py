#!/usr/bin/env python3

""" To fit the model for Thermal Performence Curves """
# Script Name: TPC_main.py
# Author: Jingkai Sun (ks3020@ic.ac.uk)

# Import packages
import sys
from lmfit import Minimizer, minimize, Parameters, report_fit
import numpy as np
import pandas as pd
import matplotlib.pylab as plt
from model_fittings import *
import os

def data_cleaning():
    d = pd.DataFrame(pd.read_csv("../Data/ThermRespData.csv"))
    d = d[d['OriginalTraitValue'] > 0]
    d['LoggedOriginalTraitValue'] = np.log(d['OriginalTraitValue'])
    d = d[["ID", "ConTemp", "OriginalTraitValue", "StandardisedTraitName", "LoggedOriginalTraitValue", "OrignalTraitName",  "OriginalTraitUnit"]]
    d["ConTemp_Kelvin"] = d["ConTemp"] + 273.15
    d["Transferred_ConTemp"] = 1/k*d['ConTemp_Kelvin']+1/283.15*k
    d["Transferred_kT"] = 1/k*d["ConTemp_Kelvin"]

    # Find the subset that has different biological indicators
    d_photosyn = d.loc[(d['StandardisedTraitName'] == d['StandardisedTraitName'].unique()[0]) | (d['StandardisedTraitName'] == d['StandardisedTraitName'].unique()[1])]
    d_respiration = d.loc[d['StandardisedTraitName'] == d['StandardisedTraitName'].unique()[2]]
    return d_photosyn, d_respiration


def model_fitting(data):
    # Arrhenius_after = pd.DataFrame(data=None)
    # Arrhenius_before = pd.DataFrame(data=None)
    quadratic = pd.DataFrame(data=None)
    cubic = pd.DataFrame(data=None)
    Briere = pd.DataFrame(data=None)
    Schoolfield = pd.DataFrame(data=None)
    data_name = data["StandardisedTraitName"].unique()[0].replace(" ", "")

    ## Model Fitting for Net Phototh data
    for id in data["ID"].unique():
        if data[data["ID"] == id].shape[0] > 5:
            print("Fitting Models of ID %d" % id)
            quad_info, cubic_info, fit_quad, fit_cubic = linear_models(id, data)
            quadratic = quadratic.append(quad_info)
            cubic = cubic.append(cubic_info)

            # print("Fitting Briere Modle of ID %d" % id)
            briere_info, p_bri = Briere_model(id, data, min_round=65, max_round=75)
            Briere = Briere.append(briere_info)

            # print("Fitting Simplified Schoolfield Modle of ID %d" % id)
            school_info, p_sch= Schoolfield_model(id, data, min_round=65, max_round=75)
            Schoolfield = Schoolfield.append(school_info)

            # Save pictures
            models_visualisation(id, data, (fit_quad, fit_cubic, p_bri, p_sch))
            try:
                plt.savefig("../Results/%s/Figures/TPC_fitting%d.pdf" % (data_name, id))
            except FileNotFoundError:
                print(os.system("mkdir -p ../Results/" + data_name + "/Figures"))
                plt.savefig("../Results/%s/Figures/TPC_fitting%d.pdf" % (data_name, id))
            plt.close()
        else:
            print("ID %d has only %d observation(s)" % (id, data[data["ID"] == id].shape[0]))

    print("-"*80, '\n')
    try:
        quadratic.to_csv("../Results/" + data_name + "/fitInfos/quad_info.csv", sep=",", index=False)
        cubic.to_csv("../Results/" + data_name + "/fitInfos/cubic_info.csv", sep=",", index=False)
        Briere.to_csv("../Results/" + data_name + "/fitInfos/Briere_info.csv", sep=",", index=False)
        Schoolfield.to_csv("../Results/" + data_name + "/fitInfos/Schoolfield_info.csv", sep=",", index=False)
        # Arrhenius_before.to_csv("../Results/" + data_name + "/fitInfos/Arrhenius_before_info.csv")
        # Arrhenius_after.to_csv("../Results/" + data_name + "/fitInfos/Arrhenius_after_info.csv")

    except FileNotFoundError:
        print(os.system("mkdir ../Results/"+ data_name +"/fitInfos"))
        quadratic.to_csv("../Results/"+ data_name +"/fitInfos/quad_info.csv", sep=",", index=False)
        cubic.to_csv("../Results/" + data_name +"/fitInfos/cubic_info.csv", sep=",", index=False)
        Briere.to_csv("../Results/" + data_name +"/fitInfos/Briere_info.csv", sep=",", index=False)
        Schoolfield.to_csv("../Results/" + data_name +"/fitInfos/Schoolfield_info.csv", sep=",", index=False)
        # Arrhenius_before.to_csv("../Results/" + data_name + "/fitInfos/Arrhenius_before_info.csv")
        # Arrhenius_after.to_csv("../Results/" + data_name + "/fitInfos/Arrhenius_after_info.csv")

    return "Model Fitting for " + data_name + " Rate Finished!"

def main(argv):
    ''' The main function '''
    global k
    k = 8.617e-5
    
    d_photosyn, d_respiration = data_cleaning()
    model_fitting(d_photosyn)
    model_fitting(d_respiration)
    
    return 0

if __name__ == "__main__":
    status = main(sys.argv)
    sys.exit(status)

