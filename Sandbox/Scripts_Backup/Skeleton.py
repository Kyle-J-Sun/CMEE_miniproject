#!/usr/bin/env python3

"""Perform model fitting with NLLS"""
__appname__ = ...
__author__ = ...
__version__ = '0.0.1'

#Imports
import sys
import pandas as pd
import numpy as np
from lmfit import minimize, Parameters

# Define Constants here

k = ...
#etc...

#Set seed here(for debugging)

...


#Import data and select columns to be used

data = ... 

#### Define models####

#Briere model

def Briere_model(id, data, min_round, max_round):
    """ Perform NLLS model fitting on TPC with Briere model"""
    #Variabes
    T = np.asarray(data.ConTemp[data.id == id]) #Temperature
    B = np.asarray(data.OriginalTraitValue[data.id == id]) #Trait value
    #Starting values
    B0 = 0.01
    T0 = min(T)
    Tm = max(T)

    #Calculate and minimize the residuals

    def Briere_residual(p, T, B):
        """ Residuals of data and models """
        B0 = p["B0"].value
        T0 = p["T0"].value 
        Tm = p["Tm"].value

        model = B0 * T * (T-T0) * (Tm - T)**0.5

        return model - B
    
    #Initiate an empty dictionary for results
    all = {"id":[id], 
            "B0":[B0],
            "T0":[T0],
            "Tm":[Tm],
            "chisqr": [np.NAN],
            "RSS": [np.NAN],
            "TSS": [np.NAN],
            "Rsq": [np.NAN],
            "aic":[np.NAN],
            "bic":[np.NAN],
            "aicc": [np.NAN],
            }
    
    round = 0
    while round < max_round:
        
        round += 1 #Try round +1
        #If TPC converged within mean rounds, or max round is exceeded, break the loop
        if all["aic"] != [np.NaN] and round > min_round: 
            break

        #Starting values defined for rounds
        #For the first round use setted starting values
        if round == 1 : 
            p = Parameters()
            p.add_many(("B0",B0),("T0",T0),("Tm",Tm))
        #For other rounds, assign B0 with random number distributed between 0 and B0*2
        else:
            p = Parameters() 
            p.add("B0", value = np.random.uniform(0,1), min = 0)
            p.add("T0", value = T0, min = -30, max = 130)
            p.add("Tm", value = Tm, min = -30, max = 130)

        try:   
            output = minimize(Briere_residual, p, args = (T,B))
            RSS = sum(Briere_residual(output.params, T, B)**2)
            TSS = sum((B - np.mean(B))**2)
            Rsq = 1 - (RSS/TSS)
            aicc = output.aic + 2*3*((3+1)/(len(T) - 3 - 1))
            #Replace the fitting results if aicc is smaller, 
            # or if the previous fit didn't converge
            if aicc < all["aicc"] or all["aicc"] == [np.NaN]:
                all = {"id":[id], 
                    "B0":[output.params["B0"].value],
                    "T0":[output.params["T0"].value],
                    "Tm":[output.params["Tm"].value],
                    "chisqr": [output.chisqr],
                    "RSS": RSS,
                    "TSS": TSS,
                    "Rsq": Rsq,
                    "aic":[output.aic],
                    "bic":[output.bic],
                    "aicc": aicc
                    }
    
        except ValueError:
            pass
        
        continue

    return pd.DataFrame(all)


#Schoolfield model simplified, appropriate for without high temperature data

def Schoolfield_model_nh(id, data, min_round, max_round):
    """ Perform NLLS model fitting on unlogged trait value with simplified schoolfield model omitting high temperature values """
    #Variabes
    T = np.asarray(data.ConTemp_F[data.id == id]) #Temperature in kelvin
    B = np.asarray(data.OriginalTraitValue[data.id == id]) #Unlogged Trait value
    B_log = np.asarray(data.OriginalTraitValue_log[data.id == id]) #Logged Trait value
    #Starting values(matters on converge)
    B0 = data.B0[data.id == id].iloc[0]
    E = abs(data.E[data.id == id].iloc[0])
    Tl = data.Tl[data.id == id].iloc[0]
    El = abs(data.El[data.id == id].iloc[0])
    #Add starting values
    p = Parameters()
    p.add("B0", value = B0, min = 0)
    p.add("E", value = E, min = 0)
    p.add("Tl",Tl, min = 250, max = 400)
    p.add("El",El, min = 0)

    def Schoolfield_residual_nh(p, T, B):
        """ Return residuals for unlogged trait values schoolfield model """

        B0 = p['B0'].value
        E  = p['E'].value
        Tl = p['Tl'].value
        El = p['El'].value

        model = (B0*e**((-E/k)*((1/T)-(1/283.15))))/(
                1+(e**((El/k)*((1/Tl)-(1/T)))))

        return model - B
    
    def Schoolfield_residual_log_nh(p, T, B_log):
        """ Return residuals for logged trait values and logged schoolfield model """

        B0 = p['B0'].value
        E  = p['E'].value
        Tl = p['Tl'].value
        El = p['El'].value

        model = (B0*e**((-E/k)*((1/T)-(1/283.15))))/(
            1+(e**((El/k)*((1/Tl)-(1/T)))))

        return np.log(model) - B_log

    #Initiate an empty dictionary for results
    all = {"id":[id], 
            "B0":[B0],
            "E":[E],
            "Tl":[Tl],
            "El": [El],
            "chisqr":[np.NAN],
            "RSS": [np.NAN],
            "TSS": [np.NAN],
            "Rsq": [np.NAN],
            "aic":[np.NAN],
            "bic":[np.NAN],
            "aicc":[np.NAN],
            "RSS_log": [np.NAN],
            "TSS_log": [np.NAN],
            "Rsq_log": [np.NAN],
            "aic_log": [np.NAN],
            "bic_log": [np.NAN],
            "aicc_log": [np.NAN]
            }
    
    round = 0
    while round < max_round:
    
        round += 1 #Try round +1
        #If TPC converged within mean rounds, or max round is exceeded, break the loop
        if all["aic"] != [np.NaN] and round > min_round:  
            break

        #Starting values defined for rounds
        #For the first round use set starting values
        if round == 1 : 
            #Add starting values
            p = Parameters()
            p.add("B0", value = B0, min = 0)
            p.add("E", value = E, min = 0)
            p.add("Tl",Tl, min = 250, max = 400)
            p.add("El",El, min = 0)
        #For other rounds, assign B0 with random number distributed between 0 and B0*2
        else:
            p = Parameters() 
            p.add("B0", value = np.random.uniform(0,B0*2), min = 0)
            p.add("E", value = np.random.uniform(0,E*2), min = 0)
            p.add("Tl",Tl, min = 250, max = 400)
            p.add("El", value = np.random.uniform(0, El*2),min = 0)

        #Results for model fitting on logged values
        try:    
            output = minimize(Schoolfield_residual_log_nh, p, args = (T, B_log))
            RSS_log = sum(Schoolfield_residual_log_nh(output.params,T,B_log)**2)
            TSS_log = sum((B_log - np.mean(B_log))**2)
            Rsq_log = 1 - (RSS_log/TSS_log)
            aicc_log = output.aic + 2*4*((4+1)/(len(T) - 4 - 1))
            #Calculate results for unlogged fitting based on parameters estimated from 
            #logged fitting, since they stay the same after log transformation.
            RSS = sum(Schoolfield_residual_nh(output.params,T,B)**2)
            TSS = sum((B - np.mean(B))**2)
            Rsq = 1 - (RSS/TSS)
            aic = len(T)*np.log(RSS/len(T)) + 2*4
            bic = len(T)*np.log(RSS/len(T)) + np.log(len(T))*4
            aicc = aic + 2*4*((4+1)/(len(T) - 4 - 1))

            if aicc < all["aicc"] or all["aicc"] == [np.NaN]:           
                all = {"id":[id], 
                    "B0":[output.params["B0"].value],
                    "E":[output.params["E"].value], 
                    "Tl":[output.params["Tl"].value],
                    "El": [output.params["El"].value],
                    "chisqr": [output.chisqr],
                    "RSS": RSS, 
                    "TSS": TSS, 
                    "Rsq": Rsq,
                    "aic":aic,
                    "bic":bic,
                    "aicc":aicc,
                    "RSS_log": RSS_log,
                    "TSS_log": TSS_log,
                    "Rsq_log": Rsq_log,
                    "aic_log": [output.aic],
                    "bic_log": [output.bic],
                    "aicc_log": aicc_log
                    }
        
        except ValueError:
            pass
        continue


    return(pd.DataFrame(all))


#Initialise empty dataframes for fitting results

Briere = pd.DataFrame(data = None)
Schoolfield_nh = pd.DataFrame(data = None)

### NLLS model fitting ######

#Briere model
i = 0
print("Briere model fitting")
for id in data["id"].unique():
    Briere = Briere.append(Briere_model(id, data,3, 25))
    ... etc

#Schoolfield model without high temperature parameters
i = 0
print("Schoolfield model without high temperature parameters fitting")
for id in data["id"].unique():
    Schoolfield_nh = Schoolfield_nh.append(Schoolfield_model_nh(id, data,3, 25))
    ... etc
