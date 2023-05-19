
# Importing the desired modules
import numpy as np
import matplotlib.pyplot as plt
import pickle
import pandas as pd
from copy import deepcopy
#import astropy.io.ascii as aioascii
import sys
from astropy import constants as const
from math import atan2
import ttvfast
from ttvfast.models import Planet
from shutil import copyfile
from TRAPPIST1_parameters_optimum import *
import pdb ## CP: au cas ou pour debug

#%%
def create_planet(planet_dict):
  # Define constants
  Mearth = const.M_earth
  Msun = const.M_sun
  pla = Planet(planet_dict["mass_earth"]*Mearth/Msun, planet_dict["period_days"],planet_dict["ecc"], planet_dict["inc_deg"], planet_dict["Omega_deg"], planet_dict["omega_deg"], planet_dict["M"])
  
  return pla

def get_pred_tts_T1(list_of_plas_dict = [], first_tts=[], Mstar_sun=0.1, t_start_integ_days=7240, t_end_integ_days=8540, dt=1e-2):
  npla = len(list_of_plas_dict)
  if 1:
    pla_list = []
    for i in range(npla):
      pla_list.append(create_planet(list_of_plas_dict[i]))
    # Create TTVFast object
  t_start_integ_days = t_start_integ_days
  timestep_integ_days = dt ## CP: modification pour plus grande precision
  t_end_integ_days = t_end_integ_days

    # calculate transit times
  rv_dict = ttvfast.ttvfast(pla_list,
                            Mstar_sun, t_start_integ_days,
                            timestep_integ_days, t_end_integ_days)

  planet = np.asarray(rv_dict['positions'][0])  # int indices for which pla
  times = np.asarray(rv_dict['positions'][2])  # transit times

    # isolate each planet's transit times

  pred_times_list = [] # list of arrays qui contiennent les tts de chaque planete
  for i in range(npla):
    pred_times_list.append(times[np.where((planet == i)*(times != -2))])

  # epochs of predicted transit times
  epochs_pred_list = [] # contain [b_epochs_pred, c_epochs_pred, d_epochs_pred...]
  for i in range(npla):
    thispla = list_of_plas_dict[i]
    epochs_pred = np.round((pred_times_list[i] - first_tts[i])/thispla["period_days"])
    epochs_pred_list.append(epochs_pred) 

  return (epochs_pred_list), (pred_times_list)

def calc_ttvs_T1(list_of_plas_dict, first_tts, Mstar_sun, 
                  t_start_integ_days, t_end_integ_days, dt=1e-2):

    
  # calculate transit times
  (epochs_pred_list), (pred_times_list) = get_pred_tts_T1(list_of_plas_dict = list_of_plas_dict, first_tts=first_tts, Mstar_sun=Mstar_sun, 
                                                          t_start_integ_days=t_start_integ_days, t_end_integ_days=t_end_integ_days, dt=dt)

  # compare predicted to observed transit times

  inMinutes = 24.*60.
  # Get TTVs by fitting to the predicted tts 

  # ajoute une for loop pour le faire pour chacune des npla

  npla = len(list_of_plas_dict)
  ttv_pred_list = []
  for i in range(npla):
      fit_para_pred = np.polyfit(epochs_pred_list[i], pred_times_list[i], 1)    
      fit_func = np.poly1d(fit_para_pred)  
      ttv_pred_pla = pred_times_list[i] - fit_func(epochs_pred_list[i]) # modifier pred_times, b_epochs_pred, fitb_pred
      ttv_pred_list.append(ttv_pred_pla*inMinutes)
      
  
  return epochs_pred_list, pred_times_list, ttv_pred_list

def plot_ttvs(epochs_pred_list, pred_times_list, ttv_pred_list, x_is_epoch = True):
    
    n_axes = len(epochs_pred_list)
    planets=["b", "c", "d", "e", "f", "g", "h"]
    fig, axes = plt.subplots(n_axes, 1, figsize=(8,10)) ## CP: change le figsize pour une plus grosse figure
    time0 = int(pred_times_list[0][0]*100)/100
    for i, ax in enumerate(axes):
        if x_is_epoch:
            ax.plot(epochs_pred_list[i], ttv_pred_list[i])
            ax.set_ylabel("TTV "+planets[i])
        else:
            ax.plot(pred_times_list[i]-time0, ttv_pred_list[i]) ## CP: modifie: j'avais une erreur pour l'endroit ou mettre le [i]
            ax.set_ylabel("TTV "+planets[i])
    
    ## CP: deplace pour seulement 1 titre x pour le dernier ax
    if x_is_epoch:
        ax.set_xlabel("Epoch")
    else:
        ax.set_xlabel("BJD -"+str(time0))    
        
        
    fig.tight_layout() ## CP: change pour un plus beau layout de la figure
    return fig, axes


#%%

# Inputs for the code:

# create list of plas dict
list_of_plas_dict = [pla_b, pla_c, pla_d, pla_e, pla_f, pla_g, pla_h]
# create list of first tts


# give as input the right Mstar_sun and t_start_integ_days

## CP: ajoute ici pour reproduire le "range" dans la figure de Agol+
t_start_integ_days = 7257.93115525 ## CP: reference for the start of the simulation (parameter values, Agol+2021)
t_end_integ_days = 8900 



###########################################################################################
#%% Calling the functions
print(first_tts)
epochs_pred_list, pred_times_list, ttv_pred_list = calc_ttvs_T1(list_of_plas_dict, first_tts, 
                                                                Mstar_sun, t_start_integ_days, t_end_integ_days, dt=1e-3)
fig, axes = plot_ttvs(epochs_pred_list, pred_times_list, ttv_pred_list, x_is_epoch = True)
plt.show()








