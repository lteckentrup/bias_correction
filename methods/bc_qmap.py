import pickle
import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
import sys
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
rpy2.robjects.numpy2ri.activate()
qmap = importr('qmap')
from rpy2.robjects import pandas2ri
from rpy2.robjects import r
import glob
sys.setrecursionlimit(10000)
import traceback

def bias_correction(x,y):
    q_map = qmap.fitQmap(x, y,method="RQUANT", qstep=0.01, wett_day=False)
    qm1 = qmap.doQmap(y, q_map)
    bias_corrected_output = {}
    bias_corrected_output['params'] = q_map
    bias_corrected_output['outputs'] = qm1
    return bias_corrected_output

def bias_correction_model(y,q_map):
    qm1 = qmap.doQmap(y, q_map)
    bias_corrected_output = {}
    bias_corrected_output['outputs'] = qm1
    return bias_corrected_output

observed = 'temp_CRUJRA_1951-1985.nc'
prcp_hist = 'temp_ACCESS-ESM1-5_1951-1985.nc'
prcp_LIG = 'temp_LIG.nc'

observed = xr.open_dataset(observed)
model_hist = xr.open_dataset(prcp_hist)
model_LIG = xr.open_dataset(prcp_LIG)

observed = observed['temp']
model_hist = model_hist['temp']
model_LIG = model_LIG['temp']

lats = observed.lat.values
lons = observed.lon.values

bias_corrected_results_hist = np.zeros([len(model_hist.time.values),
                                        len(model_hist.lat.values),
                                        len(model_hist.lon.values)])
bias_corrected_results_hist[:] = np.nan

bias_corrected_results_LIG = np.zeros([len(model_LIG.time.values),
                                       len(model_LIG.lat.values),
                                       len(model_LIG.lon.values)])
bias_corrected_results_LIG[:] = np.nan

model_hist_values = model_hist.values
hist_dict = {}
hist_dict['time'] = model_hist.time.values
hist_dict['lon'] =  model_hist.lon.values
hist_dict['lat'] =  model_hist.lat.values

modelLIG_values = model_LIG.values
LIG_dict = {}
LIG_dict['time'] = model_LIG.time.values
LIG_dict['lon'] =  model_LIG.lon.values
LIG_dict['lat'] =  model_LIG.lat.values

observation_attr_values = observed.values

correct_params = []
for i,lat in enumerate(lats):
    for j,lon in enumerate(lons):
        params_dict = {}
        if np.isnan(model_hist_values[0,i,j]) or np.isnan(observation_attr_values[0,i,j]):
            bias_corrected_results_hist[:,i,j] = np.nan
            params_dict['lat'] = lat
            params_dict['lon'] = lon
            params_dict['params'] = np.nan
        else:
            try:
                y = model_hist_values[:,i,j]
                x = observation_attr_values[:,i,j]

                y_LIG = modelLIG_values[:,i,j]

                temp = bias_correction(x,y)
                q_map = temp['params']
                temp_LIG = bias_correction_model(y_LIG,q_map)

                bias_corrected_results_LIG[:,i,j] = temp_LIG['outputs']
                bias_corrected_results_hist[:,i,j] = temp['outputs']

                if i%5==0 and j%5==0:
                    print(lat,lon)
                params_dict['lat'] = lat
                params_dict['lon'] = lon
                params_dict['params'] = temp['params']

            except:
                bias_corrected_results_hist[:,i,j] = np.nan
                bias_corrected_results_LIG[:,i,j] = np.nan

                params_dict['lat'] = lat
                params_dict['lon'] = lon
                params_dict['params'] = np.nan

        correct_params.append(params_dict)

ds_hist = xr.Dataset({'temp': (('time', 'lat','lon'),
                               bias_corrected_results_hist)},
                     coords={'lat': lats,
                             'lon': lons,
                             'time':hist_dict['time'] })
ds_sspLIG = xr.Dataset({'temp': (('time', 'lat','lon'),
                                 bias_corrected_results_LIG)},
                       coords={'lat': lats,
                               'lon': lons,
                               'time':LIG_dict['time'] })

ds_hist.to_netcdf('hist_temp_cor.nc')
ds_sspLIG.to_netcdf('LIG_temp_cor.nc')
