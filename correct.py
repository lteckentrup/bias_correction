import numpy as np
import xarray as xr
import pandas as pd
import matplotlib.pyplot as plt
from dpm import delta
from dpm import scaling_add
from dpm import scaling_multi
from dpm import mva
from dpm import dpm
from eqm import eqm
from qdm import qdm
import argparse

def Bias_Correction(var, model, method, method_long):
    observed = ('../../CRUJRA/'+var+'/crujra.v2.0.'+var+'.std-ordering.nc')
    prcp_hist = ('../../../australia_climate/'+var+'/'+var+'_'+model+
                 '_SSP245_r1i1p1f1_K_1850_2100.nc')
    prcp_COR = ('../../../australia_climate/'+var+'/'+var+'_'+model+
                '_SSP245_r1i1p1f1_K_1850_2100.nc')

    observed = xr.open_dataset(observed)
    model_hist = xr.open_dataset(prcp_hist)
    model_COR = xr.open_dataset(prcp_COR)

    observed = observed.sel(time = slice('1951-01-01','2015-12-31'))
    model_hist = model_hist.sel(time = slice('1951-01-01','2015-12-31'))
    model_COR = model_COR.sel(time = slice('1951-01-01','2015-12-31'))

    observed = observed[var]
    model_hist = model_hist[var]
    model_COR = model_COR[var]

    lats = observed.lat.values
    lons = observed.lon.values

    bias_corrected_results_hist = np.zeros([len(model_hist.time.values),
                                            len(model_hist.lat.values),
                                            len(model_hist.lon.values)])
    bias_corrected_results_hist[:] = np.nan

    bias_corrected_results_COR = np.zeros([len(model_COR.time.values),
                                           len(model_COR.lat.values),
                                           len(model_COR.lon.values)])
    bias_corrected_results_COR[:] = np.nan

    if var == 'prec':
        model_hist_values = model_hist.values*86400
    else:
        model_hist_values = model_hist.values
    hist_dict = {}
    hist_dict['time'] = model_hist.time.values
    hist_dict['lon'] =  model_hist.lon.values
    hist_dict['lat'] =  model_hist.lat.values

    if var == 'prec':
        modelCOR_values = model_COR.values*86400
    else:
        modelCOR_values = model_COR.values

    COR_dict = {}
    COR_dict['time'] = model_COR.time.values
    COR_dict['lon'] =  model_COR.lon.values
    COR_dict['lat'] =  model_COR.lat.values

    observation_attr_values = observed.values

    correct_params = []
    for i,lat in enumerate(lats):
        for j,lon in enumerate(lons):
            params_dict = {}
            if (np.isnan(model_hist_values[0,i,j]) or
                np.isnan(observation_attr_values[0,i,j])):
                bias_corrected_results_hist[:,i,j] = np.nan
                params_dict['lat'] = lat
                params_dict['lon'] = lon
                params_dict['params'] = np.nan
            else:
                try:
                    y = model_hist_values[:,i,j]
                    x = observation_attr_values[:,i,j]
                    y_COR = modelCOR_values[:,i,j]

                    if method == 'dpm':
                        temp_COR = dpm(var, x, y, y_COR, 0, None, True)
                    elif method == 'eqm':
                        temp_COR = eqm(var, x, y, y_COR, 0, None, 'constant')
                    elif method == 'qdm':
                        temp_COR = qdm(var, x, y, y_COR, 0, None)
                    elif method == 'delta':
                        temp_COR = delta(x, y, y_COR)
                    elif method == 'scaling_add':
                        temp_COR = scaling_add(x, y, y_COR)
                    elif method == 'scaling_multi':
                        temp_COR = scaling_multi(x, y, y_COR)
                    elif method == 'mva':
                        temp_COR = mva(x, y, y_COR)
                    bias_corrected_results_COR[:,i,j] = temp_COR

                    if i%5==0 and j%5==0:
                        print(lat,lon)

                    params_dict['lat'] = lat
                    params_dict['lon'] = lon
                    params_dict['params'] = temp_COR

                except:
                    bias_corrected_results_hist[:,i,j] = np.nan
                    bias_corrected_results_COR[:,i,j] = np.nan

                    params_dict['lat'] = lat
                    params_dict['lon'] = lon
                    params_dict['params'] = np.nan

            correct_params.append(params_dict)

    ds_COR = xr.Dataset({var:(('time', 'lat','lon'),
                                 bias_corrected_results_COR)},
                           coords={'lat': lats,
                                   'lon': lons,
                                   'time':COR_dict['time']})

    ds_COR['lat'].attrs={'units':'degrees', 'long_name':'latitude'}
    ds_COR['lon'].attrs={'units':'degrees', 'long_name':'longitude'}
    ds_COR.attrs={'Conventions':'CF-1.6',
                  'Model':model+' CMIP6',
                  'Experiment':'SSP245',
                  'Realisation':'r1i1p1f1',
                  'Bias correction': method_long,
                  'Date_Created':str(date_created)}
    ds_COR.to_netcdf(method+'_'+var+'_'+model+'_cor.nc',
                        encoding={'time':{'dtype': 'double'},
                                  'lat':{'dtype': 'double'},
                                  'lon':{'dtype': 'double'},
                                  var:{'dtype': 'float32'}
                                  }
                        )

parser = argparse.ArgumentParser(description='Pass variable, model and method for bias correction')
parser.add_argument('argument_names', nargs='+', help='nada')
args = parser.parse_args()

Bias_Correction(args.argument_names[0], args.argument_names[1],
                args.argument_names[2], args.argument_names[3])
