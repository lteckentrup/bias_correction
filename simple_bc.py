import xarray as xr
import matplotlib.pyplot as plt

var='prec'
model='CanESM5'

def readin(var, model):
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

    if var == 'prec':
        model_hist['prec'] = model_hist['prec']*86400
        model_COR['prec'] = model_COR['prec']*86400
    else:
        pass

    return(observed, model_hist, model_COR)

observed, model_hist, model_COR = readin('prec', 'CanESM5')

def method(var, model, method):

    if method == 'delta':
        clim_hist = model_hist.groupby('time.day').mean(dim='time')
        clim_cor = model_COR.groupby('time.day').mean(dim='time')

        diff_delta = clim_cor-clim_hist
        bc = observed.groupby('time.day') + diff_delta
    elif method == 'scaling_add':
        clim_obs = observed.groupby('time.day').mean(dim='time')
        clim_hist = model_hist.groupby('time.day').mean(dim='time')

        diff_scaling = clim_hist-clim_obs
        bc = model_COR.groupby('time.day') - diff_scaling
    elif method == 'scaling_multi':
        clim_obs = observed.groupby('time.day').mean(dim='time')
        clim_hist = model_hist.groupby('time.day').mean(dim='time')

        quotient_scaling = clim_obs/clim_hist
        bc = model_COR.groupby('time.day') * quotient_scaling
    elif method == 'mva':
        clim_obs = observed.groupby('time.day').mean(dim='time')
        clim_hist = model_hist.groupby('time.day').mean(dim='time')

        std_obs = observed.groupby('time.day').std(dim='time')
        std_cor = model_COR.groupby('time.day').std(dim='time')

        quotient_std = std_obs/std_cor
        prel = clim_obs + quotient_std - clim_hist
        bc = model_COR.groupby('time.day')+prel

    bc['prec'].attrs={'long_name':'Precipitation',
                      'standard_name':'precipitation_amount',
                      'units':'kg m-2'}

    bc['lat'].attrs={'units':'degrees',
                     'long_name':'latitude',
                     'standard_name':'latitude',
                     'axis':'Y'}
    bc['lon'].attrs={'units':'degrees',
                     'long_name':'longitude',
                     'standard_name':'longitude',
                     'axis':'X'}

    bc.attrs={'Conventions':'CF-1.6',
              'Model':model+' CMIP6',
              'Experiment':'SSP245',
              'Realisation':'r1i1p1f1'}

    bc.to_netcdf(method+'_'+var+'_'+model+'_cor1.nc',
                 encoding={'time':{'dtype': 'double'},
                           'lat':{'dtype': 'double'},
                           'lon':{'dtype': 'double'},
                           var:{'dtype': 'float32'}
                           }
                 )

# method('prec', 'CanESM5', 'delta')
# method('prec', 'CanESM5', 'scaling_add')
# method('prec', 'CanESM5', 'scaling_multi')
method('prec', 'CanESM5', 'mva')
