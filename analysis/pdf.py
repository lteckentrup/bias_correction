import matplotlib.pylab as plt
import pandas as pd
import netCDF4 as nc
import numpy as np
import glob
import xarray as xr

def readin_reanalysis(var, model):
    observed = ('../../'+model+'/'+var+'/crujra.v2.0.'+var+'.std-ordering.nc')
    obs = xr.open_dataset(observed)
    obs = obs.sel(time = slice('1989-01-01','2018-12-31'))
    obs_monthly = obs.groupby('time.month').mean('time')
    obs_annual = obs.groupby('time.year').mean('time')

def readin_cmip(var, model, method):
    if method == 'original':
        data = ('../../../../australia_climate/'+var+'/'+var+'_'+model+
                 '_SSP245_r1i1p1f1_K_1850_2100.nc')
    else:
        data = ('../'+method+'/'+method+'_'+var+'_'+model+'_cor.nc')

    data_daily = xr.open_dataset(data)
    data_daily = data_daily.sel(time = slice('1989-01-01','2018-12-31'))

    if var == 'prec':
        if method == 'original':
            data_daily['prec'] = data_daily['prec']*86400
        else:
            pass
    else:
        pass

    if var == 'prec':
        data_monthly = data_daily.resample(time="1MS").sum()
        data_annual = data_daily.groupby('time.year').sum('time')
    else:
        data_monthly = data_daily.resample(time="1MS").mean()
        data_annual = data_daily.groupby('time.year').mean('time')

    # data_daily_array = data_daily[var].values
    data_monthly_array = data_monthly[var].values
    data_annual_array = data_annual[var].values

    df_daily = pd.DataFrame()
    df_monthly = pd.DataFrame()
    df_annual = pd.DataFrame()

    # df_daily[model] = data_daily_array.flatten()
    df_monthly[model] = data_monthly_array.flatten()
    df_annual[model] = data_annual_array.flatten()

    df_monthly = df_monthly.dropna()
    df_annual = df_annual.dropna()

    return(df_monthly, df_annual)

fig = plt.figure(0,figsize=(10.2,10.2))

fig.subplots_adjust(hspace=0.2)
fig.subplots_adjust(wspace=0.36)
fig.subplots_adjust(right=0.94)
fig.subplots_adjust(left=0.09)
fig.subplots_adjust(bottom=0.1)
fig.subplots_adjust(top=0.95)

plt.rcParams['text.usetex'] = False
plt.rcParams['axes.labelsize'] = 12
plt.rcParams['font.size'] = 11
plt.rcParams['legend.fontsize'] = 12
plt.rcParams['xtick.labelsize'] = 11
plt.rcParams['ytick.labelsize'] = 11

ax1 = fig.add_subplot(3,3,1)
ax2 = fig.add_subplot(3,3,2)
ax3 = fig.add_subplot(3,3,3)
ax4 = fig.add_subplot(3,3,4)
ax5 = fig.add_subplot(3,3,5)
ax6 = fig.add_subplot(3,3,6)
ax7 = fig.add_subplot(3,3,7)
ax8 = fig.add_subplot(3,3,8)
ax9 = fig.add_subplot(3,3,9)

cmip6_prec=glob.glob('../../../../australia_climate/prec/prec_*', recursive=True)
cmip6_prec = sorted(cmip6_prec)
prel = [w.replace('_SSP245_r1i1p1f1_K_1850_2100.nc', '') for w in cmip6_prec]
prel1 = [w.replace('../../../../australia_climate/prec/prec_', '') for w in prel]
model_names = [w.replace('ensmean.nc', '') for w in prel1]

for mn in model_names:
    df_monthly, df_annual = readin_cmip('temp', mn, 'original')
    # ax1 = df_daily['CanESM5'].plot.kde(ax=ax1, lw=2.0, label='CanESM5')
    ax2 = df_monthly[mn].plot.kde(ax=ax2, lw=2.0, label=mn)
    ax3 = df_annual[mn].plot.kde(ax=ax3, lw=2.0, label=mn)

plt.show()
