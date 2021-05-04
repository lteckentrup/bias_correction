import matplotlib.pylab as plt
import pandas as pd
import netCDF4 as nc
import numpy as np
import glob
import xarray as xr
from matplotlib.pyplot import cm

def readin_reanalysis(var, model):
    observed = ('../../../'+model+'/'+var+'/crujra.v2.0.'+var+'.std-ordering.nc')
    obs = xr.open_dataset(observed)
    obs = obs.sel(time = slice('1989-01-01','2018-12-31'))
    if var == 'temp':
        obs[var] = obs[var]-273.15

    obs_monthly = obs.groupby('time.month').mean('time')
    obs_annual = obs.groupby('time.year').mean('time')

    obs_daily_array = obs[var].values
    obs_monthly_array = obs_monthly[var].values
    obs_annual_array = obs_annual[var].values

    df_daily = pd.DataFrame()
    df_monthly = pd.DataFrame()
    df_annual = pd.DataFrame()

    df_daily[model] = obs_daily_array.flatten()
    df_monthly[model] = obs_monthly_array.flatten()
    df_annual[model] = obs_annual_array.flatten()

    df_daily = df_daily.dropna()
    df_monthly = df_monthly.dropna()
    df_annual = df_annual.dropna()

    return(df_daily, df_monthly, df_annual)

def readin_cmip(var, model, method):
    if method == 'original':
        data = ('../../../../australia_climate/'+var+'/'+var+'_'+model+
                 '_SSP245_r1i1p1f1_K_1850_2100.nc')
    else:
        data = ('../'+method+'/'+method+'_'+var+'_'+model+'_cor.nc')

    data_daily = xr.open_dataset(data)
    data_daily = data_daily.sel(time = slice('1989-01-01','2018-12-31'))

    if method == 'original':
        if var == 'prec':
            data_daily[var] = data_daily[var]*86400
        elif var == 'temp':
            data_daily[var] = data_daily[var]-273.15
            pass
    else:
        pass

    if var == 'prec':
        data_monthly = data_daily.resample(time="1MS").sum()
        data_annual = data_daily.groupby('time.year').sum('time')
    else:
        data_monthly = data_daily.resample(time="1MS").mean()
        data_annual = data_daily.groupby('time.year').mean('time')

    data_daily_array = data_daily[var].values
    data_monthly_array = data_monthly[var].values
    data_annual_array = data_annual[var].values

    df_daily = pd.DataFrame()
    df_monthly = pd.DataFrame()
    df_annual = pd.DataFrame()

    df_daily[model] = data_daily_array.flatten()
    df_monthly[model] = data_monthly_array.flatten()
    df_annual[model] = data_annual_array.flatten()

    df_daily = df_daily.dropna()
    df_monthly = df_monthly.dropna()
    df_annual = df_annual.dropna()

    return(df_daily, df_monthly, df_annual)

fig = plt.figure(0,figsize=(8.2,10.2))

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

ax1 = fig.add_subplot(6,3,1)
ax2 = fig.add_subplot(6,3,2)
ax3 = fig.add_subplot(6,3,3)
ax4 = fig.add_subplot(6,3,4)
ax5 = fig.add_subplot(6,3,5)
ax6 = fig.add_subplot(6,3,6)
ax7 = fig.add_subplot(6,3,7)
ax8 = fig.add_subplot(6,3,8)
ax9 = fig.add_subplot(6,3,9)
ax10 = fig.add_subplot(6,3,10)
ax11 = fig.add_subplot(6,3,11)
ax12 = fig.add_subplot(6,3,12)
ax13 = fig.add_subplot(6,3,13)
ax14 = fig.add_subplot(6,3,14)
ax15 = fig.add_subplot(6,3,15)
ax16 = fig.add_subplot(6,3,16)
ax17 = fig.add_subplot(6,3,17)
ax18 = fig.add_subplot(6,3,18)

cmip6_prec=glob.glob('../original/*_prec_timmean.nc', recursive=True)
cmip6_prec = sorted(cmip6_prec)
names_prel = [w.replace('_prec_timmean.nc', '') for w in cmip6_prec]
model_names = [w.replace('../original/', '') for w in names_prel]

color=cm.tab20(np.linspace(0,1,16))
for mn, c in zip(model_names,color):
    df_daily, df_monthly, df_annual = readin_cmip('temp', mn, 'original')
    ax1 = df_daily[m].plot.kde(ax=ax1, lw=2.0, label='_Hidden')
    ax2 = df_monthly[mn].plot.kde(ax=ax2, lw=2.0, color=c, label='_Hidden')
    ax3 = df_annual[mn].plot.kde(ax=ax3, lw=2.0, color=c, label='_Hidden')

    df_daily, df_monthly, df_annual = readin_cmip('temp', mn, 'scaling_additive')
    ax4 = df_daily[m].plot.kde(ax=ax4, lw=2.0, label='_Hidden')
    ax5 = df_monthly[mn].plot.kde(ax=ax5, lw=2.0, color=c, label='_Hidden')
    ax6 = df_annual[mn].plot.kde(ax=ax6, lw=2.0, color=c, label='_Hidden')

    df_daily, df_monthly, df_annual = readin_cmip('temp', mn, 'mva')
    ax7 = df_daily[m].plot.kde(ax=ax7, lw=2.0, label='_Hidden')
    ax8 = df_monthly[mn].plot.kde(ax=ax8, lw=2.0, color=c, label='_Hidden')
    ax9 = df_annual[mn].plot.kde(ax=ax9, lw=2.0, color=c, label='_Hidden')

    df_daily, df_monthly, df_annual = readin_cmip('temp', mn, 'dpm')
    ax10 = df_daily[m].plot.kde(ax=ax10, lw=2.0, label='_Hidden')
    ax11 = df_monthly[mn].plot.kde(ax=ax11, lw=2.0, color=c, label=mn)
    ax12 = df_annual[mn].plot.kde(ax=ax12, lw=2.0, color=c, label='_Hidden')

    df_daily, df_monthly, df_annual = readin_cmip('temp', mn, 'eqm')
    ax13 = df_daily[m].plot.kde(ax=ax13, lw=2.0, label='_Hidden')
    ax14 = df_monthly[mn].plot.kde(ax=ax14, lw=2.0, color=c, label='_Hidden')
    ax15 = df_annual[mn].plot.kde(ax=ax15, lw=2.0, color=c, label='_Hidden')

    df_daily, df_monthly, df_annual = readin_cmip('temp', mn, 'qdm')
    ax16 = df_daily[m].plot.kde(ax=ax16, lw=2.0, label='_Hidden')
    ax17 = df_monthly[mn].plot.kde(ax=ax17, lw=2.0, color=c, label='_Hidden')
    ax18 = df_annual[mn].plot.kde(ax=ax18, lw=2.0, color=c, label='_Hidden')

df_daily, df_monthly, df_annual = readin_reanalysis('temp', 'CRUJRA')

ax1 = df_daily['CRUJRA'].plot.kde(ax=ax1, lw=4.0, color='k', label='CRUJRA')
ax2 = df_monthly['CRUJRA'].plot.kde(ax=ax2, lw=4.0, color='k', label='_Hidden')
ax3 = df_annual['CRUJRA'].plot.kde(ax=ax3, lw=4.0, color='k', label='_Hidden')

# ax4 = df_daily['CRUJRA'].plot.kde(ax=ax4, lw=4.0, color='k', label='CRUJRA')
# ax5 = df_monthly['CRUJRA'].plot.kde(ax=ax5, lw=4.0, color='k', label='CRUJRA')
# ax6 = df_annual['CRUJRA'].plot.kde(ax=ax6, lw=4.0, color='k', label='CRUJRA')
#
# ax7 = df_daily['CRUJRA'].plot.kde(ax=ax7, lw=4.0, color='k', label='CRUJRA')
# ax8 = df_monthly['CRUJRA'].plot.kde(ax=ax8, lw=4.0, color='k', label='CRUJRA')
# ax9 = df_annual['CRUJRA'].plot.kde(ax=ax9, lw=4.0, color='k', label='CRUJRA')
#
# ax10 = df_daily['CRUJRA'].plot.kde(ax=ax10, lw=4.0, color='k', label='CRUJRA')
# ax11 = df_monthly['CRUJRA'].plot.kde(ax=ax11, lw=4.0, color='k', label='CRUJRA')
# ax12 = df_annual['CRUJRA'].plot.kde(ax=ax12, lw=4.0, color='k', label='CRUJRA')
#
# ax13 = df_daily['CRUJRA'].plot.kde(ax=ax13, lw=4.0, color='k', label='CRUJRA')
# ax14 = df_monthly['CRUJRA'].plot.kde(ax=ax14, lw=4.0, color='k', label='CRUJRA')
# ax15 = df_annual['CRUJRA'].plot.kde(ax=ax15, lw=4.0, color='k', label='CRUJRA')
#
# ax16 = df_daily['CRUJRA'].plot.kde(ax=ax16, lw=4.0, color='k', label='CRUJRA')
# ax17 = df_monthly['CRUJRA'].plot.kde(ax=ax17, lw=4.0, color='k', label='CRUJRA')
# ax18 = df_annual['CRUJRA'].plot.kde(ax=ax18, lw=4.0, color='k', label='CRUJRA')

ax1.set_title('Annual PDF')
ax2.set_title('Monthly PDF')
ax3.set_title('Daily PDF')

ax11.legend(loc='upper center', bbox_to_anchor=(0.4, -2.8), ncol=4, frameon = True)
for a in (ax2,ax3,ax5,ax6,ax8,ax9,ax11,ax12,ax14,ax15,ax17,ax18):
    a.set_ylabel('')

for a in (ax16,ax17,ax18):
    a.set_xlabel('Temperature')

for a in (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15):
    a.set_xticklabels([])

for a in (ax1,ax2,ax3,ax4,ax5,ax6,ax7,ax8,ax9,ax10,ax11,ax12,ax13,ax14,ax15,ax16,ax17,ax18):
    a.spines['right'].set_visible(False)
    a.spines['top'].set_visible(False)

ax1.set_ylabel('Density \n (not corrected)')
ax4.set_ylabel('Density \n (scaling add.)')
ax7.set_ylabel('Density \n (MVA)')
ax10.set_ylabel('Density \n (QDM)')
ax13.set_ylabel('Density \n (EQM)')
ax16.set_ylabel('Density \n (QDM)')

plt.show()
