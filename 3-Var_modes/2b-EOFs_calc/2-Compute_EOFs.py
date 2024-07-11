# %% [markdown]
# # 2. Compute EOFs

# This script is used to compute the EOFs and PCs. 
# from monthly seasonal forescasts for the hindcast period.
# 
# EOFs are calculted for ERA5 data with the Eof library.
# Seasonal forecasts are then projected into these EOFs to obtain the PCs.
#
# First we have to decide a forecast system (institution and system name) and a start month. 

#%%
print("2. Compute EOFs") 

import os
import sys
from dotenv import load_dotenv
import xarray as xr
import pandas as pd
import numpy as np
from eofs.xarray import Eof
import xskillscore as xs
import locale
import calendar
from dateutil.relativedelta import relativedelta
import time
import datetime as dt
# Libraries for plotting and geospatial data visualisation
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import warnings
warnings.filterwarnings('ignore')

# If variables are introduced from the command line, read them
if len(sys.argv) > 2:
    institution = str(sys.argv[1]).replace('"', '')
    name = str(sys.argv[2]).replace('"', '')
    startmonth = int(sys.argv[3])
# If no variables were introduced, ask for them
else:
    # Which model institution
    institution = input("Usar modelo del siguiente organismo [ ECMWF , Météo France , Met Office , DWD , CMCC , NCEP , JMA , ECCC ]: ")

    # Which model system
    if institution=='ECMWF':
        name = input("Sistema del modelo [ System 4 , SEAS5 , SEAS5.1 ]: ")
    elif institution=='Météo France':
        name = input("Sistema del modelo [ System 5 , System 6 , System 7 , System 8 ]: ")
    elif institution=='Met Office':
        name = input("Sistema del modelo [ System 12 , System 13 , System 14 , System 15 , GloSea6 , GloSea6.1 , GloSea6.2 , GloSea6.3 ]: ")
    elif institution=='DWD':
        name = input("Sistema del modelo [ GCFS2.0 , GCFS2.1 ]: ")
    elif institution=='CMCC':
        name = input("Sistema del modelo [ SPSv3.0 , SPSv3.5 ]: ")
    elif institution=='NCEP':
        name = input("Sistema del modelo [ CFSv2 ]: ")
    elif institution=='JMA':
        name = input("Sistema del modelo [ CPS2 , CPS3 ]: ")
    elif institution=='ECCC':
        name = input("Sistema del modelo [ GEM-NEMO , CanCM4i , GEM5-NEMO ]: ")
    else:
        sys.exit()

    # Which start month
    startmonth = int(input("Mes de inicialización (en número): "))

# Dictionary to link full system names and simplier names
full_name = {'ECMWF-System 4': ['ecmwf','4'],
             'ECMWF-SEAS5': ['ecmwf', '5'],
             'ECMWF-SEAS5.1': ['ecmwf', '51'],
             'Météo France-System 5': [ 'meteo_france', '5'],
             'Météo France-System 6': [ 'meteo_france', '6'],
             'Météo France-System 7': [ 'meteo_france', '7'],
             'Météo France-System 8': [ 'meteo_france', '8'],
             'Met Office-System 12': ['ukmo', '12'],
             'Met Office-System 13': ['ukmo', '13'],
             'Met Office-System 14': ['ukmo', '14'],
             'Met Office-System 15': ['ukmo', '15'],
             'Met Office-GloSea6': ['ukmo', '600'],
             'Met Office-GloSea6.1': ['ukmo', '601'],
             'Met Office-GloSea6.2': ['ukmo', '602'],
             'Met Office-GloSea6.3': ['ukmo', '603'],
             'DWD-GCFS2.0': ['dwd', '2'],
             'DWD-GCFS2.1': ['dwd', '21'],
             'CMCC-SPSv3.0': ['cmcc', '3'],
             'CMCC-SPSv3.5': ['cmcc', '35'],
             'NCEP-CFSv2': ['ncep', '2'],
             'JMA-CPS2': ['jma', '2'],
             'JMA-CPS3': ['jma', '3'],
             'ECCC-GEM-NEMO': ['eccc', '1'],
             'ECCC-CanCM4i': ['eccc', '2'],
             'ECCC-GEM5-NEMO': ['eccc', '3']
            }
# Save the full model name in origin_labels
origin_labels = {'institution': institution, 'name': name}
# Save the simplier model and system name
model = full_name[institution+'-'+name][0]
system = full_name[institution+'-'+name][1]

# Here we save the configuration
config = dict(
    list_vars = 'geopotential',
    pressure_level = '500',
    hcstarty = 1993,
    hcendy = 2016,
    start_month = startmonth,
    origin = model,
    system = system,
    isLagged = False if model in ['ecmwf', 'meteo_france', 'dwd', 'cmcc', 'eccc'] else True
)

# Directory selection
load_dotenv() # Data is saved in a path defined in file .env
DATADIR = os.getenv('DATA_DIR')
MODESDIR = DATADIR + '/modes'
SCOREDIR = DATADIR + '/scores'
PLOTSDIR = f'./plots/stmonth{config["start_month"]:02d}'
PLOTSDIR_EOFpatterns = f'./plots/EOFs/patterns'
PLOTSDIR_EOFpcs = f'./plots/EOFs/pcs'
PLOTSDIR_EOFcorrelations = f'./plots/EOFs/correlations'
PLOTSDIR_EOFvariance = f'./plots/EOFs/variance'
# Directory creation
for directory in [PLOTSDIR_EOFpatterns, PLOTSDIR_EOFpcs, PLOTSDIR_EOFcorrelations, PLOTSDIR_EOFvariance]:
    # Check if the directory exists
    if not os.path.exists(directory):
        # If it doesn't exist, create it
        try:
            os.makedirs(directory)
        except FileExistsError:
            pass

# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)
# File name for hindcast
hcst_fname = f'{DATADIR}/{hcst_bname}.grib'
# Base name for observations
obs_bname = 'era5_monthly_stmonth{start_month:02d}_{hcstarty}-{hcendy}'.format(**config)
# File name for observations
obs_fname = f'{DATADIR}/{obs_bname}.grib'
# File name for climatologies
clim_fname = f'{DATADIR}/era5_monthly.grib'
clim_fname_sf = f'{DATADIR}/era5_monthly_sf.grib'

# Check if files exist
if not os.path.exists(obs_fname):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()
elif not os.path.exists(hcst_fname):
    print('No se descargaron aún los datos de este modelo y sistema')
    sys.exit()
elif not os.path.exists(clim_fname):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()
elif not os.path.exists(clim_fname_sf):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()

# %% [markdown]
# ## 2.1 Hindcast anomalies

# We calculate the monthly and 3-months anomalies for the hindcast data.

#%%
print("2.1 Hindcast anomalies")

# For the re-shaping of time coordinates in xarray.Dataset we need to select the right one 
#  -> burst mode ensembles (e.g. ECMWF SEAS5) use "time". This is the default option in this notebook
#  -> lagged start ensembles (e.g. MetOffice GloSea6) use "indexing_time" (see CDS documentation about nominal start date)
st_dim_name = 'time' if not config.get('isLagged',False) else 'indexing_time'

# Reading hindcast data from file
hcst = xr.open_dataset(hcst_fname,engine='cfgrib', backend_kwargs=dict(time_dims=('forecastMonth', st_dim_name)))
# We use dask.array with chunks on leadtime, latitude and longitude coordinate
hcst = hcst.chunk({'forecastMonth':1, 'latitude':'auto', 'longitude':'auto'})  #force dask.array using chunks on leadtime, latitude and longitude coordinate
# Reanme coordinates to match those of observations
hcst = hcst.rename({'latitude':'lat','longitude':'lon', st_dim_name:'start_date'})

# Add start_month to the xr.Dataset
start_month = pd.to_datetime(hcst.start_date.values[0]).month
hcst = hcst.assign_coords({'start_month':start_month})
# Add valid_time to the xr.Dataset
vt = xr.DataArray(dims=('start_date','forecastMonth'), coords={'forecastMonth':hcst.forecastMonth,'start_date':hcst.start_date})
vt.data = [[pd.to_datetime(std)+relativedelta(months=fcmonth-1) for fcmonth in vt.forecastMonth.values] for std in vt.start_date.values]
hcst = hcst.assign_coords(valid_time=vt)

# Calculate 3-month aggregations
hcst_3m = hcst.rolling(forecastMonth=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
hcst_3m = hcst_3m.where(hcst_3m.forecastMonth>=3,drop=True)

# Calculate 1m anomalies
hcmean = hcst.mean(['number','start_date'])
hcanom = hcst - hcmean
hcanom = hcanom.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
# Calculate 3m anomalies
hcmean_3m = hcst_3m.mean(['number','start_date'])
hcanom_3m = hcst_3m - hcmean_3m
hcanom_3m = hcanom_3m.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))

# %% [markdown]
# ## 2.2 Observations anomalies

# We calculate the monthly and 3-months anomalies for the ERA5 data.

#%%
print("2.2 Observations anomalies")  

# Reading observations from file
era5_1deg = xr.open_dataset(obs_fname, engine='cfgrib')

# Renaming to match hindcast names 
era5_1deg = era5_1deg.rename({'latitude':'lat','longitude':'lon','time':'start_date'}).swap_dims({'start_date':'valid_time'})
# Assign 'forecastMonth' coordinate values
fcmonths = [mm+1 if mm>=0 else mm+13 for mm in [t.month - config['start_month'] for t in pd.to_datetime(era5_1deg.valid_time.values)] ]
era5_1deg = era5_1deg.assign_coords(forecastMonth=('valid_time',fcmonths))
# Drop obs values not needed (earlier than first start date) - this is useful to create well shaped 3-month aggregations from obs.
era5_1deg = era5_1deg.where(era5_1deg.valid_time>=np.datetime64('{hcstarty}-{start_month:02d}-01'.format(**config)),drop=True)

# Calculate 3-month aggregations
era5_1deg_3m = era5_1deg.rolling(valid_time=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
era5_1deg_3m = era5_1deg_3m.where(era5_1deg_3m.forecastMonth>=3,drop=True)

# As we don't need it anymore at this stage, we can safely remove 'forecastMonth'
era5_1deg = era5_1deg.drop('forecastMonth')
era5_1deg_3m = era5_1deg_3m.drop('forecastMonth')

# Calculate 1m anomalies
obmean = era5_1deg.groupby('valid_time.month').mean('valid_time')
obanom = era5_1deg.groupby('valid_time.month') - obmean
# Calculate 3m anomalies
obmean_3m = era5_1deg_3m.groupby('valid_time.month').mean('valid_time')
obanom_3m = era5_1deg_3m.groupby('valid_time.month') - obmean_3m

# We can remove unneeded coordinates
obanom = obanom.drop(['number', 'step', 'start_date', 'isobaricInhPa', 'month'])
obanom_3m = obanom_3m.drop(['number', 'step', 'start_date', 'isobaricInhPa', 'month'])

# %% [markdown]
# ## 2.3 Climatology anomalies

# We calculate the monthly and 3-months anomalies for the ERA5 climatological data.

#%%
print("2.3 Climatology anomalies")  

# Reading climatology from file | This is for obtaining the eofs (just DJF)
era5_clim = xr.open_dataset(clim_fname, engine='cfgrib')
# Obtain desired months
d_months = [12, 1, 2]
# Subsetting data
era5_clim = era5_clim.where(
          #era5_clim.time.dt.year > 1978, # &
          era5_clim.time.dt.month.isin(d_months),
          drop = True).sel(time=slice('1941-{:02d}-01'.format(d_months[0]),'2020-{:02d}-01'.format(d_months[-1]))).rename({'latitude':'lat','longitude':'lon'})
# Calculate 3-month aggregations
era5_clim_3m = era5_clim.rolling(time=3).mean()
era5_clim_3m = era5_clim_3m.where(era5_clim_3m['time.month']==2, drop=True).drop(['number', 'step'])
# Calculate 3m anomalies
clmean_3m = era5_clim_3m.mean('time')
clanom_3m = era5_clim_3m - clmean_3m

# Reading climatology from file | This is for projecting the eofs (corresponding 3-months season)
era5_clim = xr.open_dataset(clim_fname, engine='cfgrib')
# Obtain desired months
d_months = [(config['start_month']+leadm)%12 if config['start_month']+leadm!=12 else 12 for leadm in range(6)]
# Subsetting data
era5_clim = era5_clim.where(
          #era5_clim.time.dt.year > 1978, # &
          era5_clim.time.dt.month.isin(d_months),
          drop = True).sel(time=slice('1941-{:02d}-01'.format(d_months[0]),'2020-{:02d}-01'.format(d_months[-1]))).rename({'latitude':'lat','longitude':'lon'})
# Assign 'forecastMonth' coordinate values
fcmonths = [mm+1 if mm>=0 else mm+13 for mm in [t.month - config['start_month'] for t in pd.to_datetime(era5_clim.valid_time.values)] ]
era5_clim = era5_clim.assign_coords(forecastMonth=('time',fcmonths))
# Calculate 3-month aggregations
era5_clim_3m = era5_clim.rolling(time=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
era5_clim_3m = era5_clim_3m.where(era5_clim_3m.forecastMonth>=3, drop=True).drop(['number', 'step'])
# Calculate 1m anomalies
clmean = era5_clim.groupby('time.month').mean('time')
clanom_pl = era5_clim.groupby('time.month') - clmean
# Calculate 3m anomalies
clmean_3m = era5_clim_3m.mean('time')
clanom_pl_3m = era5_clim_3m - clmean_3m

# Reading surface climatology from file
era5_clim_sf_notp = xr.open_dataset(clim_fname_sf, engine='cfgrib', backend_kwargs={'filter_by_keys': {'step': 0}})
era5_clim_sf_tp = xr.open_dataset(clim_fname_sf, engine='cfgrib', backend_kwargs={'filter_by_keys': {'shortName': 'tp'}})
era5_clim_sf_tp = era5_clim_sf_tp.assign_coords(time=era5_clim_sf_notp.time.values)
era5_clim_sf = xr.merge([era5_clim_sf_notp,era5_clim_sf_tp],compat='override')
del era5_clim_sf_notp, era5_clim_sf_tp
# Obtain desired months
d_months = [(config['start_month']+leadm)%12 if config['start_month']+leadm!=12 else 12 for leadm in range(6)]
# Subsetting data
era5_clim_sf = era5_clim_sf.where(
          #era5_clim.time.dt.year > 1978, # &
          era5_clim_sf.time.dt.month.isin(d_months),
          drop = True).sel(time=slice('1941-{:02d}-01'.format(d_months[0]),'2020-{:02d}-01'.format(d_months[-1]))).rename({'latitude':'lat','longitude':'lon'})
# Assign 'forecastMonth' coordinate values
fcmonths = [mm+1 if mm>=0 else mm+13 for mm in [t.month - config['start_month'] for t in pd.to_datetime(era5_clim_sf.valid_time.values)] ]
era5_clim_sf = era5_clim_sf.assign_coords(forecastMonth=('time',fcmonths))

# Calculate 3-month aggregations
era5_clim_sf_3m = era5_clim_sf.rolling(time=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
era5_clim_sf_3m = era5_clim_sf_3m.where(era5_clim_sf_3m.forecastMonth>=3, drop=True).drop(['number', 'step'])

# Calculate 1m anomalies
clmean_sf = era5_clim_sf.groupby('time.month').mean('time')
clanom_sf = era5_clim_sf.groupby('time.month') - clmean_sf
# Calculate 3m anomalies
clmean_sf_3m = era5_clim_sf_3m.groupby('time.month').mean('time')
clanom_sf_3m = era5_clim_sf_3m.groupby('time.month') - clmean_sf_3m

# %% [markdown]
# ## 2.4  Climatological EOF and PCs

# We calculate the EOFs for the climatological ERA5 data, and their corresponding PCs.
# 
# To calculate the EOFs, we have to weight each grid point by its area.

#%%
print("2.4 Climatological EOF and PCs")  

# Square-root of cosine of latitude weights are applied before the computation of EOFs
coslat = np.cos(np.deg2rad(clanom_3m.coords['lat'].values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]

list_solvers = {"forecastMonth":[],"validMonth":[],"EOFs":[]} 
# For each forecast month
for f in hcanom.forecastMonth:
    # Create an EOF solver to do the EOF analysis.
    cl_array = clanom_3m.copy() #.where(clanom.forecastMonth==f, drop=True)
    cl_array_pl = clanom_pl.where(clanom_pl.forecastMonth==f, drop=True)
    solver = Eof(cl_array.z, weights=wgts)
    list_solvers["forecastMonth"].append(int(f.values))
    list_solvers["validMonth"].append(int(cl_array_pl['valid_time.month'].values[0]))
    list_solvers["EOFs"].append(solver)

    validmonth = config['start_month']+f if config['start_month']+f<=12 else config['start_month']+f-12
    explained_var = solver.varianceFraction()
    variance = pd.DataFrame({'EOF': range(1,5), 'VAR':explained_var[0:4]})
    variance.to_csv(f'{MODESDIR}/ERA5_VAR_{validmonth:02d}.csv')


list_solvers_3m = {"forecastMonth":[],"validMonth":[],"EOFs":[]} 
# For each forecast month
for f in hcanom_3m.forecastMonth:
    # Create an EOF solver to do the EOF analysis.
    cl_array = clanom_3m.copy() #.where(clanom_3m.forecastMonth==f, drop=True)
    cl_array_pl = clanom_pl_3m.where(clanom_pl_3m.forecastMonth==f, drop=True)
    cl_array_sf = clanom_sf_3m.where(clanom_sf_3m.forecastMonth==f, drop=True)
    solver = Eof(cl_array.z, weights=wgts)
    list_solvers_3m["forecastMonth"].append(int(f.values))
    list_solvers_3m["validMonth"].append(int(cl_array_pl['valid_time.month'].values[0]))
    list_solvers_3m["EOFs"].append(solver)
    # Retrieve the leading EOF, expressed as the correlation between the leading PC time series and the input SLP anomalies at each grid point
    eofs_corr = solver.eofsAsCorrelation(neofs=4)
    explained_var = solver.varianceFraction()
    # Extracting the principal component
    #pcs = solver.pcs(npcs=4, pcscaling=1)
    pcs = solver.projectField(cl_array_pl.z, neofs=4, eofscaling=1)

    # Months string
    locale.setlocale(locale.LC_ALL, 'en_GB')
    validmonths = [vm if vm<=12 else vm-12 for vm in [config['start_month'] + (int(f.values)-1) - shift for shift in range(3)]]
    validmonths = [calendar.month_abbr[vm][0] for vm in reversed(validmonths)]
    tit_line = "".join(validmonths)

    # PLOTs EOFs
    fig = plt.figure(figsize=(15,11))
    gs = fig.add_gridspec(2,2)
    ax1 = fig.add_subplot(gs[0,0],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
    fill = ax1.contourf(eofs_corr.lon,eofs_corr.lat,eofs_corr.sel(mode=0).squeeze(),
                    levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                    transform=ccrs.PlateCarree())
    title = 'EOF1 correlation - NAO\n Explained variance: {:.2f} %'.format(explained_var[0].values*100.)
    plt.title(title, fontsize=12)
    ax1.coastlines(resolution='50m')
    ax1.gridlines()
    ax2 = fig.add_subplot(gs[0,1],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
    fill = ax2.contourf(eofs_corr.lon,eofs_corr.lat,eofs_corr.sel(mode=1).squeeze(),
                    levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                    transform=ccrs.PlateCarree())
    title = 'EOF2 correlation - EA \n Explained variance: {:.2f} %'.format(explained_var[1].values*100.)
    plt.title(title, fontsize=12)
    ax2.coastlines(resolution='50m')
    ax2.gridlines()
    ax3 = fig.add_subplot(gs[1,0],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
    fill = ax3.contourf(eofs_corr.lon,eofs_corr.lat,eofs_corr.sel(mode=2).squeeze(),
                    levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                    transform=ccrs.PlateCarree())
    title = 'EOF3 correlation - EAWR \n Explained variance: {:.2f} %'.format(explained_var[2].values*100.)
    plt.title(title, fontsize=12)
    ax3.coastlines(resolution='50m')
    ax3.gridlines()
    ax4 = fig.add_subplot(gs[1,1],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
    fill = ax4.contourf(eofs_corr.lon,eofs_corr.lat,eofs_corr.sel(mode=3).squeeze(),
                    levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                    transform=ccrs.PlateCarree())
    title = 'EOF4 correlation - SCA \n Explained variance: {:.2f} %'.format(explained_var[3].values*100.)
    plt.title(title, fontsize=12)
    ax4.coastlines(resolution='50m')
    ax4.gridlines()
    cbar_ax = fig.add_axes([0.05, 0.05, 0.9, 0.01])
    cb = fig.colorbar(fill, cax=cbar_ax, orientation='horizontal',label=' ')
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, hspace=0.05, wspace=0.07)
    fig.savefig(f'{PLOTSDIR_EOFpatterns}/ERA5_EOFs_{tit_line}.png')

    # PLOTs correlations
    names = ['NAO', 'EA', 'EAWR', 'SCA']
    for modo in range(4):
        mslp_corr = xs.spearman_r(cl_array_sf.msl, pcs.sel(mode=modo), dim='time')
        t2m_corr = xs.spearman_r(cl_array_sf.t2m, pcs.sel(mode=modo), dim='time')
        tp_corr = xs.spearman_r(cl_array_sf.tp, pcs.sel(mode=modo), dim='time')

        fig = plt.figure(figsize=(15,4))
        gs = fig.add_gridspec(1,3)
        ax1 = fig.add_subplot(gs[0,0],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
        fill = ax1.contourf(mslp_corr.lon,mslp_corr.lat,mslp_corr.squeeze(),
                        levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                        transform=ccrs.PlateCarree())
        title = 'Mean sea level pressure - {} correlation'.format(names[modo])
        plt.title(title, fontsize=12)
        ax1.coastlines(resolution='50m')
        ax1.gridlines()
        ax2 = fig.add_subplot(gs[0,1],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
        fill = ax2.contourf(t2m_corr.lon,t2m_corr.lat,t2m_corr.squeeze(),
                        levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                        transform=ccrs.PlateCarree())
        title = '2m temperature - {} correlation'.format(names[modo])
        plt.title(title, fontsize=12)
        ax2.coastlines(resolution='50m')
        ax2.gridlines()
        ax3 = fig.add_subplot(gs[0,2],projection=ccrs.Orthographic(central_longitude=-15, central_latitude=50))
        fill = ax3.contourf(tp_corr.lon,tp_corr.lat,tp_corr.squeeze(),
                        levels=np.linspace(-1.,1.,11), cmap=plt.cm.RdBu_r,
                        transform=ccrs.PlateCarree())
        title = 'Total precipitation - {} correlation'.format(names[modo])
        plt.title(title, fontsize=12)
        ax3.coastlines(resolution='50m')
        ax3.gridlines()
        cbar_ax = fig.add_axes([0.05, 0.08, 0.9, 0.04])
        cb = fig.colorbar(fill, cax=cbar_ax, orientation='horizontal',label=' ')
        #fig.suptitle('', fontsize=16)
        plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, hspace=0.05, wspace=0.07)
        fig.savefig(f'{PLOTSDIR_EOFcorrelations}/ERA5_EOF{str(modo+1)}_{tit_line}.png')

    # PLOTs PCs
    fig = plt.figure(figsize=(8,12))
    gs = fig.add_gridspec(4,1)
    years = pcs.time.dt.year
    ax1 = fig.add_subplot(gs[0,0])
    plt.fill_between(years, 0, pcs.sel(mode=0), where=(pcs.sel(mode=0) >= 0), color='r',
                    interpolate=True)
    plt.fill_between(years, pcs.sel(mode=0), 0, where=(pcs.sel(mode=0) <= 0), color='b',
                    interpolate=True)
    plt.axhline(0, color='k')
    plt.xlabel('Year')
    plt.ylabel('PC1')
    plt.title('$NAO_{EOF}$', loc='left')
    ax1.grid(True,axis='x')
    plt.ylim([-3., 3.])
    ax2 = fig.add_subplot(gs[1,0])
    plt.fill_between(years, 0, pcs.sel(mode=1), where=(pcs.sel(mode=1) >= 0), color='r',
                    interpolate=True)
    plt.fill_between(years, pcs.sel(mode=1), 0, where=(pcs.sel(mode=1) <= 0), color='b',
                    interpolate=True)
    plt.axhline(0, color='k')
    plt.xlabel('Year')
    plt.ylabel('PC2')
    plt.title('$EA_{EOF}$', loc='left')
    ax2.grid(True,axis='x')
    plt.ylim([-3., 3.])
    ax3 = fig.add_subplot(gs[2,0])
    plt.fill_between(years, 0, pcs.sel(mode=2), where=(pcs.sel(mode=2) >= 0), color='r',
                    interpolate=True)
    plt.fill_between(years, pcs.sel(mode=2), 0, where=(pcs.sel(mode=2) <= 0), color='b',
                    interpolate=True)
    plt.axhline(0, color='k')
    plt.xlabel('Year')
    plt.ylabel('PC3')
    plt.title('$EAWR_{EOF}$', loc='left')
    ax3.grid(True,axis='x')
    plt.ylim([-3., 3.])
    ax4 = fig.add_subplot(gs[3,0])
    plt.fill_between(years, 0, pcs.sel(mode=3), where=(pcs.sel(mode=3) >= 0), color='r',
                    interpolate=True)
    plt.fill_between(years, pcs.sel(mode=3), 0, where=(pcs.sel(mode=3) <= 0), color='b',
                    interpolate=True)
    plt.axhline(0, color='k')
    plt.xlabel('Year')
    plt.ylabel('PC4')
    plt.title('$SCA_{EOF}$', loc='left')
    ax4.grid(True,axis='x')
    plt.ylim([-3., 3.])
    plt.subplots_adjust(left=0.1, right=0.9, top=0.95, bottom=0.05, hspace=0.25, wspace=0.1)
    fig.savefig(f'{PLOTSDIR_EOFpcs}/ERA5_PCs_{tit_line}.png')

    plt.figure(figsize=(8, 4))
    eof_num = range(1, 16)
    plt.plot(eof_num, explained_var[0:15], linewidth=2)
    plt.plot(eof_num, explained_var[0:15], linestyle='None', marker="o", color='r', markersize=4)
    plt.axhline(0, color='k')
    plt.xticks(range(1, 16))
    plt.title('Fraction of the total variance represented by each EOF')
    plt.xlabel('EOF #')
    plt.ylabel('Variance Fraction')
    plt.xlim(1, 15)
    plt.ylim([0, 0.6])
    plt.savefig(f'{PLOTSDIR_EOFvariance}/ERA5_VAR_{tit_line}.png')
    variance = pd.DataFrame({'EOF': range(1,5), 'VAR':explained_var[0:4]})
    variance.to_csv(f'{MODESDIR}/ERA5_VAR_{tit_line}.csv')

# %% [markdown]
# ## 2.5  Hindcast PCs

# We project the seasonal forecast data into the four main ERA5 EOFs, 
# obtaining the main climate variability modes: NAO, EA, EA/WR and SCA.

#%%
print("2.5 Hindcast PCs")  

print('- Computing PCs for 1m aggregation (hindcast)')
number_of_eofs = 4
nn = 0
list1_hcpcs = list()
TIME = time.time()
hcanom = hcanom.load()
# For each model member
for n in hcanom.number:
    print('{:.2f}%'.format((float(nn)/float(hcanom.number.size)*100.)))
    list2_hcpcs = list() 
    # For each forecast month
    for f in hcanom.forecastMonth:
        list3_hcpcs = list() 
        # For each year
        for t in hcanom.start_date:
            # Project the z500hPa field in the EOF solver
            solver = list_solvers["EOFs"][list_solvers["forecastMonth"].index(int(f.values))]
            pcs = solver.projectField(hcanom.sel(number=n,forecastMonth=f,start_date=t).z, neofs=number_of_eofs, eofscaling=1)
            list3_hcpcs.append(pcs.assign_coords({'start_date':t}))
        list3 = xr.concat(list3_hcpcs,dim='start_date')                    
        list2_hcpcs.append(list3.assign_coords({'forecastMonth':f}))
    list2 = xr.concat(list2_hcpcs,dim='forecastMonth')                    
    list1_hcpcs.append(list2.assign_coords({'number':n}))
    nn+=1
hcpcs = xr.concat(list1_hcpcs,dim='number').assign_coords({'valid_time':hcanom.valid_time})
print("TIME: --- {} ---".format(dt.timedelta(seconds=(time.time() - TIME)))) 

print('- Computing PCs for 3m aggregation (hindcast)')
nn = 0
list1_hcpcs = list()
TIME = time.time()
hcanom_3m = hcanom_3m.load()
# For each model member
for n in hcanom_3m.number:
    print('{:.2f}%'.format((float(nn)/float(hcanom_3m.number.size)*100.)))
    list2_hcpcs = list()
    # For each forecast month
    for f in hcanom_3m.forecastMonth:
        list3_hcpcs = list()
        # For each year
        for t in hcanom_3m.start_date:
            # Project the z500hPa field in the EOF solver
            solver = list_solvers_3m["EOFs"][list_solvers_3m["forecastMonth"].index(int(f.values))]
            pcs = solver.projectField(hcanom_3m.sel(number=n,forecastMonth=f,start_date=t).z, neofs=number_of_eofs, eofscaling=1)
            list3_hcpcs.append(pcs.assign_coords({'start_date':t}))
        list3 = xr.concat(list3_hcpcs,dim='start_date')                    
        list2_hcpcs.append(list3.assign_coords({'forecastMonth':f}))
    list2 = xr.concat(list2_hcpcs,dim='forecastMonth')                    
    list1_hcpcs.append(list2.assign_coords({'number':n}))
    nn+=1
hcpcs_3m = xr.concat(list1_hcpcs,dim='number').assign_coords({'valid_time':hcanom_3m.valid_time})
print("TIME: --- {} ---".format(dt.timedelta(seconds=(time.time() - TIME)))) 

# %% [markdown]
# ## 2.6  Observed PCs

# The same as in the previous section but with the ERA5 data.

#%%
print("2.6 Obseved PCs")

print('- Computing PCs for 1m aggregation (observation)')
list1_obpcs = list()
obanom = obanom.load()
# For each year
for t in obanom.valid_time:
    # Project the z500hPa field in the EOF solver
    solver = list_solvers["EOFs"][list_solvers["validMonth"].index(int(t['valid_time.month'].values))]
    pcs = solver.projectField(obanom.sel(valid_time=t).z, neofs=number_of_eofs, eofscaling=1)
    list1_obpcs.append(pcs.assign_coords({'valid_time':t}))
obpcs = xr.concat(list1_obpcs,dim='valid_time')      

print('- Computing PCs for 3m aggregation (observation)')
list1_obpcs = list() 
obanom_3m = obanom_3m.load()
# For each year
for t in obanom_3m.valid_time:
    # Project the z500hPa field in the EOF solver
    solver = list_solvers_3m["EOFs"][list_solvers_3m["validMonth"].index(int(t['valid_time.month'].values))]
    pcs = solver.projectField(obanom_3m.sel(valid_time=t).z, neofs=number_of_eofs, eofscaling=1)
    list1_obpcs.append(pcs.assign_coords({'valid_time':t}))
obpcs_3m = xr.concat(list1_obpcs,dim='valid_time')                   

# Saving pcs to netCDF files
hcpcs.to_netcdf(f'{MODESDIR}/{hcst_bname}.1m.PCs.nc')
hcpcs_3m.to_netcdf(f'{MODESDIR}/{hcst_bname}.3m.PCs.nc')
obpcs.to_netcdf(f'{MODESDIR}/{obs_bname}.1m.PCs.nc')
obpcs_3m.to_netcdf(f'{MODESDIR}/{obs_bname}.3m.PCs.nc')
