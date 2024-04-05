# %% [markdown]
# # 2. Compute NAO

# This script is used to compute the EOFs and PCs .
# from daily (or subdaily) seasonal forescasts for the hindcast period.
# 
# EOFs are calculted for ERA5 data with the Eof library.
# Seasonal forecasts are then projected into these EOFs to obtain the PCs.
#
# First we have to decide a forecast system (institution and system name) and a start month. 

#%%
print("2. Compute EOFs") 

import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
from eofs.xarray import Eof
import xskillscore as xs
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
DATADIR = os.environ['MAS'] + '/Seasonal_Verification/3-Var_modes/2-EOFs_calc/daily_data_extension/data'
MODESDIR = DATADIR + '/modes'
SCOREDIR = DATADIR + '/scores'
PLOTSDIR = f'./plots/stmonth{config["start_month"]:02d}'
# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_daily'.format(**config)
# Base name for observations
obs_bname = 'era5_daily_stmonth{start_month:02d}_{hcstarty}-{hcendy}'.format(**config)
# File name for observations
obs_fname = f'{DATADIR}/{obs_bname}.grib'
# File name for climatologies
clim_fname = f'../{DATADIR}/era5_monthly.grib'
clim_fname_sf = f'../{DATADIR}/era5_monthly_sf.grib'

# Check if files exist
if not os.path.exists(obs_fname):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()
elif not os.path.exists(clim_fname):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()
elif not os.path.exists(clim_fname_sf):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()

# %% [markdown]
# ## 2.1 Convert daily/subdaily data to monthly data

# We need to reshape daily/subdaily data to monthly data by calculating the monthly average.
# 
# Sice seasonal forecasts are indexed by leadtime hour, 
# we need to create a valid_time coordinate with the correct date information before computing the monthly averages.
# 
# Then we assign to the dataset a forecastMonth coordinate.
    
#%%
print("2.1 Convert daily/subdaily data to monthly data")

# In "burst" mode models all the members are initialized on the same start date ['ecmwf', 'meteo_france', 'dwd', 'cmcc', 'eccc']
# In "lagged" mode models, members are initialized on different start dates ['ukmo', 'ncep', 'jma']
# Here we select the start dates:
days_on_year = pd.date_range('1993-01-01','1993-12-31',freq='d')
if '{origin}'.format(**config)=='ukmo':
    # Met Office forecasting system is initialized on days 1, 9, 17 and 25 of each month    
    start_days = pd.DatetimeIndex(['{hcstarty}-{start_month:02d}-01'.format(**config),
                                '{hcstarty}-{start_month:02d}-09'.format(**config),
                                '{hcstarty}-{start_month:02d}-17'.format(**config),
                                '{hcstarty}-{start_month:02d}-25'.format(**config)])
elif '{origin}'.format(**config)=='ncep':
    # NCEP forecasting system is initialized each 5 days and each 6h   
    start_days = days_on_year[::5]
    start_days = start_days[start_days.month==config['start_month']]
    start_days = pd.DatetimeIndex([d for d in start_days]+  
                                  [d+relativedelta(hours=6) for d in start_days]+
                                  [d+relativedelta(hours=12) for d in start_days]+
                                  [d+relativedelta(hours=18) for d in start_days])

elif '{origin}'.format(**config)=='jma':
    # JMA forecasting system is initialized each 15 days    
    start_days = days_on_year[::15]
    start_days = start_days[start_days.month==config['start_month']]
else:
    start_days = pd.DatetimeIndex(['{hcstarty}-{start_month:02d}-01'.format(**config)])

f = 0
l_files_hcst=list()
reset_start = pd.DatetimeIndex([str(year)+'-{start_month:02d}-01'.format(**config) for year in range(config['hcstarty'],config['hcendy']+1)])
# We read each start day from a separated file
for start_day in start_days:

    # File name for hindcast
    hcst_fname = f'{DATADIR}/{hcst_bname}_{start_day.day:02d}.grib'

    # Reading hindcast data from file
    hcst_daily = xr.open_dataset(hcst_fname,engine='cfgrib', backend_kwargs=dict(time_dims=('step', 'time')))
    # We use dask.array with chunks on leadtime, latitude and longitude coordinate
    hcst_daily = hcst_daily.chunk({'time': 1,'latitude':'auto', 'longitude':'auto'})
    # Reanme coordinates to match those of observations
    hcst_daily = hcst_daily.rename({'latitude':'lat','longitude':'lon', 'time':'start_date'})
    # Here we have to select the desired hours among the start_dates (this step is only needed for NCEP forecasting system) 
    hcst_daily = hcst_daily.where(hcst_daily['start_date.hour']==start_day.hour,drop=True)
    
    # Add start_month to the xr.Dataset
    start_month = pd.to_datetime(hcst_daily.start_date.values[0]).month
    hcst_daily = hcst_daily.assign_coords({'start_month':start_month})
    # Add valid_time to the xr.Dataset
    vt = xr.DataArray(dims=('start_date','step'), coords={'step':hcst_daily.step,'start_date':hcst_daily.start_date})
    vt.data = [[pd.to_datetime(std)+relativedelta(hours=int(step)) for step in vt.step.values.astype('timedelta64[h]')/np.timedelta64(1, 'h')] for std in vt.start_date.values]
    hcst_daily = hcst_daily.assign_coords(valid_time=vt)

    # Resample daily/subdaily data to monthly data
    l_means_hcst=list()
    # For each year
    for start in hcst_daily.start_date.values:
        # Monthly resample
        hcst_mean = hcst_daily.sel(start_date=start).resample(valid_time="1MS").mean()
        # Select desired months (start_month + forecastMonth [5,6,7])
        hcst_mean = hcst_mean.isel(valid_time=hcst_mean.valid_time.dt.month.isin([(config['start_month']+leadm)%12 if config['start_month']+leadm!=12 else 12 for leadm in range(4,7)]))
        # Rename valid_time to forecastMonth
        hcst_mean = hcst_mean.rename({'valid_time':'forecastMonth'})
        # Assign forecast month values
        hcst_mean = hcst_mean.assign_coords(forecastMonth=[leadm+1 for leadm in range(4,7)])
        # Append loop results
        l_means_hcst.append(hcst_mean)
    # Concatenate along start_date
    hcst_init = xr.concat(l_means_hcst,dim='start_date').load()                    
    # Rename ensemble numbers, to concate along them latter    
    hcst_init['number'] = np.arange(hcst_init.number.size)+hcst_init.number.size*f
    # Rename start_date
    hcst_init['start_date'] = reset_start
    # Append loop results
    l_files_hcst.append(hcst_init)

    f+=1
# Concatenate along numbers (i.e we rename new start days as ensemble members)
hcst = xr.concat(l_files_hcst,dim='number')  

# %% [markdown]
# ## 2.2 Hindcast anomalies

# We calculate the monthly and 3-months anomalies for the hindcast data.

#%%
print("2.2 Hindcast anomalies")

# We use dask.array with chunks on leadtime, latitude and longitude coordinate
hcst = hcst.chunk({'forecastMonth':1, 'lat':'auto', 'lon':'auto'})

# Add (again) valid_time to the xr.Dataset
vt = xr.DataArray(dims=('start_date','forecastMonth'), coords={'forecastMonth':hcst.forecastMonth,'start_date':hcst.start_date})
vt.data = [[pd.to_datetime(std)+relativedelta(months=fcmonth-1) for fcmonth in vt.forecastMonth.values] for std in vt.start_date.values]
hcst = hcst.assign_coords(valid_time=vt)

# Calculate 3-month aggregations
hcst_3m = hcst.rolling(forecastMonth=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
hcst_3m = hcst_3m.where(hcst_3m.forecastMonth>=7,drop=True)

# Calculate 1m anomalies
hcmean = hcst.mean(['number','start_date'])
hcanom = hcst - hcmean
hcanom = hcanom.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
# Calculate 3m anomalies
hcmean_3m = hcst_3m.mean(['number','start_date'])
hcanom_3m = hcst_3m - hcmean_3m
hcanom_3m = hcanom_3m.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))

# %% [markdown]
# ## 2.3 Observations anomalies

# We calculate the monthly and 3-months anomalies for the ERA5 data.
# 
# The conversion between daily/subdaily data and monthly data is also done here, 
# as it is simpler than with forecasts.

#%%
print("2.3 Observations anomalies")  

# Reading observations from file
era5_1deg_daily = xr.open_dataset(obs_fname, engine='cfgrib')
# Monthly resample
era5_1deg = era5_1deg_daily.resample(time="1MS").mean()

# Renaming to match hindcast names 
era5_1deg = era5_1deg.rename({'latitude':'lat','longitude':'lon','time':'valid_time'})
# Assign 'forecastMonth' coordinate values
fcmonths = [mm+1 if mm>=0 else mm+13 for mm in [t.month - config['start_month'] for t in pd.to_datetime(era5_1deg.valid_time.values)] ]
era5_1deg = era5_1deg.assign_coords(forecastMonth=('valid_time',fcmonths))
# Drop obs values not needed (forecastMonth not equal to 5, 6 or 7)
era5_1deg = era5_1deg.where(era5_1deg.forecastMonth.isin([leadm+1 for leadm in range(4,7)]),drop=True)
# Drop obs values not needed (earlier than first start date) - this is useful to create well shaped 3-month aggregations from obs.
era5_1deg = era5_1deg.where(era5_1deg.valid_time>=np.datetime64('{hcstarty}-{start_month:02d}-01'.format(**config)),drop=True)

# Calculate 3-month aggregations
era5_1deg_3m = era5_1deg.rolling(valid_time=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
era5_1deg_3m = era5_1deg_3m.where(era5_1deg_3m.forecastMonth>=7,drop=True)

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
obanom = obanom.drop(['number', 'step', 'isobaricInhPa', 'month'])
obanom_3m = obanom_3m.drop(['number', 'step', 'isobaricInhPa', 'month'])

# %% [markdown]
# ## 2.4 Climatology anomalies

# We calculate the monthly and 3-months anomalies for the ERA5 climatological data.

#%%
print("2.4 Climatology anomalies")  

# Reading climatology from file
era5_clim = xr.open_dataset(clim_fname, engine='cfgrib')

# Subsetting data
era5_clim = era5_clim.where(
          #era5_clim.time.dt.year > 1978, # &
          era5_clim.time.dt.month.isin([12, 1, 2]),
          drop = True).sel(time=slice('1941-11-01','2020-11-01')).rename({'latitude':'lat','longitude':'lon'})

# Calculate 3-month aggregations
era5_clim_3m = era5_clim.rolling(time=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
era5_clim_3m = era5_clim_3m.where(era5_clim_3m.time.dt.month.isin(2), drop=True).drop(['number', 'step', 'valid_time'])

# Calculate anomalies
clmean = era5_clim_3m.mean(dim='time')
clanom = era5_clim_3m - clmean

# Reading surface climatology from file
era5_clim_sf_notp = xr.open_dataset(clim_fname_sf, engine='cfgrib', backend_kwargs={'filter_by_keys': {'step': 0}})
era5_clim_sf_tp = xr.open_dataset(clim_fname_sf, engine='cfgrib', backend_kwargs={'filter_by_keys': {'shortName': 'tp'}})
era5_clim_sf_tp = era5_clim_sf_tp.assign_coords(time=era5_clim_sf_notp.time.values)
era5_clim_sf = xr.merge([era5_clim_sf_notp,era5_clim_sf_tp],compat='override')
del era5_clim_sf_notp, era5_clim_sf_tp

# Subsetting data
era5_clim_sf = era5_clim_sf.where(
          #era5_clim_sf.time.dt.year > 1978, # &
          era5_clim_sf.time.dt.month.isin([12, 1, 2]),
          drop = True).sel(time=slice('1941-11-01','2020-11-01')).rename({'latitude':'lat','longitude':'lon'})

# Calculate 3-month aggregations
era5_clim_sf_3m = era5_clim_sf.rolling(time=3).mean()
# rolling() assigns the label to the end of the N month period, so the first N-1 elements have NaN and can be dropped
era5_clim_sf_3m = era5_clim_sf_3m.where(era5_clim_sf_3m.time.dt.month.isin(2), drop=True).drop(['number', 'step', 'valid_time'])

# Calculate anomalies
clmean_sf = era5_clim_sf_3m.mean(dim='time')
clanom_sf = era5_clim_sf_3m - clmean_sf

# %% [markdown]
# ## 2.5  Climatological EOF and PCs

# We calculate the EOFs for the climatological ERA5 data, and their corresponding PCs.
# 
# To calculate the EOFs, we have to weight each grid point by its area.

#%%
print("2.5 Climatological EOF and PCs")  

# Square-root of cosine of latitude weights are applied before the computation of EOFs
coslat = np.cos(np.deg2rad(clanom.coords['lat'].values)).clip(0., 1.)
wgts = np.sqrt(coslat)[..., np.newaxis]
# Create an EOF solver to do the EOF analysis.
solver = Eof(clanom.z, weights=wgts)
# Retrieve the leading EOF, expressed as the correlation between the leading PC time series and the input SLP anomalies at each grid point
eofs_corr = solver.eofsAsCorrelation(neofs=4)
explained_var = solver.varianceFraction()
# Extracting the principal component
pcs = solver.pcs(npcs=4, pcscaling=1)

# %% [markdown]
# ## 2.6  Hindcast PCs

# We project the seasonal forecast data into the four main ERA5 EOFs, 
# obtaining the main climate variability modes: NAO, EA, EA/WR and SCA.

#%%
print("2.6 Hindcast PCs")  

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
# ## 2.7  Observed PCs

# The same as in the previous section but with the ERA5 data.

#%%
print("2.7 Obseved PCs")

print('- Computing PCs for 1m aggregation (observation)')
list1_obpcs = list()
obanom = obanom.load()
# For each year
for t in obanom.valid_time:
    # Project the z500hPa field in the EOF solver
    pcs = solver.projectField(obanom.sel(valid_time=t).z, neofs=number_of_eofs, eofscaling=1)
    list1_obpcs.append(pcs.assign_coords({'valid_time':t}))
obpcs = xr.concat(list1_obpcs,dim='valid_time')

print('- Computing PCs for 3m aggregation (observation)')
list1_obpcs = list() 
obanom_3m = obanom_3m.load()
# For each year
for t in obanom_3m.valid_time:
    # Project the z500hPa field in the EOF solver
    pcs = solver.projectField(obanom_3m.sel(valid_time=t).z, neofs=number_of_eofs, eofscaling=1)
    list1_obpcs.append(pcs.assign_coords({'valid_time':t}))
obpcs_3m = xr.concat(list1_obpcs,dim='valid_time')                  

# Saving NAO index to netCDF files
hcpcs.to_netcdf(f'{MODESDIR}/{hcst_bname}.1m.PCs.nc')
hcpcs_3m.to_netcdf(f'{MODESDIR}/{hcst_bname}.3m.PCs.nc')
obpcs.to_netcdf(f'{MODESDIR}/{obs_bname}.1m.PCs.nc')
obpcs_3m.to_netcdf(f'{MODESDIR}/{obs_bname}.3m.PCs.nc')
