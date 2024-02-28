# %% [markdown]
# # 2. Compute NAO

# This script is used to compute the North Atlantic Oscillation (NAO) 
# from monthly seasonal forescasts for the hindcast period.
# 
# NAO is calculated as the difference of the mean sea level pressure normalized anomaly 
# between a north and a south lat-lon box region in the North Atlantic. 
#
# First we have to decide a forecast system (institution and system name) and a start month. 

#%%
print("2. Compute NAO") 

import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
from dateutil.relativedelta import relativedelta
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
        name = input("Sistema del modelo [ System 12 , System 13 , System 14 , System 15 , GloSea6 , GloSea6.1 , GloSea6.2 ]: ")
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
    list_vars = ['mean_sea_level_pressure'],
    hcstarty = 1993,
    hcendy = 2016,
    start_month = startmonth,
    origin = model,
    system = system,
    isLagged = False if model in ['ecmwf', 'meteo_france', 'dwd', 'cmcc', 'eccc'] else True
)

# Directory selection
DATADIR = './data'
# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)
# File name for hindcast
hcst_fname = f'{DATADIR}/{hcst_bname}.grib'
# Base name for observations
obs_bname = 'era5_monthly_stmonth{start_month:02d}_{hcstarty}-{hcendy}'.format(**config)
# File name for observations
obs_fname = f'{DATADIR}/{obs_bname}.grib'

# Check if files exist
if not os.path.exists(obs_fname):
    print('No se descargaron aún los datos de ERA5')
    sys.exit()
elif not os.path.exists(hcst_fname):
    print('No se descargaron aún los datos de este modelo y sistema')
    sys.exit()

# %% [markdown]
# ## 2.1 Hindcast anomalies

# We calculate the monthly and 3-months anomalies for the hindcast data.
# 
# We also obtain the normalized anomalies dividing by the standard deviation.

#%%
print("2.1 Hindcast anomalies")

# For the re-shaping of time coordinates in xarray.Dataset we need to select the right one 
#  -> burst mode ensembles (e.g. ECMWF SEAS5) use "time". This is the default option in this notebook
#  -> lagged start ensembles (e.g. MetOffice GloSea6) use "indexing_time" (see CDS documentation about nominal start date)
st_dim_name = 'time' if not config.get('isLagged',False) else 'indexing_time'

# Reading hindcast data from file
hcst = xr.open_dataset(hcst_fname,engine='cfgrib', backend_kwargs=dict(time_dims=('forecastMonth', st_dim_name)))
# We use dask.array with chunks on leadtime, latitude and longitude coordinate
hcst = hcst.chunk({'forecastMonth':1, 'latitude':'auto', 'longitude':'auto'})
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
hcstd = hcst.std(['number','start_date'])
hcanom = hcst - hcmean
hcanom = hcanom.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
hcanom_norm = (hcst - hcmean)/hcstd
hcanom_norm = hcanom_norm.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
# Calculate 3m anomalies
hcmean_3m = hcst_3m.mean(['number','start_date'])
hcstd_3m = hcst_3m.std(['number','start_date'])
hcanom_3m = hcst_3m - hcmean_3m
hcanom_3m = hcanom_3m.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
hcanom_norm_3m = (hcst_3m - hcmean_3m)/hcstd_3m
hcanom_norm_3m = hcanom_norm_3m.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
# NOTE: Hindacast data is already splitted into years and months, with years represented in start_date and months represented in forecastMonth

# %% [markdown]
# ## 2.2 Observations anomalies

# We calculate the monthly and 3-months anomalies for the ERA5 data.
# We also obtain the normalized anomalies dividing by the standard deviation.

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
obstd = era5_1deg.groupby('valid_time.month').std('valid_time')
obanom = era5_1deg.groupby('valid_time.month') - obmean
obanom_norm = xr.apply_ufunc(
    lambda x, m, s: (x - m) / s,
    era5_1deg.groupby('valid_time.month'),
    obmean,
    obstd,
)
# Calculate 3m anomalies
obmean_3m = era5_1deg_3m.groupby('valid_time.month').mean('valid_time')
obstd_3m = era5_1deg_3m.groupby('valid_time.month').std('valid_time')
obanom_3m = era5_1deg_3m.groupby('valid_time.month') - obmean_3m
obanom_norm_3m = xr.apply_ufunc(
    lambda x, m, s: (x - m) / s,
    era5_1deg_3m.groupby('valid_time.month'),
    obmean_3m,
    obstd_3m,
)
# NOTE: Observation data is not splitted into years and months, so in order to calculate anomalies, we need to group by months

# We can remove unneeded coordinates
obanom = obanom.drop(['number', 'step', 'start_date', 'month'])
obanom_3m = obanom_3m.drop(['number', 'step', 'start_date', 'month'])
obanom_norm = obanom_norm.drop(['number', 'step', 'start_date', 'month'])
obanom_norm_3m = obanom_norm_3m.drop(['number', 'step', 'start_date', 'month'])

# %% [markdown]
# ## 2.3  Hindcast NAO index

# We calculate the NAO as the difference of the mean sea level pressure normalized anomaly 
# in a north (Iceland) and a south (Portugal) lat-lon box region.
# 
# To calculate the average of the mslp in each box, we have to weight each grid point by its area.

#%%
print("2.3 Hindcast NAO index")  

# We defined the 2 boxes for the NAO calculation
Iceland = {'lat': slice(90., 55.), 'lon': slice(-90., 60.)}
Portugal = {'lat': slice(55., 20.), 'lon': slice(-90., 60.)}

# 1m MSLP weighted-average for northern box
msl_Iceland = hcanom_norm.msl.sel(lat=Iceland['lat'], lon=Iceland['lon'])
weights = np.cos(np.deg2rad(msl_Iceland.lat))
msl_Iceland = msl_Iceland.weighted(weights).mean(("lon", "lat"))
# 1m MSLP weighted-average for southern box
msl_Portugal = hcanom_norm.msl.sel(lat=Portugal['lat'], lon=Portugal['lon'])
weights = np.cos(np.deg2rad(msl_Portugal.lat))
msl_Portugal = msl_Portugal.weighted(weights).mean(("lon", "lat"))
# 3m MSLP weighted-average for northern box
msl_Iceland_3m = hcanom_norm_3m.msl.sel(lat=Iceland['lat'], lon=Iceland['lon'])
weights = np.cos(np.deg2rad(msl_Iceland_3m.lat))
msl_Iceland_3m = msl_Iceland_3m.weighted(weights).mean(("lon", "lat"))
# 3m MSLP weighted-average for southern box
msl_Portugal_3m = hcanom_norm_3m.msl.sel(lat=Portugal['lat'], lon=Portugal['lon'])
weights = np.cos(np.deg2rad(msl_Portugal_3m.lat))
msl_Portugal_3m = msl_Portugal_3m.weighted(weights).mean(("lon", "lat"))

# 1m MSLP box difference 
hcnao_station = msl_Portugal-msl_Iceland
# 3m MSLP box difference 
hcnao_station_3m = msl_Portugal_3m-msl_Iceland_3m

# %% [markdown]
# ## 2.4  Observed NAO index

# The same as in the previous section but with the ERA5 data.

#%%
print("2.4 Obseved NAO index")

# 1m MSLP weighted-average for northern box
msl_Iceland = obanom_norm.msl.sel(lat=Iceland['lat'], lon=Iceland['lon'])
weights = np.cos(np.deg2rad(msl_Iceland.lat))
msl_Iceland = msl_Iceland.weighted(weights).mean(("lon", "lat"))
# 1m MSLP weighted-average for southern box
msl_Portugal = obanom_norm.msl.sel(lat=Portugal['lat'], lon=Portugal['lon'])
weights = np.cos(np.deg2rad(msl_Portugal.lat))
msl_Portugal = msl_Portugal.weighted(weights).mean(("lon", "lat"))
# 3m MSLP weighted-average for northern box
msl_Iceland_3m = obanom_norm_3m.msl.sel(lat=Iceland['lat'], lon=Iceland['lon'])
weights = np.cos(np.deg2rad(msl_Iceland_3m.lat))
msl_Iceland_3m = msl_Iceland_3m.weighted(weights).mean(("lon", "lat"))
# 3m MSLP weighted-average for southern box
msl_Portugal_3m = obanom_norm_3m.msl.sel(lat=Portugal['lat'], lon=Portugal['lon'])
weights = np.cos(np.deg2rad(msl_Portugal_3m.lat))
msl_Portugal_3m = msl_Portugal_3m.weighted(weights).mean(("lon", "lat"))

# 1m MSLP box difference 
obnao_station = msl_Portugal-msl_Iceland
# 3m MSLP box difference 
obnao_station_3m = msl_Portugal_3m-msl_Iceland_3m

# Saving NAO index to netCDF files
hcnao = hcnao_station.rename("nao")
hcnao_3m = hcnao_station_3m.rename("nao")
obnao = obnao_station.rename("nao")
obnao_3m = obnao_station_3m.rename("nao")
hcnao.to_netcdf(f'{DATADIR}/modes/{hcst_bname}.1m.NAO.nc')
hcnao_3m.to_netcdf(f'{DATADIR}/modes/{hcst_bname}.3m.NAO.nc')
obnao.to_netcdf(f'{DATADIR}/modes/{obs_bname}.1m.NAO.nc')
obnao_3m.to_netcdf(f'{DATADIR}/modes/{obs_bname}.3m.NAO.nc')
