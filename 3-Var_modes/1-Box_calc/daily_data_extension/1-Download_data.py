# %% [markdown]
# # 1. Download data

# This script is used to download daily (or subdaily) seasonal forescasts for the hindcast period, 
# along with ERA5 data to compute verification scores.
# 
# Since we use daily data to extend the monthly forecasts, we will download data with lead times between 4 and 7 months.
# First we have to decide a forecast system (institution and system name) and a start month. 

# %%
print("1. Download data")  

import os
import sys
import cdsapi
import numpy as np
import pandas as pd
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

# %% [markdown]
# ## 1.1 Request data from the CDS using CDS API

# Choose variables and time interval, create directories and define CDS configuration.

# %%
print("1.1 Request data from the CDS using CDS API")  

# Here we save the configuration
config = dict(
    list_vars = 'mean_sea_level_pressure',
    hcstarty = 1993,
    hcendy = 2016,
    start_month = startmonth,
    origin = model,
    system = system,
)

# Directory selection
DATADIR = './data'
MODESDIR = DATADIR + '/modes'
SCOREDIR = DATADIR + '/scores'
PLOTSDIR = DATADIR + f'/plots/stmonth{config["start_month"]:02d}'
# Directory creation
for directory in [DATADIR, MODESDIR, SCOREDIR, PLOTSDIR]:
    # Check if the directory exists
    if not os.path.exists(directory):
        # If it doesn't exist, create it
        try:
            os.makedirs(directory)
        except FileExistsError:
            pass

# CDS configuration
CDSAPI_URL = 'https://cds.climate.copernicus.eu/api/v2'
CDSAPI_KEY = input("CDS KEY: ")
c = cdsapi.Client(url=CDSAPI_URL, key=CDSAPI_KEY)

# %% [markdown]
# ## 1.2 Retrieve hindcast data

# We define a time step (24h in this example), and then we calculate the lead hours of the fourth and seventh forecast month.
# We have to download all the time steps between these two lead hours.

# %%
print("1.2 Retrieve hindcast data")

# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_daily'.format(**config)

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
    # NCEP forecasting system is initialized each 5 days    
    start_days = days_on_year[::5]
    start_days = start_days[start_days.month==config['start_month']]
elif '{origin}'.format(**config)=='jma':
    # JMA forecasting system is initialized each 15 days    
    start_days = days_on_year[::15]
    start_days = start_days[start_days.month==config['start_month']]
else:
    start_days = pd.DatetimeIndex(['{hcstarty}-{start_month:02d}-01'.format(**config)])

# We download each start day in a separated file
for start_day in start_days:
    # File name for hindcast
    hcst_fname = f'{DATADIR}/{hcst_bname}_{start_day.day:02d}.grib'

    # Compute lead_hours
    step = 24 # Values each 24 hours
    end_day1 = start_day+relativedelta(months=4)
    end_day2 = start_day+relativedelta(months=7)
    lead_hours_1 = (end_day1 - start_day).days*24
    lead_hours_2 = (end_day2 - start_day).days*24

    # Check if file exists
    if not os.path.exists(hcst_fname):
        # If it doesn't exist, download it
        c.retrieve(
            'seasonal-original-single-levels',
            {
                'format': 'grib',
                'originating_centre': config['origin'],
                'system': config['system'],
                'variable': config['list_vars'],
                'year': ['{}'.format(yy) for yy in range(config['hcstarty'],config['hcendy']+1)],
                'month': '{:02d}'.format(config['start_month']),
                'day': '{:02d}'.format(start_day.day),
                'leadtime_hour': ['{}'.format(hour) for hour in range(lead_hours_1,lead_hours_2+step,step)],
                'grid': '1/1',
                'area': [90, -90, 20, 60],
            },
            hcst_fname)

# %% [markdown]
# ## 1.3 Retrieve observational data (ERA5)
 
# Here we can simply select desired months, with no need to calculate lead hours.

# %%
print("1.3 Retrieve observational data (ERA5)")

# File name for observations
obs_fname = '{fpath}/era5_daily_stmonth{start_month:02d}_{hcstarty}-{hcendy}.grib'.format(fpath=DATADIR,**config)

# Check if file exists
if not os.path.exists(obs_fname):
    # If it doesn't exist, download it
    c.retrieve(
        'reanalysis-era5-single-levels',
        {
            'product_type': 'reanalysis',
            'variable': config['list_vars'],
            # NOTE from observations we need to go one year beyond so we have available all the right valid dates        
            'year': ['{}'.format(yy) for yy in range(config['hcstarty'],config['hcendy']+2)],
            'month': ['{:02d}'.format((config['start_month']+leadm)%12) if config['start_month']+leadm!=12 else '12' for leadm in range(4,7)],
            'day': [
                '01', '02', '03',
                '04', '05', '06',
                '07', '08', '09',
                '10', '11', '12',
                '13', '14', '15',
                '16', '17', '18',
                '19', '20', '21',
                '22', '23', '24',
                '25', '26', '27',
                '28', '29', '30',
                '31',
            ],
            'time': ['{:02d}:00'.format(hour) for hour in range(0,24,step)],
            'grid': '1/1',
            'area': [90, -90, 20, 60],
            'format': 'grib',
        },
        obs_fname)
