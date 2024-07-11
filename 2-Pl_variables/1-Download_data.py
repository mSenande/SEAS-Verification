# %% [markdown]
# # 1. Download data

# This script is used to download monthly seasonal forescasts for the hindcast period, 
# along with ERA5 data to compute verification scores.
# 
# We will download data with lead times between 1 and 6 months.
# First we have to decide a forecast system (institution and system name) and a start month. 

# %%
print("1. Download data")  

import os
import sys
from dotenv import load_dotenv
import cdsapi
import numpy as np
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

# %% [markdown]
# ## 1.1 Request data from the CDS using CDS API

# Choose variables, create directories and define CDS configuration.

# %%
print("1.1 Request data from the CDS using CDS API")  

# Here we save the configuration
config = dict(
    list_vars = [ 'geopotential', 'u_component_of_wind', 'v_component_of_wind'],
    pressure_level = '500',
    hcstarty = 1993,
    hcendy = 2016,
    start_month = startmonth,
    origin = model,
    system = system,
)

# Directory selection
load_dotenv() # Data is saved in a path defined in file .env
DATADIR = os.getenv('DATA_DIR')
SCOREDIR = DATADIR + '/scores'
PLOTSDIR = f'./plots/stmonth{config["start_month"]:02d}'
# Directory creation
for directory in [DATADIR, SCOREDIR, PLOTSDIR]:
    # Check if the directory exists
    if not os.path.exists(directory):
        # If it doesn't exist, create it
        try:
            os.makedirs(directory)
        except FileExistsError:
            pass

# CDS configuration
c = cdsapi.Client()

# %% [markdown]
# ## 1.2 Retrieve hindcast data

# We will download all the selected start month and the corresponding forecast months (from 1 to 6)
# for all the hindcast period (from 1993 to 2016).

# %%
print("1.2 Retrieve hindcast data")

# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)
# File name for hindcast
hcst_fname = f'{DATADIR}/{hcst_bname}.grib'

# Check if file exists
if not os.path.exists(hcst_fname):
    # If it doesn't exist, download it
    c.retrieve(
        'seasonal-monthly-pressure-levels',
        {
            'format': 'grib',
            'originating_centre': config['origin'],
            'system': config['system'],
            'variable': config['list_vars'],
            'pressure_level': config['pressure_level'],
            'product_type': 'monthly_mean',
            'year': ['{}'.format(yy) for yy in range(config['hcstarty'],config['hcendy']+1)],
            'month': '{:02d}'.format(config['start_month']),
            'leadtime_month': ['1', '2', '3','4', '5', '6'],
            'grid': '1/1',
            'area': [55, -35, 10, 55], #[50, -30, 25, 5],
        },
        hcst_fname)

# %% [markdown]
# ## 1.3 Retrieve observational data (ERA5)

# Here we will download the same months as for the hindcast data.

# %%
print("1.3 Retrieve observational data (ERA5)")

# File name for observations
obs_fname = '{fpath}/era5_monthly_stmonth{start_month:02d}_{hcstarty}-{hcendy}.grib'.format(fpath=DATADIR,**config)

# Check if file exists
if not os.path.exists(obs_fname):
    # If it doesn't exist, download it
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': config['list_vars'],
            'pressure_level': config['pressure_level'],
            # NOTE from observations we need to go one year beyond so we have available all the right valid dates
            # e.g. Nov.2016 start date forecast goes up to April 2017             
            'year': ['{}'.format(yy) for yy in range(config['hcstarty'],config['hcendy']+2)],
            'month': ['{:02d}'.format((config['start_month']+leadm)%12) if config['start_month']+leadm!=12 else '12' for leadm in range(6)],
            'time': '00:00',
            'grid': '1/1',
            'area': [55, -35, 10, 55], #[50, -30, 25, 5],
            'format': 'grib',
        },
        obs_fname)

