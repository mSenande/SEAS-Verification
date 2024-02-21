import os
import sys
import inquirer
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
    # Which model
    questions = [
    inquirer.List('institution',
                    message="Usar modelo del siguiente organismo",
                    choices=['ECMWF','Météo France','Met Office','DWD','CMCC','NCEP','JMA','ECCC'],
                ),
    ]
    answers = inquirer.prompt(questions)
    institution = answers["institution"]

    # Which version of each model
    if institution=='ECMWF':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['System 4','SEAS5','SEAS5.1'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='Météo France':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['System 5','System 6','System 7','System 8'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='Met Office':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['System 12','System 13','System 14','System 15','GloSea6','GloSea6.1','GloSea6.2'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='DWD':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['GCFS2.0','GCFS2.1'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='CMCC':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['SPSv3.0','SPSv3.5'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='NCEP':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['CFSv2'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='JMA':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['CPS2','CPS3'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    elif institution=='ECCC':
        questions2 = [
        inquirer.List('name',
                        message="Sistema del modelo",
                        choices=['GEM-NEMO','CanCM4i','GEM5-NEMO'],
                    ),
        ]
        answers2 = inquirer.prompt(questions2)
        name = answers2["name"]
    else:
        sys.exit()

    # Which start month
    questions3 = [
    inquirer.List('startmonth',
                    message="Mes de inicialización",
                    choices=['Enero', 'Febrero','Marzo','Abril','Mayo','Junio','Julio','Agosto','Septiembre','Octubre','Noviembre','Diciembre'],
                ),
    ]
    answers3 = inquirer.prompt(questions3)
    startmonth= np.where(np.array(['Enero', 'Febrero','Marzo','Abril','Mayo','Junio','Julio','Agosto','Septiembre','Octubre','Noviembre','Diciembre'])== answers3["startmonth"])[0][0]+1


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

### 1. Request data from the CDS using CDS API ###
##################################################
print("1. Request data from the CDS using CDS API")  

# Here we save the configuration
config = dict(
    list_vars = 'geopotential',
    pressure_level = '500',
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
CDSAPI_KEY = '7068:d94dcd65-20e0-4ee7-ab4a-15b4aa11b321'
c = cdsapi.Client(url=CDSAPI_URL, key=CDSAPI_KEY)

### 1a. Retrieve hindcast data ###
print("1a. Retrieve hindcast data")

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
            'area': [80, -90, 20, 60],
        },
        hcst_fname)

### 1b. Retrieve observational data (ERA5) ###
print("1b. Retrieve observational data (ERA5)")

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
            'area': [80, -90, 20, 60],
            'format': 'grib',
        },
        obs_fname)

### 1c. Retrieve climaological data (ERA5) ###
print("1c. Retrieve climaological data (ERA5)")

# File name for climatology
clim_fname = f'{DATADIR}/era5_monthly.grib'

# Check if file exists
if not os.path.exists(clim_fname):
    # If it doesn't exist, download it
    c.retrieve(
        'reanalysis-era5-pressure-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': config['list_vars'],
            'pressure_level': config['pressure_level'],
            'year': ['{}'.format(yy) for yy in range(1940,2021)],
            'month': ['{}'.format(mm) for mm in range(1,13)],
            'time': '00:00',
            'grid': '1/1',
            'area': [80, -90, 20, 60],
            'format': 'grib',
        },
        clim_fname)

# File name for surface climatology
clim_fname2 = f'{DATADIR}/era5_monthly_sf.grib'

# Check if file exists
if not os.path.exists(clim_fname2):
    # If it doesn't exist, download it
    c.retrieve(
        'reanalysis-era5-single-levels-monthly-means',
        {
            'product_type': 'monthly_averaged_reanalysis',
            'variable': ['2m_temperature', 'total_precipitation', 'mean_sea_level_pressure'],
            'year': ['{}'.format(yy) for yy in range(1940,2021)],
            'month': ['{}'.format(mm) for mm in range(1,13)],
            'time': '00:00',
            'grid': '1/1',
            'area': [80, -90, 20, 60],
            'format': 'grib',
        },
        clim_fname2)
