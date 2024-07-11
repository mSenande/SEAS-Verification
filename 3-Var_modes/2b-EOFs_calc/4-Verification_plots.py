# %% [markdown]
# # 4. Visualize verification plots

# This script is used to create time-series with the variability modes forecasts of each forecasting system. 
# 
# First we have to decide a forecast system (institution and system name), a start month a month aggregation and a end month. 

#%%
print("4. Visualize verification plots")  

import os
import sys
from dotenv import load_dotenv
import xarray as xr
import numpy as np
import locale
import calendar
from matplotlib import pyplot as plt
import warnings
warnings.filterwarnings('ignore')


# If variables are introduced from the command line, read them
if len(sys.argv) > 2:
    institution = str(sys.argv[1]).replace('"', '')
    name = str(sys.argv[2]).replace('"', '')
    startmonth = int(sys.argv[3])
    aggr = str(sys.argv[4])
    fcmonth = int(sys.argv[5])
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

    # Subset of plots to be produced
    aggr = input("Selecciona el tipo de agregación mensual [ 1m , 3m ]: ")

    # Forecast month
    if aggr=='1m':
        answ = input("Resultados para el mes [ Jan , Feb , Mar , Apr , May , Jun , Jul , Aug , Sep , Oct , Nov , Dec ]: ")
        endmonth = np.where(np.array(['Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']) == answ)[0][0]+1
        fcmonth = (endmonth-startmonth)+1 if (endmonth-startmonth)>=0 else (endmonth-startmonth)+13
    else:
        answ = input("Resultados para el trimestre [ NDJ , DJF , JFM , FMA , MAM , AMJ , MJJ , JJA , JAS , ASO , SON , OND ]: ")
        endmonth = np.where(np.array(['NDJ', 'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND']) == answ)[0][0]+1
        fcmonth = (endmonth-startmonth)+1 if (endmonth-startmonth)>=0 else (endmonth-startmonth)+13

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
)

# Directory selection
load_dotenv() # Data is saved in a path defined in file .env
DATADIR = os.getenv('DATA_DIR')
MODESDIR = DATADIR + '/modes'
SCOREDIR = DATADIR + '/scores'
PLOTSDIR = f'./plots/stmonth{config["start_month"]:02d}'
# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)
# Read hindcast file
hcpcs_fname = f'{MODESDIR}/{hcst_bname}.{aggr}.PCs.nc'
# Base name for observations
obs_bname = 'era5_monthly_stmonth{start_month:02d}_{hcstarty}-{hcendy}'.format(**config)
# Read observations file
obpcs_fname = f'{MODESDIR}/{obs_bname}.{aggr}.PCs.nc'

#%% [markdown]

# ## 4.1 Time series
#
# Boxplots represent the forecasted modes, whereas background solid line represents observed values.
#
# A text box is included indicating the three verification scores (Spearman's rank correlation, area under Relative Operating Characteristic (ROC) curve and Ranked Probability Skill Score (RPSS)) 

#%%
print("4.1 Time series")  

# Some predefined options to plot each score
score_options = {'bs': [np.linspace(0.,0.5,11), plt.colormaps['YlGn'], 3, 'max', 'Brier Score (BS)'],
                 'corr': [np.linspace(-1.,1.,11), plt.colormaps['RdYlBu_r'], 1, 'both', 'Spearmans Rank Correlation (stippling where significance below 95%)'],
                 'roc': [np.linspace(0.,1.,11), plt.colormaps['BrBG'], 3, 'both', 'Area under Relative Operating Characteristic (ROC) curve'],
                 'rocss': [np.linspace(-0.5,0.5,9), plt.colormaps['BrBG'], 3, 'both', 'Relative Operating Characteristic Skill Score (ROCSS)'],
                 'rps': [np.linspace(0.3,0.5,11), plt.colormaps['YlGn_r'], 1, 'max', 'Ranked Probability Score (RPS)'],
                 'rpss': [np.linspace(-0.5,0.5,11), plt.colormaps['BrBG'], 1, 'both', 'Ranked Probability Skill Score (RPSS)'],
                }

# Change titles font size
plt.rc('axes', titlesize=20)

# Prepare strings for titles
locale.setlocale(locale.LC_ALL, 'en_GB')
tit_line1 = '{institution} {name}'.format(**origin_labels)
tit_line2_base = f'Start month: {calendar.month_abbr[config["start_month"]].upper()}'
if aggr=='1m':
    validmonth = config['start_month'] + (fcmonth-1)
    validmonth = validmonth if validmonth<=12 else validmonth-12
    tit_line2 = tit_line2_base + f' - Valid month: {calendar.month_abbr[validmonth].upper()}'
elif aggr=='3m':
    validmonths = [vm if vm<=12 else vm-12 for vm in [config['start_month'] + (fcmonth-1) - shift for shift in range(3)]]
    validmonths = [calendar.month_abbr[vm][0] for vm in reversed(validmonths)]
    tit_line2 = tit_line2_base + f' - Valid months: {"".join(validmonths)}'
else:
    raise BaseException(f'Unexpected aggregation {aggr}')

# Reading HCST data from file
hcpcs = xr.open_dataset(hcpcs_fname)
pcs_hcst = hcpcs.sel(forecastMonth=fcmonth).swap_dims({'start_date':'valid_time'})
# Reading OBS data from file
obpcs = xr.open_dataset(obpcs_fname)
pcs_obs = obpcs.where(obpcs.valid_time==pcs_hcst.valid_time,drop=True)

# Mode names
name_pcs = ['PC1 (NAO)', 'PC2 (EA)', 'PC3 (EAWR)', 'PC4 (SCA)']
name_pcs2 = ['nao_eof', 'ea', 'eawr', 'sca']

# For each mode
for m in range(4):
    # Create figure
    fig = plt.figure(figsize=(12,6))
    ax = fig.add_subplot()

    # Represent
    thispc_obs = pcs_obs.pseudo_pcs.sel(mode=m)
    thispc_hcst = pcs_hcst.pseudo_pcs.sel(mode=m)
    years = thispc_obs.valid_time.dt.year
    plt.plot(years, thispc_obs, linewidth=2, color='k')
    medianprops = dict(linestyle='-', linewidth=1., color='w')
    meanprops = dict(marker='o', markeredgecolor='black', markerfacecolor='w')
    bplot = ax.boxplot([thispc_hcst.isel(valid_time=y).squeeze().values for y in range(years.size)], positions=[years[y].values for y in range(years.size)], widths=0.7, vert=True, showfliers=False, patch_artist=True, notch=False, meanprops=meanprops, medianprops=medianprops, showmeans=True)
    for y in range(years.size):
        if thispc_hcst.isel(valid_time=y).mean().values<-0.25:
            bplot['boxes'][y].set_facecolor('lightblue')
        elif thispc_hcst.isel(valid_time=y).mean().values>0.25:
            bplot['boxes'][y].set_facecolor('firebrick')
        else:
            bplot['boxes'][y].set_facecolor('lightgrey')
    plt.fill_between(years, 0, thispc_obs.squeeze(), where=(thispc_obs.squeeze() >= 0), color='firebrick',
                    interpolate=True)
    plt.fill_between(years, thispc_obs.squeeze(), 0, where=(thispc_obs.squeeze() <= 0), color='lightblue',
                    interpolate=True)
    plt.axhline(0, color='k')
    plt.xlabel('Year')
    plt.ylabel(name_pcs[m])
    plt.ylim([-4., 4.])
    ax.set_xticks(years.values[::2])    
    ax.set_xticklabels(years.values[::2])    
    ax.grid(True)
    # Scores
    corr = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.corr.nc').pseudo_pcs.sel(mode=m,forecastMonth=fcmonth).values)
    roc = [float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.roc.nc').pseudo_pcs.sel(mode=m,category=0,forecastMonth=fcmonth).values),
        float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.roc.nc').pseudo_pcs.sel(mode=m,category=1,forecastMonth=fcmonth).values),
        float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.roc.nc').pseudo_pcs.sel(mode=m,category=2,forecastMonth=fcmonth).values)]
    rpss = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.rpss.nc').pseudo_pcs.sel(mode=m,forecastMonth=fcmonth).values)
    textstr= '\n'.join((
        r'Correlation={:.2f}'.format(corr),
        r'ROC=[{:.2f},{:.2f},{:.2f}]'.format(roc[0],roc[1],roc[2]),
        r'RPSS={:.2f}'.format(rpss)))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)

    # Titles
    plt.title(tit_line1+' '+tit_line2, fontsize=12)
    ax.text(0.05, 0.95, textstr, transform=ax.transAxes, fontsize=12,
            verticalalignment='top', bbox=props)      
    fig.suptitle(name_pcs[m]+f' - Forecasts (boxplots) and observations (shading)\n', fontsize=16)
    plt.subplots_adjust(left=0.05, right=0.95, top=0.9, bottom=0.1, hspace=0.25, wspace=0.15)

    # Save figure
    figname = f'./{PLOTSDIR}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{name_pcs2[m]}.png'
    fig.savefig(figname,dpi=600)  
