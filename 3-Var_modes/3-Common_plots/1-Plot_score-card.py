# %% [markdown]
# # 1. Plot Score-cards

# This script is used to represent score-cards, 
# which show the verification score values for forecasts of a specific month or season with different systems and initializations 
# 
# The analised variables are: NAO (from the box method), NAO (from eofs), EA, EA/WR and SCA.
#
# First we have to decide a month aggregation, an end month and a verification score. 

#%%
print("1. Plot Score-cards")

import os
import sys
import xarray as xr
import numpy as np
import locale
import calendar
from matplotlib import pyplot as plt
from matplotlib.colors import BoundaryNorm
import warnings
warnings.filterwarnings('ignore')

# If variables are introduced from the command line, read them
if len(sys.argv) > 2:
    aggr = str(sys.argv[1])
    endmonth = int(sys.argv[2])
    score = str(sys.argv[3])
# If no variables were introduced, ask for them
else:
    # Subset of plots to be produced
    aggr = input("Selecciona el tipo de agregación mensual [ 1m , 3m ]: ")

    # Forecast month
    if aggr=='1m':
        answ = input("Resultados para el mes [ Jan , Feb , Mar , Apr , May , Jun , Jul , Aug , Sep , Oct , Nov , Dec ]: ")
        endmonth = np.where(np.array(['Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']) == answ)[0][0]+1
    else:
        answ = input("Resultados para el trimestre [ NDJ , DJF , JFM , FMA , MAM , AMJ , MJJ , JJA , JAS , ASO , SON , OND ]: ")
        endmonth = np.where(np.array(['NDJ', 'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND']) == answ)[0][0]+1

    # Verification score
    score = input("Usar el siguiente score [ bs , corr , roc , rocss , rps , rpss ]: ")

# Dictionary to link full system names and simplier names
full_name = {#'ECMWF-System 4': ['ecmwf','4'],
             #'ECMWF-SEAS5': ['ecmwf', '5'],
             'ECMWF-SEAS5.1': ['ecmwf', '51'],
             #'Météo France-System 5': [ 'meteo_france', '5'],
             #'Météo France-System 6': [ 'meteo_france', '6'],
             #'Météo France-System 7': [ 'meteo_france', '7'],
             'Météo France-System 8': [ 'meteo_france', '8'],
             #'Met Office-System 12': ['ukmo', '12'],
             #'Met Office-System 13': ['ukmo', '13'],
             #'Met Office-System 14': ['ukmo', '14'],
             #'Met Office-System 15': ['ukmo', '15'],
             #'Met Office-GloSea6': ['ukmo', '600'],
             #'Met Office-GloSea6.1': ['ukmo', '601'],
             'Met Office-GloSea6.2': ['ukmo', '602'],
             #'Met Office-GloSea6.3': ['ukmo', '603'],
             #'DWD-GCFS2.0': ['dwd', '2'],
             'DWD-GCFS2.1': ['dwd', '21'],
             #'CMCC-SPSv3.0': ['cmcc', '3'],
             'CMCC-SPSv3.5': ['cmcc', '35'],
             'NCEP-CFSv2': ['ncep', '2'],
             #'JMA-CPS2': ['jma', '2'],
             'JMA-CPS3': ['jma', '3'],
             #'ECCC-GEM-NEMO': ['eccc', '1'],
             #'ECCC-CanCM4i': ['eccc', '2'],
             'ECCC-GEM5-NEMO': ['eccc', '3']
            }

# Number of lead-times
lts=4
# List of initializations
if aggr=='1m':
    initialization = {calendar.month_name[endmonth-(l+1) if endmonth-(l+1)>0 else endmonth-(l+1)+12]: [endmonth-(l+1) if endmonth-(l+1)>0 else endmonth-(l+1)+12, l+2] for l in reversed(range(lts))}
elif aggr=='3m': 
    initialization = {calendar.month_name[endmonth-(l+3) if endmonth-(l+3)>0 else endmonth-(l+3)+12]: [endmonth-(l+3) if endmonth-(l+3)>0 else endmonth-(l+3)+12, l+4] for l in reversed(range(lts))}


# Common labels to be used in plot titles
VARNAMES = {
    'NAO (box)' : 'North Atlantic Oscillation (box)',
    'NAO (eof)'  : 'North Atlantic Oscillation (eof)',    
    'EA' : 'East Atlantic',
    'EAWR' : 'East Atlantic - Western Russia',
    'SCA' : 'Scandinavian',
}

# Some predefined options to plot each score
score_options = {'bs': [np.linspace(0.,0.5,11), plt.colormaps['YlGn'], 3, 'max', 'Brier Score (BS)'],
                 'corr': [np.linspace(-1.,1.,11), plt.colormaps['RdYlBu_r'], 1, 'both', 'Spearmans Rank Correlation (stippling where significance below 95%)'],
                 'roc': [np.linspace(0.,1.,11), plt.colormaps['BrBG'], 3, 'both', 'Area under Relative Operating Characteristic (ROC) curve'],
                 'rocss': [np.linspace(-0.5,0.5,9), plt.colormaps['BrBG'], 3, 'both', 'Relative Operating Characteristic Skill Score (ROCSS)'],
                 'rps': [np.linspace(0.3,0.5,11), plt.colormaps['YlGn_r'], 1, 'max', 'Ranked Probability Score (RPS)'],
                 'rpss': [np.linspace(-0.5,0.5,11), plt.colormaps['BrBG'], 1, 'both', 'Ranked Probability Skill Score (RPSS)'],
                }

# %% [markdown]
# ## 1.1 Array construction

# First we construct an array with all results for the selected score, 
# considering different initializations and forecasting systems.

#%%
print("1.1 Array construction")

# Size of solution's array
n_models = len(full_name)
n_init = len(initialization)
n_vars = len(VARNAMES)
n_terciles = score_options[score][2]
DATA = np.zeros([n_vars*n_terciles,n_models*n_init])*np.nan

m = 0
# For each model
for label in full_name:
    # Save the simplier model and system name
    model = full_name[label][0]
    system = full_name[label][1]

    i = 0
    # For each initialization
    for init in initialization:
        startmonth = initialization[init][0]
        fcmonth = initialization[init][1]
        # Here we save the configuration
        config = dict(
            hcstarty = 1993,
            hcendy = 2016,
            start_month = startmonth,
            origin = model,
            system = system,
            aggr = aggr,
            score = score,
        )
        # Read files
        if fcmonth>6:
            fname1 = '../1-Box_calc/daily_data_extension/data/scores/{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_daily.{aggr}.{score}.nc'.format(**config)
            fname2 = '../2-EOFs_calc/daily_data_extension/data/scores/{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_daily.{aggr}.{score}.nc'.format(**config)
            try:
                score1 = xr.open_dataset(fname1).sel(forecastMonth=fcmonth)
                score2 = xr.open_dataset(fname2).sel(forecastMonth=fcmonth)
            except:
                i+=1
                continue
        else:
            fname1 = '../1-Box_calc/data/scores/{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly.{aggr}.{score}.nc'.format(**config)
            fname2 = '../2-EOFs_calc/data/scores/{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly.{aggr}.{score}.nc'.format(**config)
            score1 = xr.open_dataset(fname1).sel(forecastMonth=fcmonth)
            score2 = xr.open_dataset(fname2).sel(forecastMonth=fcmonth)
        
        # For each tercile
        for t in range(n_terciles):
            if n_terciles>1:
                # Filling the array
                DATA[t+0*n_terciles,m*n_init+i] = float(score1['nao'].sel(category=t).values)
                DATA[t+1*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(category=t,mode=0).values)
                DATA[t+2*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(category=t,mode=1).values)
                DATA[t+3*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(category=t,mode=2).values)
                DATA[t+4*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(category=t,mode=3).values)
                print('Model {}, Variable NAO (box), Initialization {}, Tercile {}:'.format(model,init,t))
                print(float(score1['nao'].sel(category=t).values))
                print('Model {}, Variable NAO (eof), Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(category=t,mode=0).values))
                print('Model {}, Variable EA, Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(category=t,mode=1).values))
                print('Model {}, Variable EAWR, Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(category=t,mode=2).values))
                print('Model {}, Variable SCA, Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(category=t,mode=3).values))          
            else:
                # Filling the array
                DATA[t+0*n_terciles,m*n_init+i] = float(score1['nao'].values)
                DATA[t+1*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(mode=0).values)
                DATA[t+2*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(mode=1).values)
                DATA[t+3*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(mode=2).values)
                DATA[t+4*n_terciles,m*n_init+i] = float(score2['pseudo_pcs'].sel(mode=3).values)
                print('Model {}, Variable NAO (box), Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score1['nao'].values))
                print('Model {}, Variable NAO (eof), Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(mode=0).values))
                print('Model {}, Variable EA, Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(mode=1).values))
                print('Model {}, Variable EAWR, Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(mode=2).values))
                print('Model {}, Variable SCA, Inicializacion {}, Tercile {}:'.format(model,init,t))
                print(float(score2['pseudo_pcs'].sel(mode=3).values))      
        i+=1
    m+=1

# %% [markdown]
# ## 1.2 Score-cards

# Then we represent the results.

#%%
print("1.2 Score-cards")

# Prepare strings for titles
locale.setlocale(locale.LC_ALL, 'en_GB')
if aggr=='1m':
    validmonth = config['start_month'] + (fcmonth-1)
    validmonth = validmonth if validmonth<=12 else validmonth-12
#    tit_line2 = tit_line2_base + rf' - Valid month: $\bf{calendar.month_abbr[validmonth].upper()}$'
    tit_line2 = f'Valid month: {calendar.month_abbr[validmonth].upper()}'
elif aggr=='3m':
    validmonths = [vm if vm<=12 else vm-12 for vm in [config['start_month'] + (fcmonth-1) - shift for shift in range(3)]]
    validmonths = [calendar.month_abbr[vm][0] for vm in reversed(validmonths)]
#    tit_line2 = tit_line2_base + rf' - Valid months: $\bf{"".join(validmonths)}$'
    tit_line2 = f'Valid months: {"".join(validmonths).upper()}'
else:
    raise BaseException(f'Unexpected aggregation {aggr}')
tit_line1 = 'Euro-Atlantic teleconnections - '+tit_line2+f'\n{score_options[score][4]}'
figname = f'./Score-card_{score}_{"".join(validmonths)}_{aggr}.png'

# Create figure
fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot()

# Representation
x = np.arange(DATA.shape[0])
y = np.arange(DATA.shape[1])
levels = score_options[score][0]
cmap = score_options[score][1]
norm = BoundaryNorm(levels, ncolors=cmap.N, clip=True)
c = ax.pcolor(DATA.T, edgecolors='w', linewidths=4, cmap=cmap, norm=norm)
cb = plt.colorbar(c,
                orientation='horizontal',
                location='bottom',
                aspect=60,
                pad=0.05,
                label=score.upper())

# Labels and title
ax.set_yticks(y[::n_init]+0.5*n_init,[label for label in full_name])
ax.set_xticks(x[::n_terciles]+0.5*n_terciles,[tick for tick in VARNAMES])
init_labels = [init for init in initialization] *n_models
ax.set_yticks(y+0.499,init_labels,minor=True)
if n_terciles>1:
    terciles_labels = ['lower', 'middle', 'upper'] * n_vars
    ax.set_xticks(x+0.499,terciles_labels,minor=True)
ax.tick_params(axis='x', which="minor", 
               bottom=False, left=False, top=True, right=True,
               labelbottom=False, labelleft=False, labeltop=True, labelright=True,
               labelrotation=25.)
ax.tick_params(axis='y', which="minor", 
               bottom=False, left=False, top=True, right=True,
               labelbottom=False, labelleft=False, labeltop=True, labelright=True,
               labelrotation=0.)
plt.hlines(y[::n_init],0,len(x),color='k')
plt.vlines(x[::n_terciles],0,len(y),color='k')
plt.title(tit_line1, fontsize=14, loc='center')

# Numbers inside box
for y in range(DATA.shape[1]):
    for x in range(DATA.shape[0]):
        if (DATA[x,y]>levels[-2]) | (DATA[x,y]<=levels[1]):
            plt.text(x + 0.5, y + 0.5, '%.2f' % DATA[x,y],
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='w'
                )
        else:
            plt.text(x + 0.5, y + 0.5, '%.2f' % DATA[x,y],
                        horizontalalignment='center',
                        verticalalignment='center',
                        color='k'
                )
        
# Save figure
plt.subplots_adjust(left=0.2, right=0.9, top=0.85, bottom=0.05)
fig.savefig(figname,dpi=600,bbox_inches='tight')  

# Create csv
import pandas as pd
df = pd.DataFrame(DATA.T)
if n_terciles>1:
    df = df.set_axis(pd.MultiIndex.from_product([[var for var in VARNAMES],['lower', 'middle', 'upper']]), axis=1)
else:
    df = df.set_axis([var for var in VARNAMES], axis=1)
df = df.set_axis(pd.MultiIndex.from_product([[label for label in full_name],[init for init in initialization]]), axis=0)
csvname = f'./Score-card_{score}_{"".join(validmonths)}_{aggr}.csv'
df.to_csv(csvname)