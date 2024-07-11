# %% [markdown]
# # 7. Multi-season score cards

# This script is used to represent score-cards, 
# which show the verification score values for forecasts of a specific month or season with different systems and initializations 
#
# First we have to decide a month aggregation, an end month and a verification score. 

#%%
print("7. Multi-season score cards")

import os
import sys
from dotenv import load_dotenv
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
    score = str(sys.argv[2])
# If no variables were introduced, ask for them
else:
    # Subset of plots to be produced
    aggr = input("Selecciona el tipo de agregación mensual [ 1m , 3m ]: ")

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
lts=3

# Directory selection
load_dotenv() # Data is saved in a path defined in file .env
DATADIR = os.getenv('DATA_DIR')

var = 'NAO'
# Common labels to be used in plot titles
VARNAMES = {
    'NAO' : 'North Atlantic Oscillation (box)',
    # 'eof1' : 'EOF 1',    
    # 'eof2' : 'EOF 2',
    # 'eof3' : 'EOF 3',
    # 'eof4' : 'EOF 4',
}

if aggr=='1m':
    var_seasons = ["Jan", "Feb", "Mar", "Apr", "May", "Jun", "Jul", "Aug", "Sep", "Oct", "Nov", "Dec"]
elif aggr=='3m': 
    var_seasons = ["JFM", "FMA", "MAM", "AMJ", "MJJ", "JJA", "JAS", "ASO", "SON", "OND", "NDJ", "DJF"]

# Some predefined options to plot each score
score_options = {'bs': [np.linspace(0.,0.5,11), plt.colormaps['YlGn'], 3, 'max', 'Brier Score (BS)'],
                 'corr': [np.linspace(-1.,1.,11), plt.colormaps['RdYlBu_r'], 1, 'both', 'Spearmans Rank Correlation'],
                 'roc': [np.linspace(0.,1.,11), plt.colormaps['BrBG'], 3, 'both', 'Area under Relative Operating Characteristic (ROC) curve'],
                 'rocss': [np.linspace(-0.5,0.5,9), plt.colormaps['BrBG'], 3, 'both', 'Relative Operating Characteristic Skill Score (ROCSS)'],
                 'rps': [np.linspace(0.3,0.5,11), plt.colormaps['YlGn_r'], 1, 'max', 'Ranked Probability Score (RPS)'],
                 'rpss': [np.linspace(-0.5,0.5,11), plt.colormaps['BrBG'], 1, 'both', 'Ranked Probability Skill Score (RPSS)'],
                }

# %% [markdown]
# ## 7.1 Array construction

# First we construct an array with all results for the selected score, 
# considering different initializations and forecasting systems.

#%%
print("7.1 Array construction")

# Size of solution's array
n_models = len(full_name)
n_init = lts
n_seasons = len(var_seasons)
n_terciles = score_options[score][2]
DATA = np.zeros([n_seasons*n_terciles,n_models*n_init])*np.nan

m = 0
# For each model
for label in full_name:
    # Save the simplier model and system name
    model = full_name[label][0]
    system = full_name[label][1]

    s = 0
    # For each season
    for season in var_seasons:
        # List of initializations
        if aggr=='1m':
            endmonth = np.where(np.array(['Jan', 'Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']) == season)[0][0]+1
            initialization = {calendar.month_name[endmonth-(l+1) if endmonth-(l+1)>0 else endmonth-(l+1)+12]: [endmonth-(l+1) if endmonth-(l+1)>0 else endmonth-(l+1)+12, l+2] for l in reversed(range(lts))}
        elif aggr=='3m':
            endmonth = np.where(np.array(['NDJ', 'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND']) == season)[0][0]+1
            initialization = {calendar.month_name[endmonth-(l+3) if endmonth-(l+3)>0 else endmonth-(l+3)+12]: [endmonth-(l+3) if endmonth-(l+3)>0 else endmonth-(l+3)+12, l+4] for l in reversed(range(lts))}
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
            fname = '{datadir}/scores/{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly.{aggr}.{score}.nc'.format(datadir=DATADIR,**config)
            scorefile = xr.open_dataset(fname).sel(forecastMonth=fcmonth)
            
            # For each tercile
            for t in range(n_terciles):
                if n_terciles>1:
                    # Filling the array
                    if var=='NAO':
                        DATA[t+s*n_terciles,m*n_init+i] = float(scorefile['nao'].sel(category=t).values)
                        print('Model {}, Variable: {}, Initialization {}, Tercile {}:'.format(model,var,init,t))
                        print(float(scorefile['nao'].sel(category=t).values))
                    else:
                        DATA[t+s*n_terciles,m*n_init+i] = float(scorefile['pseudo_pcs'].sel(category=t,mode=int(var[-1])-1).values)
                        print('Model {}, Variable: {}, Initialization {}, Tercile {}:'.format(model,var,init,t))
                        print(float(scorefile['pseudo_pcs'].sel(category=t,mode=int(var[-1])-1).values))
                else:
                    # Filling the array
                    if var=='NAO':
                        DATA[t+s*n_terciles,m*n_init+i] = float(scorefile['nao'].values)
                        print('Model {}, Variable: {}, Initialization {}:'.format(model,var,init))
                        print(float(scorefile['nao'].values))
                    else:
                        DATA[t+s*n_terciles,m*n_init+i] = float(scorefile['pseudo_pcs'].sel(mode=int(var[-1])-1).values)
                        print('Model {}, Variable: {}, Initialization {}:'.format(model,var,init))
                        print(float(scorefile['pseudo_pcs'].sel(mode=int(var[-1])-1).values))
            i+=1
        s+=1
    m+=1

# %% [markdown]
# ## 7.2 Score-cards

# Then we represent the results.

#%%
print("7.2 Score-cards")

# Directory creation
PLOTSDIR = f'./plots/scorecards'
# Check if the directory exists
if not os.path.exists(PLOTSDIR):
    # If it doesn't exist, create it
    try:
        os.makedirs(PLOTSDIR)
    except FileExistsError:
        pass

# Prepare strings for titles
locale.setlocale(locale.LC_ALL, 'en_GB')
# if aggr=='1m':
#     validmonth = config['start_month'] + (fcmonth-1)
#     validmonth = validmonth if validmonth<=12 else validmonth-12
# #    tit_line2 = tit_line2_base + rf' - Valid month: $\bf{calendar.month_abbr[validmonth].upper()}$'
#     tit_line2 = f'Valid month: {calendar.month_abbr[validmonth].upper()}'
# elif aggr=='3m':
#     validmonths = [vm if vm<=12 else vm-12 for vm in [config['start_month'] + (fcmonth-1) - shift for shift in range(3)]]
#     validmonths = [calendar.month_abbr[vm][0] for vm in reversed(validmonths)]
# #    tit_line2 = tit_line2_base + rf' - Valid months: $\bf{"".join(validmonths)}$'
#     tit_line2 = f'Valid months: {"".join(validmonths).upper()}'
# else:
#     raise BaseException(f'Unexpected aggregation {aggr}')
tit_line1 = f'{VARNAMES[var]}'+f'\n{score_options[score][4]}'
figname = f'{PLOTSDIR}/Score-card_{score}_{aggr}_{var}.png'

# Create figure
fig = plt.figure(figsize=(17,8))
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
ax.set_xticks(x[::n_terciles]+0.5*n_terciles,[tick for tick in var_seasons])
init_labels = ['lead3', 'lead2', 'lead1'] *n_models
ax.set_yticks(y+0.499,init_labels,minor=True)
if n_terciles>1:
    terciles_labels = ['lower', 'middle', 'upper'] * n_seasons
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
    df = df.set_axis(pd.MultiIndex.from_product([var_seasons,['lower', 'middle', 'upper']]), axis=1)
else:
    df = df.set_axis(var_seasons, axis=1)
df = df.set_axis(pd.MultiIndex.from_product([[label for label in full_name],['lead3', 'lead2', 'lead1']]), axis=0)
csvname = f'{PLOTSDIR}/Score-card_{score}_{aggr}_{var}.csv'
df.to_csv(csvname)