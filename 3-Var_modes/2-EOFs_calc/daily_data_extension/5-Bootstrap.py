# %% [markdown]
# # 5. Compute score uncertainties with bootstrap

# This script is used to compute uncertainties of different verification scores
# for daily seasonal forescasts of climate variability modes.
# 
# The computed scores are: Spearman's rank correlation, area under Relative Operating Characteristic (ROC) curve, 
# Relative Operating Characteristic Skill Score (ROCSS), Ranked Probability Score (RPS), Ranked Probability Skill Score (RPSS) and Brier Score (BS).
#
# The uncertainty is computed through a bootstrap test. 
#
# First we have to decide a forecast system (institution and system name), a start month a month aggregation and a end month. 

#%%
print("5. Compute score uncertainties with bootstrap")

import os
import sys
import xarray as xr
import numpy as np
import pandas as pd
import xskillscore as xs
import locale
import calendar
import matplotlib.pyplot as plt
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
    isLagged = False if model in ['ecmwf', 'meteo_france', 'dwd', 'cmcc', 'eccc'] else True
)

# Directory selection
DATADIR = os.environ['MAS'] + '/Seasonal_Verification/3-Var_modes/2-EOFs_calc/daily_data_extension/data'
MODESDIR = DATADIR + '/modes'
SCOREDIR = DATADIR + '/scores'
PLOTSDIR = f'./plots/stmonth{config["start_month"]:02d}'
# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_daily'.format(**config)
# File name for hindcast
hcpcs_fname = f'{MODESDIR}/{hcst_bname}.1m.PCs.nc'
hcpcs_3m_fname = f'{MODESDIR}/{hcst_bname}.3m.PCs.nc'
# Base name for observations
obs_bname = 'era5_daily_stmonth{start_month:02d}_{hcstarty}-{hcendy}'.format(**config)
# File name for observations
obpcs_fname = f'{MODESDIR}/{obs_bname}.1m.PCs.nc'
obpcs_3m_fname = f'{MODESDIR}/{obs_bname}.3m.PCs.nc'

# Check if files exist
if not os.path.exists(obpcs_fname) & os.path.exists(obpcs_3m_fname):
    print('No se calcularon aún los índices PCs de ERA5')
    sys.exit()
elif not os.path.exists(hcpcs_fname) & os.path.exists(hcpcs_3m_fname):
    print('No se calcularon aún los índices PCs de este modelo y sistema')
    sys.exit()


# %% [markdown]
# ## 5.1 Hindcast iteration over time

# Here we resample hindcast data over time, with n iterations. 

# %% 
print("5.1 Hindcast iteration over time")

# Reading HCST and OBS data from file
if aggr=='1m':
    hcpcs = xr.open_dataset(hcpcs_fname)
    obpcs = xr.open_dataset(obpcs_fname)
elif aggr=='3m':
    hcpcs = xr.open_dataset(hcpcs_3m_fname)
    obpcs = xr.open_dataset(obpcs_3m_fname)
else:
    raise BaseException(f'Unexpected aggregation {aggr}')

# We have to re-arange the hindcast and observed datasets and merge them into one single dataset
hcpcs = hcpcs.sel(forecastMonth=fcmonth).swap_dims({'start_date':'valid_time'})
obpcs = obpcs.where(obpcs.valid_time==hcpcs.valid_time,drop=True)
pcs_hcob = xr.merge([hcpcs, obpcs.rename({"pseudo_pcs": "pseudo_pcs_obs"})])

# Resample hindcast data over time, with n iterations
n = 10000
pcs_boot = xs.resample_iterations_idx(pcs_hcob, n, 'valid_time')

h = pcs_boot[['pseudo_pcs']]
o = pcs_boot[['pseudo_pcs_obs']].rename({"pseudo_pcs_obs": "pseudo_pcs"})

# %% [markdown]
# ## 5.2 Probabilities for tercile categories

# Here we get the probabilities for tercile categories of the hindcast data, 
# by counting the number of ensemble members found in each tercile.

# %% 
print("5.2 Probabilities for tercile categories")

# We define a function to calculate the boundaries of forecast categories defined by quantiles
def get_thresh(icat,quantiles,xrds,dims=['number','valid_time']):

    if not all(elem in xrds.dims for elem in dims):           
        raise Exception('Some of the dimensions in {} is not present in the xr.Dataset {}'.format(dims,xrds)) 
    else:
        if icat == 0:
            xrds_lo = -np.inf
            xrds_hi = xrds.quantile(quantiles[icat],dim=dims,skipna=True)      
            
        elif icat == len(quantiles):
            xrds_lo = xrds.quantile(quantiles[icat-1],dim=dims,skipna=True)
            xrds_hi = np.inf
            
        else:
            xrds_lo = xrds.quantile(quantiles[icat-1],dim=dims,skipna=True)
            xrds_hi = xrds.quantile(quantiles[icat],dim=dims,skipna=True)
      
    return xrds_lo,xrds_hi

# Calculate probabilities for tercile categories by counting members within each category
quantiles = [1/3., 2/3.]
numcategories = len(quantiles)+1

l_probs_hcst=list()
# For each quantile
for icat in range(numcategories):
    # Get the lower and higher threshold
    h_lo,h_hi = get_thresh(icat, quantiles, h)
    # Count the number of member between the threshold
    probh = np.logical_and(h>h_lo, h<=h_hi).sum('number')/float(h.number.size)
    # Instead of using the coordinate 'quantile' coming from the hindcast xr.Dataset
    # we will create a new coordinate called 'category'
    if 'quantile' in list(probh.coords):
        probh = probh.drop('quantile')
    l_probs_hcst.append(probh.assign_coords({'category':icat}))

    # Concatenating tercile probs categories
    probs = xr.concat(l_probs_hcst,dim='category')                    

# %% [markdown]
# ## 5.3 Compute deterministic scores

# Here we calculate the Spearman's rank correlation and thei p-values. 
# 
# This score is based on the ensemble mean, not on the probabilities for each tercile.

# %% 
print("5.3 Compute deterministic scores")

# Check if hindcast data is ensemble
is_fullensemble = 'number' in h.dims

# Compute ensemble mean (if data is an ensemble)
thishcst_em = h if not is_fullensemble else h.mean('number')
# Calculate Spearman's rank correlation
corr = xs.spearman_r(thishcst_em, o, dim='valid_time')

# %% [markdown]
# ## 5.4 Compute probabilistic scores for tercile categories

# Here we calculate the probabilistic scores: area under Relative Operating Characteristic (ROC) curve, 
# Relative Operating Characteristic Skill Score (ROCSS), Ranked Probability Score (RPS), Ranked Probability Skill Score (RPSS) and Brier Score (BS). 

# %% 
print("5.4 Compute probabilistic scores for tercile categories")

thishcst = probs.copy()

l_probs_obs=list()
l_probs_clim=list()
# For each quantile
for icat in range(numcategories):
    # Get the lower and higher threshold
    o_lo,o_hi = get_thresh(icat, quantiles, o, dims=['valid_time'])
    # Count the number of "members" between the threshold (1 or 0)
    probo = 1. * np.logical_and(o>o_lo, o<=o_hi)
    if 'quantile' in list(probo.coords):
        probo=probo.drop('quantile')
    l_probs_obs.append(probo.assign_coords({'category':icat}))
    # Count the number of months between the threshold (1 or 0)
    probc = np.logical_and(o>o_lo, o<=o_hi).sum('valid_time')/float(o.valid_time.size)        
    try:
        probc=probc.drop('quantile')
    except:
        pass
    l_probs_clim.append(probc.assign_coords({'category':icat}))
# Concatenate observations and climatology probabilities
thisobs = xr.concat(l_probs_obs, dim='category')
thisclim = xr.concat(l_probs_clim, dim='category')
# Here we have to correct climatological probabilities.
# Bootstrap method consist on repeating values over time,
# so the computed percentiles are not entirely "correct" because some values are exactly the same.
# This results in some climatological probabilities for the tercile categories different from 1/3
thisclim = thisclim*0.+1./3.

# Calculate the probabilistic (tercile categories) scores
roc = xr.Dataset()
rps = xr.Dataset()
rpsclim = xr.Dataset()
rpss = xr.Dataset()
rocss = xr.Dataset()
bs = xr.Dataset()
# For each variable
for var in hcpcs.data_vars:
    # Compute Area ROC
    roc[var] = xs.roc(thisobs[var],thishcst[var], dim='valid_time', bin_edges=np.linspace(0,1,101))
    # Compute RPS
    rps[var] = xs.rps(thisobs[var],thishcst[var], dim='valid_time', category_edges=None, input_distributions='p')
    # Compute climatological RPS           
    rpsclim[var] = xs.rps(thisobs[var],thisclim[var], dim='valid_time', category_edges=None, input_distributions='p')
    # Compute RPSS           
    rpss[var] = 1.-rps[var]/rpsclim[var]
    # Compute ROCSS
    rocss[var] = (roc[var] - 0.5) / (1. - 0.5)
    # Compute Brier Score
    bscat = list()
    for cat in thisobs[var].category:
        thisobscat = thisobs[var].sel(category=cat)
        thishcstcat = thishcst[var].sel(category=cat)
        bscat.append(xs.brier_score(thisobscat, thishcstcat, dim='valid_time'))
    bs[var] = xr.concat(bscat,dim='category')

corr_ci = corr.quantile([0.025, 0.975], dim="iteration")
rps_ci = rps.quantile([0.025, 0.975], dim="iteration")
rpss_ci = rpss.quantile([0.025, 0.975], dim="iteration")
bs_ci = bs.quantile([0.025, 0.975], dim="iteration")
roc_ci = roc.quantile([0.025, 0.975], dim="iteration")
rocss_ci = rocss.quantile([0.025, 0.975], dim="iteration")

#%% [markdown]

# ## 5.5 Histogram

# We represent the probability density function of the bootstrap test for each score.

#%%
print("5.5 Histograms")  

# Some predefined options to plot each score
score_options = {'bs': [np.linspace(0.,0.5,11), plt.colormaps['YlGn'], 3, 'max', 'Brier Score (BS)'],
                 'corr': [np.linspace(-1.,1.,11), plt.colormaps['RdYlBu_r'], 1, 'both', 'Spearmans Rank Correlation'],
                 'roc': [np.linspace(0.,1.,11), plt.colormaps['BrBG'], 3, 'both', 'Area under Relative Operating Characteristic (ROC) curve'],
                 'rocss': [np.linspace(-0.5,0.5,9), plt.colormaps['BrBG'], 3, 'both', 'Relative Operating Characteristic Skill Score (ROCSS)'],
                 'rps': [np.linspace(0.3,0.5,11), plt.colormaps['YlGn_r'], 1, 'max', 'Ranked Probability Score (RPS)'],
                 'rpss': [np.linspace(-0.5,0.5,11), plt.colormaps['BrBG'], 1, 'both', 'Ranked Probability Skill Score (RPSS)'],
                }

# Prepare strings for titles
locale.setlocale(locale.LC_ALL, 'en_GB')
tit_line1_base = '{institution} {name}'.format(**origin_labels)
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

# Mode names
name_pcs = ['PC1 (NAO)', 'PC2 (EA)', 'PC3 (EAWR)', 'PC4 (SCA)']
name_pcs2 = ['nao_eof', 'ea', 'eawr', 'sca']

# For each mode
for m in range(4):
    # Correlation
    score = 'corr'
    tit_line = f'{name_pcs[m]} - Bootstrap: {score_options[score][4]}\n'+tit_line1_base+' - '+tit_line2
    score_values = pd.Series(corr.pseudo_pcs.sel(mode=m).values)
    score_ref = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.{score}.nc').sel(forecastMonth=fcmonth,mode=m).pseudo_pcs.values)
    fig = plt.figure(figsize=(6,6))
    ax1 = fig.add_subplot()
    # Pdfs 
    bins = np.linspace(score_options[score][0][0],score_options[score][0][-1],50)
    n, bins, patch = plt.hist(score_values, bins=bins, density=True, color='lightgrey', alpha=0.7, rwidth=0.85, zorder=-1)
    score_values.plot.density(ax=ax1, color='k', lw=2., label='Bootstrap') # identical to s.plot.kde(...)
    # Boxplot
    height = round(n.max()*1.8)
    medianprops = dict(linestyle='-', linewidth=2., color='k')
    meanprops = dict(marker='o', markeredgecolor='black', markerfacecolor='w')
    bplot = ax1.boxplot(score_values, positions=[height*0.75], widths=[height*0.12], vert=False, showfliers=False, patch_artist=True, notch=True, showmeans=True, meanprops=meanprops, medianprops=medianprops)
    bplot['boxes'][0].set_facecolor('lightgrey')
    ax1.text(score_values.mean(), height*0.85, '{:.3f}'.format(score_values.mean()),family='sans-serif',weight='bold',size=12, horizontalalignment='center', verticalalignment='top')
    # Other things
    ax1.vlines(0., 0., height, color='k', linestyle='-', lw=1)
    ax1.vlines(score_ref, 0., height, color='k', linestyle='--', lw=1.5, alpha=0.6, zorder=0)
    ax1.set_ylim([0.,height])
    ax1.set_xlim([score_options[score][0][0],score_options[score][0][-1]])
    #ax1.set_title('Present - Pre-industial (~1850)  |  Thermodynamic',loc='left',weight='bold',size=8,pad=4.)
    ax1.set_yticks(np.linspace(0.,height,height+1),labels=np.linspace(0.,height,height+1))
    plt.title(tit_line, fontsize=12)
    ax1.set_xlabel(score.upper())
    textstr= '\n'.join((
        r'Reference={:.2f}'.format(score_ref),
        r'Conf. interval=({:.2f},{:.2f})'.format(float(corr_ci.pseudo_pcs.sel(quantile=0.025,mode=m).values),float(corr_ci.pseudo_pcs.sel(quantile=0.975,mode=m).values))))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.01, 0.99, textstr, transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)      
    # Save figure
    figname = f'./{PLOTSDIR}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{name_pcs2[m]}.{score}-bootstrap.png'
    fig.savefig(figname,dpi=600)  

    # RPSS
    score = 'rpss'
    tit_line = f'{name_pcs[m]} - Bootstrap: {score_options[score][4]}\n'+tit_line1_base+' - '+tit_line2
    score_values = pd.Series(rpss.pseudo_pcs.sel(mode=m).values)
    score_ref = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.{score}.nc').sel(forecastMonth=fcmonth,mode=m).pseudo_pcs.values)
    fig = plt.figure(figsize=(6,6))
    ax1 = fig.add_subplot()
    # Pdfs 
    bins = np.linspace(score_options[score][0][0],score_options[score][0][-1],50)
    n, bins, patch = plt.hist(score_values, bins=bins, density=True, color='lightgrey', alpha=0.7, rwidth=0.85, zorder=-1)
    score_values.plot.density(ax=ax1, color='k', lw=2., label='Bootstrap') # identical to s.plot.kde(...)
    # Boxplot
    height = round(n.max()*1.8)
    medianprops = dict(linestyle='-', linewidth=2., color='k')
    meanprops = dict(marker='o', markeredgecolor='black', markerfacecolor='w')
    bplot = ax1.boxplot(score_values, positions=[height*0.75], widths=[height*0.12], vert=False, showfliers=False, patch_artist=True, notch=True, showmeans=True, meanprops=meanprops, medianprops=medianprops)
    bplot['boxes'][0].set_facecolor('lightgrey')
    ax1.text(score_values.mean(), height*0.85, '{:.3f}'.format(score_values.mean()),family='sans-serif',weight='bold',size=12, horizontalalignment='center', verticalalignment='top')
    # Other things
    ax1.vlines(0., 0., height, color='k', linestyle='-', lw=1)
    ax1.vlines(score_ref, 0., height, color='k', linestyle='--', lw=1.5, alpha=0.6, zorder=0)
    ax1.set_ylim([0.,height])
    ax1.set_xlim([score_options[score][0][0],score_options[score][0][-1]])
    ax1.set_yticks(np.linspace(0.,height,height+1),labels=np.linspace(0.,height,height+1))
    plt.title(tit_line, fontsize=12)
    ax1.set_xlabel(score.upper())
    textstr= '\n'.join((
        r'Reference={:.2f}'.format(score_ref),
        r'Conf. interval=({:.2f},{:.2f})'.format(float(rpss_ci.pseudo_pcs.sel(quantile=0.025,mode=m).values),float(rpss_ci.pseudo_pcs.sel(quantile=0.975,mode=m).values))))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.01, 0.99, textstr, transform=ax1.transAxes, fontsize=10,
            verticalalignment='top', bbox=props)      
    # Save figure
    figname = f'./{PLOTSDIR}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{name_pcs2[m]}.{score}-bootstrap.png'
    fig.savefig(figname,dpi=600)  

    # ROC
    score = 'roc'
    tit_line = f'{name_pcs[m]} - Bootstrap: \n{score_options[score][4]}\n'+tit_line1_base+' - '+tit_line2
    score_values0 = pd.Series(roc.sel(category=0,mode=m).pseudo_pcs.values)
    score_values1 = pd.Series(roc.sel(category=1,mode=m).pseudo_pcs.values)
    score_values2 = pd.Series(roc.sel(category=2,mode=m).pseudo_pcs.values)
    score_ref0 = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.{score}.nc').sel(category=0,forecastMonth=fcmonth,mode=m).pseudo_pcs.values)
    score_ref1 = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.{score}.nc').sel(category=1,forecastMonth=fcmonth,mode=m).pseudo_pcs.values)
    score_ref2 = float(xr.open_dataset(f'{SCOREDIR}/{hcst_bname}.{aggr}.{score}.nc').sel(category=2,forecastMonth=fcmonth,mode=m).pseudo_pcs.values)
    fig = plt.figure(figsize=(6,6))
    ax1 = fig.add_subplot()
    # Pdfs 
    bins = np.linspace(score_options[score][0][0],score_options[score][0][-1],50)
    n, bins, patch = plt.hist(score_values0, bins=bins, density=True, color='dodgerblue', alpha=0.4, rwidth=0.85,zorder=-1)
    n, bins, patch = plt.hist(score_values1, bins=bins, density=True, color='lightgrey', alpha=0.4, rwidth=0.85,zorder=-1)
    n, bins, patch = plt.hist(score_values2, bins=bins, density=True, color='firebrick', alpha=0.4, rwidth=0.85,zorder=-1)
    score_values0.plot.density(ax=ax1, color='dodgerblue', lw=2., label='Lower') # identical to s.plot.kde(...)
    score_values1.plot.density(ax=ax1, color='k', lw=2., label='Middle') # identical to s.plot.kde(...)
    score_values2.plot.density(ax=ax1, color='firebrick', lw=2., label='Upper') # identical to s.plot.kde(...)
    # Boxplot
    height = round(n.max()*1.8)
    medianprops = dict(linestyle='-', linewidth=2., color='k')
    meanprops = dict(marker='o', markeredgecolor='black', markerfacecolor='w')
    bplot = ax1.boxplot(score_values0, positions=[height*0.75-height*0.07], widths=[height*0.04], vert=False, showfliers=False, patch_artist=True, notch=True, showmeans=True, meanprops=meanprops, medianprops=medianprops)
    bplot['boxes'][0].set_facecolor('dodgerblue')
    bplot = ax1.boxplot(score_values1, positions=[height*0.75], widths=[height*0.04], vert=False, showfliers=False, patch_artist=True, notch=True, showmeans=True, meanprops=meanprops, medianprops=medianprops)
    bplot['boxes'][0].set_facecolor('lightgrey')
    bplot = ax1.boxplot(score_values2, positions=[height*0.75+height*0.07], widths=[height*0.04], vert=False, showfliers=False, patch_artist=True, notch=True, showmeans=True, meanprops=meanprops, medianprops=medianprops)
    bplot['boxes'][0].set_facecolor('firebrick')
    ax1.text(score_values0.mean(), height*0.8-height*0.07, '{:.3f}'.format(score_values0.mean()),family='sans-serif',weight='bold',size=11, horizontalalignment='center', verticalalignment='top')
    ax1.text(score_values1.mean(), height*0.8, '{:.3f}'.format(score_values1.mean()),family='sans-serif',weight='bold',size=11, horizontalalignment='center', verticalalignment='top')
    ax1.text(score_values2.mean(), height*0.8+height*0.07, '{:.3f}'.format(score_values2.mean()),family='sans-serif',weight='bold',size=11, horizontalalignment='center', verticalalignment='top')
    # Other things
    ax1.vlines(0.5, 0., height, color='k', linestyle='-', lw=1)
    ax1.vlines(score_ref0, 0., height, color='dodgerblue', linestyle='--', lw=1.5, alpha=0.6, zorder=0)
    ax1.vlines(score_ref1, 0., height, color='k', linestyle='--', lw=1.5, alpha=0.6, zorder=0)
    ax1.vlines(score_ref2, 0., height, color='firebrick', linestyle='--', lw=1.5, alpha=0.6, zorder=0)
    ax1.set_ylim([0.,height])
    ax1.set_xlim([score_options[score][0][0],score_options[score][0][-1]])
    ax1.set_yticks(np.linspace(0.,height,height+1),labels=np.linspace(0.,height,height+1))
    plt.title(tit_line, fontsize=12)
    ax1.set_xlabel(score.upper())
    textstr= '\n'.join((
        'Lower:',
        r'Reference={:.2f}'.format(score_ref0),
        r'Conf. int.=({:.2f},{:.2f})'.format(float(roc_ci.pseudo_pcs.sel(quantile=0.025,category=0,mode=m).values),float(roc_ci.pseudo_pcs.sel(quantile=0.975,category=0,mode=m).values))))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.01, 0.99, textstr, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='left', bbox=props)
    textstr= '\n'.join((
        'Middle:',
        r'Reference={:.2f}'.format(score_ref1),
        r'Conf. int.=({:.2f},{:.2f})'.format(float(roc_ci.pseudo_pcs.sel(quantile=0.025,category=1,mode=m).values),float(roc_ci.pseudo_pcs.sel(quantile=0.975,category=1,mode=m).values))))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.5, 0.99, textstr, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='center', bbox=props)      
    textstr= '\n'.join((
        'Upper:',
        r'Reference={:.2f}'.format(score_ref2),
        r'Conf. int.=({:.2f},{:.2f})'.format(float(roc_ci.pseudo_pcs.sel(quantile=0.025,category=2,mode=m).values),float(roc_ci.pseudo_pcs.sel(quantile=0.975,category=2,mode=m).values))))
    props = dict(boxstyle='round', facecolor='wheat', alpha=0.5)
    ax1.text(0.99, 0.99, textstr, transform=ax1.transAxes, fontsize=9,
            verticalalignment='top', horizontalalignment='right', bbox=props)      
    ax1.legend(fontsize=10, loc='center left')
    # Save figure
    figname = f'./{PLOTSDIR}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{name_pcs2[m]}.{score}-bootstrap.png'
    fig.savefig(figname,dpi=600)  
