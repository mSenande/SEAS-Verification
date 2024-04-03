# %% [markdown]
# # 2. Compute deterministic and probabilistic scores

# This script is used to compute different verification scores 
# for monthly seasonal forescasts of pressure level variables.
# 
# The computed scores are: Spearman's rank correlation, area under Relative Operating Characteristic (ROC) curve, 
# Relative Operating Characteristic Skill Score (ROCSS), Ranked Probability Score (RPS), Ranked Probability Skill Score (RPSS) and Brier Score (BS).
#
# First we have to decide a forecast system (institution and system name) and a start month. 

#%%
print("2. Compute deterministic and probabilistic scores")

import os
import sys
import xarray as xr
import pandas as pd
import numpy as np
import xskillscore as xs
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
    list_vars = [ 'geopotential', 'u_component_of_wind', 'v_component_of_wind'],
    pressure_level = '500',
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
anom = hcst - hcmean
anom = anom.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))
# Calculate 3m anomalies
hcmean_3m = hcst_3m.mean(['number','start_date'])
anom_3m = hcst_3m - hcmean_3m
anom_3m = anom_3m.assign_attrs(reference_period='{hcstarty}-{hcendy}'.format(**config))

# %% [markdown]
# ## 2.2 Probabilities for tercile categories

# Here we get the probabilities for tercile categories of the hindcast data, 
# by counting the number of ensemble members found in each tercile.

# %% 
print("2.2 Probabilities for tercile categories")

# We define a function to calculate the boundaries of forecast categories defined by quantiles
def get_thresh(icat,quantiles,xrds,dims=['number','start_date']):
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
# For each aggregation
for aggr,h in [("1m",hcst), ("3m",hcst_3m)]:
    l_probs_hcst=list()
    # For each quantile
    for icat in range(numcategories):
        # Get the lower and higher threshold
        h_lo,h_hi = get_thresh(icat, quantiles, h)
        # Count the number of member between the threshold
        probh = np.logical_and(h>h_lo, h<=h_hi).sum('number')/float(h.dims['number'])
        # Instead of using the coordinate 'quantile' coming from the hindcast xr.Dataset
        # we will create a new coordinate called 'category'
        if 'quantile' in probh:
            probh = probh.drop('quantile')
        l_probs_hcst.append(probh.assign_coords({'category':icat}))

    # Concatenating tercile probs categories
    if aggr=='1m':
        probs_1m = xr.concat(l_probs_hcst,dim='category')                    
    elif aggr=='3m':
        probs_3m = xr.concat(l_probs_hcst,dim='category')                    

# %% [markdown]
# ## 2.3 Read observation data

# We read the monthly ERA5 data and obtain 3-months means.

#%%
print("2.3 Read observation data")  

# Read observations data into a xr.Dataset
era5_1deg = xr.open_dataset(obs_fname, engine='cfgrib')

# Renaming to match hindcast names 
era5_1deg = era5_1deg.rename({'latitude':'lat','longitude':'lon','time':'start_date'}).swap_dims({'start_date':'valid_time'})

# Assign 'forecastMonth' coordinate values
fcmonths = [mm+1 if mm>=0 else mm+13 for mm in [t.month - config['start_month'] for t in pd.to_datetime(era5_1deg.valid_time.values)] ]
era5_1deg = era5_1deg.assign_coords(forecastMonth=('valid_time',fcmonths))
# Drop obs values not needed (earlier than first start date) - this is useful to create well shaped 3-month aggregations from obs.
era5_1deg = era5_1deg.where(era5_1deg.valid_time>=np.datetime64('{hcstarty}-{start_month:02d}-01'.format(**config)),drop=True)

# Calculate 3-month AGGREGATIONS
era5_1deg_3m = era5_1deg.rolling(valid_time=3).mean()
era5_1deg_3m = era5_1deg_3m.where(era5_1deg_3m.forecastMonth>=3)

# As we don't need it anymore at this stage, we can safely remove 'forecastMonth'
era5_1deg = era5_1deg.drop('forecastMonth')
era5_1deg_3m = era5_1deg_3m.drop('forecastMonth')

# %% [markdown]
# ## 2.4 Compute deterministic scores
        
# Here we calculate the Spearman's rank correlation and their p-values. 
# 
# This score is based on the ensemble mean, not on the probabilities for each tercile.

# %% 
print("2.4 Compute deterministic scores")

# Loop over aggregations
for aggr in ['1m','3m']:
    if aggr=='1m':
        o = era5_1deg
        h = anom
    elif aggr=='3m':
        o = era5_1deg_3m
        h = anom_3m
    else:
        raise BaseException(f'Unknown aggregation {aggr}')

    # Check if hindcast data is ensemble
    is_fullensemble = 'number' in h.dims

    l_corr=list()
    l_corr_pval=list()
    # For each forecast month
    for this_fcmonth in h.forecastMonth.values:
        # Select hindcast values
        thishcst = h.sel(forecastMonth=this_fcmonth).swap_dims({'start_date':'valid_time'})
        # Select observation values for this hindcast
        thisobs = o.where(o.valid_time==thishcst.valid_time,drop=True)
        # Compute ensemble mean (if data is an ensemble)
        thishcst_em = thishcst if not is_fullensemble else thishcst.mean('number')
        # Calculate Spearman's rank correlation
        l_corr.append( xs.spearman_r(thishcst_em, thisobs, dim='valid_time') )
        # Calculate p-value
        l_corr_pval.append ( xs.spearman_r_p_value(thishcst_em, thisobs, dim='valid_time') )

    # Concatenating (by fcmonth) correlation
    corr=xr.concat(l_corr,dim='forecastMonth')
    corr_pval=xr.concat(l_corr_pval,dim='forecastMonth')
    
    # Saving to netCDF file correlation   
    corr.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr.nc')
    corr_pval.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr_pval.nc')

# %% [markdown]
# ## 2.5 Compute probabilistic scores for tercile categories
        
# Here we calculate the probabilistic scores: area under Relative Operating Characteristic (ROC) curve, 
# Relative Operating Characteristic Skill Score (ROCSS), Ranked Probability Score (RPS), Ranked Probability Skill Score (RPSS) and Brier Score (BS). 

# %% 
print("2.5 Compute probabilistic scores for tercile categories")

# Loop over aggregations
for aggr in ['1m','3m']:
    if aggr=='1m':
        o = era5_1deg
        probs_hcst = probs_1m
    elif aggr=='3m':
        o = era5_1deg_3m
        probs_hcst = probs_3m
    else:
        raise BaseException(f'Unknown aggregation {aggr}')
   
    l_roc=list()
    l_rps=list()
    l_rpss=list()
    l_rpssd=list()
    l_rocss=list()
    l_bs=list()
    # For each forecast month
    for this_fcmonth in probs_hcst.forecastMonth.values:
        # Select hindcast values
        thishcst = probs_hcst.sel(forecastMonth=this_fcmonth).swap_dims({'start_date':'valid_time'})
        # Select observation values for this hindcast
        thiso = o.where(o.valid_time==thishcst.valid_time,drop=True)

        # Calculate probabilities from observations and climatology
        l_probs_obs=list()
        l_probs_clim=list()
        # For each quantile
        for icat in range(numcategories):
            # Get the lower and higher threshold
            o_lo,o_hi = get_thresh(icat, quantiles, thiso, dims=['valid_time'])
            # Count the number of "members" between the threshold (1 or 0)
            probo = 1. * np.logical_and(thiso>o_lo, thiso<=o_hi)
            if 'quantile' in probo:
                probo=probo.drop('quantile')
            l_probs_obs.append(probo.assign_coords({'category':icat}))
            # Count the number of months between the threshold (1 or 0)
            probc = np.logical_and(thiso>o_lo, thiso<=o_hi).sum('valid_time')/float(thiso.dims['valid_time'])        
            if 'quantile' in probc:
                probc=probc.drop('quantile')
            l_probs_clim.append(probc.assign_coords({'category':icat}))
        # Concatenate observations and climatology probabilities
        thisobs = xr.concat(l_probs_obs, dim='category')
        thisclim = xr.concat(l_probs_clim, dim='category')

        # Calculate the probabilistic (tercile categories) scores
        thisroc = xr.Dataset()
        thisrps = xr.Dataset()
        rpsclim = xr.Dataset()
        thisrpss = xr.Dataset()
        thisrpssd = xr.Dataset()
        thisrocss = xr.Dataset()
        thisbs = xr.Dataset()
        # For each variable
        for var in thishcst.data_vars:
            # Compute Area ROC
            thisroc[var] = xs.roc(thisobs[var],thishcst[var], dim='valid_time', bin_edges=np.linspace(0,1,101))
            # Compute RPS
            thisrps[var] = xs.rps(thisobs[var],thishcst[var], dim='valid_time', category_edges=None, input_distributions='p')
            # Compute climatological RPS
            rpsclim[var] = xs.rps(thisobs[var],thisclim[var], dim='valid_time', category_edges=None, input_distributions='p')
            # Compute RPSS
            thisrpss[var] = 1.-thisrps[var]/rpsclim[var]
            # Compute ROCSS
            thisrocss[var] = (thisroc[var] - 0.5) / (1. - 0.5)
            # Compute Brier Score
            bscat = list()
            for cat in thisobs[var].category:
                thisobscat = thisobs[var].sel(category=cat)
                thishcstcat = thishcst[var].sel(category=cat)
                bscat.append(xs.brier_score(thisobscat, thishcstcat, dim='valid_time'))
            thisbs[var] = xr.concat(bscat,dim='category')
        l_roc.append(thisroc)
        l_rps.append(thisrps)
        l_rpss.append(thisrpss)
        l_rocss.append(thisrocss)
        l_bs.append(thisbs)

    # Concatenate along forecast month
    roc=xr.concat(l_roc,dim='forecastMonth')
    rps=xr.concat(l_rps,dim='forecastMonth')
    rpss=xr.concat(l_rpss,dim='forecastMonth')
    rocss=xr.concat(l_rocss,dim='forecastMonth')
    bs=xr.concat(l_bs,dim='forecastMonth')

    # Save scores to netcdf
    rps.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rps.nc')
    rpss.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rpss.nc')
    bs.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.bs.nc')
    roc.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.roc.nc')
    rocss.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rocss.nc')

