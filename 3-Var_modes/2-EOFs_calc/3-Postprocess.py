import os
import sys
import inquirer
import xarray as xr
import numpy as np
import xskillscore as xs
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
DATADIR = './data'
MODESDIR = './data/modes'
# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)
# File name for hindcast
hcpcs_fname = f'{MODESDIR}/{hcst_bname}.1m.PCs.nc'
hcpcs_3m_fname = f'{MODESDIR}/{hcst_bname}.3m.PCs.nc'
# Base name for observations
obs_bname = 'era5_monthly_stmonth{start_month:02d}_{hcstarty}-{hcendy}'.format(**config)
# File name for observations
obpcs_fname = f'{MODESDIR}/{obs_bname}.1m.PCs.nc'
obpcs_3m_fname = f'{MODESDIR}/{obs_bname}.3m.PCs.nc'


# Check if files exist
if not os.path.exists(obpcs_fname) & os.path.exists(obpcs_3m_fname):
    print('No se calcularon aún las PCs de ERA5')
    sys.exit()
elif not os.path.exists(hcpcs_fname) & os.path.exists(hcpcs_3m_fname):
    print('No se calcularon aún las PCs de este modelo y sistema')
    sys.exit()

### 3. Compute deterministic and probabilistic scores ###
#########################################################
print("3. Compute deterministic and probabilistic scores")  

# Reading HCST data from file
# Reading OBS data from file
hcpcs = xr.open_dataset(hcpcs_fname)
hcpcs_3m = xr.open_dataset(hcpcs_3m_fname)

# Reading OBS data from file
obpcs = xr.open_dataset(obpcs_fname)
obpcs_3m = xr.open_dataset(obpcs_3m_fname)

### 3a. Probabilities for tercile categories ###
print("3a. Probabilities for tercile categories")

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
for aggr,h in [("1m",hcpcs), ("3m",hcpcs_3m)]:
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
    if aggr=='1m':
        probs_1m = xr.concat(l_probs_hcst,dim='category')                    
    elif aggr=='3m':
        probs_3m = xr.concat(l_probs_hcst,dim='category')                    

### 3b. Compute deterministic scores ###
print("3b. Compute deterministic scores")

# Loop over aggregations
for aggr in ['1m','3m']:

    if aggr=='1m':
        o = obpcs
        h = hcpcs
    elif aggr=='3m':
        o = obpcs_3m
        h = hcpcs_3m
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
    corr = xr.concat(l_corr,dim='forecastMonth')                    
    corr_pval = xr.concat(l_corr_pval,dim='forecastMonth')

    # Saving to netCDF file correlation   
    corr.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr.nc')
    corr_pval.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr_pval.nc')

### 3c. Compute probabilistic scores for tercile categories ###
print("3c. Compute probabilistic scores for tercile categories")

# Loop over aggregations
for aggr in ['1m','3m']:

    if aggr=='1m':
        o = obpcs
        probs_hcst = probs_1m
    elif aggr=='3m':
        o = obpcs_3m
        probs_hcst = probs_3m
    else:
        raise BaseException(f'Unknown aggregation {aggr}')
  
    l_roc=list()
    l_rps=list()
    l_rpss=list()
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
            if 'quantile' in list(probo.coords):
                probo=probo.drop('quantile')
            l_probs_obs.append(probo.assign_coords({'category':icat}))
            # Count the number of months between the threshold (1 or 0)
            probc = np.logical_and(thiso>o_lo, thiso<=o_hi).sum('valid_time')/float(thiso.valid_time.size)        
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
    roc=xr.concat(l_roc,dim='forecastMonth').squeeze()
    rps=xr.concat(l_rps,dim='forecastMonth').squeeze()
    rpss=xr.concat(l_rpss,dim='forecastMonth').squeeze()
    rocss=xr.concat(l_rocss,dim='forecastMonth').squeeze()
    bs=xr.concat(l_bs,dim='forecastMonth').squeeze()

    # Save scores to netcdf
    rps.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rps.nc')
    rpss.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rpss.nc')
    bs.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.bs.nc')
    roc.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.roc.nc')
    rocss.to_netcdf(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rocss.nc')

