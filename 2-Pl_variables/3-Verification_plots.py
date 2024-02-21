import os
import sys
import inquirer
import xarray as xr
import numpy as np
import locale
import calendar
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
import cartopy.feature as cfeature
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

    # Subset of plots to be produced
    questions4 = [
    inquirer.List('aggregation',
                    message="Selecciona el tipo de agregación mensual",
                    choices=['1m','3m'],
                ),
    ]
    answers4 = inquirer.prompt(questions4)
    aggr = answers4["aggregation"]

    # Forecast month
    if aggr=='1m':
        questions5 = [
        inquirer.List('endmonth',
                        message="Mes de forecast",
                        choices=['Enero', 'Febrero','Marzo','Abril','Mayo','Junio','Julio','Agosto','Septiembre','Octubre','Noviembre','Diciembre'],
                    ),
        ]
        answers5 = inquirer.prompt(questions5)
        endmonth = np.where(np.array(['Enero', 'Febrero','Marzo','Abril','Mayo','Junio','Julio','Agosto','Septiembre','Octubre','Noviembre','Diciembre'])== answers5["endmonth"])[0][0]+1
        fcmonth = (endmonth-startmonth)+1 if (endmonth-startmonth)>=0 else (endmonth-startmonth)+13
    else:
        questions5 = [
        inquirer.List('endmonth',
                        message="Mes de forecast",
                        choices=['NDJ', 'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND'],
                    ),
        ]
        answers5 = inquirer.prompt(questions5)
        endmonth = np.where(np.array(['NDJ', 'DJF','JFM','FMA','MAM','AMJ','MJJ','JJA','JAS','ASO','SON','OND'])== answers5["endmonth"])[0][0]+1
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
)

# Directory selection
DATADIR = './data'
# Base name for hindcast
hcst_bname = '{origin}_s{system}_stmonth{start_month:02d}_hindcast{hcstarty}-{hcendy}_monthly'.format(**config)

### 3. Visualize verification plots ###
#######################################
print("3. Visualize verification plots")  

# Common labels to be used in plot titles
VARNAMES = {
    'z' : 'gepotential height',
    'u' : 'u wind component',
    'v' : 'v wind component',
}
CATNAMES=['lower tercile', 'middle tercile', 'upper tercile']
PRESSURES = {
    '500hPa' : 500.,
}

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

### 3a. Correlation ###
print("3a. Correlation")  

# Check if file exists
if not os.path.exists(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr.nc'):
    print('No se crearon los coeficientes de correlación')
    sys.exit()
elif not os.path.exists(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr_pval.nc'):
    print('No se crearon los p-valores')
    sys.exit()

# Read the data files
corr = xr.open_dataset(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr.nc')
corr_pval = xr.open_dataset(f'{DATADIR}/scores/{hcst_bname}.{aggr}.corr_pval.nc')

# For each pressure level
for press in PRESSURES:

    # Select forecast month and pressure level
    try:
        thiscorr = corr.sel(forecastMonth=fcmonth,level=PRESSURES[press])
        thiscorrpval = corr_pval.sel(forecastMonth=fcmonth,level=PRESSURES[press])
    except:
        thiscorr = corr.sel(forecastMonth=fcmonth)
        thiscorrpval = corr_pval.sel(forecastMonth=fcmonth)

    # For each var
    for var in thiscorr.data_vars:
        # Create figure
        fig = plt.figure(figsize=(9,8))
        ax = fig.add_subplot(projection=ccrs.PlateCarree())
        ax.set_extent([-30, 5, 25, 50], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2.)
        ax.gridlines(draw_labels=True)

        # Check if data values matrices need to be transposed
        corrvalues = thiscorr[var].values 
        corrpvalvalues = thiscorrpval[var].values
        if corrvalues.T.shape == (thiscorr[var].lat.size, thiscorr[var].lon.size):
            corrvalues = corrvalues.T
            corrpvalvalues = corrpvalvalues.T
        elif corrvalues.shape == (thiscorr[var].lat.size, thiscorr[var].lon.size):
            pass                           
        else:
            raise BaseException(f'Unexpected data value matrix shape: {corrvalues.shape}' )

        # Represent
        cs = ax.contourf(thiscorr[var].lon,thiscorr[var].lat,corrvalues,levels=score_options['corr'][0],cmap=score_options['corr'][1],extend=score_options['corr'][3])
        plt.title(tit_line1 + f' {press} {VARNAMES[var]}' +' (stippling where significance below 95%)\n' + tit_line2,loc='left',fontsize=12)    
        cbar_ax = fig.add_axes([0.07, 0.07, 0.9, 0.02])
        cb = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal',label='Correlation')
        ax.contourf(thiscorrpval[var].lon,thiscorrpval[var].lat,corrpvalvalues,levels=[0.05,np.inf],hatches=['...',None],colors='none')
    
        # Save figure
        figname = f'{DATADIR}/plots/stmonth{config["start_month"]:02d}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{var}{press}.corr.png'
        fig.savefig(figname,dpi=600,bbox_inches='tight')    

### 3b. Area under Relative Operating Characteristic (ROC) curve ###
print("3b. Area under Relative Operating Characteristic (ROC) curve")

# Check if file exists
if not os.path.exists(f'{DATADIR}/scores/{hcst_bname}.{aggr}.roc.nc'):
    print('No se creó el archivo con el Area under Relative Operating Characteristic (ROC) curve')
    sys.exit()

# Read the data file
roc = xr.open_dataset(f'{DATADIR}/scores/{hcst_bname}.{aggr}.roc.nc')

# For each pressure level
for press in PRESSURES:
    # Select forecast month
    try:
        thisroc = roc.sel(forecastMonth=fcmonth,level=PRESSURES[press])
    except:
        thisroc = roc.sel(forecastMonth=fcmonth)

    # For each variable
    for var in thisroc.data_vars:
        # Create figure
        fig = plt.figure(figsize=(20,6))
        gs = fig.add_gridspec(1,3)
        # For each quantile
        for icat in thisroc.category.values:
            # Create subplot
            ax = fig.add_subplot(gs[0,icat],projection=ccrs.PlateCarree())
            ax.set_extent([-30, 5, 25, 50], crs=ccrs.PlateCarree())
            ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
            ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2.)
            ax.gridlines(draw_labels=True)

            # Check if data values matrices need to be transposed
            avalues = thisroc.sel(category=icat)[var].values
            if avalues.T.shape == (thisroc.sel(category=icat)[var].lat.size, thisroc.sel(category=icat)[var].lon.size):
                avalues = avalues.T
            elif avalues.shape == (thisroc.sel(category=icat)[var].lat.size, thisroc.sel(category=icat)[var].lon.size):
                pass                           
            else:
                raise BaseException(f'Unexpected data value matrix shape: {avalues.shape}' )

            # Represent
            cs = ax.contourf(thisroc[var].lon,thisroc[var].lat,avalues,levels=score_options['roc'][0],cmap=score_options['roc'][1],extend=score_options['roc'][3])
            plt.title(f' ({CATNAMES[icat]})',loc='center',fontsize=12)  
        # Save figure    
        fig.suptitle(tit_line1 + f' {press}  {VARNAMES[var]}\n' + tit_line2 ,fontsize=12)
        cbar_ax = fig.add_axes([0.07, 0.07, 0.9, 0.02])
        cb = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal',label='ROC')
        figname = f'{DATADIR}/plots/stmonth{config["start_month"]:02d}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{var}{press}.roc.png'
        fig.savefig(figname,dpi=600,bbox_inches='tight')    

### 3c. Ranked Probability Skill Score (RPSS) ###
print("3c. Ranked Probability Skill Score (RPSS)")

# Check if file exists
if not os.path.exists(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rpss.nc'):
    print('No se creó el archivo con el Ranked Probability Skill Score (RPSS)')
    sys.exit()

# Read the data file
rpss = xr.open_dataset(f'{DATADIR}/scores/{hcst_bname}.{aggr}.rpss.nc')

# For each pressure level
for press in PRESSURES:
    # Select forecast month
    try:
        thisrpss = rpss.sel(forecastMonth=fcmonth,level=PRESSURES[press])
    except:
        thisrpss = rpss.sel(forecastMonth=fcmonth)
     
    # For each variable
    for var in thisrpss.data_vars:
        # Create figure
        fig = plt.figure(figsize=(9,8))
        ax = fig.add_subplot(projection=ccrs.PlateCarree())
        ax.set_extent([-30, 5, 25, 50], crs=ccrs.PlateCarree())
        ax.add_feature(cfeature.BORDERS, edgecolor='black', linewidth=0.5)
        ax.add_feature(cfeature.COASTLINE, edgecolor='black', linewidth=2.)
        ax.gridlines(draw_labels=True)

        # Check if data values matrices need to be transposed
        avalues = thisrpss[var].values
        if avalues.T.shape == (thisrpss[var].lat.size, thisrpss[var].lon.size):
            avalues = avalues.T
        elif avalues.shape == (thisrpss[var].lat.size, thisrpss[var].lon.size):
            pass                           
        else:
            raise BaseException(f'Unexpected data value matrix shape: {avalues.shape}' )

        # Represent
        cs = ax.contourf(thisrpss[var].lon,thisrpss[var].lat,avalues,levels=score_options['rpss'][0],cmap=score_options['rpss'][1],extend=score_options['rpss'][3])
        plt.title(tit_line1 + f' {press} {VARNAMES[var]}\n' + tit_line2,loc='left',fontsize=12)
        cbar_ax = fig.add_axes([0.07, 0.07, 0.9, 0.02])
        cb = fig.colorbar(cs, cax=cbar_ax, orientation='horizontal',label='RPSS')

        # Save figure
        figname = f'{DATADIR}/plots/stmonth{config["start_month"]:02d}/{hcst_bname}.{aggr}.{"".join(validmonths)}.{var}{press}.rpss.png'
        fig.savefig(figname,dpi=600,bbox_inches='tight')    

