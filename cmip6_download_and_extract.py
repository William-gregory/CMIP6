import glob
import numpy as np
from netCDF4 import Dataset
import os
import pickle
import shutil
from mpl_toolkits.basemap import Basemap
from scipy.interpolate import griddata
from datetime import datetime,timedelta

#models = ['ACCESS-CM2','ACCESS-ESM1-5','AWI-ESM-1-1-LR','BCC-CSM2-MR','BCC-ESM1','CAMS-CSM1-0','CanESM5','CESM2','CESM2-WACCM','CESM2-WACCM-FV2','CNRM-CM6-1','CNRM-CM6-1-HR','CNRM-ESM2-1','E3SM-1-0','EC-Earth3','EC-Earth3-Veg','FGOALS-g3','FIO-ESM-2-0','GFDL-CM4','GFDL-ESM4','GISS-E2-1-G','GISS-E2-1-G-CC','GISS-E2-1-H','HadGEM3-GC31-LL','HadGEM3-GC31-MM','INM-CM4-8','INM-CM5-0','IPSL-CM6A-LR','MIROC-ES2L','MIROC6','MPI-ESM-1-2-HAM','MPI-ESM1-2-HR','MPI-ESM1-2-LR','MRI-ESM2-0','NESM3','NorCPM1','NorESM2-LM','NorESM2-MM','SAM0-UNICON','UKESM1-0-LL']

#PARAMETER SETUP
model = input('Please specify which model you would like to analyse. E.g., HadGEM3-GC31-LL:\n')
yrange = input('Please specify time window to perform analysis. E.g., 1850-2100:\n')
y1 = int(yrange.split('-')[0]) ; y2 = int(yrange.split('-')[1])
if y2 > 2014:
    ssp = input('Please specify scenario to include observations beyond 2014. E.g., ssp585:\n')
else:
    ssp = False
res = input('Please specify frequency. E.g., mon or day:\n')
variables = input('Please specify which variables. E.g., siconc,sithick:\n')
variables = variables.split(',')
nvars = len(variables)
base = '/badc/cmip6/data/CMIP6/'

def common(lists):
    return sorted(list(set.intersection(*map(set, lists))))

def get_common_ripfs(var):
    histf = sorted(glob.glob(base+'CMIP/*/'+model+'/historical/*/SI'+res+'/'+var+'/gn/latest/*.nc'))
    if len(histf) == 0:
        os.system('python cmip6_downloader.py --variable_id='+var+' --frequency='+res+' --experiment_id=historical --source_id='+model)
        histf = sorted(glob.glob('./data/s_'+model+'_e_historical_f_'+res+'_v_'+var+'/*'))
        hist_ripfs = [key.split('/')[-1].split('_')[4] for key in histf]
    else:
        hist_ripfs = [key.split('/')[9] for key in histf]
    if ssp:
        sspf = sorted(glob.glob(base+'ScenarioMIP/*/'+model+'/'+ssp+'/*/SI'+res+'/'+var+'/gn/latest/*.nc'))
        if len(sspf) == 0:
            os.system('python cmip6_downloader.py --variable_id='+var+' --frequency='+res+' --experiment_id='+ssp+' --source_id='+model)
            sspf = sorted(glob.glob('./data/s_'+model+'_e_'+ssp+'_f_'+res+'_v_'+var+'/*'))
            ssp_ripfs = [key.split('/')[-1].split('_')[4] for key in sspf]
        else:
            ssp_ripfs = [key.split('/')[9] for key in sspf]
        return common([hist_ripfs,ssp_ripfs])
    else:
        return hist_ripfs

def trim_files(files):
    if res == 'mon':
        scale = 100
    elif res == 'day':
        scale = 10000
    yrange = [key.split('/')[-1].split('_')[-1].split('.')[0].split('-') for key in files]
    ymins = [round(int(a[0])/scale) for a in yrange]
    ymaxs = [round(int(a[1])/scale) for a in yrange]
    dist_y1 = np.array([y-y1 for y in ymins])
    dist_y2 = np.array([y-y2 for y in ymaxs])
    min_dist = max([d for d in dist_y1 if d<=0])
    max_dist = min([d for d in dist_y2  if d>=0])
    f1 = np.array(yrange)[np.where(dist_y1==min_dist)]
    f2 = np.array(yrange)[np.where(dist_y2==max_dist)]
    f1 = str(f1[0][0])+'-'+str(f1[0][1])
    f2 = str(f2[0][0])+'-'+str(f2[0][1])
    yrange_str = np.array([key.split('/')[-1].split('_')[-1].split('.')[0] for key in files])
    lID = np.where(yrange_str==f1)[0][0]
    uID = np.where(yrange_str==f2)[0][0]
    return list(np.array(files)[lID:uID+1])

def read(files,ripf,var):
    read = Dataset(files[0])
    nc = [] ; times = []
    for file in files:
        nc.append(np.array(Dataset(file)[var]))
        times.append(np.array(Dataset(file)['time']))
    times = np.concatenate(times,0)
    data = np.concatenate(nc,0)
    data[data>100] = np.nan
    data[data<0] = np.nan
    if (var == 'siconc') & (np.nanmax(data)>1):
        data = data/100
    try:
        try:
            try:
                lat = np.array(read['lat'])
                lon = np.array(read['lon'])
            except IndexError as e1:
                lat = np.array(read['latitude'])
                lon = np.array(read['longitude'])
        except IndexError as e2:
            lat = np.array(read['nav_lat'])
            lon = np.array(read['nav_lon'])
    except IndexError as e3:
        lat = np.array(read['nav_latitude'])
        lon = np.array(read['nav_longitude'])
    if (lat.ndim == 1) & (lon.ndim == 1):
        lon,lat = np.meshgrid(lon,lat)
    x,y = m(lon,lat)
    NH = np.where(lat>=40)
    yrange = [key.split('/')[-1].split('_')[-1].split('.')[0].split('-') for key in files]
    years = [int(a) for b in yrange for a in b]
    dT,dX,dY = data.shape
    if res == 'mon':
        ymin = round(np.min(years)/100) ; ymax = round(np.max(years)/100)
        data = data.transpose(1,2,0).reshape(dX,dY,12,ymax-ymin+1,order='F')
        data = data[:,:,:,y1-ymin:y2-ymin+1].reshape(dX,dY,(y2-y1+1)*12,order='F')
    elif res == 'day':
        t0 = datetime.strptime("1/1/1850", "%m/%d/%Y")
        tplus = np.array([str(t0 + timedelta(days=t)).split(' ')[0] for t in times])
        startID = np.where(tplus==str(y1)+'-01-01')[0][0]
        endID = np.where(tplus==str(y2)+'-12-31')[0][0]
        data = data[startID:endID+1,:,:]
    dT = data.shape[2]
    regrid = np.zeros((dXR,dYR,dT))*np.nan
    for t in range(dT):
        regrid[:,:,t] = griddata((x[NH],y[NH]),data[:,:,t][NH],(xr,yr),'nearest')
    return regrid
   
#DEFINE A NEW 50km GRID
m = Basemap(projection='npstere',boundinglat=50,lon_0=360,resolution='l')
regrid_res = 50 #km
lonr,latr = m.makegrid(int((m.xmax-m.xmin)/(regrid_res*1000))+1, int((m.ymax-m.ymin)/(regrid_res*1000))+1)
xr,yr = m(lonr,latr)
dXR,dYR = xr.shape
psa = np.fromfile('./nsidc_nt_siconc/north/psn25area_v3.dat',dtype='<i4').reshape(448,304)/1000
lats = np.fromfile('./nsidc_nt_siconc/north/psn25lats_v3.dat',dtype='<i4').reshape(448,304)/100000
lons = np.fromfile('./nsidc_nt_siconc/north/psn25lats_v3.dat',dtype='<i4').reshape(448,304)/100000
xs,ys = m(lons,lats)
psar = ((regrid_res/25)**2) * griddata((xs.ravel(),ys.ravel()),psa.ravel(),(xr,yr),'nearest')

data = {}
#GET COMMON REALISATIONS FOR ALL VARIABLES, HISTORICAL AND FUTURE
common_internal = []
for variable in range(nvars):
    common_internal.append(get_common_ripfs(variables[variable]))
check_data = [len(n) for n in common_internal]
if 0 in check_data:
    which_data = np.where(np.array(check_data)==0)
    missing_var = np.array(variables)[which_data][0]
    print('Unable to find any consistent data files that span the range '+str(y1)+'-'+str(y2)+' for '+missing_var+'. Please try an alternative range, e.g., dates prior to or after 2014, or try a different model')
else:
    all_common = common(common_internal)
    if len(all_common) == 0:
        print('No common realisations found between selected variables')
        if (os.path.exists('./logs.txt')) & (os.path.exists('./data')) & (len(glob.glob('./s_*.txt'))>0):
            shutil.rmtree('./data',ignore_errors=True)
            os.remove('./logs.txt')
            [os.remove(f) for f in glob.glob('./s_*.txt')]
    else:
        #READ AND PROCESS DATA
        for ripf in all_common:
            for variable in range(nvars):
                print('Reading '+model+' '+ripf+' '+variables[variable]+' data...')
                histf = sorted(glob.glob(base+'CMIP/*/'+model+'/historical/'+ripf+'/SI'+res+'/'+variables[variable]+'/gn/latest/*.nc'))
                if len(histf) == 0:
                    histf = sorted(glob.glob('./data/s_'+model+'_e_historical_f_'+res+'_v_'+variables[variable]+'/*'+ripf+'*'))
                if ssp:
                    sspf = sorted(glob.glob(base+'ScenarioMIP/*/'+model+'/'+ssp+'/'+ripf+'/SI'+res+'/'+variables[variable]+'/gn/latest/*.nc'))
                    if len(sspf) == 0:
                        sspf = sorted(glob.glob('./data/s_'+model+'_e_'+ssp+'_f_'+res+'_v_'+variables[variable]+'/*'+ripf+'*'))
                    files = trim_files(histf+sspf)
                else:
                    files = trim_files(histf)
                data[model+'_'+variables[variable]+'_'+ripf] = read(files,ripf,variables[variable])

        clean = 'y'#input('Loading complete. Remove downloaded files on local disk to save space? y or n:\n')
        if clean == 'y':
            if (os.path.exists('./logs.txt')) & (os.path.exists('./data')) & (len(glob.glob('./s_*.txt'))>0):
                shutil.rmtree('./data',ignore_errors=True)
                os.remove('./logs.txt')
                [os.remove(f) for f in glob.glob('./s_*.txt')]

if nvars == 1:
    fname = './'+model+'_'+variables[0]+'_'+res+'_data_'+ssp+'_'+str(y1)+'-'+str(y2)+'.pkl'
else:
    fname = './'+model+'_multivars_'+res+'_data_'+ssp+'_'+str(y1)+'-'+str(y2)+'.pkl'
f = open(fname,'wb')
pickle.dump(data,f)
f.close()
    




