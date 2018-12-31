##--------------------------------------------------------------------------
##
##			Script to read in WRF output files, extract necessary data,
##			then save into new NetCDF file (reduces file size for archiving)
##					-- GYoung
##
##--------------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np
import time 
from datetime import datetime, timedelta 
from netCDF4 import num2date, date2num 
import constants

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				IN
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------

###################################
# Pick file
###################################
filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'

###################################
# LOAD NETCDF FILE
###################################

nc1 = Dataset(filename1, 'r')

###################################
# PROCESS WRF DATA FOR USE
###################################

data1 = {}

## Domain information
data1['dx'] = float(nc1.DX)
data1['dy'] = float(nc1.DY)
data1['x_dim'] = len(nc1.dimensions['west_east'])
data1['y_dim'] = len(nc1.dimensions['south_north'])
data1['width_meters']  = data1['dx'] * (data1['x_dim'] - 1)
data1['height_meters'] = data1['dy'] * (data1['y_dim'] - 1)
data1['cen_lat']  = float(nc1.CEN_LAT)
data1['cen_lon']  = float(nc1.CEN_LON)
data1['truelat1'] = float(nc1.TRUELAT1)
data1['truelat2'] = float(nc1.TRUELAT2)
data1['standlon'] = float(nc1.STAND_LON)
data1['xlat'] = nc1.variables['XLAT'][:]
data1['xlon'] = nc1.variables['XLONG'][:]

## Thermodynamic data
data1['theta'] = nc1.variables['T'][:]+300 # potential temperature in K
data1['p'] = (nc1.variables['P'][:]+nc1.variables['PB'][:])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data1['p']/100000)**tempvar       
data1['Tk'] = tempvar0*data1['theta']	# temperature in K
data1['rho'] = data1['p']/(constants.R*data1['Tk'])		# air density in kg/m3

## Interpolated fields
data1['Z'] = np.zeros([np.size(data1['xlat'],0),69,np.size(data1['xlat'],1),np.size(data1['xlat'],2)])
tempvar1 = (data1['p'])/9.81
data1['Z'][:,:-1,:,:] = 0.5*(tempvar1[:,0:-1,:,:] + tempvar1[:,1:,:,:]) # Z in m at theta mid-point
data1['Z'][:,68,:,:] = np.nan	# populate model top elements with nans
data1['w'][:,:-1,:,:] = 0.5*(nc1.variables['W'][:,0:-1,:,:] + nc1.variables['W'][:,1:,:,:])


## Cloud microphysics variables
data1['qcloud'] = nc1.variables['QCLOUD'][:]  # LW mixing ratio in kg/kg
data1['qnisg'] = (nc1.variables['QNICE'][:]+nc1.variables['QNSNOW'][:]+nc1.variables['QNGRAUPEL'][:]) # total ice number concentration in kg-1
data1['nisg80'] = nc1.variables['NISG80'][:]*(data1['rho']) 	# Nisg>80 in kg-1
data1['nisg50'] = data1['qnisg'] - (nc1.variables['NI50'][:] - nc1.variables['NG50'][:])*(data1['rho']) # small ice number concentration in kg-1
data1['qrain'] = nc1.variables['QRAIN'][:]

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				OUT
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------

dataset =  Dataset('OUT/test1.nc', 'w', format ='NETCDF4_CLASSIC') 

print dataset.file_format 

###################################
## Global Attributes 
###################################
dataset.description = 'CNTRL simulation from Young et al., 2016 (GRL) -- Nest (1 km grid size).'  
# dataset.history = 'Created ' + time.ctime(time.time())  
dataset.source = 'netCDF4 python' 

###################################
## Data dimensions
###################################
time = dataset.createDimension('time', np.size(data1['xlat'],0))
level = dataset.createDimension('level', np.size(data1['theta'],1)) 
lat = dataset.createDimension('lat', data1['y_dim'])
lon = dataset.createDimension('lon', data1['x_dim']) 

###################################
## Dimensions variables
###################################
times = dataset.createVariable('time', np.float32, ('time',)) 
levels = dataset.createVariable('level', np.int32, ('level',)) 
latitudes = dataset.createVariable('latitude', np.float32, ('time','lat', 'lon',))
longitudes = dataset.createVariable('longitude', np.float32, ('time','lat','lon',)) 

dx = dataset.createVariable('dx',np.int32)

# data1['dx'] = float(nc1.DX)
# data1['dy'] = float(nc1.DY)
# data1['x_dim'] = len(nc1.dimensions['west_east'])
# data1['y_dim'] = len(nc1.dimensions['south_north'])
# data1['width_meters']  = data1['dx'] * (data1['x_dim'] - 1)
# data1['height_meters'] = data1['dy'] * (data1['y_dim'] - 1)

###################################
## Create 4-d variables
###################################
temperature = dataset.createVariable('air_temperature', np.float32, ('time','level','lat','lon')) 
theta = dataset.createVariable('air_potential_temperature', np.float32, ('time','level','lat','lon')) 
Z = dataset.createVariable('height', np.float32, ('time','level','lat','lon')) 
P = dataset.createVariable('air_pressure', np.float32, ('time','level','lat','lon')) 
rho = dataset.createVariable('air_density', np.float32, ('time','level','lat','lon')) 

W = dataset.createVariable('vertical_wind_speed', np.float32, ('time','level','lat','lon')) 
qcloud = dataset.createVariable('cloud_liquid_water_mixing_ratio', np.float32, ('time','level','lat','lon')) 
qrain = dataset.createVariable('rain_water_mixing_ratio', np.float32, ('time','level','lat','lon')) 

# data1['qnisg'] = (nc1.variables['QNICE'][:]+nc1.variables['QNSNOW'][:]+nc1.variables['QNGRAUPEL'][:]) # total ice number concentration in kg-1
# data1['nisg80'] = nc1.variables['NISG80'][:]*(data1['rho']) 	# Nisg>80 in kg-1
# data1['nisg50'] = data1['qnisg'] - (nc1.variables['NI50'][:] - nc1.variables['NG50'][:])*(data1['rho']) # small ice number concentration in kg-1

###################################
## Variable Attributes  
###################################
times.units = 'hours since 2015-11-27 00:00:00'  
times.calendar = 'gregorian' 
levels.units = 'm' 
latitudes.units = 'degree_north'  
longitudes.units = 'degree_east'  

dx.units = 'km'

temperature.units = 'K' 
theta.units = 'K' 
Z.units = 'm'
P.units = 'Pa'
rho.units = 'kg m-3'

W.units = 'm s-1'
qcloud.units = 'kg kg-1'
qrain.units = 'kg kg-1'

###################################
## Fill in times
###################################
# dates = [] 
# for n in range(temp.shape[0]): 
#      dates.append(datetime(2015, 11, 27) + n * timedelta(hours=0)) 
# times[:] = date2num(dates, units = times.units, calendar = times.calendar) 
# print 'time values (in units %s): ' % times.units + '\n', times[:] 

###################################
## Fill arrays
###################################
levels[:] = np.arange(0,np.size(data1['Z'],1))
latitudes[:,:,:] = data1['xlat'][:,:,:]
longitudes[:,:,:] = data1['xlon'][:,:,:]
temperature[:,:,:,:] = data1['Tk'][:,:,:,:]
theta[:,:,:,:] = data1['theta'][:,:,:,:]
# Z[:,:,:,:] = data1['Z'][:,:,:,:]
P[:,:,:,:] = data1['p'][:,:,:,:]
rho[:,:,:,:] = data1['rho'][:,:,:,:]

W[:,:,:,:] = data1['w'][:,:,:,:]
qcloud[:,:,:,:] = data1['qcloud'][:,:,:,:]
qrain[:,:,:,:] = data1['qrain'][:,:,:,:]

###################################
## Write out file
###################################
dataset.close()