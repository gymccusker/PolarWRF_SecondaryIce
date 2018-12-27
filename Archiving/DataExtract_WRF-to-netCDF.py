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
##---------------				IN
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
tempvar1 = (ph+phb)/9.81
data1['Z'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0) # Z in m at theta mid-point
data1['w'] = 0.5*(nc1.variables['W'][:,0:-1,:,:] + nc1.variables['W'][:,1:,:,:])

## Cloud microphysics variables
data1['qcloud'] = nc1.variables['QCLOUD'][:]  # LW mixing ratio in kg/kg
data1['qnisg'] = (nc1.variables['QNICE'][:]+nc1.variables['QNSNOW'][:]+nc1.variables['QNGRAUPEL'][:]) # total ice number concentration in kg-1
data1['nisg80'] = nc1.variables['NISG80'][:]*(data1['rho']) 	# Nisg>80 in kg-1
data1['nisg50'] = data1['qnisg'] - (nc1.variables['NI50'][:] - nc1.variables['NG50'][:])*(data1['rho']) # small ice number concentration in kg-1
data1['qrain'] = nc1.variables['QRAIN'][:]
data1['qliq'] = data1['qcloud'] + data1['qrain']


##--------------------------------------------------------------------------
##---------------				OUT
##--------------------------------------------------------------------------

dataset =  Dataset('OUT/test.nc', 'w', format ='NETCDF4_CLASSIC') 

print dataset.file_format 

# Global Attributes 
dataset.description = 'test script'  
dataset.history = 'Created ' + time.ctime(time.time())  
dataset.source = 'netCDF4 python module tutorial' 


level = dataset.createDimension('level', 10) 
lat = dataset.createDimension('lat', 73)
lon = dataset.createDimension('lon', 144) 
time = dataset.createDimension('time', None)


times = dataset.createVariable('time', np.float64, ('time',)) 
levels = dataset.createVariable('level', np.int32, ('level',)) 
latitudes = dataset.createVariable('latitude', np.float32, ('lat',))
longitudes = dataset.createVariable('longitude', np.float32, ('lon',)) 
# Create the actual 4-d variable
temp = dataset.createVariable('temp', np.float32, ('time','level','lat','lon')) 

# Variable Attributes  
latitudes.units = 'degree_north'  
longitudes.units = 'degree_east'  
levels.units = 'hPa' 
temp.units = 'K' 
times.units = 'hours since 0001-01-01 00:00:00'  
times.calendar = 'gregorian' 

# Fill in times. 
dates = [] 
for n in range(temp.shape[0]): 
     dates.append(datetime(2015, 11, 27) + n * timedelta(hours=0)) 
times[:] = date2num(dates, units = times.units, calendar = times.calendar) 
print 'time values (in units %s): ' % times.units + '\n', times[:] 

# # Fill arrays
lats = np.arange(-90,91,2.5) 
lons = np.arange(-180,180,2.5) 
latitudes[:] = lats  
longitudes[:] = lons 

dataset.close()