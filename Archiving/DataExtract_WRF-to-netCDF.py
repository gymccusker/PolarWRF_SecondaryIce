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

dataset =  Dataset('OUT/test.nc', 'w', format ='NETCDF4_CLASSIC') 

print dataset.file_format 

# Global Attributes 
dataset.description = 'bogus example script'  
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
     dates.append(datetime(2015, 11, 27) + n *  
                       timedelta(hours=0)) 
times[:] = date2num(dates, units = times.units,   
                          calendar = times.calendar) 
print 'time values (in units %s): ' % times.units + 
                          '\n', times[:] 

# Fill arrays
lats = np.arange(-90,91,2.5) 
lons = np.arange(-180,180,2.5) 
latitudes[:] = lats  
longitudes[:] = lons 

dataset.close()