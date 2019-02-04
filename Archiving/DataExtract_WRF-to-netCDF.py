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
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'

runlabel_start = filename1.find('/MAC_WRF/') + 9
runlabel_end = filename1.find('/wrfout',runlabel_start)
runlabel = filename1[runlabel_start:runlabel_end]

if runlabel == '31_DeMott_WATSAT_eta70_MYNN':
	runlab = 'CNTRL'
if runlabel == '30_DeMott_WATSAT_HM_noThresh_eta70_MYNN':
	runlab = 'NoThresh'
if runlabel == '36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN':
	runlab = '2xHM'
if runlabel == '57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN':
	runlab = '5xHM'
if runlabel == '56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN':
	runlab = '10xHM'

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
data1['w'] = 0.5*(nc1.variables['W'][:,0:-1,:,:] + nc1.variables['W'][:,1:,:,:])

## Cloud microphysics variables
data1['qvap'] = nc1.variables['QVAPOR'][:]   # water vapour mixing ratio, kg kg-1
data1['qcloud'] = nc1.variables['QCLOUD'][:]  # LW mixing ratio in kg/kg
data1['qrain'] = nc1.variables['QRAIN'][:]
data1['qisg'] = nc1.variables['QICE'][:]+nc1.variables['QSNOW'][:]+nc1.variables['QGRAUP'][:] # total ice mass mixing ratio, kg kg-1
data1['qnisg'] = (nc1.variables['QNICE'][:]+nc1.variables['QNSNOW'][:]+nc1.variables['QNGRAUPEL'][:]) # total ice number concentration in kg-1
data1['nisg80'] = nc1.variables['NISG80'][:]*(data1['rho']) 	# Nisg>80 in kg-1
data1['nisg50'] = data1['qnisg'] - ((nc1.variables['NI50'][:] + nc1.variables['NG50'][:])*data1['rho']) # small ice number concentration in kg-1

## Force small negative values to zero
data1['qcloud'][data1['qcloud']<0] = 0
data1['qrain'][data1['qrain']<0] = 0
data1['qisg'][data1['qisg']<0] = 0
data1['qnisg'][data1['qnisg']<0] = 0
data1['nisg80'][data1['nisg80']<0] = 0
data1['nisg50'][data1['nisg50']<0] = 0

## Radiation fields
data1['swdnb'] = nc1.variables['SWDNB'][:,:,:] # instantaneous downwelling shortwave flux at surface, W m-2
data1['swdnbc'] = nc1.variables['SWDNBC'][:,:,:] # instantaneous clear sky downwelling shortwave flux at surface, W m-2
data1['lwdnb'] = nc1.variables['LWDNB'][:,:,:] # instantaneous downwelling longwave flux at surface, W m-2
data1['lwdnbc'] = nc1.variables['LWDNBC'][:,:,:] # instantaneous clear sky upwelling longwave flux at surface, W m-2
data1['swupb'] = nc1.variables['SWUPB'][:,:,:] # instantaneous upwelling shortwave flux at surface, W m-2
data1['swupbc'] = nc1.variables['SWUPBC'][:,:,:] # instantaneous clear sky upwelling shortwave flux at surface, W m-2
data1['lwupb'] = nc1.variables['LWUPB'][:,:,:] # instantaneous upwelling longwave flux at surface, W m-2
data1['lwupbc'] = nc1.variables['LWUPBC'][:,:,:] # instantaneous clear sky upwelling longwave flux at surface, W m-2

## Surface fields
data1['seaice'] = nc1.variables['SEAICE'][:] # sea ice concentration

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				OUT
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
###################################
## Global Attributes
###################################
str_dx = "%.1f" % data1['dx']	# x/y resolution in m
if data1['dx'] == 1000.0: str_domain = 'Nest'	# domain option
if data1['dx'] == 5000.0: str_domain = 'Parent' # domain option

###################################
## Open File
###################################
outfile = "".join(['OUT/',runlab,'_',str_domain,'_r1.nc'])
dataset =  Dataset(outfile, 'w', format ='NETCDF4_CLASSIC') 

print dataset.file_format 

###################################
## Global Attributes (Cont.)
###################################
str_levels = "%1i" % np.size(data1['theta'],1) # number of vertical levels
str_xdim = "%1i" % data1['x_dim'] # number of grid points in x
str_ydim = "%1i" % data1['y_dim'] # number of grid points in y
str_width = "%.1f" % data1['width_meters'] # domain width (x) in m
str_height = "%.1f" % data1['height_meters'] # domain height (y) in m
str_latmin = "%.4f" % np.nanmin(data1['xlat']*-1)
str_latmax = "%.4f" % np.nanmax(data1['xlat']*-1)
str_lonmin = "%.4f" % np.nanmin(data1['xlon']*-1)
str_lonmax = "%.4f" % np.nanmax(data1['xlon']*-1)
desc = runlab + ' simulation from Young et al., 2019 (GRL) -- ' + str_domain + ' domain. x/y grid size = ' + str_dx + ' m with ' + str_levels + ' vertical levels. Domain size = ' + str_xdim + ' x ' + str_ydim + ' grid points, equalling ' + str_width + ' x ' + str_height + ' m, from ' + str_latmax + ' degS to ' + str_latmin + ' degS and ' + str_lonmax + ' degW to ' + str_lonmin + ' degW.'
dataset.description = desc
dataset.history = 'Created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
dataset.source = 'Weather Research and Forecasting (WRF) model, version 3.6.1, with polar updates from the Byrd Polar Research Center (http://polarmet.osu.edu/PWRF/).' 
dataset.references = 'First published in Young et al., 2019 (GRL): Radiative effects of secondary ice enhancement in coastal Antarctic clouds.'
dataset.project = 'Microphysics of Antarctic Clouds (MAC), funded by the UK Natural Environment Research Council (Grant no. NE/K01305X/1).'
dataset.comment = 'Other WRF variables from this simulation are archived locally at the British Antarctic Survey. Contact Gillian Young (G.Young1@leeds.ac.uk) for details.'
dataset.institution = 'British Antarctic Survey.'
dataset.conventions = 'CF-1.6'

###################################
## Switch off automatic filling 
###################################
dataset.set_fill_off()

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
times = dataset.createVariable('time', np.float32, ('time',),fill_value='-9999') 
levels = dataset.createVariable('level', np.int32, ('level',),fill_value='-9999') 
latitudes = dataset.createVariable('latitude', np.float32, ('time','lat', 'lon',),fill_value='-9999')
longitudes = dataset.createVariable('longitude', np.float32, ('time','lat','lon',),fill_value='-9999') 

###################################
## Create 3-d variables
###################################
swdnb = dataset.createVariable('swdnb', np.float32, ('time','lat', 'lon',),fill_value='-9999')
swdnbc = dataset.createVariable('swdnbc', np.float32, ('time','lat', 'lon',),fill_value='-9999')
lwdnb = dataset.createVariable('lwdnb', np.float32, ('time','lat', 'lon',),fill_value='-9999')
lwdnbc = dataset.createVariable('lwdnbc', np.float32, ('time','lat', 'lon',),fill_value='-9999')
swupb = dataset.createVariable('swupb', np.float32, ('time','lat', 'lon',),fill_value='-9999')
swupbc = dataset.createVariable('swupbc', np.float32, ('time','lat', 'lon',),fill_value='-9999')
lwupb = dataset.createVariable('lwupb', np.float32, ('time','lat', 'lon',),fill_value='-9999')
lwupbc = dataset.createVariable('lwupbc', np.float32, ('time','lat', 'lon',),fill_value='-9999')
seaice = dataset.createVariable('seaice', np.float32, ('time','lat', 'lon',),fill_value='-9999')

###################################
## Create 4-d variables
###################################
temperature = dataset.createVariable('temperature', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
theta = dataset.createVariable('theta', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
height = dataset.createVariable('height', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
pressure = dataset.createVariable('pressure', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
rho = dataset.createVariable('rho', np.float32, ('time','level','lat','lon'),fill_value='-9999') 

W = dataset.createVariable('W', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
qvapor = dataset.createVariable('qvapor', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
qcloud = dataset.createVariable('qcloud', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
qrain = dataset.createVariable('qrain', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
qisg = dataset.createVariable('qisg', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
nisg =  dataset.createVariable('nisg', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
nisg80 =  dataset.createVariable('nisg80', np.float32, ('time','level','lat','lon'),fill_value='-9999') 
nisg50 =  dataset.createVariable('nisg50', np.float32, ('time','level','lat','lon'),fill_value='-9999') 

###################################
## Dimensions variables: standard names
###################################
times.standard_name = 'time'
levels.long_name = 'model_level'
latitudes.standard_name = 'latitude'
longitudes.standard_name = 'longitude'

###################################
## 3-d variables: standard names
###################################
swdnb.standard_name = 'surface_downwelling_shortwave_flux_in_air'
swdnbc.standard_name = 'surface_downwelling_shortwave_flux_in_air_assuming_clear_sky'
lwdnb.standard_name = 'surface_downwelling_longwave_flux_in_air'
lwdnbc.standard_name = 'surface_downwelling_longwave_flux_in_air_assuming_clear_sky'
swupb.standard_name = 'surface_upwelling_shortwave_flux_in_air'
swupbc.standard_name = 'surface_upwelling_shortwave_flux_in_air_assuming_clear_sky'
lwupb.standard_name = 'surface_upwelling_longwave_flux_in_air'
lwupbc.standard_name = 'surface_upwelling_longwave_flux_in_air_assuming_clear_sky'
seaice.standard_name = 'sea_ice_area_fraction'

# ###################################
# ## Create 4-d variables
# ###################################
temperature.standard_name = 'air_temperature'
theta.standard_name = 'air_potential_temperature'
height.standard_name = 'height'
pressure.standard_name = 'air_pressure'
rho.standard_name = 'air_density'

W.long_name = 'vertical_wind_speed'
qvapor.standard_name = 'humidity_mixing_ratio'
qcloud.standard_name = 'cloud_liquid_water_mixing_ratio'
qrain.long_name = 'rain_water_mixing_ratio'
qisg.standard_name = 'cloud_ice_mixing_ratio'
nisg.long_name = 'number_concentration_of_ice_crystals_in_air'
nisg80.long_name = 'number_concentration_of_ice_crystals_larger_than_80micron_in_air'
nisg50.long_name = 'number_concentration_of_ice_crystals_smaller_than_50micron_in_air'

###################################
## Variable Attributes  
###################################
times.units = 'hours since 2015-11-27 00:00:00'  
times.calendar = 'gregorian' 
levels.units = 'm' 
latitudes.units = 'degree_north'  
longitudes.units = 'degree_east'  

swdnb.units = 'W m-2'
swdnbc.units = 'W m-2'
lwdnb.units = 'W m-2'
lwdnbc.units = 'W m-2'
swupb.units = 'W m-2'
swupbc.units = 'W m-2'
lwupb.units = 'W m-2'
lwupbc.units = 'W m-2'
seaice.units = 'dimensionless'

temperature.units = 'K' 
theta.units = 'K' 
height.units = 'm'
pressure.units = 'Pa'
rho.units = 'kg m-3'

W.units = 'm s-1'
qvapor.units = 'kg kg-1'
qcloud.units = 'kg kg-1'
qrain.units = 'kg kg-1'
qisg.units = 'kg kg-1'
nisg.units = 'kg-1'
nisg80.units = 'kg-1'
nisg50.units = 'kg-1'


###################################
## Ice comments
###################################
qisg.comment = 'Sum of ice + snow + graupel particle categories'
nisg.comment = 'Sum of ice + snow + graupel particle categories'
nisg80.comment = 'Sum of ice + snow + graupel particle categories. Particle sizes calculated online following the assumption that each hydrometeor class is represented by a Gamma distribution.'
nisg50.comment = 'Sum of ice + snow + graupel particle categories. Particle sizes calculated online following the assumption that each hydrometeor class is represented by a Gamma distribution.'

###################################
## Fill in times
###################################
# dates = [] 
# for n in range(temp.shape[0]): 
# 	dates.append(datetime(2015, 11, 27) + n * timedelta(hours=0)) 
# 	times[:] = date2num(dates, units = times.units, calendar = times.calendar) 
# print 'time values (in units %s): ' % times.units + '\n', times[:] 

wrftime = nc1.variables['Times']
tim = np.zeros(np.size(data1['Tk'],0))
for i in range(np.size(data1['Tk'],0)):
	str_times = wrftime[i][11:]
	tim[i] = (np.int(str_times[0])*600 + np.int(str_times[1])*60 + np.int(str_times[3])*10)/float(60)

###################################
## Fill arrays
###################################
times[:] = tim[:]
levels[:] = np.arange(0,np.size(data1['Z'],1))
latitudes[:,:,:] = data1['xlat'][:,:,:]
longitudes[:,:,:] = data1['xlon'][:,:,:]

swdnb[:,:,:] = data1['swdnb'][:,:,:]
swdnbc[:,:,:] = data1['swdnbc'][:,:,:]
lwdnb[:,:,:] = data1['lwdnb'][:,:,:]
lwdnbc[:,:,:] = data1['lwdnbc'][:,:,:]
swupb[:,:,:] = data1['swupb'][:,:,:]
swupbc[:,:,:] = data1['swupbc'][:,:,:]
lwupb[:,:,:] = data1['lwupb'][:,:,:]
lwupbc[:,:,:] = data1['lwupbc'][:,:,:]
seaice[:,:,:] = data1['seaice'][:,:,:]

temperature[:,:,:,:] = data1['Tk'][:,:,:,:]
theta[:,:,:,:] = data1['theta'][:,:,:,:]
height[:,:,:,:] = data1['Z'][:,:,:,:]
pressure[:,:,:,:] = data1['p'][:,:,:,:]
rho[:,:,:,:] = data1['rho'][:,:,:,:]

W[:,:,:,:] = data1['w'][:,:,:,:]
qvapor[:,:,:,:] = data1['qvap'][:,:,:,:]
qcloud[:,:,:,:] = data1['qcloud'][:,:,:,:]
qrain[:,:,:,:] = data1['qrain'][:,:,:,:]
qisg[:,:,:,:] = data1['qisg'][:,:,:,:]
nisg[:,:,:,:] = data1['qnisg'][:,:,:,:]
nisg80[:,:,:,:] = data1['nisg80'][:,:,:,:]
nisg50[:,:,:,:] = data1['nisg50'][:,:,:,:]

###################################
## Write out file
###################################
dataset.close()