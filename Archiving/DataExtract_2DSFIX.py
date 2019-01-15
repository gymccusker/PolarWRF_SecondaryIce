##--------------------------------------------------------------------------
##
##			Script to read in 2DS output files and fix for CF-compliance
##					-- GYoung
##
##--------------------------------------------------------------------------

from netCDF4 import Dataset
import numpy as np
import time 
from datetime import datetime, timedelta 
from netCDF4 import num2date, date2num 
import matplotlib.pyplot as plt
# import constants

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				IN
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------

###################################
# Pick file
###################################
# filename1 = '/data/scihub-users/giyoung/MAC/FlightData/Processed2DS/UMAN_2DS_20151127_r0_Flight218.nc'
filename1 = '/data/scihub-users/giyoung/MAC/FlightData/Processed2DS/UMAN_2DS_20151127_r1_Flight219.nc'


# UMAN_2DS_20151127_r0_Flight218.nc  UMAN_2DS_20151208_r0_Flight227.nc
# UMAN_2DS_20151127_r1_Flight219.nc  UMAN_2DS_20151209_r0_Flight228.nc
# UMAN_2DS_20151128_r0_Flight220.nc  UMAN_2DS_20151209_r0_Flight229.nc
# UMAN_2DS_20151129_r0_Flight221.nc  UMAN_2DS_20151210_r0_Flight230.nc
# UMAN_2DS_20151130_r0_Flight222.nc  UMAN_2DS_20151211_r0_Flight231.nc
# UMAN_2DS_20151203_r0_Flight223.nc  UMAN_2DS_20151211_r0_Flight232.nc
# UMAN_2DS_20151206_r0_Flight224.nc  UMAN_2DS_20151212_r0_Flight233.nc
# UMAN_2DS_20151207_r0_Flight226.nc  UMAN_2DS_20151213_r0_Flight234.nc
# UMAN_2DS_20151207_r1_Flight225.nc  UMAN_2DS_20151214_r0_Flight235.nc

###################################
# LOAD NETCDF FILE
###################################

nc1 = Dataset(filename1, 'r')

revis_start = filename1.find('/UMAN_2DS_') + 20
revis_end = filename1.find('_Flight',revis_start)
revis = filename1[revis_start:revis_end]
revis = str(int(revis)+1)
 
flightno_start = filename1.find('/UMAN_2DS_') + 28
flightno_end = filename1.find('.nc',flightno_start)
flightno = filename1[flightno_start:flightno_end]

date_start = filename1.find('/UMAN_2DS_') + 10
date_end = filename1.find('_r',date_start)
date = filename1[date_start:date_end]

year = date[0:4]
month = date[4:6]
day = date[6:]

##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
##---------------				OUT
##--------------------------------------------------------------------------
##--------------------------------------------------------------------------
###################################
## Open File
###################################
outfile = "".join(['OUT/UMAN_2DS_20151127_r',revis,'_Flight',flightno,'.nc'])
dataset =  Dataset(outfile, 'w', format ='NETCDF4_CLASSIC') 

print dataset.file_format 

###################################
## Global Attributes
###################################
dataset.title = 'University of Manchester 2DS observations on board the MASIN research aircraft, flight number Flight' + flightno + ' (' + year + '-' + month + '-' + day + ').'
dataset.history = 'Revision number ' + revis + ', created ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S')
dataset.source = nc1.source
dataset.references = nc1.references
dataset.project = 'Microphysics of Antarctic Clouds (MAC), funded by the UK Natural Environment Research Council (Grant no. NE/K01305X/1).'
dataset.comment = nc1.comment
dataset.institution = nc1.institution

###################################
## Switch off automatic filling 
###################################
dataset.set_fill_off()

###################################
## Data dimensions
###################################
time_mid = dataset.createDimension('Time_mid', np.size(nc1.variables['Time_mid']))
time_edge = dataset.createDimension('Time_edge', np.size(nc1.variables['Time_edge']))
size_mid = dataset.createDimension('Size_mid', np.size(nc1.variables['Size_mid']))
size_edge = dataset.createDimension('Size_edge', np.size(nc1.variables['Size_edge']))

###################################
## Dimensions variables
###################################
#### Time_mid
time_mid = dataset.createVariable('Time_mid', np.float64, ('Time_mid',),fill_value='-9999') 
time_mid.scale_factor = float(1)
time_mid.add_offset = float(0)
time_mid.comment = 'None'
time_mid.units = ['seconds since ' + year + '-' + month + '-' + day + ' 00:00:00']
time_mid.long_name = 'Mid_point_of_time_bin' 
time_mid[:] = nc1.variables['Time_mid'][:]

#### Time_edge
time_edge = dataset.createVariable('Time_edge', np.float64, ('Time_edge',),fill_value='-9999') 
time_edge.scale_factor = float(1)
time_edge.add_offset = float(0)
time_edge.comment = 'None'
time_edge.units = ['seconds since ' + year + '-' + month + '-' + day + ' 00:00:00']
time_edge.long_name = 'Edge_point_of_time_bin' 
time_edge[:] = nc1.variables['Time_edge'][:]

#### Size_mid
size_mid = dataset.createVariable('Size_mid', np.float64, ('Size_mid',),fill_value='-9999') 
size_mid.scale_factor = float(1)
size_mid.add_offset = float(0)
size_mid.comment = 'None'
size_mid.units = 'micron' 
size_mid.long_name = 'Mid_point_of_size_bin' 
size_mid[:] = nc1.variables['Size_mid'][:]

#### Size_edge
size_edge = dataset.createVariable('Size_edge', np.float64, ('Size_edge',),fill_value='-9999') 
size_edge.scale_factor = float(1)
size_edge.add_offset = float(0)
size_edge.comment = 'None'
size_edge.units = 'micron' 
size_edge.long_name = 'Edge_point_of_size_bin' 
size_edge[:] = nc1.variables['Size_edge'][:]

###################################
## Create number concentrations
###################################
#### NC_All
nc_all = dataset.createVariable('NC_All', np.float64, ('Time_mid',),fill_value='-9999') 
nc_all.scale_factor = float(1)
nc_all.add_offset = float(0)
nc_all.comment = 'Particles in contact with the edge of the sample array have been rejected. Sum of small and low, medium, and high irregularity particle categories.'
nc_all.units = 'L-1' 
nc_all.long_name = 'Total_number_concentration_of_particles' 
nc_all[:] = nc1.variables['NC_S'][:] + nc1.variables['NC_LI'][:] + nc1.variables['NC_MI'][:] + nc1.variables['NC_HI'][:]

#### NC_S
nc_s = dataset.createVariable('NC_S', np.float64, ('Time_mid',),fill_value='-9999') 
nc_s.scale_factor = float(1)
nc_s.add_offset = float(0)
nc_s.comment = 'This category contains particles with areas between 0 and 30 pixels. Particles in contact with the edge of the sample array have been rejected.'
nc_s.units = 'L-1' 
nc_s.long_name = 'Number_concentration_of_particles_classed_as_small' 
nc_s[:] = nc1.variables['NC_S'][:]

#### NC_LI
nc_li = dataset.createVariable('NC_LI', np.float64, ('Time_mid',),fill_value='-9999') 
nc_li.scale_factor = float(1)
nc_li.add_offset = float(0)
nc_li.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 0.9 and 1.2. Particles in contact with the edge of the sample array have been rejected.'
nc_li.units = 'L-1' 
nc_li.long_name = 'Number_concentration_of_particles_classed_as_low_irregularity' 
nc_li[:] = nc1.variables['NC_LI'][:]

#### NC_MI
nc_mi = dataset.createVariable('NC_MI', np.float64, ('Time_mid',),fill_value='-9999') 
nc_mi.scale_factor = float(1)
nc_mi.add_offset = float(0)
nc_mi.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 1.2 and 1.4. Particles in contact with the edge of the sample array have been rejected.'
nc_mi.units = 'L-1' 
nc_mi.long_name = 'Number_concentration_of_particles_classed_as_medium_irregularity' 
nc_mi[:] = nc1.variables['NC_MI'][:]

#### NC_HI
nc_hi = dataset.createVariable('NC_HI', np.float64, ('Time_mid',),fill_value='-9999') 
nc_hi.scale_factor = float(1)
nc_hi.add_offset = float(0)
nc_hi.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 1.4 and inf. Particles in contact with the edge of the sample array have been rejected.'
nc_hi.units = 'L-1' 
nc_hi.long_name = 'Number_concentration_of_particles_classed_as_high_irregularity' 
nc_hi[:] = nc1.variables['NC_HI'][:]

###################################
## Create particle size distributions - number
###################################
#### PSD_Num_All
psd_num_all = dataset.createVariable('PSD_Num_All', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_num_all.scale_factor = float(1)
psd_num_all.add_offset = float(0)
psd_num_all.comment = 'Particles in contact with the edge of the sample array have been rejected. Sum of small and low, medium, and high irregularity particle categories.'
psd_num_all.units = 'L-1' 
psd_num_all.long_name = 'Number_size_distribution_of_all_particles' 
psd_num_all[:] = nc1.variables['PSD_Num_S'][:] + nc1.variables['PSD_Num_LI'][:] + nc1.variables['PSD_Num_MI'][:] + nc1.variables['PSD_Num_HI'][:]

#### PSD_Num_S
psd_num_s = dataset.createVariable('PSD_Num_S', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_num_s.scale_factor = float(1)
psd_num_s.add_offset = float(0)
psd_num_s.comment = 'This category contains particles with areas between 0 and 30 pixels. Particles in contact with the edge of the sample array have been rejected.'
psd_num_s.units = 'L-1' 
psd_num_s.long_name = 'Number_size_distribution_of_particles_classed_as_small' 
psd_num_s[:] = nc1.variables['PSD_Num_S'][:]

#### PSD_Num_LI
psd_num_li = dataset.createVariable('PSD_Num_LI', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_num_li.scale_factor = float(1)
psd_num_li.add_offset = float(0)
psd_num_li.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 0.9 and 1.2. Particles in contact with the edge of the sample array have been rejected.'
psd_num_li.units = 'L-1' 
psd_num_li.long_name = 'Number_size_distribution_of_particles_classed_as_low_irregularity' 
psd_num_li[:] = nc1.variables['PSD_Num_LI'][:]

#### PSD_Num_MI
psd_num_mi = dataset.createVariable('PSD_Num_MI', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_num_mi.scale_factor = float(1)
psd_num_mi.add_offset = float(0)
psd_num_mi.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 1.2 and 1.4. Particles in contact with the edge of the sample array have been rejected.'
psd_num_mi.units = 'L-1' 
psd_num_mi.long_name = 'Number_size_distribution_of_particles_classed_as_medium_irregularity' 
psd_num_mi[:] = nc1.variables['PSD_Num_MI'][:]

#### PSD_Num_HI
psd_num_hi = dataset.createVariable('PSD_Num_HI', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_num_hi.scale_factor = float(1)
psd_num_hi.add_offset = float(0)
psd_num_hi.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 1.4 and inf. Particles in contact with the edge of the sample array have been rejected.'
psd_num_hi.units = 'L-1' 
psd_num_hi.long_name = 'Number_size_distribution_of_particles_classed_as_high_irregularity' 
psd_num_hi[:] = nc1.variables['PSD_Num_HI'][:]

###################################
## Create particle size distributions - number
###################################
#### PSD_Mass_All
psd_mass_all = dataset.createVariable('PSD_Mass_All', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_mass_all.scale_factor = float(1)
psd_mass_all.add_offset = float(0)
psd_mass_all.comment = 'Particles in contact with the edge of the sample array have been rejected. Sum of small and low, medium, and high irregularity particle categories.'
psd_mass_all.units = 'g m-3' 
psd_mass_all.long_name = 'Mass_size_distribution_of_all_particles' 
psd_mass_all[:] = nc1.variables['PSD_Mass_S'][:] + nc1.variables['PSD_Mass_LI'][:] + nc1.variables['PSD_Mass_MI'][:] + nc1.variables['PSD_Mass_HI'][:]

#### PSD_Mass_S
psd_mass_s = dataset.createVariable('PSD_Mass_S', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_mass_s.scale_factor = float(1)
psd_mass_s.add_offset = float(0)
psd_mass_s.comment = 'This category contains particles with areas between 0 and 30 pixels. Particles in contact with the edge of the sample array have been rejected.'
psd_mass_s.units = 'g m-3' 
psd_mass_s.long_name = 'Mass_size_distribution_of_particles_classed_as_small' 
psd_mass_s[:] = nc1.variables['PSD_Mass_S'][:]

#### PSD_Mass_LI
psd_mass_li = dataset.createVariable('PSD_Mass_LI', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_mass_li.scale_factor = float(1)
psd_mass_li.add_offset = float(0)
psd_mass_li.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 0.9 and 1.2. Particles in contact with the edge of the sample array have been rejected.'
psd_mass_li.units = 'g m-3' 
psd_mass_li.long_name = 'Mass_size_distribution_of_particles_classed_as_low_irregularity' 
psd_mass_li[:] = nc1.variables['PSD_Mass_LI'][:]

#### PSD_Mass_MI
psd_mass_mi = dataset.createVariable('PSD_Mass_MI', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_mass_mi.scale_factor = float(1)
psd_mass_mi.add_offset = float(0)
psd_mass_mi.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 1.2 and 1.4. Particles in contact with the edge of the sample array have been rejected.'
psd_mass_mi.units = 'g m-3' 
psd_mass_mi.long_name = 'Mass_size_distribution_of_particles_classed_as_medium_irregularity' 
psd_mass_mi[:] = nc1.variables['PSD_Mass_MI'][:]

#### PSD_Mass_HI
psd_mass_hi = dataset.createVariable('PSD_Mass_HI', np.float64, ('Time_mid','Size_mid'),fill_value='-9999') 
psd_mass_hi.scale_factor = float(1)
psd_mass_hi.add_offset = float(0)
psd_mass_hi.comment = 'This category contains particles with areas between 30 and inf pixels. A circularity between 1.4 and inf. Particles in contact with the edge of the sample array have been rejected.'
psd_mass_hi.units = 'g m-3' 
psd_mass_hi.long_name = 'Mass_size_distribution_of_particles_classed_as_high_irregularity' 
psd_mass_hi[:] = nc1.variables['PSD_Mass_HI'][:]


###################################
## Write out file
###################################
dataset.close()