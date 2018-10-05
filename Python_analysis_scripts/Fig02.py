from netCDF4 import Dataset as NetCDFFile
import numpy as np
from datetime import datetime
import constants
from wrf_functions import wrf_load
from wrf_functions import params
from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from mpl_toolkits.basemap import Basemap, cm

###################################
# Pick file
###################################
filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename1a = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
filename2 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename2a = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
filename3 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename3a = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
filename4 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename4a = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'
filename5 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename5a = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN/wrfout_d01_2015-11-27_00:00:00'

###################################
# Extract domain number: 
###################################

domainno_start = filename1.find('/wrfout_') + 8
domainno_end = filename1.find('_2015',domainno_start)
domainno1 = filename1[domainno_start:domainno_end]
del domainno_end, domainno_start

###################################
# LOAD NETCDF FILE
###################################

nc1 = NetCDFFile(filename1, 'r')
nc2 = NetCDFFile(filename2, 'r')
nc3 = NetCDFFile(filename3, 'r')
nc4 = NetCDFFile(filename4, 'r')
nc5 = NetCDFFile(filename5, 'r')

nc1a = NetCDFFile(filename1a, 'r')
nc2a = NetCDFFile(filename2a, 'r')
nc3a = NetCDFFile(filename3a, 'r')
nc4a = NetCDFFile(filename4a, 'r')
nc5a = NetCDFFile(filename5a, 'r')

###################################
# DEFINE TEMPERATURE BINNING
###################################

min218 = 264
max218 = 271
T3D = np.arange(np.round(min218),np.round(max218),0.2)

prams = params(T3D)


###################################
###################################
### STANDARD DOMAIN VARIABLES 
###################################
###################################

times = np.arange(0,24,0.5)

time_sci = np.array((31,32,33,40,41,42,43,44,45))
data1 = {}

# x_dim and y_dim are the x and y dimensions 
# of the model domain in gridpoints
data1['x_dim'] = len(nc1.dimensions['west_east'])
data1['y_dim'] = len(nc1.dimensions['south_north'])

# Get the grid spacing
data1['dx'] = float(nc1.DX)
data1['dy'] = float(nc1.DY)

data1['width_meters']  = data1['dx'] * (data1['x_dim'] - 1)
data1['height_meters'] = data1['dy'] * (data1['y_dim'] - 1)

data1['cen_lat']  = float(nc1.CEN_LAT)
data1['cen_lon']  = float(nc1.CEN_LON)
data1['truelat1'] = float(nc1.TRUELAT1)
data1['truelat2'] = float(nc1.TRUELAT2)
data1['standlon'] = float(nc1.STAND_LON)

data1['xlat'] = nc1.variables['XLAT'][time_sci[0]]
data1['xlon'] = nc1.variables['XLONG'][time_sci[0]]

# define data to plot
####  dat1_27/2 => Basemap, lat vs lon, Qcloud @ ~350m [8]
zind1 = 11
# zind1 = 11 # 625m
zind2 = 16 # 1095m

# # define plot title
strg1 = '$CRF_{surf}$, $W m^{-2}$' 
strg2 = '$\Delta CRF_{surf}$, $W m^{-2}$'

###################################
# FILE #1
###################################
swdnb1 = nc1.variables['SWDNB'][:,:,:] # INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM"          "W m-2"
swdnbc1 = nc1.variables['SWDNBC'][:,:,:] # INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM" "W m-2"
swdown1 = np.nanmean((swdnb1 - swdnbc1),0)
swdownts1_1km = np.nanmean(np.nanmean((swdnb1 - swdnbc1),2),1)
swdownts1_5km = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:]),2),1)
swdownts1_5ov1km = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145]),2),1)

###################################
##  CRF
###################################
# swsurf1 = np.nanmean(np.nanmean((swdnb1 - swdnbc1 - nc1.variables['SWUPB'][:,:,:] + nc1.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf1 = np.nanmean(np.nanmean((swdnb1 - swdnbc1), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1 = nc1.variables['LWDNB'][:,:,:] - nc1.variables['LWDNBC'][:,:,:] - nc1.variables['LWUPB'][:,:,:] + nc1.variables['LWUPBC'][:,:,:]
crfsurfts1_1km = swsurf1 + np.nanmean(np.nanmean(lwsurf1,2),1)             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# swsurf1_map = np.nanmean((swdnb1 - swdnbc1 - nc1.variables['SWUPB'][:,:,:] + nc1.variables['SWUPBC'][:,:,:]), 0)        ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1_map = np.nanmean((nc1.variables['LWDNB'][:,:,:] - nc1.variables['LWDNBC'][:,:,:] - nc1.variables['LWUPB'][:,:,:] + nc1.variables['LWUPBC'][:,:,:]), 0)
crfsurf1_1km = swdown1 + lwsurf1_map             ### CRF AT SURFACE  - incident SW only  

# swsurf1a = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:] - nc1a.variables['SWUPB'][:,:,:] + nc1a.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf1a = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1a = nc1a.variables['LWDNB'][:,:,:] - nc1a.variables['LWDNBC'][:,:,:] - nc1a.variables['LWUPB'][:,:,:] + nc1a.variables['LWUPBC'][:,:,:]
crfsurfts1_5km = swsurf1a + np.nanmean(np.nanmean(lwsurf1a,2),1)             ### CRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# swsurf1a_5ov1 = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145] - 
#         nc1a.variables['SWUPB'][:,43:123,81:145] + nc1a.variables['SWUPBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf1a_5ov1 = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1a_5ov1 = nc1a.variables['LWDNB'][:,43:123,81:145] - nc1a.variables['LWDNBC'][:,43:123,81:145] - nc1a.variables['LWUPB'][:,43:123,81:145] + nc1a.variables['LWUPBC'][:,43:123,81:145]
crfsurfts1_5ov1km = swsurf1a_5ov1 + np.nanmean(np.nanmean(lwsurf1a_5ov1,2),1)             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swlab1 = "%2.3f" % np.nanmean(swsurf1)
crflab1 = "%2.3f" % np.nanmean(crfsurfts1_1km)
crflab1_5km = "%2.3f" % np.nanmean(crfsurfts1_5km)
crflab1_5ov1km = "%2.3f" % np.nanmean(crfsurfts1_5ov1km)

# del nc1, nc1a

runlab1 = 'CNTRL'


###################################
# FILE #2
###################################
swdnb2 = nc2.variables['SWDNB'][:,:,:] # INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM"          "W m-2"
swdnbc2 = nc2.variables['SWDNBC'][:,:,:] # INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM" "W m-2"
swdown2 = np.nanmean((swdnb2 - swdnbc2),0)

swdownts2_1km = np.nanmean(np.nanmean((swdnb2 - swdnbc2),2),1)
swdownts2_5km = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,:,:] - nc2a.variables['SWDNBC'][:,:,:]),2),1)
swdownts2_5ov1km = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,43:123,81:145] - nc2a.variables['SWDNBC'][:,43:123,81:145]),2),1)

swdown_anom2 = np.nanmean(((swdnb2 - swdnbc2) - (swdnb1 - swdnbc1)),0)

swdownts_anom2_1km = np.nanmean(np.nanmean(((swdnb2 - swdnbc2) - (swdnb1 - swdnbc1)),2),1)
swdownts_anom2_5km = np.nanmean(np.nanmean(((nc2a.variables['SWDNB'][:,:,:] - nc2a.variables['SWDNBC'][:,:,:]) -
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])),2),1)
swdownts_anom2_5ov1km = np.nanmean(np.nanmean(((nc2a.variables['SWDNB'][:,43:123,81:145] - nc2a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])),2),1)

###################################
##  CRF
###################################

swsurf2 = np.nanmean(np.nanmean((swdnb2 - swdnbc2), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf2 = nc2.variables['LWDNB'][:,:,:] - nc2.variables['LWDNBC'][:,:,:] - nc2.variables['LWUPB'][:,:,:] + nc2.variables['LWUPBC'][:,:,:]
lwsurf2a = nc2a.variables['LWDNB'][:,:,:] - nc2a.variables['LWDNBC'][:,:,:] - nc2a.variables['LWUPB'][:,:,:] + nc2a.variables['LWUPBC'][:,:,:]
lwsurf2a_5ov1 = nc2a.variables['LWDNB'][:,43:123,81:145] - nc2a.variables['LWDNBC'][:,43:123,81:145] - nc2a.variables['LWUPB'][:,43:123,81:145] + nc2a.variables['LWUPBC'][:,43:123,81:145]
# lwsurf2 = np.nanmean(np.nanmean((nc2.variables['LWDNB'][:,:,:] - nc2.variables['LWDNBC'][:,:,:] - nc2.variables['LWUPB'][:,:,:] + nc2.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts2_1km = swsurf2 + lwsurf2             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swsurf_anom2 = np.nanmean(np.nanmean(((swdnb2 - swdnbc2) - (swdnb1 - swdnbc1)), 2),1)
lwsurf_anom2 = np.nanmean(np.nanmean((lwsurf2 - lwsurf1), 2),1)
crfsurfts2_1km = swsurf_anom2 + lwsurf_anom2

swsurf_anommap2 = np.nanmean(((swdnb2 - swdnbc2) - (swdnb1 - swdnbc1)), 0)
lwsurf_anommap2 = np.nanmean((lwsurf2 - lwsurf1), 0)
crfsurf2_map1km = swsurf_anommap2 + lwsurf_anommap2

swsurf_anom2a = np.nanmean(np.nanmean(((nc2a.variables['SWDNB'][:,:,:] - nc2a.variables['SWDNBC'][:,:,:]) - 
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])), 2),1)
lwsurf_anom2a = np.nanmean(np.nanmean((lwsurf2a - lwsurf1a), 2),1)
crfsurfts2_5km = swsurf_anom2a + lwsurf_anom2a

swsurf_anom2a_5ov1 = np.nanmean(np.nanmean(((nc2a.variables['SWDNB'][:,43:123,81:145] - nc2a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])), 2),1)
lwsurf_anom2a_5ov1 = np.nanmean(np.nanmean((lwsurf2a_5ov1 - lwsurf1a_5ov1), 2),1)
crfsurfts2_5ov1km = swsurf_anom2a_5ov1 + lwsurf_anom2a_5ov1


swlab2 = "%2.3f" % np.nanmean(swsurf2)
crflab2 = "%2.3f" % np.nanmean(crfsurfts2_1km)
crflab2_5km = "%2.3f" % np.nanmean(crfsurfts2_5km)
crflab2_5ov1km = "%2.3f" % np.nanmean(crfsurfts2_5ov1km)

del nc2, nc2a

runlab2 = 'NoThresh'

# swsurf2 = np.nanmean(np.nanmean((swdnb2 - swdnbc2 - nc2.variables['SWUPB'][:,:,:] + nc2.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# swsurf1_map = np.nanmean((swdnb1 - swdnbc1 - nc1.variables['SWUPB'][:,:,:] + nc1.variables['SWUPBC'][:,:,:]), 0)        ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf2_map = np.nanmean((nc2.variables['LWDNB'][:,:,:] - nc2.variables['LWDNBC'][:,:,:] - nc2.variables['LWUPB'][:,:,:] + nc2.variables['LWUPBC'][:,:,:]), 0)
# crfsurf2_1km = swdown2 + lwsurf2_map               ### CRF AT SURFACE    - incident SW only  

# swsurf2a = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,:,:] - nc2a.variables['SWDNBC'][:,:,:] - nc2a.variables['SWUPB'][:,:,:] + nc2a.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf2a = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,:,:] - nc2a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf2a = np.nanmean(np.nanmean((nc2a.variables['LWDNB'][:,:,:] - nc2a.variables['LWDNBC'][:,:,:] - nc2a.variables['LWUPB'][:,:,:] + nc2a.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts2_5km = swsurf2a + lwsurf2a             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# swsurf2a_5ov1 = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,43:123,81:145] - nc2a.variables['SWDNBC'][:,43:123,81:145] - 
#         nc2a.variables['SWUPB'][:,43:123,81:145] + nc2a.variables['SWUPBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf2a_5ov1 = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,43:123,81:145] - nc2a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf2a_5ov1 = np.nanmean(np.nanmean((nc2a.variables['LWDNB'][:,43:123,81:145] - nc2a.variables['LWDNBC'][:,43:123,81:145] - 
#         nc2a.variables['LWUPB'][:,43:123,81:145] + nc2a.variables['LWUPBC'][:,43:123,81:145]), 2), 1)
# crfsurfts2_5ov1km = swsurf2a_5ov1 + lwsurf2a_5ov1             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
###################################
# FILE #3
###################################
swdnb3 = nc3.variables['SWDNB'][:,:,:] # INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM"          "W m-2"
swdnbc3 = nc3.variables['SWDNBC'][:,:,:] # INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM" "W m-2"
swdown3 = np.nanmean((swdnb3 - swdnbc3),0)

swdown_anom3 = np.nanmean(((swdnb3 - swdnbc3) - (swdnb1 - swdnbc1)),0)

swdownts3_1km = np.nanmean(np.nanmean((swdnb3 - swdnbc3),2),1)
swdownts3_5km = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:]),2),1)
swdownts3_5ov1km = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145]),2),1)

swdownts_anom3_1km = np.nanmean(np.nanmean(((swdnb3 - swdnbc3) - (swdnb1 - swdnbc1)),2),1)
swdownts_anom3_5km = np.nanmean(np.nanmean(((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:]) -
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])),2),1)
swdownts_anom3_5ov1km = np.nanmean(np.nanmean(((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])),2),1)

###################################
##  CRF
###################################

swsurf3 = np.nanmean(np.nanmean((swdnb3 - swdnbc3), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf3 = nc3.variables['LWDNB'][:,:,:] - nc3.variables['LWDNBC'][:,:,:] - nc3.variables['LWUPB'][:,:,:] + nc3.variables['LWUPBC'][:,:,:]
lwsurf3a = nc3a.variables['LWDNB'][:,:,:] - nc3a.variables['LWDNBC'][:,:,:] - nc3a.variables['LWUPB'][:,:,:] + nc3a.variables['LWUPBC'][:,:,:]
lwsurf3a_5ov1 = nc3a.variables['LWDNB'][:,43:123,81:145] - nc3a.variables['LWDNBC'][:,43:123,81:145] - nc3a.variables['LWUPB'][:,43:123,81:145] + nc3a.variables['LWUPBC'][:,43:123,81:145]
# lwsurf3 = np.nanmean(np.nanmean((nc3.variables['LWDNB'][:,:,:] - nc3.variables['LWDNBC'][:,:,:] - nc3.variables['LWUPB'][:,:,:] + nc3.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurftsts3_1km = swsurf3 + lwsurf3             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swsurf_anom3 = np.nanmean(np.nanmean(((swdnb3 - swdnbc3) - (swdnb1 - swdnbc1)), 2),1)
lwsurf_anom3 = np.nanmean(np.nanmean((lwsurf3 - lwsurf1), 2),1)
crfsurfts3_1km = swsurf_anom3 + lwsurf_anom3

swsurf_anommap3 = np.nanmean(((swdnb3 - swdnbc3) - (swdnb1 - swdnbc1)), 0)
lwsurf_anommap3 = np.nanmean((lwsurf3 - lwsurf1), 0)
crfsurf3_map1km = swsurf_anommap3 + lwsurf_anommap3

swsurf_anom3a = np.nanmean(np.nanmean(((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:]) - 
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])), 2),1)
lwsurf_anom3a = np.nanmean(np.nanmean((lwsurf3a - lwsurf1a), 2),1)
crfsurfts3_5km = swsurf_anom3a + lwsurf_anom3a

swsurf_anom3a_5ov1 = np.nanmean(np.nanmean(((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])), 2),1)
lwsurf_anom3a_5ov1 = np.nanmean(np.nanmean((lwsurf3a_5ov1 - lwsurf1a_5ov1), 2),1)
crfsurfts3_5ov1km = swsurf_anom3a_5ov1 + lwsurf_anom3a_5ov1


# swdownts3_1km = np.nanmean(np.nanmean((swdnb3 - swdnbc3),2),1)
# swdownts3_5km = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:]),2),1)
# swdownts3_5ov1km = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145]),2),1)

# ###################################
# ##  CRF
# ###################################
# # swsurf3 = np.nanmean(np.nanmean((swdnb3 - swdnbc3 - nc3.variables['SWUPB'][:,:,:] + nc3.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf3 = np.nanmean(np.nanmean((swdnb3 - swdnbc3), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf3 = np.nanmean(np.nanmean((nc3.variables['LWDNB'][:,:,:] - nc3.variables['LWDNBC'][:,:,:] - nc3.variables['LWUPB'][:,:,:] + nc3.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts3_1km = swsurf3 + lwsurf3             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# # swsurf1_map = np.nanmean((swdnb1 - swdnbc1 - nc1.variables['SWUPB'][:,:,:] + nc1.variables['SWUPBC'][:,:,:]), 0)        ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf3_map = np.nanmean((nc3.variables['LWDNB'][:,:,:] - nc3.variables['LWDNBC'][:,:,:] - nc3.variables['LWUPB'][:,:,:] + nc3.variables['LWUPBC'][:,:,:]), 0)
# crfsurf3_1km = swdown3 + lwsurf3_map               ### CRF AT SURFACE    - incident SW only  

# # swsurf3a = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:] - nc3a.variables['SWUPB'][:,:,:] + nc3a.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf3a = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf3a = np.nanmean(np.nanmean((nc3a.variables['LWDNB'][:,:,:] - nc3a.variables['LWDNBC'][:,:,:] - nc3a.variables['LWUPB'][:,:,:] + nc3a.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts3_5km = swsurf3a + lwsurf3a             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# # swsurf3a_5ov1 = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145] - 
# #         nc3a.variables['SWUPB'][:,43:123,81:145] + nc3a.variables['SWUPBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf3a_5ov1 = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf3a_5ov1 = np.nanmean(np.nanmean((nc3a.variables['LWDNB'][:,43:123,81:145] - nc3a.variables['LWDNBC'][:,43:123,81:145] - 
#         nc3a.variables['LWUPB'][:,43:123,81:145] + nc3a.variables['LWUPBC'][:,43:123,81:145]), 2), 1)
# crfsurfts3_5ov1km = swsurf3a_5ov1 + lwsurf3a_5ov1             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swlab3 = "%2.3f" % np.nanmean(swsurf3)
crflab3 = "%2.3f" % np.nanmean(crfsurfts3_1km)
crflab3_5km = "%2.3f" % np.nanmean(crfsurfts3_5km)
crflab3_5ov1km = "%2.3f" % np.nanmean(crfsurfts3_5ov1km)

del nc3, nc3a

runlab3 = '2xHM'

###################################
# FILE #4
###################################
swdnb4 = nc4.variables['SWDNB'][:,:,:] # INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM"          "W m-2"
swdnbc4 = nc4.variables['SWDNBC'][:,:,:] # INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM" "W m-2"
swdown4 = np.nanmean((swdnb4 - swdnbc4),0)

swdownts4_1km = np.nanmean(np.nanmean((swdnb4 - swdnbc4),2),1)
swdownts4_5km = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:]),2),1)
swdownts4_5ov1km = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145]),2),1)

swdown_anom4 = np.nanmean(((swdnb4 - swdnbc4) - (swdnb1 - swdnbc1)),0)

swdownts_anom4_1km = np.nanmean(np.nanmean(((swdnb4 - swdnbc4) - (swdnb1 - swdnbc1)),2),1)
swdownts_anom4_5km = np.nanmean(np.nanmean(((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:]) -
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])),2),1)
swdownts_anom4_5ov1km = np.nanmean(np.nanmean(((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])),2),1)

###################################
##  CRF
###################################

swsurf4 = np.nanmean(np.nanmean((swdnb4 - swdnbc4), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf4 = nc4.variables['LWDNB'][:,:,:] - nc4.variables['LWDNBC'][:,:,:] - nc4.variables['LWUPB'][:,:,:] + nc4.variables['LWUPBC'][:,:,:]
lwsurf4a = nc4a.variables['LWDNB'][:,:,:] - nc4a.variables['LWDNBC'][:,:,:] - nc4a.variables['LWUPB'][:,:,:] + nc4a.variables['LWUPBC'][:,:,:]
lwsurf4a_5ov1 = nc4a.variables['LWDNB'][:,43:123,81:145] - nc4a.variables['LWDNBC'][:,43:123,81:145] - nc4a.variables['LWUPB'][:,43:123,81:145] + nc4a.variables['LWUPBC'][:,43:123,81:145]
# lwsurf4 = np.nanmean(np.nanmean((nc4.variables['LWDNB'][:,:,:] - nc4.variables['LWDNBC'][:,:,:] - nc4.variables['LWUPB'][:,:,:] + nc4.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts4_1km = swsurf4 + lwsurf4             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swsurf_anom4 = np.nanmean(np.nanmean(((swdnb4 - swdnbc4) - (swdnb1 - swdnbc1)), 2),1)
lwsurf_anom4 = np.nanmean(np.nanmean((lwsurf4 - lwsurf1), 2),1)
crfsurfts4_1km = swsurf_anom4 + lwsurf_anom4

swsurf_anommap4 = np.nanmean(((swdnb4 - swdnbc4) - (swdnb1 - swdnbc1)), 0)
lwsurf_anommap4 = np.nanmean((lwsurf4 - lwsurf1), 0)
crfsurf4_map1km = swsurf_anommap4 + lwsurf_anommap4

swsurf_anom4a = np.nanmean(np.nanmean(((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:]) - 
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])), 2),1)
lwsurf_anom4a = np.nanmean(np.nanmean((lwsurf4a - lwsurf1a), 2),1)
crfsurfts4_5km = swsurf_anom4a + lwsurf_anom4a

swsurf_anom4a_5ov1 = np.nanmean(np.nanmean(((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])), 2),1)

lwsurf_anom4a_5ov1 = np.nanmean(np.nanmean((lwsurf4a_5ov1 - lwsurf1a_5ov1), 2),1)
crfsurfts4_5ov1km = swsurf_anom4a_5ov1 + lwsurf_anom4a_5ov1

# swdownts4_1km = np.nanmean(np.nanmean((swdnb4 - swdnbc4),2),1)
# swdownts4_5km = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:]),2),1)
# swdownts4_5ov1km = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145]),2),1)

# ###################################
# ##  CRF
# ###################################
# # swsurf4 = np.nanmean(np.nanmean((swdnb4 - swdnbc4 - nc4.variables['SWUPB'][:,:,:] + nc4.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf4 = np.nanmean(np.nanmean((swdnb4 - swdnbc4), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf4 = np.nanmean(np.nanmean((nc4.variables['LWDNB'][:,:,:] - nc4.variables['LWDNBC'][:,:,:] - nc4.variables['LWUPB'][:,:,:] + nc4.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts4_1km = swsurf4 + lwsurf4             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# # swsurf1_map = np.nanmean((swdnb1 - swdnbc1 - nc1.variables['SWUPB'][:,:,:] + nc1.variables['SWUPBC'][:,:,:]), 0)        ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf4_map = np.nanmean((nc4.variables['LWDNB'][:,:,:] - nc4.variables['LWDNBC'][:,:,:] - nc4.variables['LWUPB'][:,:,:] + nc4.variables['LWUPBC'][:,:,:]), 0)
# crfsurf4_1km = swdown4 + lwsurf4_map               ### CRF AT SURFACE    - incident SW only  

# # swsurf4a = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:] - nc4a.variables['SWUPB'][:,:,:] + nc4a.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf4a = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf4a = np.nanmean(np.nanmean((nc4a.variables['LWDNB'][:,:,:] - nc4a.variables['LWDNBC'][:,:,:] - nc4a.variables['LWUPB'][:,:,:] + nc4a.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts4_5km = swsurf4a + lwsurf4a             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# # swsurf4a_5ov1 = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145] - 
# #         nc4a.variables['SWUPB'][:,43:123,81:145] + nc4a.variables['SWUPBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf4a_5ov1 = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf4a_5ov1 = np.nanmean(np.nanmean((nc4a.variables['LWDNB'][:,43:123,81:145] - nc4a.variables['LWDNBC'][:,43:123,81:145] - 
#         nc4a.variables['LWUPB'][:,43:123,81:145] + nc4a.variables['LWUPBC'][:,43:123,81:145]), 2), 1)
# crfsurfts4_5ov1km = swsurf4a_5ov1 + lwsurf4a_5ov1             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swlab4 = "%2.3f" % np.nanmean(swsurf4)
crflab4 = "%2.3f" % np.nanmean(crfsurfts4_1km)
crflab4_5km = "%2.3f" % np.nanmean(crfsurfts4_5km)
crflab4_5ov1km = "%2.3f" % np.nanmean(crfsurfts4_5ov1km)

del nc4, nc4a

runlab4 = '5xHM'

###################################
# FILE #5
###################################
swdnb5 = nc5.variables['SWDNB'][:,:,:] # INSTANTANEOUS DOWNWELLING SHORTWAVE FLUX AT BOTTOM"          "W m-2"
swdnbc5 = nc5.variables['SWDNBC'][:,:,:] # INSTANTANEOUS DOWNWELLING CLEAR SKY SHORTWAVE FLUX AT BOTTOM" "W m-2"
swdown5 = np.nanmean((swdnb5 - swdnbc5),0)

swdownts5_1km = np.nanmean(np.nanmean((swdnb5 - swdnbc5),2),1)
swdownts5_5km = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:]),2),1)
swdownts5_5ov1km = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145]),2),1)

swdown_anom5 = np.nanmean(((swdnb5 - swdnbc5) - (swdnb1 - swdnbc1)),0)

swdownts_anom5_1km = np.nanmean(np.nanmean(((swdnb5 - swdnbc5) - (swdnb1 - swdnbc1)),2),1)
swdownts_anom5_5km = np.nanmean(np.nanmean(((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:]) -
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])),2),1)
swdownts_anom5_5ov1km = np.nanmean(np.nanmean(((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])),2),1)

###################################
##  CRF
###################################

swsurf5 = np.nanmean(np.nanmean((swdnb5 - swdnbc5), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf5 = nc5.variables['LWDNB'][:,:,:] - nc5.variables['LWDNBC'][:,:,:] - nc5.variables['LWUPB'][:,:,:] + nc5.variables['LWUPBC'][:,:,:]
lwsurf5a = nc5a.variables['LWDNB'][:,:,:] - nc5a.variables['LWDNBC'][:,:,:] - nc5a.variables['LWUPB'][:,:,:] + nc5a.variables['LWUPBC'][:,:,:]
lwsurf5a_5ov1 = nc5a.variables['LWDNB'][:,43:123,81:145] - nc5a.variables['LWDNBC'][:,43:123,81:145] - nc5a.variables['LWUPB'][:,43:123,81:145] + nc5a.variables['LWUPBC'][:,43:123,81:145]
# lwsurf5 = np.nanmean(np.nanmean((nc5.variables['LWDNB'][:,:,:] - nc5.variables['LWDNBC'][:,:,:] - nc5.variables['LWUPB'][:,:,:] + nc5.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts5_1km = swsurf5 + lwsurf5             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swsurf_anom5 = np.nanmean(np.nanmean(((swdnb5 - swdnbc5) - (swdnb1 - swdnbc1)), 2),1)
lwsurf_anom5 = np.nanmean(np.nanmean((lwsurf5 - lwsurf1), 2),1)
crfsurfts5_1km = swsurf_anom5 + lwsurf_anom5

swsurf_anommap5 = np.nanmean(((swdnb5 - swdnbc5) - (swdnb1 - swdnbc1)), 0)
lwsurf_anommap5 = np.nanmean((lwsurf5 - lwsurf1), 0)
crfsurf5_map1km = swsurf_anommap5 + lwsurf_anommap5

swsurf_anom5a = np.nanmean(np.nanmean(((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:]) - 
	(nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:])), 2),1)
lwsurf_anom5a = np.nanmean(np.nanmean((lwsurf5a - lwsurf1a), 2),1)
crfsurfts5_5km = swsurf_anom5a + lwsurf_anom5a

swsurf_anom5a_5ov1 = np.nanmean(np.nanmean(((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145]) - 
	(nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145])), 2),1)
lwsurf_anom5a_5ov1 = np.nanmean(np.nanmean((lwsurf5a_5ov1 - lwsurf1a_5ov1), 2),1)
crfsurfts5_5ov1km = swsurf_anom5a_5ov1 + lwsurf_anom5a_5ov1

# swdownts5_1km = np.nanmean(np.nanmean((swdnb5 - swdnbc5),2),1)
# swdownts5_5km = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:]),2),1)
# swdownts5_5ov1km = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145]),2),1)

# ###################################
# ##  CRF
# ###################################
# # swsurf5 = np.nanmean(np.nanmean((swdnb5 - swdnbc5 - nc5.variables['SWUPB'][:,:,:] + nc5.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf5 = np.nanmean(np.nanmean((swdnb5 - swdnbc5), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf5 = np.nanmean(np.nanmean((nc5.variables['LWDNB'][:,:,:] - nc5.variables['LWDNBC'][:,:,:] - nc5.variables['LWUPB'][:,:,:] + nc5.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts5_1km = swsurf5 + lwsurf5             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# # swsurf1_map = np.nanmean((swdnb1 - swdnbc1 - nc1.variables['SWUPB'][:,:,:] + nc1.variables['SWUPBC'][:,:,:]), 0)        ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf5_map = np.nanmean((nc5.variables['LWDNB'][:,:,:] - nc5.variables['LWDNBC'][:,:,:] - nc5.variables['LWUPB'][:,:,:] + nc5.variables['LWUPBC'][:,:,:]), 0)
# crfsurf5_1km = swdown5 + lwsurf5_map               ### CRF AT SURFACE    - incident SW only  

# # swsurf5a = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:] - nc5a.variables['SWUPB'][:,:,:] + nc5a.variables['SWUPBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf5a = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf5a = np.nanmean(np.nanmean((nc5a.variables['LWDNB'][:,:,:] - nc5a.variables['LWDNBC'][:,:,:] - nc5a.variables['LWUPB'][:,:,:] + nc5a.variables['LWUPBC'][:,:,:]), 2), 1)
# crfsurfts5_5km = swsurf5a + lwsurf5a             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

# # swsurf5a_5ov1 = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145] - 
# #         nc5a.variables['SWUPB'][:,43:123,81:145] + nc5a.variables['SWUPBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# swsurf5a_5ov1 = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
# lwsurf5a_5ov1 = np.nanmean(np.nanmean((nc5a.variables['LWDNB'][:,43:123,81:145] - nc5a.variables['LWDNBC'][:,43:123,81:145] - 
#         nc5a.variables['LWUPB'][:,43:123,81:145] + nc5a.variables['LWUPBC'][:,43:123,81:145]), 2), 1)
# crfsurfts5_5ov1km = swsurf5a_5ov1 + lwsurf5a_5ov1             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swlab5 = "%2.3f" % np.nanmean(swsurf5)
crflab5 = "%2.3f" % np.nanmean(crfsurfts5_1km)
crflab5_5km = "%2.3f" % np.nanmean(crfsurfts5_5km)
crflab5_5ov1km = "%2.3f" % np.nanmean(crfsurfts5_5ov1km)

del nc5, nc5a

runlab5 = '10xHM'

###################################
# LOAD FLIGHT DATA
###################################

# data218 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight218/M218_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
# data219 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight219/M219_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()

# ###################################################
# ###################################################
# ##### 	OBSERVATIONS
# ###################################################
# ###################################################
# ## science period for flight M218 - all at 27degW
# science27 = np.where(np.logical_and(data218['CORE']['Intp_time']>=15.3, data218['CORE']['Intp_time']<=16.7))
# newlat27 = data218['CORE']['Intp_lat'][science27]
# newlon27 = data218['CORE']['Intp_lon'][science27]

# science28 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.45, data219['CORE']['Intp_time']<=21.2))
# newlat28 = data219['CORE']['Intp_lat'][science28]
# newlon28 = data219['CORE']['Intp_lon'][science28]

# science29 = np.where(np.logical_and(data219['CORE']['Intp_time']>=21.3, data219['CORE']['Intp_time']<=22.5))
# newlat29 = data219['CORE']['Intp_lat'][science29]
# newlon29 = data219['CORE']['Intp_lon'][science29]

##################################################
##################################################
#### 	MODELLED + OBS
##################################################
##################################################

SMALL_SIZE = 9
MED_SIZE = 12
LARGE_SIZE = 14

plt.rc('font',size=SMALL_SIZE)
plt.rc('axes',titlesize=SMALL_SIZE)
plt.rc('axes',labelsize=SMALL_SIZE)
plt.rc('xtick',labelsize=SMALL_SIZE)
plt.rc('ytick',labelsize=SMALL_SIZE)
plt.rc('legend',fontsize=SMALL_SIZE)
# plt.rc('figure',titlesize=LARGE_SIZE)


## create figure and axes instances
fig = plt.figure(figsize=(9,9))

###################################
## 	CNTRL
###################################
ax  = fig.add_axes([0.3,0.7,0.2,0.3])   # left, bottom, width, height

m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.,color='w')

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
maxdat = 20
mindat = -100
clevs = np.arange(mindat,maxdat + 0.1,20)

data = np.copy(crfsurf1_1km)

cs = m.contourf(x,y,data,clevs,cmap=mpl_cm.Blues)

plt.annotate(runlab1,xy=(-78,-28),xytext=(-78,-28),fontsize=10,color='w')

cbaxes = fig.add_axes([0.15,0.74,0.02, 0.2])  # This is the position for the colorbar
cb = plt.colorbar(cs, ticks=clevs, cax = cbaxes)
# tcks = np.power(10,clevs)
# cb.ax.set_yticklabels(np.round(tcks,1))
cb.ax.axes.set_xlabel(strg1,color='k',fontsize=10)
cb.ax.xaxis.set_label_position('top')

###################################

# ax  = fig.add_axes([0.5,0.75,0.15,0.1])   # left, bottom, width, height
# # ax.xticks([], [])
# ax.yaxis.tick_right()
# plt.plot(times,swdownts1_5km,color='red',label='5km')
# plt.plot(times,swdownts1_5ov1km,color='darkorange',label='5km[1km]')
# plt.plot(times,swdownts1_1km,color='blue',label='1km')
# plt.ylabel('$SW^{\downarrow}_{surf}$, \n $W m^{-2}$',fontsize=8)
# plt.xlabel('Time, h',fontsize=8)
# ax.yaxis.set_label_position("right")
# plt.xlim([0,24])
# plt.ylim([-170,-10])
# plt.xticks([0, 6, 12, 18, 24])
# plt.yticks([-170, -130, -90, -50])
# plt.grid('on')
# plt.legend(bbox_to_anchor=(1.7, 0.5, 1., .102), loc=3, ncol=1)
# plt.annotate(swlab1,xy=(16,-164),xytext=(16,-164),fontsize=8,color='b')

###################################

# ax  = fig.add_axes([0.5,0.85,0.15,0.1])   # left, bottom, width, height
ax  = fig.add_axes([0.5,0.77,0.15,0.15])   # left, bottom, width, height
ax.yaxis.tick_right()
plt.plot(times,crfsurfts1_5km,color='red',label='5km')
plt.plot(times,crfsurfts1_5ov1km,color='darkorange',label='5km[1km]')
plt.plot(times,crfsurfts1_1km,color='blue',label='1km')
plt.ylabel('$CRF_{surf}$, $W m^{-2}$',fontsize=8)
ax.yaxis.set_label_position("right")
plt.xlim([0,24])
plt.ylim([-100,40])
plt.xticks([6, 12, 18, 24])
plt.yticks([-80, -40, 0, 40])
plt.grid('on')
plt.xlabel('Time, h',fontsize=8)
# ax.xaxis.set_ticklabels([])
plt.legend(bbox_to_anchor=(1.7, 0.2, 1., .102), loc=3, ncol=1)
plt.annotate(crflab1_5km,xy=(15,25),xytext=(15,25),fontsize=8,color='red')
plt.annotate(crflab1_5ov1km,xy=(15,13),xytext=(15,13),fontsize=8,color='darkorange')
plt.annotate(crflab1,xy=(15,1),xytext=(15,1),fontsize=8,color='b')


###################################
##  NoThresh
###################################
ax  = fig.add_axes([0.08,0.39,0.2,0.3])   # left, bottom, width, height
m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# clevs1 = np.arange(-30,30.1,20)
clevs1 = [-40.,-30.,-20.,-10.,10.,20.,30.,40.]
# clevs1 = [-30.,-20.,-10.,10.,20.,30.]

data = crfsurf2_map1km

cs = m.contourf(x,y,data,clevs1,cmap=mpl_cm.RdBu_r)

plt.annotate(runlab2,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

cbaxes = fig.add_axes([0.35,0.03,0.3, 0.02])  # This is the position for the colorbar
tcks = clevs1
cb2 = plt.colorbar(cs, ticks=tcks, cax = cbaxes,orientation='horizontal')
# cb2.ax.set_yticks(tcks)
cb2.ax.set_yticklabels(tcks)
# cb2.ax.xaxis.set_label_position('top')
# cb2.ax.axes.set_xlabel(strg2,color='k',fontsize=10)
cb2.ax.axes.set_xlabel(strg2,color='k',fontsize=10)
cb2.ax.xaxis.set_label_position('top')

###################################

# ax  = fig.add_axes([0.28,0.44,0.15,0.1])   # left, bottom, width, height
# # ax.xticks([], [])
# ax.yaxis.tick_right()
# plt.plot(times,swdownts2_5km,color='red',label='5km')
# plt.plot(times,swdownts2_5ov1km,color='darkorange',label='5km[1km]')
# plt.plot(times,swdownts2_1km,color='blue',label='1km')
# plt.ylabel('$SW^{\downarrow}_{surf}$',fontsize=8)
# plt.xlabel('Time, h',fontsize=8)
# ax.yaxis.set_label_position("right")
# plt.xlim([0,24])
# plt.ylim([-170,-10])
# plt.xticks([0, 6, 12, 18, 24])
# plt.yticks([-170, -130, -90, -50])
# plt.grid('on')
# plt.annotate(swlab2,xy=(16,-164),xytext=(16,-164),fontsize=8,color='b')

###################################

# ax  = fig.add_axes([0.28,0.54,0.15,0.1])   # left, bottom, width, height
ax  = fig.add_axes([0.28,0.46,0.15,0.15])   # left, bottom, width, height
ax.yaxis.tick_right()
plt.plot(times,crfsurfts2_5km,color='red',label='5km')
plt.plot(times,crfsurfts2_5ov1km,color='darkorange',label='5km[1km]')
plt.plot(times,crfsurfts2_1km,color='blue',label='1km')
plt.ylabel('$\Delta CRF_{surf}$, $W m^{-2}$',fontsize=8)
ax.yaxis.set_label_position("right")
plt.xlim([0,24])
plt.ylim([-4,8])
plt.xticks([6, 12, 18, 24])
plt.yticks([-4,0,4,8])
plt.grid('on')
plt.xlabel('Time, h',fontsize=8)
# ax.xaxis.set_ticklabels([])
plt.annotate(crflab2_5km,xy=(18,6.5),xytext=(18,6.5),fontsize=8,color='red')
plt.annotate(crflab2_5ov1km,xy=(17.5,5.5),xytext=(17.5,5.5),fontsize=8,color='darkorange')
plt.annotate(crflab2,xy=(18,4.5),xytext=(18,4.5),fontsize=8,color='b')



###################################
## 	2xHM
###################################
ax  = fig.add_axes([0.55,0.39,0.2,0.3])    # left, bottom, width, height
m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

## make lat/lon grid
# lons28, lats28 = m3.makegrid(150, len(np.unique(data1['lonindex_udom'][1]))) # get lat/lons of ny by nx evenly space grid.

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

data = crfsurf3_map1km

cs = m.contourf(x,y,data,clevs1,cmap=mpl_cm.RdBu_r)

# plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab3,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

# ax  = fig.add_axes([0.75,0.44,0.15,0.1])   # left, bottom, width, height
# # ax.xticks([], [])
# ax.yaxis.tick_right()
# plt.plot(times,swdownts3_5km,color='red',label='5km')
# plt.plot(times,swdownts3_5ov1km,color='darkorange',label='5km[1km]')
# plt.plot(times,swdownts3_1km,color='blue',label='1km')
# plt.ylabel('$SW^{\downarrow}_{surf}$',fontsize=8)
# plt.xlabel('Time, h',fontsize=8)
# ax.yaxis.set_label_position("right")
# plt.xlim([0,24])
# plt.ylim([-170,-10])
# plt.xticks([0, 6, 12, 18, 24])
# plt.yticks([-170, -130, -90, -50])
# plt.grid('on')
# plt.annotate(swlab3,xy=(16,-164),xytext=(16,-164),fontsize=8,color='b')

###################################

# ax  = fig.add_axes([0.75,0.54,0.15,0.1])   # left, bottom, width, height
ax  = fig.add_axes([0.75,0.46,0.15,0.15])   # left, bottom, width, height
ax.yaxis.tick_right()
plt.plot(times,crfsurfts3_5km,color='red',label='5km')
plt.plot(times,crfsurfts3_5ov1km,color='darkorange',label='5km[1km]')
plt.plot(times,crfsurfts3_1km,color='blue',label='1km')
plt.ylabel('$\Delta CRF_{surf}$, $W m^{-2}$',fontsize=8)
ax.yaxis.set_label_position("right")
plt.xlim([0,24])
plt.ylim([-4,8])
plt.xticks([6, 12, 18, 24])
plt.yticks([-4,0,4,8])
plt.grid('on')
plt.xlabel('Time, h',fontsize=8)
# ax.xaxis.set_ticklabels([])
plt.annotate(crflab3_5km,xy=(18,6.5),xytext=(18,6.5),fontsize=8,color='red')
plt.annotate(crflab3_5ov1km,xy=(17.5,5.5),xytext=(17.5,5.5),fontsize=8,color='darkorange')
plt.annotate(crflab3,xy=(18,4.5),xytext=(18,4.5),fontsize=8,color='b')



###################################
## 	5xHM
###################################

ax  = fig.add_axes([0.08,0.08,0.2,0.3])    # left, bottom, width, height
m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

## make lat/lon grid
# lons29, lats29 = m.makegrid(150, len(np.unique(data1['lonindex_udom'][1]))) # get lat/lons of ny by nx evenly space grid.

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

data = crfsurf4_map1km

# contour levels
# clevs = np.arange(0.0,1.1,0.1) 
cs = m.contourf(x,y,data,clevs1,cmap=mpl_cm.RdBu_r)

# add colorbar.
# plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab4,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

# ax  = fig.add_axes([0.28,0.13,0.15,0.1])   # left, bottom, width, height
# # ax.xticks([], [])
# ax.yaxis.tick_right()
# plt.plot(times,swdownts4_5km,color='red',label='5km')
# plt.plot(times,swdownts4_5ov1km,color='darkorange',label='5km[1km]')
# plt.plot(times,swdownts4_1km,color='blue',label='1km')
# plt.ylabel('$SW^{\downarrow}_{surf}$',fontsize=8)
# plt.xlabel('Time, h',fontsize=8)
# ax.yaxis.set_label_position("right")
# plt.xlim([0,24])
# plt.ylim([-170,-10])
# plt.xticks([0, 6, 12, 18, 24])
# plt.yticks([-170, -130, -90, -50])
# plt.grid('on')
# plt.annotate(swlab4,xy=(16,-164),xytext=(16,-164),fontsize=8,color='b')

###################################

# ax  = fig.add_axes([0.28,0.23,0.15,0.1])   # left, bottom, width, height
ax  = fig.add_axes([0.28,0.15,0.15,0.15])   # left, bottom, width, height
ax.yaxis.tick_right()
plt.plot(times,crfsurfts4_5km,color='red',label='5km')
plt.plot(times,crfsurfts4_5ov1km,color='darkorange',label='5km[1km]')
plt.plot(times,crfsurfts4_1km,color='blue',label='1km')
plt.ylabel('$\Delta CRF_{surf}$, $W m^{-2}$',fontsize=8)
ax.yaxis.set_label_position("right")
plt.xlim([0,24])
plt.ylim([-4,8])
plt.xticks([6, 12, 18, 24])
plt.yticks([-4,0,4,8])
plt.grid('on')
plt.xlabel('Time, h',fontsize=8)
# ax.xaxis.set_ticklabels([])
plt.annotate(crflab4_5km,xy=(17.5,6.5),xytext=(17.5,6.5),fontsize=8,color='red')
plt.annotate(crflab4_5ov1km,xy=(17.5,5.5),xytext=(17.5,5.5),fontsize=8,color='darkorange')
plt.annotate(crflab4,xy=(18,4.5),xytext=(18,4.5),fontsize=8,color='b')


###################################
##  10xHM
###################################
ax  = fig.add_axes([0.55,0.08,0.2,0.3])   # left, bottom, width, height
m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

data = crfsurf5_map1km

cs = m.contourf(x,y,data,clevs1,cmap=mpl_cm.RdBu_r)

# x29,y29 = m(newlon29, newlat29)
# plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab5,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

# ax  = fig.add_axes([0.75,0.13,0.15,0.1])   # left, bottom, width, height
# # ax.xticks([], [])
# ax.yaxis.tick_right()
# plt.plot(times,swdownts5_5km,color='red',label='5km')
# plt.plot(times,swdownts5_5ov1km,color='darkorange',label='5km[1km]')
# plt.plot(times,swdownts5_1km,color='blue',label='1km')
# plt.ylabel('$SW^{\downarrow}_{surf}$',fontsize=8)
# plt.xlabel('Time, h',fontsize=8)
# ax.yaxis.set_label_position("right")
# plt.xlim([0,24])
# plt.ylim([-170,-10])
# plt.xticks([0, 6, 12, 18, 24])
# plt.yticks([-170, -130, -90, -50])
# plt.grid('on')
# plt.annotate(swlab5,xy=(16,-164),xytext=(16,-164),fontsize=8,color='b')

###################################

# ax  = fig.add_axes([0.75,0.23,0.15,0.1])   # left, bottom, width, height
ax  = fig.add_axes([0.75,0.15,0.15,0.15])   # left, bottom, width, height
ax.yaxis.tick_right()
plt.plot(times,crfsurfts5_5km,color='red',label='5km')
plt.plot(times,crfsurfts5_5ov1km,color='darkorange',label='5km[1km]')
plt.plot(times,crfsurfts5_1km,color='blue',label='1km')
plt.ylabel('$\Delta CRF_{surf}$, $W m^{-2}$',fontsize=8)
ax.yaxis.set_label_position("right")
plt.xlim([0,24])
plt.ylim([-4,8])
plt.xticks([6, 12, 18, 24])
plt.yticks([-4,0,4,8])
plt.grid('on')
plt.xlabel('Time, h',fontsize=8)
# ax.xaxis.set_ticklabels([])
plt.annotate(crflab5_5km,xy=(18,6.5),xytext=(18,6.5),fontsize=8,color='red')
plt.annotate(crflab5_5ov1km,xy=(17.5,5.5),xytext=(17.5,5.5),fontsize=8,color='darkorange')
plt.annotate(crflab5,xy=(18,4.5),xytext=(18,4.5),fontsize=8,color='b')

# plt.savefig('/data/scihub-users/giyoung/PYTHON/WRF/FIGS/Misc/CRF_v3.svg')
plt.show()
