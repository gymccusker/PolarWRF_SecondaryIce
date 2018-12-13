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
swsurf1 = np.nanmean(np.nanmean((swdnb1 - swdnbc1), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1 = nc1.variables['LWDNB'][:,:,:] - nc1.variables['LWDNBC'][:,:,:] - nc1.variables['LWUPB'][:,:,:] + nc1.variables['LWUPBC'][:,:,:]
crfsurfts1_1km = swsurf1 + np.nanmean(np.nanmean(lwsurf1,2),1)             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

lwsurf1_map = np.nanmean((nc1.variables['LWDNB'][:,:,:] - nc1.variables['LWDNBC'][:,:,:] - nc1.variables['LWUPB'][:,:,:] + nc1.variables['LWUPBC'][:,:,:]), 0)
crfsurf1_1km = swdown1 + lwsurf1_map             ### CRF AT SURFACE  - incident SW only  

swsurf1a = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,:,:] - nc1a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1a = nc1a.variables['LWDNB'][:,:,:] - nc1a.variables['LWDNBC'][:,:,:] - nc1a.variables['LWUPB'][:,:,:] + nc1a.variables['LWUPBC'][:,:,:]
crfsurfts1_5km = swsurf1a + np.nanmean(np.nanmean(lwsurf1a,2),1)             ### CRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

swsurf1a_5ov1 = np.nanmean(np.nanmean((nc1a.variables['SWDNB'][:,43:123,81:145] - nc1a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
lwsurf1a_5ov1 = nc1a.variables['LWDNB'][:,43:123,81:145] - nc1a.variables['LWDNBC'][:,43:123,81:145] - nc1a.variables['LWUPB'][:,43:123,81:145] + nc1a.variables['LWUPBC'][:,43:123,81:145]

crfsurfts1_5ov1km = swsurf1a_5ov1 + np.nanmean(np.nanmean(lwsurf1a_5ov1,2),1)             ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

lwdw1 = nc1.variables['LWDNB'][:,:,:] - nc1.variables['LWDNBC'][:,:,:]
lwup1 = nc1.variables['LWUPBC'][:,:,:] - nc1.variables['LWUPB'][:,:,:] 
lwdw1a = nc1a.variables['LWDNB'][:,:,:] - nc1a.variables['LWDNBC'][:,:,:]
lwup1a = nc1a.variables['LWUPBC'][:,:,:] - nc1a.variables['LWUPB'][:,:,:] 
lwdw1a_5ov1 = nc1a.variables['LWDNB'][:,43:123,81:145] - nc1a.variables['LWDNBC'][:,43:123,81:145]
lwup1a_5ov1 = nc1a.variables['LWUPBC'][:,43:123,81:145] - nc1a.variables['LWUPB'][:,43:123,81:145]

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

swsurf2a = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,:,:] - nc2a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf2a_5ov1 = np.nanmean(np.nanmean((nc2a.variables['SWDNB'][:,43:123,81:145] - nc2a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

lwdw2 = nc2.variables['LWDNB'][:,:,:] - nc2.variables['LWDNBC'][:,:,:]
lwup2 = nc2.variables['LWUPBC'][:,:,:] - nc2.variables['LWUPB'][:,:,:] 
lwdw2a = nc2a.variables['LWDNB'][:,:,:] - nc2a.variables['LWDNBC'][:,:,:]
lwup2a = nc2a.variables['LWUPBC'][:,:,:] - nc2a.variables['LWUPB'][:,:,:] 
lwdw2a_5ov1 = nc2a.variables['LWDNB'][:,43:123,81:145] - nc2a.variables['LWDNBC'][:,43:123,81:145]
lwup2a_5ov1 = nc2a.variables['LWUPBC'][:,43:123,81:145] - nc2a.variables['LWUPB'][:,43:123,81:145]

swlab2 = "%2.3f" % np.nanmean(swsurf2)
crflab2 = "%2.3f" % np.nanmean(crfsurfts2_1km)
crflab2_5km = "%2.3f" % np.nanmean(crfsurfts2_5km)
crflab2_5ov1km = "%2.3f" % np.nanmean(crfsurfts2_5ov1km)

del nc2, nc2a

runlab2 = 'NoThresh'

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

swsurf3a = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,:,:] - nc3a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf3a_5ov1 = np.nanmean(np.nanmean((nc3a.variables['SWDNB'][:,43:123,81:145] - nc3a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

lwdw3 = nc3.variables['LWDNB'][:,:,:] - nc3.variables['LWDNBC'][:,:,:]
lwup3 = nc3.variables['LWUPBC'][:,:,:] - nc3.variables['LWUPB'][:,:,:] 
lwdw3a = nc3a.variables['LWDNB'][:,:,:] - nc3a.variables['LWDNBC'][:,:,:]
lwup3a = nc3a.variables['LWUPBC'][:,:,:] - nc3a.variables['LWUPB'][:,:,:] 
lwdw3a_5ov1 = nc3a.variables['LWDNB'][:,43:123,81:145] - nc3a.variables['LWDNBC'][:,43:123,81:145]
lwup3a_5ov1 = nc3a.variables['LWUPBC'][:,43:123,81:145] - nc3a.variables['LWUPB'][:,43:123,81:145]

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

swsurf4a = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,:,:] - nc4a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf4a_5ov1 = np.nanmean(np.nanmean((nc4a.variables['SWDNB'][:,43:123,81:145] - nc4a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

lwdw4 = nc4.variables['LWDNB'][:,:,:] - nc4.variables['LWDNBC'][:,:,:]
lwup4 = nc4.variables['LWUPBC'][:,:,:] - nc4.variables['LWUPB'][:,:,:] 
lwdw4a = nc4a.variables['LWDNB'][:,:,:] - nc4a.variables['LWDNBC'][:,:,:]
lwup4a = nc4a.variables['LWUPBC'][:,:,:] - nc4a.variables['LWUPB'][:,:,:] 
lwdw4a_5ov1 = nc4a.variables['LWDNB'][:,43:123,81:145] - nc4a.variables['LWDNBC'][:,43:123,81:145]
lwup4a_5ov1 = nc4a.variables['LWUPBC'][:,43:123,81:145] - nc4a.variables['LWUPB'][:,43:123,81:145]

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

swsurf5a = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,:,:] - nc5a.variables['SWDNBC'][:,:,:]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)
swsurf5a_5ov1 = np.nanmean(np.nanmean((nc5a.variables['SWDNB'][:,43:123,81:145] - nc5a.variables['SWDNBC'][:,43:123,81:145]), 2),1)          ### SWCRF AT SURFACE    (i.e. how much clouds cool the surface => -ve!)

lwdw5 = nc5.variables['LWDNB'][:,:,:] - nc5.variables['LWDNBC'][:,:,:]
lwup5 = nc5.variables['LWUPBC'][:,:,:] - nc5.variables['LWUPB'][:,:,:] 
lwdw5a = nc5a.variables['LWDNB'][:,:,:] - nc5a.variables['LWDNBC'][:,:,:]
lwup5a = nc5a.variables['LWUPBC'][:,:,:] - nc5a.variables['LWUPB'][:,:,:] 
lwdw5a_5ov1 = nc5a.variables['LWDNB'][:,43:123,81:145] - nc5a.variables['LWDNBC'][:,43:123,81:145]
lwup5a_5ov1 = nc5a.variables['LWUPBC'][:,43:123,81:145] - nc5a.variables['LWUPB'][:,43:123,81:145]

swlab5 = "%2.3f" % np.nanmean(swsurf5)
crflab5 = "%2.3f" % np.nanmean(crfsurfts5_1km)
crflab5_5km = "%2.3f" % np.nanmean(crfsurfts5_5km)
crflab5_5ov1km = "%2.3f" % np.nanmean(crfsurfts5_5ov1km)

del nc5, nc5a

runlab5 = '10xHM'


###################################
# Print stats to terminal
###################################


print runlab1,":"
print " CRF_5km = ", np.nanmean(crfsurfts1_5km), "; CRF_5km[1km] = ",np.nanmean(crfsurfts1_5ov1km),"; CRF_1km = ",np.nanmean(crfsurfts1_1km)
print " SW_DW_5km = ", np.nanmean(swsurf1a), "; SW_DW_5km[1km] = ",np.nanmean(swsurf1a_5ov1),"; SW_DW_1km = ",np.nanmean(swsurf1)
print " LW_DW_5km = ", np.nanmean(lwdw1a), "; LW_DW_5km[1km] = ",np.nanmean(lwdw1a_5ov1),"; LW_DW_1km = ",np.nanmean(lwdw1)
print " LW_UP_5km = ", np.nanmean(lwup1a), "; LW_UP_5km[1km] = ",np.nanmean(lwup1a_5ov1)," LW_UP_1km = ",np.nanmean(lwup1) 

print runlab2,":"
print " CRF_5km = ", np.nanmean(crfsurfts2_5km), "; CRF_5km[1km] = ",np.nanmean(crfsurfts2_5ov1km),"; CRF_1km = ",np.nanmean(crfsurfts2_1km)
print " SW_DW_5km = ", np.nanmean(swsurf2a), "; SW_DW_5km[1km] = ",np.nanmean(swsurf2a_5ov1),"; SW_DW_1km = ",np.nanmean(swsurf2)
print " LW_DW_5km = ", np.nanmean(lwdw2a), "; LW_DW_5km[1km] = ",np.nanmean(lwdw2a_5ov1),"; LW_DW_1km = ",np.nanmean(lwdw2)
print " LW_UP_5km = ", np.nanmean(lwup2a), "; LW_UP_5km[1km] = ",np.nanmean(lwup2a_5ov1)," LW_UP_1km = ",np.nanmean(lwup2) 

print runlab3,":"
print " CRF_5km = ", np.nanmean(crfsurfts3_5km), "; CRF_5km[1km] = ",np.nanmean(crfsurfts3_5ov1km),"; CRF_1km = ",np.nanmean(crfsurfts3_1km)
print " SW_DW_5km = ", np.nanmean(swsurf3a), "; SW_DW_5km[1km] = ",np.nanmean(swsurf3a_5ov1),"; SW_DW_1km = ",np.nanmean(swsurf3)
print " LW_DW_5km = ", np.nanmean(lwdw3a), "; LW_DW_5km[1km] = ",np.nanmean(lwdw3a_5ov1),"; LW_DW_1km = ",np.nanmean(lwdw3)
print " LW_UP_5km = ", np.nanmean(lwup3a), "; LW_UP_5km[1km] = ",np.nanmean(lwup3a_5ov1)," LW_UP_1km = ",np.nanmean(lwup3) 

print runlab4,":"
print " CRF_5km = ", np.nanmean(crfsurfts4_5km), "; CRF_5km[1km] = ",np.nanmean(crfsurfts4_5ov1km),"; CRF_1km = ",np.nanmean(crfsurfts4_1km)
print " SW_DW_5km = ", np.nanmean(swsurf4a), "; SW_DW_5km[1km] = ",np.nanmean(swsurf4a_5ov1),"; SW_DW_1km = ",np.nanmean(swsurf4)
print " LW_DW_5km = ", np.nanmean(lwdw4a), "; LW_DW_5km[1km] = ",np.nanmean(lwdw4a_5ov1),"; LW_DW_1km = ",np.nanmean(lwdw4)
print " LW_UP_5km = ", np.nanmean(lwup4a), "; LW_UP_5km[1km] = ",np.nanmean(lwup4a_5ov1)," LW_UP_1km = ",np.nanmean(lwup4) 

print runlab5,":"
print " CRF_5km = ", np.nanmean(crfsurfts5_5km), "; CRF_5km[1km] = ",np.nanmean(crfsurfts5_5ov1km),"; CRF_1km = ",np.nanmean(crfsurfts5_1km)
print " SW_DW_5km = ", np.nanmean(swsurf5a), "; SW_DW_5km[1km] = ",np.nanmean(swsurf5a_5ov1),"; SW_DW_1km = ",np.nanmean(swsurf5)
print " LW_DW_5km = ", np.nanmean(lwdw5a), "; LW_DW_5km[1km] = ",np.nanmean(lwdw5a_5ov1),"; LW_DW_1km = ",np.nanmean(lwdw5)
print " LW_UP_5km = ", np.nanmean(lwup5a), "; LW_UP_5km[1km] = ",np.nanmean(lwup5a_5ov1)," LW_UP_1km = ",np.nanmean(lwup5) 