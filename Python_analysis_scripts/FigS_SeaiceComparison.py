from netCDF4 import Dataset as NetCDFFile
import numpy as np
from   mpl_toolkits.basemap import Basemap
from datetime import datetime
import constants
from wrf_functions import wrf_load
from wrf_functions import params
from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from matplotlib.patches import Polygon

###################################
# Pick file
###################################
filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename2 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename3 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename4 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename5 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'

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
# nc2 = NetCDFFile(filename2, 'r')
# nc3 = NetCDFFile(filename3, 'r')
# nc4 = NetCDFFile(filename4, 'r')
# nc5 = NetCDFFile(filename5, 'r')

#########################################################################################################
#########################################################################################################


###################################
# DOMAIN
###################################

#===============================  DEFINE GEO FILES and specify resolution

nc_dom1 = "/data/scihub-users/giyoung/WRF_V3.6/WPS/geo_em.d01.nc"
nc_dom2 = "/data/scihub-users/giyoung/WRF_V3.6/WPS/geo_em.d02.nc"
nc_dom3 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'

#===============================  OPEN GEO FILES

geo1= NetCDFFile(nc_dom1,'r')
geo2= NetCDFFile(nc_dom2,'r')
geo3= NetCDFFile(nc_dom3,'r')

#=============================== get geo 1 variables

lat1      = geo1.variables['XLAT_M']
lon1      = geo1.variables['XLONG_M']
topo1     = geo1.variables['HGT_M']

lat1  = lat1[0,:,:]
lon1  = lon1[0,:,:]
topo1 = topo1[0,:,:]

n1x = len(geo1.dimensions['west_east']) # dimensions of domain
n1y = len(geo1.dimensions['south_north'])

dx1 = float(geo1.DX) # resolution 
dy1 = float(geo1.DY)

width_meters  = dx1 * (n1x - 1) # dimensions of domain in meters
height_meters = dy1 * (n1y - 1)
                      
cen_lat  = float(geo1.CEN_LAT)
cen_lon  = float(geo1.CEN_LON)
truelat1 = float(geo1.TRUELAT1)
truelat2 = float(geo1.TRUELAT2)
standlon = float(geo1.STAND_LON)

#=============================== get geo 2 variables

lat2      = geo2.variables['XLAT_M']
lon2      = geo2.variables['XLONG_M']
topo2     = geo2.variables['HGT_M']

lat2  = lat2[0,:,:]
lon2  = lon2[0,:,:]
topo2 = topo2[0,:,:]

dx2 = float(geo2.DX)
dy2 = float(geo2.DY)

n2x = len(geo2.dimensions['west_east'])
n2y = len(geo2.dimensions['south_north'])

#=============================== MICRO DOMAIN

lonindex_udom = np.where(np.logical_and(lon2>=-29.5, lon2<=-26.5))

lat3  = lat2[190:340,np.unique(lonindex_udom[1])]
lon3  = lon2[190:340,np.unique(lonindex_udom[1])]

n3x = np.size(lon3,1)
n3y = np.size(lat3,0)

#=============================== get geo 2 variables

dom1 = 'd01'
dom2 = 'd02'

#==============================  Additionnal points
lat_max      = -72.
lat_min_cont = -78.
lat_min_sea  = -76.
lon_max      = -25.
lon_min      = -39.

latH         = -75.583 #Halley
lonH         = -26.65
   
#============================== SEAICE


filename_seaice = '/data/scihub-users/giyoung/MAC/Seaice/27-Nov-2015_Seaice_ComisoBootstrapV2.nc'

nc_seaice = NetCDFFile(filename_seaice, 'r')

lat_cice = nc_seaice.variables['lat'][:]
lon_cice = nc_seaice.variables['lon'][:]
data_cice = np.squeeze(nc_seaice.variables['sic'][:])
   
#============================== SEAICE


filename_mslp = '/data/scihub-users/giyoung/MAC/Reanalyses/DATA/MSLP_27NOV2015_1800_EI_MSLP_CLIVARM.nc'

nc_mslp = NetCDFFile(filename_mslp, 'r')

lat_mslp = nc_mslp.variables['lat'][:]
lon_mslp = nc_mslp.variables['lon'][:]
mslp = np.squeeze(nc_mslp.variables['psl'][:,:])

# Make levels which span the dataset
contourf_levels = np.arange(980,1029,1)
contour_levels = np.arange(940,1040,2)


#########################################################################################################
#########################################################################################################
data1 = {}
data1['x_dim'] = len(nc1.dimensions['west_east'])
data1['y_dim'] = len(nc1.dimensions['south_north'])
data1['seaice'] = nc1.variables['SEAICE'][36] 	# sea ice at 1800 UTC for comparison with NSIDC
data1['xlat'] = nc1.variables['XLAT'][36]
data1['xlon'] = nc1.variables['XLONG'][36]


##################################################
##################################################
#### 	MODELLED + OBS
##################################################
##################################################

SMALL_SIZE = 10
MED_SIZE = 12
LARGE_SIZE = 14

plt.rc('font',size=MED_SIZE)
plt.rc('axes',titlesize=MED_SIZE)
plt.rc('axes',labelsize=MED_SIZE)
plt.rc('xtick',labelsize=SMALL_SIZE)
plt.rc('ytick',labelsize=SMALL_SIZE)
plt.rc('legend',fontsize=SMALL_SIZE)
# plt.rc('figure',titlesize=LARGE_SIZE)


## create figure and axes instances
fig = plt.figure(figsize=(6,4))

#########################################################################################################

ax  = fig.add_axes([0.12,0.12,0.35,0.7])	# left, bottom, width, height
m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=width_meters,height=height_meters,\
        lat_0=cen_lat,lon_0=cen_lon,lat_1=truelat1)#,lat_2=truelat2) truelat2 not used in polar proj.

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),labels=[1,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),labels=[0,0,1,0],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)  

m.drawmapboundary(fill_color='lightgrey')
m.fillcontinents(color='white')

#============================== PLOT sub-domaines 2, 3,...  boundaries in domain 1

xd1_1,yd1_1 = m(lon1[n1y-1,0],lat1[n1y-1,0]) 
xd1_2,yd1_2 = m(lon1[0,0],lat1[0,0]) 
xd1_3,yd1_3 = m(lon1[0,n1x-1],lat1[0,n1x-1]) 
xd1_4,yd1_4 = m(lon1[n1y-1,n1x-1],lat1[n1y-1,n1x-1]) 

xd2_1,yd2_1 = m(lon2[n2y-1,0],lat2[n2y-1,0]) 
xd2_2,yd2_2 = m(lon2[0,0],lat2[0,0]) 
xd2_3,yd2_3 = m(lon2[0,n2x-1],lat2[0,n2x-1]) 
xd2_4,yd2_4 = m(lon2[n2y-1,n2x-1],lat2[n2y-1,n2x-1]) 

xd3_1,yd3_1 = m(lon3[n3y-1,0],lat3[n3y-1,0]) 
xd3_2,yd3_2 = m(lon3[0,0],lat3[0,0]) 
xd3_3,yd3_3 = m(lon3[0,n3x-1],lat3[0,n3x-1]) 
xd3_4,yd3_4 = m(lon3[n3y-1,n3x-1],lat3[n3y-1,n3x-1])

# draw nests
p2 =  Polygon([(xd2_1,yd2_1),(xd2_2,yd2_2),(xd2_3,yd2_3),(xd2_4,yd2_4)],\
              facecolor='none',edgecolor='k',linewidth=1)
plt.gca().add_patch(p2)

p3 =  Polygon([(xd3_1,yd3_1),(xd3_2,yd3_2),(xd3_3,yd3_3),(xd3_4,yd3_4)],\
              facecolor='none',linestyle='--',edgecolor='k',linewidth=1)
plt.gca().add_patch(p3)

plt.annotate(str(int(dx1/1000.))+' km',xy=(0.85,0.95),xycoords='axes fraction',fontsize=10,fontweight='bold')
plt.annotate(str(int(dx2/1000.))+' km',xy=(xd2_2,yd2_1*1.02),xycoords='data',fontsize=10,fontweight='bold')


#============================== SEAICE
x_cice, y_cice = m(lon_cice, lat_cice)

# plot data
data_cice[data_cice >= 1.0] = 1.0
clevs = np.arange(0.4,1.06,0.05)

csf = m.contourf(x_cice,y_cice,data_cice,clevs,cmap=mpl_cm.Blues_r)

#============================== MSLP
 
# make lat/lon grid
lons_mslp, lats_mslp = np.meshgrid(lon_mslp,lat_mslp)
x_mslp, y_mslp = m(lons_mslp, lats_mslp)

cs = m.contour(x_mslp,y_mslp,mslp,contour_levels,colors='grey',linewidth=0.1)
plt.clabel(cs,fontsize=8, inline=True, fmt='%1.f')	
# inline=10,
#============================== ADDITIONAL LOCATIONS

xH,yH = m(lonH, latH)

plt.legend(bbox_to_anchor=(0.3, 0.76, 1., .102), loc=3, ncol=1)

#============================== COLOURBAR

# add colorbar.
cbaxes = fig.add_axes([0.42,0.68, 0.2, 0.02])  # This is the position for the colorbar
cb = plt.colorbar(csf, cax = cbaxes, orientation = 'horizontal')
cb.ax.xaxis.set_ticks_position('top')
cb.ax.xaxis.set_label_position('top')
cb.ax.axes.set_ylabel('Sea ice fraction')

#########################################################################################################

#########################################################################################################

ax  = fig.add_axes([0.12,0.12,0.35,0.7])	# left, bottom, width, height
m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

data = data1['seaice']

cs = m.contourf(x,y,data,clevs,cmap=mpl_cm.Blues_r)

# plt.savefig('../Figures/S_SeaiceComparison.svg')
plt.show()


