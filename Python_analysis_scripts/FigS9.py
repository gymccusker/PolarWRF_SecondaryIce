from netCDF4 import Dataset as NetCDFFile
import numpy as np
from datetime import datetime
import constants
from wrf_functions import wrf_load
from wrf_functions import params
from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
from mpl_toolkits.basemap import Basemap, cm

###################################
# Pick file
###################################
filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename2 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename3 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename4 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/57_DeMott_WATSAT_5xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
filename5 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/56_DeMott_WATSAT_10xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'

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
zind3 = 21 
timeindex = 1

# # define plot title
strg1 = '$N_{isg>80}$, $L^{-1}$' 
strg2 = 'W, $ms^{-1}$'
# strg2 = '$Q_{liq}$, $g kg^{-1}$'

###################################
# FILE #1
###################################
data1['theta'] = nc1.variables['T'][time_sci]+300 # potential temperature in K
data1['p'] = (nc1.variables['P'][time_sci]+nc1.variables['PB'][time_sci])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data1['p']/100000)**tempvar       
data1['Tk'] = tempvar0*data1['theta']
data1['rho'] = data1['p']/(constants.R*data1['Tk'])
data1['nisg80'] = nc1.variables['NISG80'][time_sci,:,:,:]*(data1['rho'])*(data1['rho'])
# data1['nisg80'][data1['nisg80']<=0] = np.nan
ph = nc1.variables['PH'][time_sci]
phb = nc1.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data1['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)
# data1['qliq'] = nc1.variables['QCLOUD'][time_sci,:,:,:] + nc1.variables['QRAIN'][time_sci,:,:,:]
# data1['qliq'][data1['qliq']<0]=0

# ni1 = np.nanmean(np.nanmean(data1['nisg80'][0:3,0:16,:,:],0),0)/float(1e3)
# ql1 = np.nanmean(np.nanmean(data1['qliq'][0:3,0:16,:,:],0),0)*float(1e3)

ni1 = data1['nisg80'][timeindex,zind1,:,:]/float(1e3)
t1 = data1['Tk'][timeindex,zind3,:,:]

data1['qke'] = nc1.variables['QKE'][time_sci]
pbl1 = np.where(data1['qke'][timeindex]/2 >= 1e-6)

data1['qcloud'] = nc1.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data1['qcloud'][data1['qcloud']<0]=0
dz1 = data1['Zsci'][23:30,:,:] - data1['Zsci'][22:29,:,:]
tempvar = np.squeeze(data1['qcloud'][timeindex,22:29,:,:])*1e3*dz1*data1['rho'][timeindex,22:29,:,:]
lwp1 = np.nansum(tempvar,0)

data1['qisg'] = (nc1.variables['QICE'][time_sci,:,:,:]+
        nc1.variables['QSNOW'][time_sci,:,:,:]+
        nc1.variables['QGRAUP'][time_sci,:,:,:])
data1['qisg'][data1['qisg']<0]=0

dz1 = data1['Zsci'][23:30,:,:] - data1['Zsci'][22:29,:,:]
tempvar = np.squeeze(data1['qisg'][timeindex,22:29,:,:])*1e3*dz1*data1['rho'][timeindex,22:29,:,:]
iwp1 = np.nansum(tempvar,0)

w1 = nc1.variables['W'][time_sci[timeindex],zind3,:,:]

#del nc1

runlab1 = 'CNTRL'

# if data1['Zsci'][zind1,0,0] < 1000: strg2 = "%3.f" % data1['Zsci'][zind1,0,0]
# if data1['Zsci'][zind1,0,0] > 1000: strg2 = "%4.f" % data1['Zsci'][zind1,0,0]

# if data1['Zsci'][zind2,0,0] < 1000: strg4 = "%3.f" % data1['Zsci'][zind2,0,0]
# if data1['Zsci'][zind2,0,0] > 1000: strg4 = "%4.f" % data1['Zsci'][zind2,0,0]

# title1 = ''.join(['Z',strg2,'m'])
# title2 = ''.join(['Z',strg4,'m'])


###################################
# FILE #2
###################################
data2 = {}
data2['theta'] = nc2.variables['T'][time_sci]+300 # potential temperature in K
data2['p'] = (nc2.variables['P'][time_sci]+nc2.variables['PB'][time_sci])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data2['p']/100000)**tempvar       
data2['Tk'] = tempvar0*data2['theta']
data2['rho'] = data2['p']/(constants.R*data2['Tk'])
data2['nisg80'] = nc2.variables['NISG80'][time_sci,:,:,:]*(data2['rho'])*(data2['rho'])
# data2['nisg80'][data2['nisg80']<=0] = np.nan
data2['qliq'] = nc2.variables['QCLOUD'][time_sci,:,:,:] + nc2.variables['QRAIN'][time_sci,:,:,:]
data2['qliq'][data2['qliq']<0]=0
ph = nc2.variables['PH'][time_sci]
phb = nc2.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data2['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# ni2 = np.nanmean(np.nanmean(data2['nisg80'][0:3,0:16,:,:],0),0)/float(1e3)
# ql2 = np.nanmean(np.nanmean(data2['qliq'][0:3,0:16,:,:],0),0)*float(1e3)

ni2 = data2['nisg80'][1,zind1,:,:]/float(1e3)
t2 = data2['Tk'][1,zind3,:,:]

data2['qcloud'] = nc2.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data2['qcloud'][data2['qcloud']<0]=0
dz2 = data2['Zsci'][23:30,:,:] - data2['Zsci'][22:29,:,:]
tempvar = np.squeeze(data2['qcloud'][timeindex,22:29,:,:])*1e3*dz2*data2['rho'][timeindex,22:29,:,:]
lwp2 = np.nansum(tempvar,0)

data2['qisg'] = (nc2.variables['QICE'][time_sci,:,:,:]+
        nc2.variables['QSNOW'][time_sci,:,:,:]+
        nc2.variables['QGRAUP'][time_sci,:,:,:])
data2['qisg'][data2['qisg']<0]=0

dz2 = data2['Zsci'][23:30,:,:] - data2['Zsci'][22:29,:,:]
tempvar = np.squeeze(data2['qisg'][timeindex,22:29,:,:])*1e3*dz2*data2['rho'][timeindex,22:29,:,:]
iwp2 = np.nansum(tempvar,0)

w2 = nc2.variables['W'][time_sci[timeindex],zind3,:,:]

del nc2
del data2

runlab2 = 'NoThresh'

###################################
# FILE #3
###################################
data3 = {}
data3['theta'] = nc3.variables['T'][time_sci]+300 # potential temperature in K
data3['p'] = (nc3.variables['P'][time_sci]+nc3.variables['PB'][time_sci])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data3['p']/100000)**tempvar       
data3['Tk'] = tempvar0*data3['theta']
data3['rho'] = data3['p']/(constants.R*data3['Tk'])
data3['nisg80'] = nc3.variables['NISG80'][time_sci,:,:,:]*(data3['rho'])*(data3['rho'])
# data3['nisg80'][data3['nisg80']<=0] = np.nan
# data3['qliq'] = nc3.variables['QCLOUD'][time_sci,:,:,:] + nc3.variables['QRAIN'][time_sci,:,:,:]
# data3['qliq'][data3['qliq']<0]=0
ph = nc3.variables['PH'][time_sci]
phb = nc3.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data3['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# ni3 = np.nanmean(np.nanmean(data3['nisg80'][0:3,0:16,:,:],0),0)/float(1e3)
# ql3 = np.nanmean(np.nanmean(data3['qliq'][0:3,0:16,:,:],0),0)*float(1e3)
ni3 = data3['nisg80'][1,zind1,:,:]/float(1e3)
t3 = data3['Tk'][1,zind3,:,:]

data3['qcloud'] = nc3.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data3['qcloud'][data3['qcloud']<0]=0
dz3 = data3['Zsci'][23:30,:,:] - data3['Zsci'][22:29,:,:]
tempvar = np.squeeze(data3['qcloud'][timeindex,22:29,:,:])*1e3*dz3*data3['rho'][timeindex,22:29,:,:]
lwp3 = np.nansum(tempvar,0)

data3['qisg'] = (nc3.variables['QICE'][time_sci,:,:,:]+
        nc3.variables['QSNOW'][time_sci,:,:,:]+
        nc3.variables['QGRAUP'][time_sci,:,:,:])
data3['qisg'][data3['qisg']<0]=0

dz3 = data3['Zsci'][23:30,:,:] - data3['Zsci'][22:29,:,:]
tempvar = np.squeeze(data3['qisg'][timeindex,22:29,:,:])*1e3*dz3*data3['rho'][timeindex,22:29,:,:]
iwp3 = np.nansum(tempvar,0)

w3 = nc3.variables['W'][time_sci[timeindex],zind3,:,:]

del nc3
del data3

runlab3 = '2xHM'

###################################
# FILE #4
###################################
data4 = {}
data4['theta'] = nc4.variables['T'][time_sci]+300 # potential temperature in K
data4['p'] = (nc4.variables['P'][time_sci]+nc4.variables['PB'][time_sci])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data4['p']/100000)**tempvar       
data4['Tk'] = tempvar0*data4['theta']
data4['rho'] = data4['p']/(constants.R*data4['Tk'])
data4['nisg80'] = nc4.variables['NISG80'][time_sci,:,:,:]*(data4['rho'])*(data4['rho'])
# data4['nisg80'][data4['nisg80']<=0] = np.nan
# data4['qliq'] = nc4.variables['QCLOUD'][time_sci,:,:,:] + nc4.variables['QRAIN'][time_sci,:,:,:]
# data4['qliq'][data4['qliq']<0]=0
ph = nc4.variables['PH'][time_sci]
phb = nc4.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data4['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# ni4 = np.nanmean(np.nanmean(data4['nisg80'][0:3,0:16,:,:],0),0)/float(1e3)
# ql4 = np.nanmean(np.nanmean(data4['qliq'][0:3,0:16,:,:],0),0)*float(1e3)

ni4 = data4['nisg80'][1,zind1,:,:]/float(1e3)
t4 = data4['Tk'][1,zind3,:,:]

data4['qcloud'] = nc4.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data4['qcloud'][data4['qcloud']<0]=0
dz4 = data4['Zsci'][23:30,:,:] - data4['Zsci'][22:29,:,:]
tempvar = np.squeeze(data4['qcloud'][timeindex,22:29,:,:])*1e3*dz4*data4['rho'][timeindex,22:29,:,:]
lwp4 = np.nansum(tempvar,0)

data4['qisg'] = (nc4.variables['QICE'][time_sci,:,:,:]+
        nc4.variables['QSNOW'][time_sci,:,:,:]+
        nc4.variables['QGRAUP'][time_sci,:,:,:])
data4['qisg'][data4['qisg']<0]=0

dz4 = data4['Zsci'][23:30,:,:] - data4['Zsci'][22:29,:,:]
tempvar = np.squeeze(data4['qisg'][timeindex,22:29,:,:])*1e3*dz4*data4['rho'][timeindex,22:29,:,:]
iwp4 = np.nansum(tempvar,0)

w4 = nc4.variables['W'][time_sci[timeindex],zind3,:,:]

del nc4
del data4

runlab4 = '5xHM'

###################################
# FILE #5
###################################
data5 = {}
data5['theta'] = nc5.variables['T'][time_sci]+300 # potential temperature in K
data5['p'] = (nc5.variables['P'][time_sci]+nc5.variables['PB'][time_sci])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data5['p']/100000)**tempvar       
data5['Tk'] = tempvar0*data5['theta']
data5['rho'] = data5['p']/(constants.R*data5['Tk'])
data5['nisg80'] = nc5.variables['NISG80'][time_sci,:,:,:]*(data5['rho'])*(data5['rho'])
# data5['nisg80'][data5['nisg80']<=0] = np.nan
# data5['qliq'] = nc5.variables['QCLOUD'][time_sci,:,:,:] + nc5.variables['QRAIN'][time_sci,:,:,:]
# data5['qliq'][data5['qliq']<0]=0
ph = nc5.variables['PH'][time_sci]
phb = nc5.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data5['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# ni5 = np.nanmean(np.nanmean(data5['nisg80'][0:3,0:16,:,:],0),0)/float(1e3)
# ql5 = np.nanmean(np.nanmean(data5['qliq'][0:3,0:16,:,:],0),0)*float(1e3)

ni5 = data5['nisg80'][1,zind1,:,:]/float(1e3)
t5 = data5['Tk'][1,zind3,:,:]

data5['qcloud'] = nc5.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data5['qcloud'][data5['qcloud']<0]=0
dz5 = data5['Zsci'][23:30,:,:] - data5['Zsci'][22:29,:,:]
tempvar = np.squeeze(data5['qcloud'][timeindex,22:29,:,:])*1e3*dz5*data5['rho'][timeindex,22:29,:,:]
lwp5 = np.nansum(tempvar,0)

data5['qisg'] = (nc5.variables['QICE'][time_sci,:,:,:]+
        nc5.variables['QSNOW'][time_sci,:,:,:]+
        nc5.variables['QGRAUP'][time_sci,:,:,:])
data5['qisg'][data5['qisg']<0]=0

dz5 = data5['Zsci'][23:30,:,:] - data5['Zsci'][22:29,:,:]
tempvar = np.squeeze(data5['qisg'][timeindex,22:29,:,:])*1e3*dz5*data5['rho'][timeindex,22:29,:,:]
iwp5 = np.nansum(tempvar,0)

w5 = nc5.variables['W'][time_sci[timeindex],zind3,:,:]

del nc5
del data5

runlab5 = '10xHM'

###################################
# LOAD FLIGHT DATA
###################################

data218 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight218/M218_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
data219 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight219/M219_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()

###################################################
###################################################
##### 	OBSERVATIONS
###################################################
###################################################
## science period for flight M218 - all at 27degW
science27 = np.where(np.logical_and(data218['CORE']['Intp_time']>=15.3, data218['CORE']['Intp_time']<=16.7))
newlat27 = data218['CORE']['Intp_lat'][science27]
newlon27 = data218['CORE']['Intp_lon'][science27]

science28 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.45, data219['CORE']['Intp_time']<=21.2))
newlat28 = data219['CORE']['Intp_lat'][science28]
newlon28 = data219['CORE']['Intp_lon'][science28]

science29 = np.where(np.logical_and(data219['CORE']['Intp_time']>=21.3, data219['CORE']['Intp_time']<=22.5))
newlat29 = data219['CORE']['Intp_lat'][science29]
newlon29 = data219['CORE']['Intp_lon'][science29]

##################################################
##################################################
#### 	MODELLED + OBS
##################################################
##################################################

SMALL_SIZE = 9
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
fig = plt.figure(figsize=(8,9))

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
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
maxdat = 5
clevs = np.arange(0.0,np.log10(maxdat + 0.1),0.1)

# data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
data = ni1
data[data == 0] = np.nan
data[data > maxdat] = maxdat

cs = m.contourf(x,y,np.log10(data),clevs,cmap=mpl_cm.binary)

x27,y27 = m(newlon27, newlat27)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab1,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

cbaxes = fig.add_axes([0.15,0.74,0.02, 0.2])  # This is the position for the colorbar
cb = plt.colorbar(cs, ticks=clevs, cax = cbaxes)
tcks = np.power(10,clevs)
cb.ax.set_yticklabels(np.round(tcks,1))
cb.ax.xaxis.set_label_position('top')
cb.ax.axes.set_xlabel(strg1,color='k',fontsize=10)

###################################

ax  = fig.add_axes([0.5,0.7,0.2,0.3])   # left, bottom, width, height

m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
#maxdat2 = 1
#clevs2 = np.arange(-1,maxdat2 + 0.01,0.2)
maxdat2 = 0.2
mindat2 = -0.2
clevs2 = np.arange(mindat2,maxdat2 + 0.01,0.05)

# data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
data = w1
data[data < mindat2] = mindat2
data[data > maxdat2] = maxdat2

cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)

x27,y27 = m(newlon27, newlat27)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab1,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

cbaxes = fig.add_axes([0.8,0.74,0.02, 0.2])  # This is the position for the colorbar
cb = plt.colorbar(cs, ticks=clevs2, cax = cbaxes)
# tcks = np.power(10,clevs)
# cb.ax.set_yticklabels(clevs2)
cb.ax.xaxis.set_label_position('top')
cb.ax.axes.set_xlabel(strg2,color='k',fontsize=10)

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

data = ni2
data[data > maxdat] = maxdat

cs = m.contourf(x,y,np.log10(data),clevs,cmap=mpl_cm.binary)

plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab2,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

ax  = fig.add_axes([0.28,0.39,0.2,0.3])   # left, bottom, width, height

m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
# maxdat = 0.08
# clevs = np.arange(0.0,maxdat2 + 0.01,0.01)

# data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
data = w2
data[data < mindat2] = mindat2
data[data > maxdat2] = maxdat2

cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)

x27,y27 = m(newlon27, newlat27)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab2,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


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

data = ni3
data[data > maxdat] = maxdat

cs = m.contourf(x,y,np.log10(data),clevs,cmap=mpl_cm.binary)

plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab3,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

ax  = fig.add_axes([0.75,0.39,0.2,0.3])   # left, bottom, width, height

m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
# maxdat = 0.08
# clevs = np.arange(0.0,maxdat2 + 0.01,0.01)

# data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
data = w3
data[data < mindat2] = mindat2
data[data > maxdat2] = maxdat2

cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)

x27,y27 = m(newlon27, newlat27)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab3,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


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

# data = np.nanmean(data1['nisg80'][6:9,zind1,:,:],0)
data = ni4
data[data > maxdat] = maxdat

# contour levels
# clevs = np.arange(0.0,1.1,0.1) 
cs = m.contourf(x,y,np.log10(data),clevs,cmap=mpl_cm.binary)

# add colorbar.
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab4,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

ax  = fig.add_axes([0.28,0.08,0.2,0.3])   # left, bottom, width, height

m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
# maxdat = 0.08
# clevs = np.arange(0.0,maxdat2 + 0.01,0.01)

# data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
data = w4
data[data < mindat2] = mindat2
data[data > maxdat2] = maxdat2

cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)

x27,y27 = m(newlon27, newlat27)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab4,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


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

data = ni5
data[data > maxdat] = maxdat

cs = m.contourf(x,y,np.log10(data),clevs,cmap=mpl_cm.binary)

x29,y29 = m(newlon29, newlat29)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab5,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

###################################

ax  = fig.add_axes([0.75,0.08,0.2,0.3])   # left, bottom, width, height

m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
        width=data1['width_meters'],height=data1['height_meters'],\
        lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# define parallels/meridians
m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
m.drawcoastlines(linewidth=1.)

lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
x, y = m(lons, lats) # compute map proj coordinates.

# contour levels
# maxdat = 0.08
# clevs = np.arange(0.0,maxdat + 0.01,0.01)

# data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
data = w5
data[data < mindat2] = mindat2
data[data > maxdat2] = maxdat2

cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)

x27,y27 = m(newlon27, newlat27)
plt.plot(x27,y27,'r',linewidth=1)
plt.annotate(runlab5,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


###################################
##  CRF
###################################
# ax  = fig.add_axes([0.68,0.52,0.25,0.34])   # left, bottom, width, height
# ax.set_xlim([0,1])
# ax.set_ylim([0,1])
# plt.axis('off')

# swupb1 = nc1.variables['SWUPB'][:,:,:]
# swupbc1 = nc1.variables['SWUPBC'][:,:,:]
# swdnb1 = nc1.variables['SWDNB'][:,:,:]
# swdnbc1 = nc1.variables['SWDNBC'][:,:,:]
# swsurf1 = np.nanmean((swdnb1 - swdnbc1 - swupb1 + swupbc1), 0)		### SWCRF AT SURFACE 	(i.e. how much clouds cool the surface => -ve!)
# swsurf1_dailyav = np.nanmean(swsurf1)
# swsurf1_std = np.nanstd(swdnb1 - swdnbc1 - swupb1 + swupbc1)

# swupb2 = nc2.variables['SWUPB'][:,:,:]
# swupbc2 = nc2.variables['SWUPBC'][:,:,:]
# swdnb2 = nc2.variables['SWDNB'][:,:,:]
# swdnbc2 = nc2.variables['SWDNBC'][:,:,:]
# swsurf2 = np.nanmean((swdnb2 - swdnbc2 - swupb2 + swupbc2), 0)
# swsurf2_dailyav = np.nanmean(swsurf2)
# swsurf2_std = np.nanstd(swdnb2 - swdnbc2 - swupb2 + swupbc2)

# swupb3 = nc3.variables['SWUPB'][:,:,:]
# swupbc3 = nc3.variables['SWUPBC'][:,:,:]
# swdnb3 = nc3.variables['SWDNB'][:,:,:]
# swdnbc3 = nc3.variables['SWDNBC'][:,:,:]
# swsurf3 = np.nanmean((swdnb3 - swdnbc3 - swupb3 + swupbc3), 0)
# swsurf3_dailyav = np.nanmean(swsurf3)
# swsurf3_std = np.nanstd(swdnb3 - swdnbc3 - swupb3 + swupbc3)

# swupb4 = nc4.variables['SWUPB'][:,:,:]
# swupbc4 = nc4.variables['SWUPBC'][:,:,:]
# swdnb4 = nc4.variables['SWDNB'][:,:,:]
# swdnbc4 = nc4.variables['SWDNBC'][:,:,:]
# swsurf4 = np.nanmean((swdnb4 - swdnbc4 - swupb4 + swupbc4), 0)
# swsurf4_dailyav = np.nanmean(swsurf4)
# swsurf4_std = np.nanstd(swdnb4 - swdnbc4 - swupb4 + swupbc4)

# swupb5 = nc5.variables['SWUPB'][:,:,:]
# swupbc5 = nc5.variables['SWUPBC'][:,:,:]
# swdnb5 = nc5.variables['SWDNB'][:,:,:]
# swdnbc5 = nc5.variables['SWDNBC'][:,:,:]
# swsurf5 = np.nanmean((swdnb5 - swdnbc5 - swupb5 + swupbc5), 0)
# swsurf5_dailyav = np.nanmean(swsurf5)
# swsurf5_std = np.nanstd(+swdnb5 - swdnbc5 - swupb5 + swupbc5)

# plt.annotate("CRF_{SW,Surf}=",swsurf1_dailyav," +/- ",swsurf1_std,xy=(0.1,0.9),xytext=(0.11,0.91),fontsize=SMALL_SIZE)

# plt.savefig('/data/scihub-users/giyoung/PYTHON/WRF/FIGS/Misc/05_GRL_IcePatches_W1500m_z21.svg')
plt.savefig('/data/scihub-users/giyoung/PYTHON/WRF/FIGS/Misc/05_GRL_IcePatches_W1500m_z21.png',dpi=300)
plt.show()


