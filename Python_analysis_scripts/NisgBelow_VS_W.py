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
strg1 = '$N_{isg}$, $L^{-1}$ \n Below BL' 
strg2 = 'W, $ms^{-1}$ \n At BL'
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
data1['nisg80'][data1['nisg80'] < 0.005] = np.nan
ph = nc1.variables['PH'][time_sci]
phb = nc1.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data1['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0) # Z at theta mid-point

data1['qcloud'] = nc1.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data1['qcloud'][data1['qcloud']<0]=0

data1['qnisg'] = (nc1.variables['QNICE'][time_sci,:,:,:]+
        nc1.variables['QNSNOW'][time_sci,:,:,:]+
        nc1.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

data1['nisg50'] = data1['qnisg'] - (nc1.variables['NI50'][time_sci,:,:,:] - 
        nc1.variables['NG50'][time_sci,:,:,:])*(data1['rho'])*(data1['rho'])

data1['qrain'] = nc1.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data1['qrain'][data1['qrain']<0]=0
data1['qliq'] = data1['qcloud'] + data1['qrain']

# W_THETA(i,j,k) = 0.5*(W(i,j,k) + W(i,j,k+1))
w_theta1 = 0.5*(nc1.variables['W'][time_sci[timeindex],0:-1,:,:] + nc1.variables['W'][time_sci[timeindex],1:,:,:])

ind = {}
theta = data1['theta'][timeindex,:,:,:]
Z = data1['Zsci'][:,:,:]
bl1_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
bl1_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
temp1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
w1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
blindex1 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
allicebelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
smallicebelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
largeicebelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
liqbelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
iceabove1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# qnisg_above1 = np.zeros(shape=(np.size(Z,0),np.size(Z,1),np.size(Z,2)))
# Z_above1 = np.zeros(shape=(np.size(Z,0),np.size(Z,1),np.size(Z,2)))
for i in range(0,np.size(Z,2)):
        strgi = "%1.f" % (i+1) # string of longitude
        for j in range(0,np.size(Z,1)):
                strgj = "%1.f" % (j+1) # string of latitude
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
                                bl1_1[j,i] = Z[k,j,i]
                                break
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
                                bl1_2[j,i] = Z[k,j,i]
				w1[j,i] = w_theta1[k,j,i]
                                allicebelow1[j,i] = np.nanpercentile(data1['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
                                break
#del nc1
# qnisg_above1 = data1['qnisg'][timeindex,blindex]/float(1e3)

runlab1 = 'CNTRL'


## data1['qnisg'][data1['qnisg']<1e-1] = np.nan
## plt.plot(np.ndarray.flatten(nc1.variables['W'][time_sci[timeindex],0:22,:,:]),
##         np.ndarray.flatten(data1['qnisg'][timeindex,0:22,:,:]/float(1e3)),'.');
# data1['qnisg'][data1['qnisg']<1e-1] = np.nan
# plt.plot(np.ndarray.flatten(iceabove1),
#         np.ndarray.flatten(largeicebelow1),'.');
# ax = plt.gca();
# # ax.set_yscale("log", nonposy='clip');
# plt.show()


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

data2['qcloud'] = nc2.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data2['qcloud'][data2['qcloud']<0]=0

data2['qnisg'] = (nc2.variables['QNICE'][time_sci,:,:,:]+
        nc2.variables['QNSNOW'][time_sci,:,:,:]+
        nc2.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

data2['nisg50'] = data2['qnisg'] - (nc2.variables['NI50'][time_sci,:,:,:] - 
        nc2.variables['NG50'][time_sci,:,:,:])*(data2['rho'])*(data2['rho'])

data2['qrain'] = nc2.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data2['qrain'][data2['qrain']<0]=0
data2['qliq'] = data2['qcloud'] + data2['qrain']

w_theta2 = 0.5*(nc2.variables['W'][time_sci[timeindex],0:-1,:,:] + nc2.variables['W'][time_sci[timeindex],1:,:,:])

ind = {}
theta = data2['theta'][timeindex,:,:,:]
Z = data2['Zsci'][:,:,:]
bl2_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
bl2_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
temp2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
w2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
blindex2 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
allicebelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
smallicebelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
largeicebelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
liqbelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
iceabove2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
for i in range(0,np.size(Z,2)):
        strgi = "%1.f" % (i+1) # string of longitude
        for j in range(0,np.size(Z,1)):
                strgj = "%1.f" % (j+1) # string of latitude
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
                                bl2_1[j,i] = Z[k,j,i]
                                break
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
                                bl2_2[j,i] = Z[k,j,i]
                                w2[j,i] = w_theta2[k,j,i]
                                allicebelow2[j,i] = np.nanpercentile(data2['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
                                break

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

data3['qcloud'] = nc3.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data3['qcloud'][data3['qcloud']<0]=0

data3['qnisg'] = (nc3.variables['QNICE'][time_sci,:,:,:]+
        nc3.variables['QNSNOW'][time_sci,:,:,:]+
        nc3.variables['QNGRAUPEL'][time_sci,:,:,:])*(data3['rho'])

data3['nisg50'] = data3['qnisg'] - (nc3.variables['NI50'][time_sci,:,:,:] - 
        nc3.variables['NG50'][time_sci,:,:,:])*(data3['rho'])*(data3['rho'])

data3['qrain'] = nc3.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data3['qrain'][data3['qrain']<0]=0
data3['qliq'] = data3['qcloud'] + data3['qrain']

w_theta3 = 0.5*(nc3.variables['W'][time_sci[timeindex],0:-1,:,:] + nc3.variables['W'][time_sci[timeindex],1:,:,:])

ind = {}
theta = data3['theta'][timeindex,:,:,:]
Z = data3['Zsci'][:,:,:]
bl3_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
bl3_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
temp3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
w3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
blindex3 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
allicebelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
smallicebelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
largeicebelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
liqbelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
iceabove3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
for i in range(0,np.size(Z,2)):
        strgi = "%1.f" % (i+1) # string of longitude
        for j in range(0,np.size(Z,1)):
                strgj = "%1.f" % (j+1) # string of latitude
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
                                bl3_1[j,i] = Z[k,j,i]
                                break
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
                                bl3_2[j,i] = Z[k,j,i]
                                w3[j,i] = w_theta3[k,j,i]
                                allicebelow3[j,i] = np.nanpercentile(data3['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
                                break


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

data4['qcloud'] = nc4.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data4['qcloud'][data4['qcloud']<0]=0

data4['qnisg'] = (nc4.variables['QNICE'][time_sci,:,:,:]+
        nc4.variables['QNSNOW'][time_sci,:,:,:]+
        nc4.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

data4['nisg50'] = data4['qnisg'] - (nc4.variables['NI50'][time_sci,:,:,:] - 
        nc4.variables['NG50'][time_sci,:,:,:])*(data4['rho'])*(data4['rho'])

data4['qrain'] = nc4.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data4['qrain'][data4['qrain']<0]=0
data4['qliq'] = data4['qcloud'] + data4['qrain']

w_theta4 = 0.5*(nc4.variables['W'][time_sci[timeindex],0:-1,:,:] + nc4.variables['W'][time_sci[timeindex],1:,:,:])

ind = {}
theta = data4['theta'][timeindex,:,:,:]
Z = data4['Zsci'][:,:,:]
bl4_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
bl4_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
temp4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
w4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
blindex4 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
allicebelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
smallicebelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
largeicebelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
liqbelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
iceabove4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
for i in range(0,np.size(Z,2)):
        strgi = "%1.f" % (i+1) # string of longitude
        for j in range(0,np.size(Z,1)):
                strgj = "%1.f" % (j+1) # string of latitude
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
                                bl4_1[j,i] = Z[k,j,i]
                                break
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
                                bl4_2[j,i] = Z[k,j,i]
                                w4[j,i] = w_theta4[k,j,i]
                                allicebelow4[j,i] = np.nanpercentile(data4['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
                                break

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

data5['qcloud'] = nc5.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data5['qcloud'][data5['qcloud']<0]=0

data5['qnisg'] = (nc5.variables['QNICE'][time_sci,:,:,:]+
        nc5.variables['QNSNOW'][time_sci,:,:,:]+
        nc5.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

data5['nisg50'] = data5['qnisg'] - (nc5.variables['NI50'][time_sci,:,:,:] - 
        nc5.variables['NG50'][time_sci,:,:,:])*(data5['rho'])*(data5['rho'])

data5['qrain'] = nc5.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
data5['qrain'][data5['qrain']<0]=0
data5['qliq'] = data5['qcloud'] + data5['qrain']

w_theta5 = 0.5*(nc5.variables['W'][time_sci[timeindex],0:-1,:,:] + nc5.variables['W'][time_sci[timeindex],1:,:,:])

ind = {}
theta = data5['theta'][timeindex,:,:,:]
Z = data5['Zsci'][:,:,:]
bl5_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
bl5_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
temp5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
w5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
blindex5 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
allicebelow5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
smallicebelow5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
largeicebelow5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
liqbelow5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
iceabove5 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
for i in range(0,np.size(Z,2)):
        strgi = "%1.f" % (i+1) # string of longitude
        for j in range(0,np.size(Z,1)):
                strgj = "%1.f" % (j+1) # string of latitude
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
                                bl5_1[j,i] = Z[k,j,i]
                                break
                for k in range(2,np.size(Z,0)-2):
                        strgk = "%1.f" % (k+1) # string of altitude
                        if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
                                bl5_2[j,i] = Z[k,j,i]
				w5[j,i] = w_theta5[k,j,i]
                                allicebelow5[j,i] = np.nanpercentile(data5['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
                                break

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


# ## create figure and axes instances
# fig = plt.figure(figsize=(8,9))

# ###################################
# ## 	CNTRL
# ###################################
# ax  = fig.add_axes([0.3,0.7,0.2,0.3])   # left, bottom, width, height

# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # contour levels
# maxdat1 = 10
# mindat1 = 0
# #clevs = np.arange(0.0,maxdat,200)

# # data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
# data = allicebelow1 # w1 # bl1_1
# # data[data == 0] = np.nan
# # data[data > maxdat] = maxdat

# cs = m.pcolor(x,y,data,vmin=mindat1,vmax=maxdat1,cmap=mpl_cm.Reds)

# x27,y27 = m(newlon27, newlat27)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab1,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# #??? from here until ???END lines may have been inserted/deleted
# cbaxes = fig.add_axes([0.15,0.74,0.02, 0.2])  # This is the position for the colorbar
# cb = plt.colorbar(cs, cax = cbaxes)
# # cb = plt.colorbar(cs, ticks=clevs, cax = cbaxes)
# # tcks = np.power(10,clevs)
# # cb.ax.set_yticklabels(np.round(tcks,1))
# cb.ax.xaxis.set_label_position('top')
# cb.ax.axes.set_xlabel(strg1,color='k',fontsize=10)

# ###################################

# ax  = fig.add_axes([0.5,0.7,0.2,0.3])   # left, bottom, width, height

# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # contour levels
# mindat2 = -0.15
# #clevs2 = np.arange(-1,maxdat2 + 0.01,0.2)
# maxdat2 = 0.15
# #clevs2 = np.arange(0,2500.01,500)

# # data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
# data = w1 # bl1_2 #iwp1 #w1
# #data[data == 0] = np.nan
# # data[data > maxdat2] = maxdat2

# cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)
# #cs = m.pcolor(x,y,data,vmin=0,vmax=2500,cmap=mpl_cm.viridis)

# x27,y27 = m(newlon27, newlat27)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab1,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# cbaxes = fig.add_axes([0.8,0.74,0.02, 0.2])  # This is the position for the colorbar
# cb = plt.colorbar(cs, cax = cbaxes)
# # # tcks = np.power(10,clevs)
# # # cb.ax.set_yticklabels(clevs2)
# cb.ax.xaxis.set_label_position('top')
# cb.ax.axes.set_xlabel(strg2,color='k',fontsize=10)

# ###################################
# ##  NoThresh
# ###################################
# ax  = fig.add_axes([0.08,0.39,0.2,0.3])   # left, bottom, width, height
# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# data = allicebelow2 #w2 # bl2_1
# # data[data > maxdat] = maxdat

# cs = m.pcolor(x,y,data,vmin=mindat1,vmax=maxdat1,cmap=mpl_cm.Reds)

# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab2,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# ###################################

# ax  = fig.add_axes([0.28,0.39,0.2,0.3])   # left, bottom, width, height

# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # contour levels
# # maxdat = 0.08
# # clevs = np.arange(0.0,maxdat2 + 0.01,0.01)

# # data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
# data = w2 # bl2_2 # iwp2 #w2
# #data[data == 0] = np.nan
# # data[data > maxdat2] = maxdat2

# cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)
# #cs = m.pcolor(x,y,data,vmin=0,vmax=2500,cmap=mpl_cm.viridis)

# x27,y27 = m(newlon27, newlat27)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab2,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


# ###################################
# ## 	2xHM
# ###################################
# ax  = fig.add_axes([0.55,0.39,0.2,0.3])    # left, bottom, width, height
# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# ## make lat/lon grid
# # lons28, lats28 = m3.makegrid(150, len(np.unique(data1['lonindex_udom'][1]))) # get lat/lons of ny by nx evenly space grid.

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# data = allicebelow3 # w3 # bl3_1
# # data[data > maxdat] = maxdat

# cs = m.pcolor(x,y,data,vmin=mindat1,vmax=maxdat1,cmap=mpl_cm.Reds)

# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab3,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# ###################################

# ax  = fig.add_axes([0.75,0.39,0.2,0.3])   # left, bottom, width, height

# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # contour levels
# # maxdat = 0.08
# # clevs = np.arange(0.0,maxdat2 + 0.01,0.01)

# # data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
# data = w3 #bl3_2 # iwp3 #w3
# #data[data == 0] = np.nan
# # data[data > maxdat2] = maxdat2

# cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)
# # cs = m.pcolor(x,y,data,vmin=0,vmax=2500,cmap=mpl_cm.viridis)

# x27,y27 = m(newlon27, newlat27)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab3,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


# ###################################
# ## 	5xHM
# ###################################

# ax  = fig.add_axes([0.08,0.08,0.2,0.3])    # left, bottom, width, height
# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# ## make lat/lon grid
# # lons29, lats29 = m.makegrid(150, len(np.unique(data1['lonindex_udom'][1]))) # get lat/lons of ny by nx evenly space grid.

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # data = np.nanmean(data1['nisg80'][6:9,zind1,:,:],0)
# data = allicebelow4 # w4 # bl4_1
# # data[data > maxdat] = maxdat

# # contour levels
# # clevs = np.arange(0.0,1.1,0.1) 
# cs = m.pcolor(x,y,data,vmin=mindat1,vmax=maxdat1,cmap=mpl_cm.Reds)

# # add colorbar.
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab4,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# ###################################

# ax  = fig.add_axes([0.28,0.08,0.2,0.3])   # left, bottom, width, height

# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # contour levels
# # maxdat = 0.08
# # clevs = np.arange(0.0,maxdat2 + 0.01,0.01)

# # data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
# data = w4 # bl4_2 # iwp4 #w4
# #data[data == 0] = np.nan
# # data[data > maxdat2] = maxdat2

# cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)
# #cs = m.pcolor(x,y,data,vmin=0,vmax=2500,cmap=mpl_cm.viridis)

# x27,y27 = m(newlon27, newlat27)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab4,xy=(-78,-28),xytext=(-78,-28),fontsize=10)


# ###################################
# ##  10xHM
# ###################################
# ax  = fig.add_axes([0.55,0.08,0.2,0.3])   # left, bottom, width, height
# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[1,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# data = allicebelow5 # w5 # bl5_1
# # data[data > maxdat] = maxdat

# cs = m.pcolor(x,y,data,vmin=mindat1,vmax=maxdat1,cmap=mpl_cm.Reds)

# x29,y29 = m(newlon29, newlat29)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab5,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# ###################################

# ax  = fig.add_axes([0.75,0.08,0.2,0.3])   # left, bottom, width, height

# m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
#         width=data1['width_meters'],height=data1['height_meters'],\
#         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # define parallels/meridians
# m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# m.drawcoastlines(linewidth=1.)

# lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# x, y = m(lons, lats) # compute map proj coordinates.

# # contour levels
# # maxdat = 0.08
# # clevs = np.arange(0.0,maxdat + 0.01,0.01)

# # data = np.nanmean(data1['nisg80'][0:3,zind1,:,:],0)
# data = w5 # bl5_2 #iwp5 #w5
# #data[data == 0] = np.nan
# # data[data > maxdat2] = maxdat2

# cs = m.pcolor(x,y,data,vmin=mindat2,vmax=maxdat2,cmap=mpl_cm.RdBu_r)
# #cs = m.pcolor(x,y,data,vmin=0,vmax=2500,cmap=mpl_cm.viridis)

# x27,y27 = m(newlon27, newlat27)
# plt.plot(x27,y27,'r',linewidth=1)
# plt.annotate(runlab5,xy=(-78,-28),xytext=(-78,-28),fontsize=10)

# plt.savefig('/data/scihub-users/giyoung/PYTHON/WRF/FIGS/Misc/05_GRL_NisgBelow_W_v2.png',dpi=300)
# plt.show()


# # plt.subplot(1,2,1); plt.plot(np.nanmean(data1['theta'][timeindex,0:40,190:340,107:202],2),np.nanmean(data1['Zsci'][0:40,190:340,107:202],2));
# # plt.grid('on'); plt.subplot(1,2,2); 
# # plt.plot(np.nanmean(nc1.variables['W'][time_sci[timeindex],0:40,190:340,107:202],2),np.nanmean(data1['Zsci'][0:40,190:340,107:202],2),'.') ;
# # plt.grid('on'); plt.show()


# ##### Picking out BL top
# # ind[strgi] = np.where(np.logical_and(np.logical_and(W[i,:]>-0.07,W[i,:]<0.07),theta[i+1,:]<theta[i,:]+0.05))


# # ind = {}
# # # theta = data1['theta'][timeindex,0:35,190:340,107:202]
# # # Z = data1['Zsci'][0:35,190:340,107:202]
# # theta = data1['theta'][timeindex,:,:,:]
# # Z = data1['Zsci'][:,:,:]
# # bl_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# # bl_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# # for i in range(0,np.size(Z,2)):
# #         strgi = "%1.f" % (i+1) # string of longitude
# #         for j in range(0,np.size(Z,1)):
# #                 strgj = "%1.f" % (j+1) # string of latitude
# #                 for k in range(0,np.size(Z,0)-2):
# #                         strgk = "%1.f" % (k+1) # string of altitude
# #                         if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
# #                                 bl_1[j,i] = Z[k,j,i]
# #                                 break
# #                 for k in range(0,np.size(Z,0)-2):
# #                         strgk = "%1.f" % (k+1) # string of altitude
# #                         if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
# #                                 bl_2[j,i] = Z[k,j,i]
# #                                 break

# #                         # elif theta[k,j,i] < theta[k+1,j,i]-0.4:         # big jump inversion
# #                         #         bl_2[j,i] = Z[k,j,i]        
# #                                 # break
# #                         # OR: a larger theta jump than ?? deg?

# #                         # ind[strgk] = np.where(theta[k,j,i] < theta[k+1,j,i]-0.5)
# #                         # if np.size(ind[strgk]) > 0: 
# #                         #         index = np.array([ind[strgk][0][0],ind[strgk][1][0]])
# #                         #         bl[index[0],index[1]] = Z[k,index[0],index[1]]

# # plt.subplot(231)
# # m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
# #         width=data1['width_meters'],height=data1['height_meters'],\
# #         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # ##define parallels/meridians
# # m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# # m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# # m.drawcoastlines(linewidth=1.)

# # lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# # x, y = m(lons, lats) # compute map proj coordinates.

# # cs = m.pcolor(x,y,largeicebelow1,vmin=0,vmax=3);
# # # plt.colorbar(cs);
# # plt.title('Large ice below')

# # plt.subplot(232)
# # m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
# #         width=data1['width_meters'],height=data1['height_meters'],\
# #         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # ##define parallels/meridians
# # m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# # m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# # m.drawcoastlines(linewidth=1.)

# # lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# # x, y = m(lons, lats) # compute map proj coordinates.

# # cs = m.pcolor(x,y,liqbelow1,vmin=0,vmax=0.1);
# # # plt.colorbar(cs);
# # plt.title('Mean Liquid below')

# # plt.subplot(233)
# # m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
# #         width=data1['width_meters'],height=data1['height_meters'],\
# #         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # ##define parallels/meridians
# # m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# # m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# # m.drawcoastlines(linewidth=1.)

# # lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# # x, y = m(lons, lats) # compute map proj coordinates.

# # cs = m.pcolor(x,y,smallicebelow1,vmin=0,vmax=3);
# # # plt.colorbar(cs);
# # plt.title('Small ice below')


# # plt.subplot(234)
# # m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
# #         width=data1['width_meters'],height=data1['height_meters'],\
# #         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # ##define parallels/meridians
# # m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# # m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# # m.drawcoastlines(linewidth=1.)

# # lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# # x, y = m(lons, lats) # compute map proj coordinates.

# # cs = m.pcolor(x,y,iceabove1,vmin=0,vmax=3);
# # # plt.colorbar(cs);
# # plt.title('All ice above')


# # plt.subplot(235)
# # m = Basemap(resolution='i',projection='stere', rsphere=6370000.0, \
# #         width=data1['width_meters'],height=data1['height_meters'],\
# #         lat_0=data1['cen_lat'],lon_0=data1['cen_lon'],lat_1=data1['truelat1'])

# # ##define parallels/meridians
# # m.drawparallels(np.arange(-90.,-60.,2.),color='k',labels=[0,0,0,0],linewidth=0.8,fontsize=10)
# # m.drawmeridians(np.arange(-180.,181.,5.),color='k',labels=[0,0,0,1],linewidth=0.8,fontsize=10)
# # m.drawcoastlines(linewidth=1.)

# # lons, lats = m.makegrid(data1['x_dim'], data1['y_dim']) # get lat/lons of ny by nx evenly space grid.
# # x, y = m(lons, lats) # compute map proj coordinates.

# # cs = m.pcolor(x,y,bl1_2,vmin=0,vmax=2500);
# # # plt.colorbar(cs);
# # plt.title('BL Height')

# # plt.show()


# # plt.subplot(1,2,1);plt.plot(theta[0:42,40:50,25],Z[0:42,40:50,25],'o--');plt.grid('on');plt.ylim([0,5000]);
# # plt.subplot(1,2,2);plt.plot(bl5_1[40:50,25]);plt.plot(bl5_2[40:50,25]); plt.grid('on');plt.ylim([0,5000]);
# # plt.show()



bins = np.arange(-1.0,2.0,0.1)

ni1 = {}
ni1_med = 0.
ni1_nanmean = 0.
ni2 = {}
ni2_med = 0.
ni2_nanmean = 0.
ni3 = {}
ni3_med = 0.
ni3_nanmean = 0.
ni4 = {}
ni4_med = 0.
ni4_nanmean = 0.
ni5 = {}
ni5_med = 0.
ni5_nanmean = 0.
for i in range(0,len(bins)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(w1>=bins[i]-0.05, w1<bins[i]+0.05))
    ni1[strgi] = allicebelow1[ind[strgi]];
    if i==0:
        ni1_nanmean = np.nanmean(ni1[strgi])
    if i>0:
        ni1_med = np.nanmean(ni1[strgi])
        ni1_nanmean = np.append(ni1_nanmean,ni1_med)

for i in range(0,len(bins)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(w2>=bins[i]-0.05, w2<bins[i]+0.05))
    ni2[strgi] = allicebelow2[ind[strgi]];
    if i==0:
        ni2_nanmean = np.nanmean(ni2[strgi])
    if i>0:
        ni2_med = np.nanmean(ni2[strgi])
        ni2_nanmean = np.append(ni2_nanmean,ni2_med)

for i in range(0,len(bins)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(w3>=bins[i]-0.05, w3<bins[i]+0.05))
    ni3[strgi] = allicebelow3[ind[strgi]];
    if i==0:
        ni3_nanmean = np.nanmean(ni3[strgi])
    if i>0:
        ni3_med = np.nanmean(ni3[strgi])
        ni3_nanmean = np.append(ni3_nanmean,ni3_med)

for i in range(0,len(bins)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(w4>=bins[i]-0.05, w4<bins[i]+0.05))
    ni4[strgi] = allicebelow4[ind[strgi]];
    if i==0:
        ni4_nanmean = np.nanmean(ni4[strgi])
    if i>0:
        ni4_med = np.nanmean(ni4[strgi])
        ni4_nanmean = np.append(ni4_nanmean,ni4_med)


for i in range(0,len(bins)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(w5>=bins[i]-0.05, w5<bins[i]+0.05))
    ni5[strgi] = allicebelow5[ind[strgi]];
    if i==0:
        ni5_nanmean = np.nanmean(ni5[strgi])
    if i>0:
        ni5_med = np.nanmean(ni5[strgi])
        ni5_nanmean = np.append(ni5_nanmean,ni5_med)                

plt.subplot(231)
plt.plot(np.ndarray.flatten(w1),np.ndarray.flatten(allicebelow1)/np.nanmax(np.ndarray.flatten(allicebelow1)),'.',markersize=2)
plt.plot(bins,ni1_nanmean)
plt.grid('on')
plt.xlim([-1.0,2.0])
plt.title('CNTRL')
plt.ylabel('Normalised total $N_{isg}$, \n $L^{-1}$')

plt.subplot(232)
plt.plot(np.ndarray.flatten(w2),np.ndarray.flatten(allicebelow2)/np.nanmax(np.ndarray.flatten(allicebelow2)),'.',markersize=2)
plt.plot(bins,ni2_nanmean)
plt.grid('on')
plt.xlim([-1.0,2.0])
plt.title('NoThresh')

plt.subplot(234)
plt.plot(np.ndarray.flatten(w3),np.ndarray.flatten(allicebelow3)/np.nanmax(np.ndarray.flatten(allicebelow3)),'.',markersize=2)
plt.plot(bins,ni3_nanmean)
plt.grid('on')
plt.xlim([-1.0,2.0])
plt.title('2xHM')
plt.ylabel('Normalised total $N_{isg}$, \n $L^{-1}$')
plt.xlabel('W, $ms^{-1}$')

plt.subplot(235)
plt.plot(np.ndarray.flatten(w4),np.ndarray.flatten(allicebelow4)/np.nanmax(np.ndarray.flatten(allicebelow4)),'.',markersize=2)
plt.plot(bins,ni4_nanmean)
plt.grid('on')
plt.xlim([-1.0,2.0])
plt.title('5xHM')
plt.xlabel('W, $ms^{-1}$')

plt.subplot(236)
plt.plot(np.ndarray.flatten(w5),np.ndarray.flatten(allicebelow5)/np.nanmax(np.ndarray.flatten(allicebelow5)),'.',markersize=2)
plt.plot(bins,ni5_nanmean/np.nanmax(ni5_nanmean))
plt.grid('on')
plt.xlim([-1.0,2.0])
plt.title('10xHM')
plt.xlabel('W, $ms^{-1}$')

plt.show()
