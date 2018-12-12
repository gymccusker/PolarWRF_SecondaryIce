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
from sklearn.metrics import r2_score
from sklearn import linear_model
from scipy import stats

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

# ###################################
# # FILE #1
# ###################################
# data1['theta'] = nc1.variables['T'][time_sci]+300 # potential temperature in K
# data1['p'] = (nc1.variables['P'][time_sci]+nc1.variables['PB'][time_sci])   # pressure in Pa
# tempvar = constants.R/float(1005)
# tempvar0 = (data1['p']/100000)**tempvar       
# data1['Tk'] = tempvar0*data1['theta']
# data1['rho'] = data1['p']/(constants.R*data1['Tk'])
# data1['nisg80'] = nc1.variables['NISG80'][time_sci,:,:,:]*(data1['rho'])*(data1['rho'])
# data1['nisg80'][data1['nisg80'] < 0.005] = np.nan
# ph = nc1.variables['PH'][time_sci]
# phb = nc1.variables['PHB'][time_sci]
# tempvar1 = (ph+phb)/9.81
# data1['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0) # Z at theta mid-point

# data1['qcloud'] = nc1.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data1['qcloud'][data1['qcloud']<0]=0

# data1['qnisg'] = (nc1.variables['QNICE'][time_sci,:,:,:]+
#         nc1.variables['QNSNOW'][time_sci,:,:,:]+
#         nc1.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

# data1['nisg50'] = data1['qnisg'] - (nc1.variables['NI50'][time_sci,:,:,:] - 
#         nc1.variables['NG50'][time_sci,:,:,:])*(data1['rho'])*(data1['rho'])

# data1['qrain'] = nc1.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data1['qrain'][data1['qrain']<0]=0
# data1['qliq'] = data1['qcloud'] + data1['qrain']

# # W_THETA(i,j,k) = 0.5*(W(i,j,k) + W(i,j,k+1))
# w_theta1 = 0.5*(nc1.variables['W'][time_sci[timeindex],0:-1,:,:] + nc1.variables['W'][time_sci[timeindex],1:,:,:])

# ind = {}
# theta = data1['theta'][timeindex,:,:,:]
# Z = data1['Zsci'][:,:,:]
# bl1_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# bl1_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# temp1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# w1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# blindex1 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
# allicebelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# smallicebelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# largeicebelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# liqbelow1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# iceabove1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# # qnisg_above1 = np.zeros(shape=(np.size(Z,0),np.size(Z,1),np.size(Z,2)))
# # Z_above1 = np.zeros(shape=(np.size(Z,0),np.size(Z,1),np.size(Z,2)))
# for i in range(0,np.size(Z,2)):
#         strgi = "%1.f" % (i+1) # string of longitude
#         for j in range(0,np.size(Z,1)):
#                 strgj = "%1.f" % (j+1) # string of latitude
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
#                                 bl1_1[j,i] = Z[k,j,i]
#                                 break
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
#                                 bl1_2[j,i] = Z[k,j,i]
# 				w1[j,i] = w_theta1[k,j,i]
#                                 allicebelow1[j,i] = np.nanpercentile(data1['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
#                                 break
# #del nc1
# # qnisg_above1 = data1['qnisg'][timeindex,blindex]/float(1e3)

# runlab1 = 'CNTRL'


# ## data1['qnisg'][data1['qnisg']<1e-1] = np.nan
# ## plt.plot(np.ndarray.flatten(nc1.variables['W'][time_sci[timeindex],0:22,:,:]),
# ##         np.ndarray.flatten(data1['qnisg'][timeindex,0:22,:,:]/float(1e3)),'.');
# # data1['qnisg'][data1['qnisg']<1e-1] = np.nan
# # plt.plot(np.ndarray.flatten(iceabove1),
# #         np.ndarray.flatten(largeicebelow1),'.');
# # ax = plt.gca();
# # # ax.set_yscale("log", nonposy='clip');
# # plt.show()


# ###################################
# # FILE #2
# ###################################
# data2 = {}
# data2['theta'] = nc2.variables['T'][time_sci]+300 # potential temperature in K
# data2['p'] = (nc2.variables['P'][time_sci]+nc2.variables['PB'][time_sci])   # pressure in Pa
# tempvar = constants.R/float(1005)
# tempvar0 = (data2['p']/100000)**tempvar       
# data2['Tk'] = tempvar0*data2['theta']
# data2['rho'] = data2['p']/(constants.R*data2['Tk'])
# data2['nisg80'] = nc2.variables['NISG80'][time_sci,:,:,:]*(data2['rho'])*(data2['rho'])
# # data2['nisg80'][data2['nisg80']<=0] = np.nan
# data2['qliq'] = nc2.variables['QCLOUD'][time_sci,:,:,:] + nc2.variables['QRAIN'][time_sci,:,:,:]
# data2['qliq'][data2['qliq']<0]=0
# ph = nc2.variables['PH'][time_sci]
# phb = nc2.variables['PHB'][time_sci]
# tempvar1 = (ph+phb)/9.81
# data2['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# data2['qcloud'] = nc2.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data2['qcloud'][data2['qcloud']<0]=0

# data2['qnisg'] = (nc2.variables['QNICE'][time_sci,:,:,:]+
#         nc2.variables['QNSNOW'][time_sci,:,:,:]+
#         nc2.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

# data2['nisg50'] = data2['qnisg'] - (nc2.variables['NI50'][time_sci,:,:,:] - 
#         nc2.variables['NG50'][time_sci,:,:,:])*(data2['rho'])*(data2['rho'])

# data2['qrain'] = nc2.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data2['qrain'][data2['qrain']<0]=0
# data2['qliq'] = data2['qcloud'] + data2['qrain']

# w_theta2 = 0.5*(nc2.variables['W'][time_sci[timeindex],0:-1,:,:] + nc2.variables['W'][time_sci[timeindex],1:,:,:])

# ind = {}
# theta = data2['theta'][timeindex,:,:,:]
# Z = data2['Zsci'][:,:,:]
# bl2_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# bl2_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# temp2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# w2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# blindex2 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
# allicebelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# smallicebelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# largeicebelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# liqbelow2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# iceabove2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# for i in range(0,np.size(Z,2)):
#         strgi = "%1.f" % (i+1) # string of longitude
#         for j in range(0,np.size(Z,1)):
#                 strgj = "%1.f" % (j+1) # string of latitude
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
#                                 bl2_1[j,i] = Z[k,j,i]
#                                 break
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
#                                 bl2_2[j,i] = Z[k,j,i]
#                                 w2[j,i] = w_theta2[k,j,i]
#                                 allicebelow2[j,i] = np.nanpercentile(data2['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
#                                 break

# del nc2
# del data2

# runlab2 = 'NoThresh'

# ###################################
# # FILE #3
# ###################################
# data3 = {}
# data3['theta'] = nc3.variables['T'][time_sci]+300 # potential temperature in K
# data3['p'] = (nc3.variables['P'][time_sci]+nc3.variables['PB'][time_sci])   # pressure in Pa
# tempvar = constants.R/float(1005)
# tempvar0 = (data3['p']/100000)**tempvar       
# data3['Tk'] = tempvar0*data3['theta']
# data3['rho'] = data3['p']/(constants.R*data3['Tk'])
# data3['nisg80'] = nc3.variables['NISG80'][time_sci,:,:,:]*(data3['rho'])*(data3['rho'])
# # data3['nisg80'][data3['nisg80']<=0] = np.nan
# # data3['qliq'] = nc3.variables['QCLOUD'][time_sci,:,:,:] + nc3.variables['QRAIN'][time_sci,:,:,:]
# # data3['qliq'][data3['qliq']<0]=0
# ph = nc3.variables['PH'][time_sci]
# phb = nc3.variables['PHB'][time_sci]
# tempvar1 = (ph+phb)/9.81
# data3['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# data3['qcloud'] = nc3.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data3['qcloud'][data3['qcloud']<0]=0

# data3['qnisg'] = (nc3.variables['QNICE'][time_sci,:,:,:]+
#         nc3.variables['QNSNOW'][time_sci,:,:,:]+
#         nc3.variables['QNGRAUPEL'][time_sci,:,:,:])*(data3['rho'])

# data3['nisg50'] = data3['qnisg'] - (nc3.variables['NI50'][time_sci,:,:,:] - 
#         nc3.variables['NG50'][time_sci,:,:,:])*(data3['rho'])*(data3['rho'])

# data3['qrain'] = nc3.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data3['qrain'][data3['qrain']<0]=0
# data3['qliq'] = data3['qcloud'] + data3['qrain']

# w_theta3 = 0.5*(nc3.variables['W'][time_sci[timeindex],0:-1,:,:] + nc3.variables['W'][time_sci[timeindex],1:,:,:])

# ind = {}
# theta = data3['theta'][timeindex,:,:,:]
# Z = data3['Zsci'][:,:,:]
# bl3_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# bl3_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# temp3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# w3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# blindex3 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
# allicebelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# smallicebelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# largeicebelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# liqbelow3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# iceabove3 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# for i in range(0,np.size(Z,2)):
#         strgi = "%1.f" % (i+1) # string of longitude
#         for j in range(0,np.size(Z,1)):
#                 strgj = "%1.f" % (j+1) # string of latitude
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
#                                 bl3_1[j,i] = Z[k,j,i]
#                                 break
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
#                                 bl3_2[j,i] = Z[k,j,i]
#                                 w3[j,i] = w_theta3[k,j,i]
#                                 allicebelow3[j,i] = np.nanpercentile(data3['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
#                                 break


# del nc3
# del data3

# runlab3 = '2xHM'

# ###################################
# # FILE #4
# ###################################
# data4 = {}
# data4['theta'] = nc4.variables['T'][time_sci]+300 # potential temperature in K
# data4['p'] = (nc4.variables['P'][time_sci]+nc4.variables['PB'][time_sci])   # pressure in Pa
# tempvar = constants.R/float(1005)
# tempvar0 = (data4['p']/100000)**tempvar       
# data4['Tk'] = tempvar0*data4['theta']
# data4['rho'] = data4['p']/(constants.R*data4['Tk'])
# data4['nisg80'] = nc4.variables['NISG80'][time_sci,:,:,:]*(data4['rho'])*(data4['rho'])
# # data4['nisg80'][data4['nisg80']<=0] = np.nan
# # data4['qliq'] = nc4.variables['QCLOUD'][time_sci,:,:,:] + nc4.variables['QRAIN'][time_sci,:,:,:]
# # data4['qliq'][data4['qliq']<0]=0
# ph = nc4.variables['PH'][time_sci]
# phb = nc4.variables['PHB'][time_sci]
# tempvar1 = (ph+phb)/9.81
# data4['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

# data4['qcloud'] = nc4.variables['QCLOUD'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data4['qcloud'][data4['qcloud']<0]=0

# data4['qnisg'] = (nc4.variables['QNICE'][time_sci,:,:,:]+
#         nc4.variables['QNSNOW'][time_sci,:,:,:]+
#         nc4.variables['QNGRAUPEL'][time_sci,:,:,:])*(data1['rho'])

# data4['nisg50'] = data4['qnisg'] - (nc4.variables['NI50'][time_sci,:,:,:] - 
#         nc4.variables['NG50'][time_sci,:,:,:])*(data4['rho'])*(data4['rho'])

# data4['qrain'] = nc4.variables['QRAIN'][time_sci]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
# data4['qrain'][data4['qrain']<0]=0
# data4['qliq'] = data4['qcloud'] + data4['qrain']

# w_theta4 = 0.5*(nc4.variables['W'][time_sci[timeindex],0:-1,:,:] + nc4.variables['W'][time_sci[timeindex],1:,:,:])

# ind = {}
# theta = data4['theta'][timeindex,:,:,:]
# Z = data4['Zsci'][:,:,:]
# bl4_1 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# bl4_2 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# temp4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# w4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# blindex4 = np.zeros(shape=(1,np.size(Z,1),np.size(Z,2)))
# allicebelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# smallicebelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# largeicebelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# liqbelow4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# iceabove4 = np.zeros(shape=(np.size(Z,1),np.size(Z,2)))
# for i in range(0,np.size(Z,2)):
#         strgi = "%1.f" % (i+1) # string of longitude
#         for j in range(0,np.size(Z,1)):
#                 strgj = "%1.f" % (j+1) # string of latitude
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.2:           # small inversion - typically ~500m
#                                 bl4_1[j,i] = Z[k,j,i]
#                                 break
#                 for k in range(2,np.size(Z,0)-2):
#                         strgk = "%1.f" % (k+1) # string of altitude
#                         if theta[k,j,i] < theta[k+1,j,i]-0.4:           # large inversion - typically ~1500m
#                                 bl4_2[j,i] = Z[k,j,i]
#                                 w4[j,i] = w_theta4[k,j,i]
#                                 allicebelow4[j,i] = np.nanpercentile(data4['qnisg'][timeindex,0:k,j,i],99.7)/float(1e3)
#                                 break

# del nc4
# del data4

# runlab4 = '5xHM'

###################################
# FILE #1
###################################
data1 = {}
data1['theta'] = nc1.variables['T'][time_sci]+300 # potential temperature in K
data1['p'] = (nc1.variables['P'][time_sci]+nc1.variables['PB'][time_sci])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data1['p']/100000)**tempvar       
data1['Tk'] = tempvar0*data1['theta']
data1['rho'] = data1['p']/(constants.R*data1['Tk'])
data1['nisg80'] = nc1.variables['NISG80'][time_sci,:,:,:]*(data1['rho'])*(data1['rho'])
# data1['nisg80'][data1['nisg80']<=0] = np.nan
# data1['qliq'] = nc1.variables['QCLOUD'][time_sci,:,:,:] + nc1.variables['QRAIN'][time_sci,:,:,:]
# data1['qliq'][data1['qliq']<0]=0
ph = nc1.variables['PH'][time_sci]
phb = nc1.variables['PHB'][time_sci]
tempvar1 = (ph+phb)/9.81
data1['Zsci'] = 0.5*(tempvar1[:,0:-1,:,:]+tempvar1[:,1:,:,:])

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

w_theta1 = 0.5*(nc1.variables['W'][time_sci,0:-1,:,:] + nc1.variables['W'][time_sci,1:,:,:])

ind = {}
theta = data1['theta'][:,:,:,:]
Z = data1['Zsci'][:,:,:,:]
bl1_1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
bl1_2 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
temp1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
w1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
blindex1 = np.zeros(shape=(np.size(Z,0),1,np.size(Z,2),np.size(Z,3)))
allicebelow1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
smallicebelow1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
largeicebelow1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
liqbelow1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
iceabove1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
blice1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
tropice1 = np.zeros(shape=(np.size(Z,0),np.size(Z,2),np.size(Z,3)))
for t in range(0,np.size(time_sci)):
        for i in range(0,np.size(Z,3)):
                strgi = "%1.f" % (i+1) # string of longitude
                for j in range(0,np.size(Z,2)):
                        strgj = "%1.f" % (j+1) # string of latitude
        		heightindex = np.where(Z[t,:,j,i]<=2500.0)
                        for k in range(2,np.size(Z,1)-3):
                                strgk = "%1.f" % (k+1) # string of altitude
                                if theta[t,k,j,i] < theta[t,k+1,j,i]-0.2:           # small inversion - typically ~500m
                                        bl1_1[t,j,i] = Z[t,k,j,i]
                                        break
                        for k in range(2,np.size(Z,1)-3):
                                strgk = "%1.f" % (k+1) # string of altitude
                                if theta[t,k,j,i] < theta[t,k+1,j,i]-0.4:           # large inversion - typically ~1500m
                                        bl1_2[t,j,i] = np.squeeze(Z[t,k,j,i])
        				w1[t,j,i] = w_theta1[t,k-1,j,i]
                                        allicebelow1[t,j,i] = np.nanmax(data1['qnisg'][t,0:k,j,i])/float(1e3)   # /L
                                        if np.size(data1['qnisg'][t,k:heightindex[0][-1],j,i])>0:
                                                iceabove1[t,j,i] = np.nanmax(data1['qnisg'][t,k:heightindex[0][-1],j,i])/float(1e3) # k to height index (3000m)
                                                # iceabove1[t,j,i] = data1['qnisg'][t,k+2,j,i]/float(1e3) # k+1 only
                                        if np.nanpercentile(data1['qnisg'][t,0:k,j,i],99.7)>0.5:
                                                blice1[t,j,i] = 1.0
                                       	if np.nanpercentile(data1['qnisg'][t,k:heightindex[0][-1],j,i],99.7)>0.5:
                                            	    tropice1[t,j,i] = 1.0
                                        	# tropice1[t,j,i] = 0.0
                                        # if iceabove1[j,i]==np.nan:
                                        #         allicebelow1[j,i] = []
                                        #         iceabove1[j,i] = []
                                        break

del nc1
del data1

runlab1 = '10xHM'

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

binwidth = 0.4
bins = np.arange(-2.0,3.0,binwidth)

ni1 = {}
ni1_med = 0.
ni1_nanmedian = 0.
              


icebelow = np.ndarray.flatten(allicebelow1)
iceabove = np.ndarray.flatten(iceabove1)
watBL = np.ndarray.flatten(w1)

icebelow[icebelow<=0.001] = np.nan
iceabove1[iceabove1<=0.001] = np.nan

# icebelow[icebelow<0.005] = np.nan
# iceabove[iceabove<0.005] = np.nan

iceindex = np.where(np.logical_and(np.ndarray.flatten(blice1)==1,np.ndarray.flatten(tropice1)==1))
print("Percentage ice above+below BL:", np.float(np.size(iceindex[0]))/np.float(np.size(np.ndarray.flatten(blice1)))*100.0)
perc_ice = np.float(np.size(iceindex[0]))/np.float(np.size(np.ndarray.flatten(blice1)))*100.0
perc_ice_lab = "%.1f" % perc_ice
strg1 = '% ice above+below = ' + perc_ice_lab + '%'
# np.float((np.size(bliceindex[0]))/np.float(np.size(np.ndarray.flatten(tropice1))))*100.0


mask1 = ~np.isnan(icebelow) & ~np.isnan(iceabove)
slope1, intercept1, r_value1, p_value1, std_err1 = stats.linregress(iceabove[mask1], icebelow[mask1])
line1 = slope1*iceabove+intercept1
print("r-squared1:", r_value1**2)
r2_lab = "%.2f" % r_value1**2
strg2 = '$R^2$ = '+ r2_lab

mask2 = ~np.isnan(icebelow) & ~np.isnan(watBL)
slope2, intercept2, r_value2, p_value2, std_err2 = stats.linregress(watBL[mask2], icebelow[mask2])
line2 = slope2*watBL+intercept2
print("r-squared2:", r_value2**2)

for i in range(0,len(bins)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(watBL>=bins[i]-binwidth/2.0, watBL<bins[i]+binwidth/2.0))
    ni1[strgi] = icebelow[ind[strgi]];
    ni1[strgi] = ni1[strgi][~np.isnan(ni1[strgi])]
    if i==0:
        ni1_nanpercentile = np.nanpercentile(ni1[strgi],99.7)
        # ni1_array = [[ni1[strgi]]]
    if i>0:
        ni1_med = np.nanpercentile(ni1[strgi],99.7)
        ni1_nanpercentile = np.append(ni1_nanpercentile,ni1_med)  
        # ni1_array = np.append(ni1_array,ni1[strgi])

ni1_array = [[ni1['1'],ni1['2'],ni1['3'],ni1['4'],ni1['5'],ni1['6'],ni1['7'],ni1['8'],ni1['9'],ni1['10'],ni1['11'],ni1['12'],ni1['13']]]
# ,ni1['14'],ni1['15'],ni1['16'],ni1['17'],ni1['18'],ni1['19'],ni1['20']]]

fig = plt.figure(figsize=(8,5))

# Also manually adjust the spacings which are used when creating subplots
plt.gcf().subplots_adjust(top=0.9,bottom=0.1)

plt.subplot(121)
plt.plot(iceabove,icebelow,'.',markersize=2)
plt.plot(iceabove,line1)
plt.grid('on')
#plt.ylim([0.0,1.0])
#plt.xlim([0.0,1.0])
plt.title('CNTRL')
ax = plt.gca();
# ax.set_yscale("log", nonposy='clip'); plt.ylim([1e-10,4e2])
# ax.set_xscale("log", nonposy='clip'); plt.xlim([1e-10,4e2])
plt.ylabel('Max $N_{isg}$ within BL, $L^{-1}$')
plt.xlabel('$N_{isg}$ above BL, $L^{-1}$')
plt.annotate(strg1,xy=(10,320),xytext=(10,321),fontsize=8)
plt.annotate(strg2,xy=(15,300),xytext=(15,301),fontsize=8)

plt.subplot(122)
plt.plot(watBL,icebelow,'.',markersize=2)
# plt.boxplot(ni1_array[0],whis=[5, 95]) # showfliers=False
# plt.plot(watBL,line2,'r-')
plt.plot(bins,ni1_nanpercentile,'-o')
plt.grid('on')
# plt.ylim([0.0,1.0])
plt.title('CNTRL')
plt.xlabel('W, $ms^{-1}$')
# plt.xlim([-1.0,1.5])
# ax = plt.gca();
# a = ax.get_xticks().tolist()
# for m in range(0,len(bins)):
#         a[m] = "%.1f" % bins[m] 
# ax.set_xticklabels(a)

# ax.set_yscale("log", nonposy='clip');
plt.savefig('../Figures/FigS_IceCorrelations_CNTRL.png',dpi=300)
plt.show()
