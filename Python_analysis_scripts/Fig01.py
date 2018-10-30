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

###################################
# DEFINE TEMPERATURE BINNING
###################################

min218 = 264
max218 = 271
T3D = np.arange(np.round(min218),np.round(max218),0.5)

prams = params(T3D)

MAC = {}

######  1 degree binning
# MAC = {'Ni': np.array([         np.nan,          np.nan,          np.nan,          np.nan,
#                  np.nan,          np.nan,          np.nan,          np.nan,
#                  np.nan,          np.nan,          np.nan,  21.96826689,
#           3.13097291,   4.5253612 ,   6.228519  ,   6.68025217,
#           1.80675494,   0.59784199,   4.58989635,   0.19594505,
#           0.67525392,   1.10383036,   0.83145725,   1.02527887,
#           1.38177539,   2.66830338,   2.32073501,   1.9944053 ,
#           1.69453773,   3.34051168,   6.61572263,   1.29867026,
#                  np.nan,          np.nan]),
#  'Ni_997': np.array([          np.nan,           np.nan,           np.nan,           np.nan,
#                   np.nan,           np.nan,           np.nan,           np.nan,
#                   np.nan,           np.nan,           np.nan,   85.63879356,
#            7.54329109,   11.61947521,   37.70058998,   12.83772923,
#           10.8749094 ,   13.61566114,  111.86622199,    0.60738587,
#            3.36478262,    6.23009254,    4.59521796,   14.46106642,
#           24.01726235,   77.13600325,   29.81895953,   16.52320508,
#           40.97757256,   94.92269009,  202.14030993,    2.2221433 ,
#                   np.nan,           np.nan]),
#  'Temp': np.array([ 242.,  243.,  244.,  245.,  246.,  247.,  248.,  249.,  250.,
#          251.,  252.,  253.,  254.,  255.,  256.,  257.,  258.,  259.,
#          260.,  261.,  262.,  263.,  264.,  265.,  266.,  267.,  268.,
#          269.,  270.,  271.,  272.,  273.,  274.,  275.])}

# ######  0.5 degree binning, <0.0088 nan
# MAC = {'Ni': np.array([        np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
#                 np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
#                 np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
#                 np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
#                 np.nan,  4.7939509 , 19.58003293, 26.98484124,  3.41381331,
#          3.51164931,  4.1458934 ,  6.34079878,  4.30670299,  6.44044787,
#          7.38321449,  2.27236623,  2.850047  ,  0.50701689,  0.20072794,
#          0.76649342,  6.85228396,  1.08666838,  0.17847381,  0.35541903,
#          0.75171746,  0.83739474,  1.33952273,  0.94481887,  0.43285198,
#          0.85847229,  1.45952731,  0.80730945,  1.89670776,  1.1952711 ,
#          4.05061299,  2.20088065,  2.78611867,  1.31969669,  2.53736152,
#          1.71804551,  1.51567957,  1.53332914,  4.0773457 ,  6.46200846,
#          2.58872069,  0.82919983,  0.04520275,  2.22574742,         np.nan,
#                 np.nan,         np.nan,         np.nan]),
#  'Ni_997': np.array([           np.nan,            np.nan,            np.nan,            np.nan,
#                    np.nan,            np.nan,            np.nan,            np.nan,
#                    np.nan,            np.nan,            np.nan,            np.nan,
#                    np.nan,            np.nan,            np.nan,            np.nan,
#                    np.nan,            np.nan,            np.nan,            np.nan,
#                    np.nan, 8.93480415e+00, 8.55340363e+01, 8.30454327e+01,
#         7.45882405e+00, 7.45581660e+00, 8.46213101e+00, 3.77041017e+01,
#         1.40231851e+01, 1.20206070e+01, 1.28455528e+01, 1.02491204e+01,
#         1.12594099e+01, 2.27551557e+00, 5.26622545e-01, 1.41257602e+01,
#         1.15824678e+02, 1.16380342e+01, 6.06799656e-01, 1.55293335e+00,
#         3.51515404e+00, 3.88164228e+00, 6.26402925e+00, 4.64161066e+00,
#         1.68631294e+00, 8.93367193e+00, 1.60083057e+01, 5.28493720e+00,
#         1.11617828e+02, 2.28222971e+01, 8.00638273e+01, 3.49230070e+01,
#         3.13747212e+01, 1.16956804e+01, 1.94984672e+01, 3.49473605e+01,
#         1.48864422e+01, 1.91703088e+01, 1.02415623e+02, 1.90974381e+02,
#         8.75802607e+01, 1.62028543e+00, 4.52027529e-02, 2.22574742e+00,
#                    np.nan,            np.nan,            np.nan,            np.nan]),
#  'Temp': np.array([242. , 242.5, 243. , 243.5, 244. , 244.5, 245. , 245.5, 246. ,
#         246.5, 247. , 247.5, 248. , 248.5, 249. , 249.5, 250. , 250.5,
#         251. , 251.5, 252. , 252.5, 253. , 253.5, 254. , 254.5, 255. ,
#         255.5, 256. , 256.5, 257. , 257.5, 258. , 258.5, 259. , 259.5,
#         260. , 260.5, 261. , 261.5, 262. , 262.5, 263. , 263.5, 264. ,
#         264.5, 265. , 265.5, 266. , 266.5, 267. , 267.5, 268. , 268.5,
#         269. , 269.5, 270. , 270.5, 271. , 271.5, 272. , 272.5, 273. ,
#         273.5, 274. , 274.5, 275. , 275.5])}

######  0.5 degree binning, <0.005 = Nan
MAC = {'Ni': np.array([        np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
                np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
                np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
                np.nan,         np.nan,         np.nan,         np.nan,         np.nan,
                np.nan,  4.7939509 , 19.58003293, 26.98484124,  3.41381331,
         3.51164931,  4.1458934 ,  6.34079878,  4.30670299,  6.44044787,
         7.38321449,  2.27236623,  2.850047  ,  0.50701689,  0.20072794,
         0.76649342,  6.85228396,  1.08666838,  0.17847381,  0.35541903,
         0.74483711,  0.82702152,  1.33952273,  0.94481887,  0.43285198,
         0.85565033,  1.45952731,  0.80730945,  1.89670776,  1.18961869,
         4.05061299,  2.20088065,  2.78611867,  1.31969669,  2.53736152,
         1.71636388,  1.51567957,  1.53332914,  4.0773457 ,  6.46200846,
         2.58872069,  0.82919983,  0.04520275,  2.22574742,         np.nan,
                np.nan,         np.nan,         np.nan]),
 'Ni_997': np.array([           np.nan,            np.nan,            np.nan,            np.nan,
                   np.nan,            np.nan,            np.nan,            np.nan,
                   np.nan,            np.nan,            np.nan,            np.nan,
                   np.nan,            np.nan,            np.nan,            np.nan,
                   np.nan,            np.nan,            np.nan,            np.nan,
                   np.nan, 8.93480415e+00, 8.55340363e+01, 8.30454327e+01,
        7.45882405e+00, 7.45581660e+00, 8.46213101e+00, 3.77041017e+01,
        1.40231851e+01, 1.20206070e+01, 1.28455528e+01, 1.02491204e+01,
        1.12594099e+01, 2.27551557e+00, 5.26622545e-01, 1.41257602e+01,
        1.15824678e+02, 1.16380342e+01, 6.06799656e-01, 1.55293335e+00,
        3.51325060e+00, 3.88083247e+00, 6.26402925e+00, 4.64161066e+00,
        1.68631294e+00, 8.92482084e+00, 1.60083057e+01, 5.28493720e+00,
        1.11617828e+02, 2.26609651e+01, 8.00638273e+01, 3.49230070e+01,
        3.13747212e+01, 1.16956804e+01, 1.94984672e+01, 3.49328818e+01,
        1.48864422e+01, 1.91703088e+01, 1.02415623e+02, 1.90974381e+02,
        8.75802607e+01, 1.62028543e+00, 4.52027529e-02, 2.22574742e+00,
                   np.nan,            np.nan,            np.nan,            np.nan]),
 'Temp': np.array([242. , 242.5, 243. , 243.5, 244. , 244.5, 245. , 245.5, 246. ,
        246.5, 247. , 247.5, 248. , 248.5, 249. , 249.5, 250. , 250.5,
        251. , 251.5, 252. , 252.5, 253. , 253.5, 254. , 254.5, 255. ,
        255.5, 256. , 256.5, 257. , 257.5, 258. , 258.5, 259. , 259.5,
        260. , 260.5, 261. , 261.5, 262. , 262.5, 263. , 263.5, 264. ,
        264.5, 265. , 265.5, 266. , 266.5, 267. , 267.5, 268. , 268.5,
        269. , 269.5, 270. , 270.5, 271. , 271.5, 272. , 272.5, 273. ,
        273.5, 274. , 274.5, 275. , 275.5])}


###################################
# PROCESS WRF DATA FOR USE
###################################

bins = np.arange(0,91,0.2)
time_sci = np.array((31,32,33,40,41,42,43,44,45))
data1 = {}
data1['xlat'] = nc1.variables['XLAT'][time_sci[0]]
data1['xlon'] = nc1.variables['XLONG'][time_sci[0]]
lonindex_udom = np.where(np.logical_and(data1['xlon']>=-29.5, data1['xlon']<=-26.5))

###################################
# FILE #1
###################################
data1['theta'] = nc1.variables['T'][time_sci]+300 # potential temperature in K
data1['p'] = (nc1.variables['P'][time_sci[0]]+nc1.variables['PB'][time_sci[0]])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data1['p']/100000)**tempvar       
data1['Tk'] = tempvar0*data1['theta']
data1['rho'] = data1['p']/(constants.R*data1['Tk'])
data1['nisg80'] = nc1.variables['NISG80'][time_sci,:,:,:]*(data1['rho'])*(data1['rho'])
data1['nisg80'][data1['nisg80'] < 0.005] = np.nan
data1['udom_nisg80'] = data1['nisg80'][:,:,190:340,np.unique(lonindex_udom[1])]
data1['T'] = data1['Tk'][:,:,190:340,np.unique(lonindex_udom[1])]
del nc1

##########
###     #1: CNTRL
########## 	MIN 2DS = 0.008992138822380407
udom1_ni = data1['udom_nisg80'][:,0:22,:,:]/float(1e3)
udom1_ni[udom1_ni < 0.005] = np.nan
udom1_t = data1['T'][:,0:22,:,:]
Nisg_1 = data1['udom_nisg80'][:,0:22,:,:]/float(1e3)
Nisg_1 = Nisg_1[~np.isnan(Nisg_1)]
NisgFreq_1 = np.histogram(Nisg_1,bins)
del data1

runlab1 = 'CNTRL'

ni1_nanmean_tot = np.nanmean(udom1_ni)

ni1 = {}
ind = {}
ni1_med = 0.
ni1_nanmean = 0.
ni1_perc = 0.
ni1_997 = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(udom1_t>=T3D[i]-0.25, udom1_t<T3D[i]+0.25))
    ni1[strgi] = udom1_ni[ind[strgi]];
    if i==0:
        ni1_nanmean = np.nanmean(ni1[strgi])
        ni1_997 = np.nanpercentile(ni1[strgi],99.7)     
    if i>0:
        ni1_med = np.nanmean(ni1[strgi])
        ni1_nanmean = np.append(ni1_nanmean,ni1_med)
        ni1_perc = np.nanpercentile(ni1[strgi],99.7)
        ni1_997 = np.append(ni1_997,ni1_perc)

###################################
# FILE #2
###################################
data2 = {}
data2['theta'] = nc2.variables['T'][time_sci]+300 # potential temperature in K
data2['p'] = (nc2.variables['P'][time_sci[0]]+nc2.variables['PB'][time_sci[0]])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data2['p']/100000)**tempvar       
data2['Tk'] = tempvar0*data2['theta']
data2['rho'] = data2['p']/(constants.R*data2['Tk'])
data2['nisg80'] = nc2.variables['NISG80'][time_sci,:,:,:]*(data2['rho'])*(data2['rho'])
data2['nisg80'][data2['nisg80'] < 0.005] = np.nan
data2['udom_nisg80'] = data2['nisg80'][:,:,190:340,np.unique(lonindex_udom[1])]
data2['T'] = data2['Tk'][:,:,190:340,np.unique(lonindex_udom[1])]
del nc2

##########
###     #2: NoThresh
##########
udom2_ni = data2['udom_nisg80'][:,0:22,:,:]/float(1e3)
udom2_ni[udom2_ni < 0.005] = np.nan
udom2_t = data2['T'][:,0:22,:,:]
Nisg_2 = data2['udom_nisg80'][:,0:22,:,:]/float(1e3)
Nisg_2 = Nisg_2[~np.isnan(Nisg_2)]
NisgFreq_2 = np.histogram(Nisg_2,bins)
del data2

runlab2 = 'NoThresh'

ni2_nanmean_tot = np.nanmean(udom2_ni)

ni2 = {}
ni2_med = 0.
ni2_nanmean = 0.
ni2_perc = 0.
ni2_997 = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(udom2_t>=T3D[i]-0.25, udom2_t<T3D[i]+0.25))
    ni2[strgi] = udom2_ni[ind[strgi]];
    if i==0:
        ni2_nanmean = np.nanmean(ni2[strgi])
        ni2_997 = np.nanpercentile(ni2[strgi],99.7)     
    if i>0:
        ni2_med = np.nanmean(ni2[strgi])
        ni2_nanmean = np.append(ni2_nanmean,ni2_med)
        ni2_perc = np.nanpercentile(ni2[strgi],99.7)
        ni2_997 = np.append(ni2_997,ni2_perc)

###################################
# FILE #3
###################################
data3 = {}
data3['theta'] = nc3.variables['T'][time_sci]+300 # potential temperature in K
data3['p'] = (nc3.variables['P'][time_sci[0]]+nc3.variables['PB'][time_sci[0]])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data3['p']/100000)**tempvar       
data3['Tk'] = tempvar0*data3['theta']
data3['rho'] = data3['p']/(constants.R*data3['Tk'])
data3['nisg80'] = nc3.variables['NISG80'][time_sci,:,:,:]*(data3['rho'])*(data3['rho'])
data3['nisg80'][data3['nisg80'] < 0.005] = np.nan
data3['udom_nisg80'] = data3['nisg80'][:,:,190:340,np.unique(lonindex_udom[1])]
data3['T'] = data3['Tk'][:,:,190:340,np.unique(lonindex_udom[1])]
del nc3

##########
###     #3: 2XHM
##########
udom3_ni = data3['udom_nisg80'][:,0:22,:,:]/float(1e3)
udom3_ni[udom3_ni < 0.005] = np.nan
udom3_t = data3['T'][:,0:22,:,:]
Nisg_3 = data3['udom_nisg80'][:,0:22,:,:]/float(1e3)
Nisg_3 = Nisg_3[~np.isnan(Nisg_3)]
NisgFreq_3 = np.histogram(Nisg_3,bins)
del data3

runlab3 = '2xHM'

ni3_nanmean_tot = np.nanmean(udom3_ni)

ni3 = {}
ni3_med = 0.
ni3_nanmean = 0.
ni3_perc = 0.
ni3_997 = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(udom3_t>=T3D[i]-0.25, udom3_t<T3D[i]+0.25))
    ni3[strgi] = udom3_ni[ind[strgi]];
    if i==0:
        ni3_nanmean = np.nanmean(ni3[strgi])
        ni3_997 = np.nanpercentile(ni3[strgi],99.7)     
    if i>0:
        ni3_med = np.nanmean(ni3[strgi])
        ni3_nanmean = np.append(ni3_nanmean,ni3_med)
        ni3_perc = np.nanpercentile(ni3[strgi],99.7)
        ni3_997 = np.append(ni3_997,ni3_perc)

###################################
# FILE #4
###################################
data4 = {}
data4['theta'] = nc4.variables['T'][time_sci]+300 # potential temperature in K
data4['p'] = (nc4.variables['P'][time_sci[0]]+nc4.variables['PB'][time_sci[0]])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data4['p']/100000)**tempvar       
data4['Tk'] = tempvar0*data4['theta']
data4['rho'] = data4['p']/(constants.R*data4['Tk'])
data4['nisg80'] = nc4.variables['NISG80'][time_sci,:,:,:]*(data4['rho'])*(data4['rho'])
data4['nisg80'][data4['nisg80'] < 0.005] = np.nan
data4['udom_nisg80'] = data4['nisg80'][:,:,190:340,np.unique(lonindex_udom[1])]
data4['T'] = data4['Tk'][:,:,190:340,np.unique(lonindex_udom[1])]
del nc4

##########
###     #4: 5xHM
##########
udom4_ni = data4['udom_nisg80'][:,0:22,:,:]/float(1e3)
udom4_ni[udom4_ni < 0.005] = np.nan
udom4_t = data4['T'][:,0:22,:,:]
Nisg_4 = data4['udom_nisg80'][:,0:22,:,:]/float(1e3)
Nisg_4 = Nisg_4[~np.isnan(Nisg_4)]
NisgFreq_4 = np.histogram(Nisg_4,bins)
del data4

runlab4 = '5xHM'

ni4_nanmean_tot = np.nanmean(udom4_ni)

ni4 = {}
ni4_med = 0.
ni4_nanmean = 0.
ni4_perc = 0.
ni4_997 = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(udom4_t>=T3D[i]-0.25, udom4_t<T3D[i]+0.25))
    ni4[strgi] = udom4_ni[ind[strgi]];
    if i==0:
        ni4_nanmean = np.nanmean(ni4[strgi])
        ni4_997 = np.nanpercentile(ni4[strgi],99.7)     
    if i>0:
        ni4_med = np.nanmean(ni4[strgi])
        ni4_nanmean = np.append(ni4_nanmean,ni4_med)
        ni4_perc = np.nanpercentile(ni4[strgi],99.7)
        ni4_997 = np.append(ni4_997,ni4_perc)

###################################
# FILE #5
###################################
data5 = {}
data5['theta'] = nc5.variables['T'][time_sci]+300 # potential temperature in K
data5['p'] = (nc5.variables['P'][time_sci[0]]+nc5.variables['PB'][time_sci[0]])   # pressure in Pa
tempvar = constants.R/float(1005)
tempvar0 = (data5['p']/100000)**tempvar       
data5['Tk'] = tempvar0*data5['theta']
data5['rho'] = data5['p']/(constants.R*data5['Tk'])
data5['nisg80'] = nc5.variables['NISG80'][time_sci,:,:,:]*(data5['rho'])*(data5['rho'])
data5['nisg80'][data5['nisg80'] < 0.005] = np.nan
data5['udom_nisg80'] = data5['nisg80'][:,:,190:340,np.unique(lonindex_udom[1])]
data5['T'] = data5['Tk'][:,:,190:340,np.unique(lonindex_udom[1])]
del nc5

##########
###     #5: 10xHM
##########
udom5_ni = data5['udom_nisg80'][:,0:22,:,:]/float(1e3)
udom5_ni[udom5_ni < 0.005] = np.nan	# min from M218/M219
udom5_t = data5['T'][:,0:22,:,:]
Nisg_5 = data5['udom_nisg80'][:,0:22,:,:]/float(1e3)
Nisg_5 = Nisg_5[~np.isnan(Nisg_5)]
NisgFreq_5 = np.histogram(Nisg_5,bins)
del data5

runlab5 = '10xHM'

ni5_nanmean_tot = np.nanmean(udom5_ni)

ni5 = {}
ni5_med = 0.
ni5_nanmean = 0.
ni5_perc = 0.
ni5_997 = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(udom5_t>=T3D[i]-0.25, udom5_t<T3D[i]+0.25))
    ni5[strgi] = udom5_ni[ind[strgi]];
    if i==0:
        ni5_nanmean = np.nanmean(ni5[strgi])
        ni5_997 = np.nanpercentile(ni5[strgi],99.7)     
    if i>0:
        ni5_med = np.nanmean(ni5[strgi])
        ni5_nanmean = np.append(ni5_nanmean,ni5_med)
        ni5_perc = np.nanpercentile(ni5[strgi],99.7)
        ni5_997 = np.append(ni5_997,ni5_perc)

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

###################################################
###################################################
##### 	27degW
###################################################
###################################################

## science period for flight M218 - all at 27degW
science27 = np.where(np.logical_and(data218['CORE']['Intp_time']>=15.3, data218['CORE']['Intp_time']<=16.7))
newalt27 = data218['CORE']['Intp_alt'][science27]
newnisg27 = data218['TWODS']['Intp_Nice'][science27]
newnisg27[newnisg27 < 0.005] = np.nan
newcip27 = data218['CIP25']['Intp_Nice'][science27]
newcip27[newcip27 == 0] = np.nan
newtemp27 = data218['CORE']['Intp_tempdi'][science27]

## temp ranges from 260 to ~274K
ind = {}
ni27 = {}
ni27_nanmean = 0.
ncip27 = {}
ncipi27_nanmean = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(newtemp27>=T3D[i]-0.25, newtemp27<T3D[i]+0.25)); 
    ni27[strgi] = newnisg27[ind[strgi]];
    ncip27[strgi] = newcip27[ind[strgi]];
    if i==0:
    	ni27_nanmean = np.nanmean(newnisg27[ind[strgi]])
    	ncip27_nanmean = np.nanmean(newcip27[ind[strgi]])
    if i>0:
    	ni27_med = np.nanmean(newnisg27[ind[strgi]])
    	ni27_nanmean = np.append(ni27_nanmean,ni27_med)
    	ncip27_med = np.nanmean(newcip27[ind[strgi]])
    	ncip27_nanmean = np.append(ncip27_nanmean,ncip27_med)

###################################################
###################################################
##### 	28degW
###################################################
###################################################
science28 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.45, data219['CORE']['Intp_time']<=21.2))
newalt28 = data219['CORE']['Intp_alt'][science28]
newnisg28 = data219['TWODS']['Intp_Nice'][science28]
newnisg28[newnisg28 < 0.005] = np.nan
newcip28 = data219['CIP25']['Intp_Nice'][science28]
newcip28[newcip28 == 0] = np.nan
newtemp28 = data219['CORE']['Intp_tempdi'][science28]

## temp ranges from 260 to ~274K
ind = {}
ni28 = {}
ni28_nanmean = 0.
ncip28 = {}
ncipi28_nanmean = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(newtemp28>=T3D[i]-0.25, newtemp28<T3D[i]+0.25)); 
    ni28[strgi] = newnisg28[ind[strgi]];
    ncip28[strgi] = newcip28[ind[strgi]];
    if i==0:
    	ni28_nanmean = np.nanmean(newnisg28[ind[strgi]])
     	ncip28_nanmean = np.nanmean(newcip28[ind[strgi]])   	
    if i>0:
    	ni28_med = np.nanmean(newnisg28[ind[strgi]])
    	ni28_nanmean = np.append(ni28_nanmean,ni28_med)
    	ncip28_med = np.nanmean(newcip28[ind[strgi]])
    	ncip28_nanmean = np.append(ncip28_nanmean,ncip28_med)    	

###################################################
###################################################
##### 	29degW
###################################################
###################################################

## science period for M219 @ 29degW
science29 = np.where(np.logical_and(data219['CORE']['Intp_time']>=21.3, data219['CORE']['Intp_time']<=22.5))
newalt29 = data219['CORE']['Intp_alt'][science29]
newnisg29 = data219['TWODS']['Intp_Nice'][science29]
newnisg29[newnisg29 < 0.005] = np.nan
newcip29 = data219['CIP25']['Intp_Nice'][science29]
newcip29[newcip29 == 0] = np.nan
newtemp29 = data219['CORE']['Intp_tempdi'][science29]

## temp ranges from 260 to ~274K
ind = {}
ni29 = {}
ni29_nanmean = 0.
ncip29 = {}
ncipi29_nanmean = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    ind[strgi] = np.where(np.logical_and(newtemp29>=T3D[i]-0.25, newtemp29<T3D[i]+0.25)); 
    ni29[strgi] = newnisg29[ind[strgi]];
    ncip29[strgi] = newcip29[ind[strgi]];    
    if i==0:
    	ni29_nanmean = np.nanmean(newnisg29[ind[strgi]])
     	ncip29_nanmean = np.nanmean(newcip29[ind[strgi]])      	
    if i>0:
    	ni29_med = np.nanmean(newnisg29[ind[strgi]])
    	ni29_nanmean = np.append(ni29_nanmean,ni29_med)
    	ncip29_med = np.nanmean(newcip29[ind[strgi]])
    	ncip29_nanmean = np.append(ncip29_nanmean,ncip29_med)     	


###################################################
###################################################
##### 	TOTAL OBSERVED
###################################################
###################################################

niobs_2728 = np.append(newnisg27,newnisg28)
niobs_all = np.append(niobs_2728,newnisg29)
niobs_nanmean_tot = np.nanmean(niobs_all)

niobs = {}
ni_med = 0.
niobs_nanmean = 0.
niobs_perc = 0.
niobs_997 = 0.
ncipobs = {}
ncip_med = 0.
ncipobs_nanmean = 0.
for i in range(0,len(T3D)):
    strgi = "%1.f" % (i+1) # string of index number
    nia = ni27[strgi]
    nib = np.append(nia,ni28[strgi])
    niobs[strgi] = np.append(nib,ni29[strgi])
    ncipa = ncip27[strgi]
    ncipb = np.append(ncipa,ncip28[strgi])
    ncipobs[strgi] = np.append(ncipb,ncip29[strgi])    
    if i==0:
    	niobs_nanmean = np.nanmean(niobs[strgi])
    	niobs_997 = np.nanpercentile(niobs[strgi],99.7)
    	ncipobs_nanmean = np.nanmean(ncipobs[strgi])
    if i>0:
    	niobs_med = np.nanmean(niobs[strgi])
    	niobs_nanmean = np.append(niobs_nanmean,niobs_med)
    	niobs_perc = np.nanpercentile(niobs[strgi],99.7)
    	niobs_997 = np.append(niobs_997,niobs_perc)
    	ncipobs_med = np.nanmean(ncipobs[strgi])
    	ncipobs_nanmean = np.append(ncipobs_nanmean,ncipobs_med)

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
fig = plt.figure(figsize=(9,9))

#########################################################################################################

ax  = fig.add_axes([0.12,0.66,0.28,0.28])	# left, bottom, width, height
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
              facecolor='none',edgecolor='k',linewidth=2)
plt.gca().add_patch(p2)

p3 =  Polygon([(xd3_1,yd3_1),(xd3_2,yd3_2),(xd3_3,yd3_3),(xd3_4,yd3_4)],\
              facecolor='none',linestyle='--',edgecolor='k',linewidth=2)
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

plt.scatter(xH,yH,marker='x',color='k',linewidth=4,s=80)


data219 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight219/M219_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
x_m219,y_m219 = m(data219['CORE']['lon'], data219['CORE']['lat'])
plt.scatter(x_m219,y_m219,marker='o',color='mediumorchid',linewidth=1,s=1.,label='M219')


data218 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight218/M218_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
x_m218,y_m218 = m(data218['CORE']['lon'], data218['CORE']['lat'])
plt.scatter(x_m218,y_m218,marker='o',color='darkorange',linewidth=1,s=1.,label='M218')

plt.legend(bbox_to_anchor=(0.3, 0.76, 1., .102), loc=3, ncol=1)
# plt.annotate('(a)',xy=(-70,-40),xytext=(-70,-40),fontsize=10)

#============================== COLOURBAR

# add colorbar.
cbaxes = fig.add_axes([0.42,0.68, 0.02, 0.24])  # This is the position for the colorbar
cb = plt.colorbar(csf, cax = cbaxes)
# cb.ax.xaxis.set_ticks_position('top')
# cb.ax.xaxis.set_label_position('top')
cb.ax.axes.set_ylabel('NSIDC sea ice fraction')

#########################################################################################################

###################################
##  Freq distribution
###################################
ax  = fig.add_axes([0.6,0.667,0.25,0.27])   # left, bottom, width, height


data218['TWODS']['Intp_Nice'] = data218['TWODS']['Intp_Nice'][~np.isnan(data218['TWODS']['Intp_Nice'])]

lon2829index = np.where(~np.isnan(data219['TWODS']['Intp_Nice']))
data219['TWODS']['Intp_Nice'] = data219['TWODS']['Intp_Nice'][lon2829index]
data219['CORE']['Intp_time'] = data219['CORE']['Intp_time'][lon2829index]
science2829 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.45, data219['CORE']['Intp_time']<=22.5))

data218['NiceFreq'] = np.histogram(data218['TWODS']['Intp_Nice'][science27],bins)
data219['NiceFreq'] = np.histogram(data219['TWODS']['Intp_Nice'][science2829],bins)

NiceFreq_tot = np.concatenate([data218['TWODS']['Intp_Nice'][science27],data219['TWODS']['Intp_Nice'][science2829]])
NiceFreq = np.histogram(NiceFreq_tot,bins)
NiceFreqNorm = NiceFreq[0]/float(np.nansum(NiceFreq[0]))
NiceFreqNorm[NiceFreqNorm==0] = np.nan

plt.step(NiceFreq[1][:-1],NiceFreqNorm,color='k',linewidth=3,label='27-Nov')
plt.step(NisgFreq_1[1][:-1],NisgFreq_1[0]/float(np.nansum(NisgFreq_1[0])),color='slategray',label=runlab1)
plt.step(NisgFreq_2[1][:-1],NisgFreq_2[0]/float(np.nansum(NisgFreq_2[0])),color='indigo',label=runlab2)
plt.step(NisgFreq_3[1][:-1],NisgFreq_3[0]/float(np.nansum(NisgFreq_3[0])),color='darkorchid',label=runlab3)
plt.step(NisgFreq_4[1][:-1],NisgFreq_4[0]/float(np.nansum(NisgFreq_4[0])),color='blue',label=runlab4)
plt.step(NisgFreq_5[1][:-1],NisgFreq_5[0]/float(np.nansum(NisgFreq_5[0])),color='dodgerblue',label=runlab5)
plt.xlabel('$N_{isg>80}$ [$L^{-1}$]')
plt.ylabel('Normalised Frequency')
ax = plt.gca()
# ax.patch.set_facecolor('lightgrey')
plt.legend(bbox_to_anchor=(1.01, 0.2, 1., .102), loc=3, ncol=1)
ax.set_yscale("log", nonposy='clip')
ax.set_xscale("log", nonposy='clip')
ax.set_ylim([0,1])
ax.set_xlim([bins[0],bins[-1]])
yl = ax.get_ylim()
plt.plot((np.nanpercentile(NiceFreq_tot,99.7),np.nanpercentile(NiceFreq_tot,99.7)),yl,'kd-.',linewidth=3)
plt.plot((np.nanpercentile(Nisg_1,99.7),np.nanpercentile(Nisg_1,99.7)),yl,'d-.',color='slategray')
plt.plot((np.nanpercentile(Nisg_2,99.7),np.nanpercentile(Nisg_2,99.7)),yl,'d-.',color='indigo')
plt.plot((np.nanpercentile(Nisg_3,99.7),np.nanpercentile(Nisg_3,99.7)),yl,'d-.',color='darkorchid')
plt.plot((np.nanpercentile(Nisg_4,99.7),np.nanpercentile(Nisg_4,99.7)),yl,'d-.',color='blue')
plt.plot((np.nanpercentile(Nisg_5,99.7),np.nanpercentile(Nisg_5,99.7)),yl,'d-.',color='dodgerblue')
# plt.grid('on')
plt.annotate('(b)',xy=(0.25,0.3),xytext=(0.26,0.31),fontsize=14)

#########################################################################################################


###################################
## 	CNTRL
###################################
ax  = fig.add_axes([0.1,0.36,0.23,0.23])   # left, bottom, width, height
ax.plot(newtemp27-273.15,newnisg27,'.',color='lightgrey',markersize=4,label='27-Nov-2015')
ax.plot(newtemp28-273.15,newnisg28,'.',color='lightgrey',markersize=4)
ax.plot(newtemp29-273.15,newnisg29,'.',color='lightgrey',markersize=4)
ax.plot(MAC['Temp']-273.15, MAC['Ni'],color='k',label='Mean (entire MAC)',linewidth=3)
ax.plot(MAC['Temp']-273.15, MAC['Ni_997'],'k--',label='99.7 (entire MAC)')
ax.plot(T3D-273.15, niobs_nanmean,color='dimgrey',label='Mean (27-Nov)',linewidth=3)
ax.plot(T3D-273.15, niobs_997,'--',color='dimgrey',label='99.7 (27-Nov)')
ax.plot(T3D-273.15, ni1_nanmean,'r',label='Mean (PWRF)',linewidth=3)
ax.plot(T3D-273.15, ni1_997,'r--',label='99.7 (PWRF)')
ax.set_xlim([-9,-3])
ax.set_ylim([5e-3,1e2])
ax.set_yscale("log", nonposy='clip')
plt.ylabel('$N_{isg>80}$ [$L^{-1}$]')
# plt.title(runlab1)
plt.annotate('(c)',xy=(-9.0,0.5e2),xytext=(-8.9,0.5e2),fontsize=14)
plt.annotate(runlab1,xy=(-5.3,6e-3),xytext=(-5.4,6.1e-3),fontsize=14)



###################################
##  NoThresh
###################################
ax  = fig.add_axes([0.4,0.36,0.23,0.23])   # llighteft, bottom, width, height
ax.plot(newtemp27-273.15,newnisg27,'.',color='lightgrey',markersize=4,label='27-Nov-2015')
ax.plot(newtemp28-273.15,newnisg28,'.',color='lightgrey',markersize=4)
ax.plot(newtemp29-273.15,newnisg29,'.',color='lightgrey',markersize=4)
ax.plot(MAC['Temp']-273.15, MAC['Ni'],color='k',label='Mean (entire MAC)',linewidth=3)
ax.plot(MAC['Temp']-273.15, MAC['Ni_997'],'k--',label='99.7 (entire MAC)')
ax.plot(T3D-273.15, niobs_nanmean,color='dimgrey',label='Mean (27-Nov)',linewidth=3)
ax.plot(T3D-273.15, niobs_997,'--',color='dimgrey',label='99.7 (27-Nov)')
ax.plot(T3D-273.15, ni2_nanmean,'r',label='Mean (PWRF)',linewidth=3)
ax.plot(T3D-273.15, ni2_997,'r--',label='99.7 (PWRF)')
ax.set_xlim([-9,-3])
ax.set_ylim([5e-3,1e2])
ax.set_yscale("log", nonposy='clip')
# plt.title(runlab2)
plt.annotate('(d)',xy=(-9.0,0.5e2),xytext=(-8.9,0.5e2),fontsize=14)
plt.annotate(runlab2,xy=(-6.2,6e-3),xytext=(-6.1,6.1e-3),fontsize=14)
plt.legend(bbox_to_anchor=(1.3, 0.11, 1., .102), loc=3, ncol=1)

###################################
## 	2xHM
###################################
ax  = fig.add_axes([0.1,0.1,0.23,0.23])    # left, bottom, width, height
ax.plot(newtemp27-273.15,newnisg27,'.',color='lightgrey',markersize=4,label='27-Nov-2015')
ax.plot(newtemp28-273.15,newnisg28,'.',color='lightgrey',markersize=4)
ax.plot(newtemp29-273.15,newnisg29,'.',color='lightgrey',markersize=4)
ax.plot(MAC['Temp']-273.15, MAC['Ni'],color='k',label='Mean (entire MAC)',linewidth=3)
ax.plot(MAC['Temp']-273.15, MAC['Ni_997'],'k--',label='99.7 (entire MAC)')
ax.plot(T3D-273.15, niobs_nanmean,color='dimgrey',label='Mean (27-Nov)',linewidth=3)
ax.plot(T3D-273.15, niobs_997,'--',color='dimgrey',label='99.7 (27-Nov)')
ax.plot(T3D-273.15, ni3_nanmean,'r',label='Mean (PWRF)',linewidth=3)
ax.plot(T3D-273.15, ni3_997,'r--',label='99.7 (PWRF)')
ax.set_xlim([-9,-3])
ax.set_ylim([5e-3,1e2])
plt.ylabel('$N_{isg>80}$ [$L^{-1}$]')
plt.xlabel('Temperature [$^{o}C$]')
ax.set_yscale("log", nonposy='clip')
# plt.title(runlab3)
plt.annotate('(e)',xy=(-9.0,0.5e2),xytext=(-8.9,0.5e2),fontsize=14)
plt.annotate(runlab3,xy=(-5.0,6e-3),xytext=(-5.1,6.1e-3),fontsize=14)


###################################
## 	5xHM
###################################

ax  = fig.add_axes([0.4,0.1,0.23,0.23])    # left, bottom, width, height
ax.plot(newtemp27-273.15,newnisg27,'.',color='lightgrey',markersize=4,label='27-Nov-2015')
ax.plot(newtemp28-273.15,newnisg28,'.',color='lightgrey',markersize=4)
ax.plot(newtemp29-273.15,newnisg29,'.',color='lightgrey',markersize=4)
ax.plot(MAC['Temp']-273.15, MAC['Ni'],color='k',label='Mean (entire MAC)',linewidth=3)
ax.plot(MAC['Temp']-273.15, MAC['Ni_997'],'k--',label='99.7 (entire MAC)')
ax.plot(T3D-273.15, niobs_nanmean,color='dimgrey',label='Mean (27-Nov)',linewidth=3)
ax.plot(T3D-273.15, niobs_997,'--',color='dimgrey',label='99.7 (27-Nov)')
ax.plot(T3D-273.15, ni4_nanmean,'r',label='Mean (PWRF)',linewidth=3)
ax.plot(T3D-273.15, ni4_997,'r--',label='99.7 (PWRF)')
ax.set_xlim([-9,-3])
ax.set_ylim([5e-3,1e2])
plt.xlabel('Temperature [$^{o}C$]')
ax.set_yscale("log", nonposy='clip')
# plt.title(runlab4)
plt.annotate('(f)',xy=(-9.0,0.5e2),xytext=(-8.9,0.5e2),fontsize=14)
plt.annotate(runlab4,xy=(-5.0,6e-3),xytext=(-5.1,6.1e-3),fontsize=14)

###################################
##  10xHM
###################################
ax  = fig.add_axes([0.7,0.1,0.23,0.23])   # left, bottom, width, height
ax.plot(newtemp27-273.15,newnisg27,'.',color='lightgrey',markersize=4,label='27-Nov-2015')
ax.plot(newtemp28-273.15,newnisg28,'.',color='lightgrey',markersize=4)
ax.plot(newtemp29-273.15,newnisg29,'.',color='lightgrey',markersize=4)
ax.plot(MAC['Temp']-273.15, MAC['Ni'],color='k',label='Mean (entire MAC)',linewidth=3)
ax.plot(MAC['Temp']-273.15, MAC['Ni_997'],'k--',label='99.7 (entire MAC)')
ax.plot(T3D-273.15, niobs_nanmean,color='dimgrey',label='Mean (27-Nov)',linewidth=3)
ax.plot(T3D-273.15, niobs_997,'--',color='dimgrey',label='99.7 (27-Nov)')
ax.plot(T3D-273.15, ni5_nanmean,'r',label='Mean (PWRF)',linewidth=3)
ax.plot(T3D-273.15, ni5_997,'r--',label='99.7 (PWRF)')
ax.set_xlim([-9,-3])
ax.set_ylim([5e-3,1e2])
ax.set_yscale("log", nonposy='clip')
# plt.ylabel('$N_{isg>80}$ [$L^{-1}$]')
plt.xlabel('Temperature [$^{o}C$]')
# plt.title(runlab5)
plt.annotate('(g)',xy=(-9.0,0.5e2),xytext=(-8.9,0.5e2),fontsize=14)
plt.annotate(runlab5,xy=(-5.0,6e-3),xytext=(-5.1,6.1e-3),fontsize=14)


##################################################
##################################################
####    Stats
##################################################
##################################################

# print "Obs = Mean:",niobs_nanmean_tot," 99.7:",np.nanpercentile(NiceFreq_tot,99.7)
# print runlab1," = Mean:",ni1_nanmean_tot," 99.7:",np.nanpercentile(Nisg_1,99.7)
# print runlab2," = Mean:",ni2_nanmean_tot," 99.7:",np.nanpercentile(Nisg_2,99.7)
# print runlab3," = Mean:",ni3_nanmean_tot," 99.7:",np.nanpercentile(Nisg_3,99.7)
# print runlab4," = Mean:",ni4_nanmean_tot," 99.7:",np.nanpercentile(Nisg_4,99.7)
# print runlab5," = Mean:",ni5_nanmean_tot," 99.7:",np.nanpercentile(Nisg_5,99.7)

plt.savefig('/data/scihub-users/giyoung/PYTHON/WRF/FIGS/Misc/01_Domain_FreqDist_NisgvsT_v2.svg')
plt.show()


