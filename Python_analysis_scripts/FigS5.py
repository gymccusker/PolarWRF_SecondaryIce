from netCDF4 import Dataset as NetCDFFile
import numpy as np
from datetime import datetime
import constants
from wrf_functions import wrf_load
from wrf_functions import params
from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib
# from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm



###################################
# LOAD FLIGHT DATA
###################################

data218 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight218/M218_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
data219 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight219/M219_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()

###################################
# IMF - WHOLE MAC CAMPAIGN. O'SHEA ET AL., 2017
###################################

IMF = np.array([0.025,0.075,0.125,0.175,0.225,0.275,
	0.325,0.375,0.425,0.475,0.525,0.575,0.625,0.675,0.725,0.775,0.825,0.875,0.925,0.975])

IMFlower = IMF - 0.025
IMFupper = IMF + 0.025

IMF_STRATUS = np.array([0.89568,0.00871094,0.00481951,0.00309361,0.0031913,0.00267027,0.00205155,
	0.00185616,0.00216552,0.00162821,0.00183988,0.00149796,0.00172591,0.00159565,0.00174219,
	0.00183988,0.00211668,0.0024586,0.00402169,0.0552941])

###################################
# IMF - M218 and M219
###################################

data218['CAS']['Intp_LWC+2DSLI'] = data218['CAS']['Intp_LWC'] + data218['TWODS']['Intp_QLI']
data219['CAS']['Intp_LWC+2DSLI'] = data219['CAS']['Intp_LWC'] + data219['TWODS']['Intp_QLI']

science27 = np.where(np.logical_and(data218['CORE']['Intp_time']>=15.3, data218['CORE']['Intp_time']<=16.7))
data218['IMF_TS'] = data218['TWODS']['Intp_Qice'][science27] / (data218['CAS']['Intp_LWC+2DSLI'][science27] + data218['TWODS']['Intp_Qice'][science27])
science2829 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.45, data219['CORE']['Intp_time']<=22.5))
data219['IMF_TS'] = data219['TWODS']['Intp_Qice'][science2829] / (data219['CAS']['Intp_LWC+2DSLI'][science2829] + data219['TWODS']['Intp_Qice'][science2829])

data218['IMF_binned'] = np.histogram(data218['IMF_TS'],IMFupper)
data219['IMF_binned'] = np.histogram(data219['IMF_TS'],IMFupper)

IMF_TS = np.concatenate([data218['IMF_TS'],data219['IMF_TS']])
IMF_binned = np.histogram(IMF_TS,IMFupper)


###################################
# MAKE IMF PLOT
###################################

plt.plot(IMF_binned[1][:-1],IMF_binned[0]/float(np.nansum(IMF_binned[0])),'bo-',label='M218+M219')
plt.plot(IMF,IMF_STRATUS,'ko-',linewidth=3,label='O''Shea et al., 2017')
plt.xlabel('Ice mass fraction')
plt.ylabel('Frequency')
# plt.title('29$^{\circ}W$')
axes = plt.gca()
axes.set_xlim([0,1])
# axes.set_ylim([0.001,1])
plt.legend()
axes.set_yscale("log", nonposy='clip')
plt.grid('on')
plt.show()


