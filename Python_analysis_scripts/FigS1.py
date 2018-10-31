###################################
###################################
###################################
#####	SCRIPT TO PLOT PARAMETERISATIONS
###################################
###################################
###################################

####################################
#
# Load modules
#
####################################

import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm
import numpy as np
from mpl_toolkits.axes_grid.inset_locator import inset_axes

T3D_D10=np.arange(238,270,1)
T3D_T13=np.arange(238,265,1)
T3D_Y17=np.arange(252,266,1)
T3D_C86=np.arange(238,266,1)
T3D_M92=np.arange(238,270,1)
T3D_HM=np.arange(265,271,1)
T3D_FL=np.arange(254.6,271.3,0.6)
rn1fac = 0.01
rn2fac = -0.6
rn3fac = -25.0
demott_a = 5.94e-5
demott_b = 3.33
demott_c = 0.0264
demott_d = 0.0033

## constants for Tobo et al., 2013 calculation
tobo_alpha = -0.074
tobo_beta = 3.8
tobo_gamma = 0.414
tobo_delta = -9.671

###################################
# M218 - Average Grimm out of cloud 0.5<d<1.6, CAS>/3micron, 2DS<=0.01/L
###################################

## np.nanmean(np.sum(np.sum(data218['GRIMM']['Intp_conc_binned'][outofcloud218,6:14],0),1))/1000
M218_GRIMM = 0.56429092673704939	# scm-3
M218_GRIMM_STD = 0.38987378926270483 	# scm-3
M218_GRIMM_UPPER = M218_GRIMM + M218_GRIMM_STD
M218_GRIMM_LOWER = M218_GRIMM - M218_GRIMM_STD

###################################
# M219 - Average Grimm out of cloud 0.5<d<1.6, CAS>/3micron, 2DS<=0.01/L
###################################

## np.nanmean(np.sum(np.sum(data218['GRIMM']['Intp_conc_binned'][outofcloud218,6:14],0),1))/1000
M219_GRIMM = 0.4081844790333829	# scm-3
M219_GRIMM_STD = 0.21779681522472463 # scm-3
M219_GRIMM_UPPER = M219_GRIMM + M219_GRIMM_STD
M219_GRIMM_LOWER = M219_GRIMM - M219_GRIMM_STD


###################################
# Inter-flight average
###################################

GRIMM_MEAN = np.mean((M218_GRIMM,M219_GRIMM))
GRIMM_STD = GRIMM_MEAN*np.sqrt(np.square(M218_GRIMM_STD/M218_GRIMM)+np.square(M219_GRIMM_STD/M219_GRIMM))

GRIMM_UPPER = GRIMM_MEAN + GRIMM_STD
GRIMM_LOWER = GRIMM_MEAN - GRIMM_STD

D10_GRIMM = demott_a*(np.power(273.15-T3D_D10,demott_b))*(np.power(GRIMM_MEAN,(demott_c*(273.15-T3D_D10)+demott_d)))*1000 # Num of ice crystals M-3
D10_GRIMM_UPPER = demott_a*(np.power(273.15-T3D_D10,demott_b))*(np.power(GRIMM_UPPER,(demott_c*(273.15-T3D_D10)+demott_d)))*1000 # Num of ice crystals M-3
D10_GRIMM_LOWER = demott_a*(np.power(273.15-T3D_D10,demott_b))*(np.power(GRIMM_LOWER,(demott_c*(273.15-T3D_D10)+demott_d)))*1000 # Num of ice crystals M-3

###################################
# M92 - CONTACT
###################################

M92 = np.exp(-2.80+0.262*(273.15-T3D_M92))*1000.

### min/max temps for 218
min218 = 260.0015351657724
max218 = 274.70261130222667

### min/max temps for 218
min219 = 261.83667254355038
max219 = 273.81650316240268

T3D = np.arange(np.round(min218),np.round(max218),1)

## create figure and axes instances
fig = plt.figure(figsize=(5.5,4))

plt.plot(T3D_D10-273.16,D10_GRIMM,'k',label='D10')
ax = plt.gca()
ax.fill_between(T3D_D10-273.16,D10_GRIMM_LOWER,D10_GRIMM_UPPER,facecolor='black',alpha=0.4)
plt.plot(T3D_M92-273.15,M92,'m',label='M92')
plt.grid('on')
plt.xlabel('Temperature [$^{\circ}C$]')
plt.ylabel('$N_{ice}$ [$m^{-3}$]')
ax.set_yscale("log", nonposy='clip')
ax.set_ylim([10,1e5])
ax.set_xlim([-35,1])
yl = ax.get_ylim()
plt.plot((T3D_M92[0]-273.16,T3D_M92[-1]-273.16),(30.0,30.0),color='darkorange',label='B53 range')
plt.plot(T3D_M92[0]-273.0,30.0,'<',color='darkorange')
plt.plot(T3D_M92[-1]-273.16,30.0,'>',color='darkorange')

plt.plot((T3D_HM[0]-273.16,T3D_HM[-1]-273.16),(70.0,70.0),color='blue',label='H-M range')
plt.plot(T3D_HM[0]-273.0,70.0,'<',color='blue')
plt.plot(T3D_HM[-1]-273.16,70.0,'>',color='blue')

plt.plot((T3D[0]-273.16,T3D[-1]-273.16),(1000.0,1000.0),'--',color='red',label='27-Nov range')
plt.plot(T3D[0]-273.0,1000.0,'<',color='red')
plt.plot(T3D[-1]-273.16,1000.0,'>',color='red')

plt.legend()


plt.savefig('/data/scihub-users/giyoung/PYTHON/WRF/FIGS/ice_nucleation_landscape_v2.svg',dpi=600)
plt.show()
