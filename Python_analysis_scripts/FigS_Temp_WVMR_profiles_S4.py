from netCDF4 import Dataset as NetCDFFile
import numpy as np
from datetime import datetime
import constants
from wrf_functions import wrf_load
from Tkinter import Tk
from tkFileDialog import askopenfilename
import matplotlib
# from mpl_toolkits.basemap import Basemap, cm
import matplotlib.pyplot as plt
import matplotlib.cm as mpl_cm

###################################
# Pick file
###################################


#Tk().withdraw() # we don't want a full GUI, so keep the root window from appearing
#print('Choose wrfout (netcdf) file:')
#filename = askopenfilename() # show an "Open" dialog box and return the path to the selected file
	
#print(filename)

### MORRISON CNTRL
filename1 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/31_DeMott_WATSAT_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'



### THOMPSON CNTRL
# filename = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/42_ThompsonMP28/wrfout_d02_2015-11-27_00:00:00'

# filename2 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/29_DeMott_WATSAT_HM_LowThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename3 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/30_DeMott_WATSAT_HM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename4 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/35_DeMott_WATSAT_2xHM_LowThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'
# filename5 = '/data/scihub-users/giyoung/PWRF_V3.6.1/RUNS/MAC_WRF/36_DeMott_WATSAT_2xHM_noThresh_eta70_MYNN/wrfout_d02_2015-11-27_00:00:00'



###################################
# LOAD NETCDF FILE
###################################

nc1 = NetCDFFile(filename1, 'r')

##################################
##  Extract domain number: 
##################################

domainno_start = filename1.find('/wrfout_') + 8
domainno_end = filename1.find('_2015',domainno_start)
domainno = filename1[domainno_start:domainno_end]


###################################
# PROCESS WRF DATA FOR USE
###################################

data = wrf_load(nc1,domainno)

###################################
# LOAD FLIGHT DATA
###################################

data218 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight218/M218_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
data219 = np.load('/data/scihub-users/giyoung/MAC/FlightData/flight219/M219_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()

# data218 = np.load('../FlightData/M218_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()
# data219 = np.load('../FlightData/M219_CORE_CAPS_CIP25_2DS_GRIMM.npy').item()

##################################
##	Extract run label: 
##################################

runlab_start = filename1.find('/RUNS/') + 6 
runlab_end = filename1.find('/wrfout',runlab_start)
data['runlab'] = filename1[runlab_start:runlab_end]

###################################################
###################################################
##### 	27degW
###################################################
###################################################

## science period for flight M218 - all at 27degW
science27 = np.where(np.logical_and(data218['CORE']['Intp_time']>=15.3, data218['CORE']['Intp_time']<=16.7))
newalt27 = data218['CORE']['Intp_alt'][science27]
newnisg27 = data218['TWODS']['Intp_Nice'][science27]
newtemp27 = data218['CORE']['Intp_tempdi'][science27]
newq01_dp27 = data218['CORE']['Intp_q01_dp'][science27]*float(1e3)
newq01_fp27 = data218['CORE']['Intp_q01_fp'][science27]*float(1e3)

# ind100 = np.where(np.logical_and(newalt27>=50, newalt27<150)); ind200 = np.where(np.logical_and(newalt27>=150, newalt27<250))
# ind300 = np.where(np.logical_and(newalt27>=250, newalt27<350)); ind400 = np.where(np.logical_and(newalt27>=350, newalt27<450))
# ind500 = np.where(np.logical_and(newalt27>=450, newalt27<550)); ind600 = np.where(np.logical_and(newalt27>=550, newalt27<650))
# ind700 = np.where(np.logical_and(newalt27>=650, newalt27<750)); ind800 = np.where(np.logical_and(newalt27>=750, newalt27<850))
# ind900 = np.where(np.logical_and(newalt27>=850, newalt27<950)); ind1000 = np.where(np.logical_and(newalt27>=950, newalt27<1050))


####  50m vertical resolution
setres = 50
for i in range(0,22):
  tempvar = (i+1)*setres
  if i==0: strg0 = "%2.f" % tempvar 
  elif i>20: strg0 = "%4.f" % tempvar
  elif i>0: strg0 = "%3.f" % tempvar
  strg1 = ''.join(['ind',strg0])
  exec(strg1+" = np.where(np.logical_and(newalt27>=tempvar-(setres/2), newalt27<tempvar+(setres/2)))")
  altbinned = np.arange(setres,1101,setres)
  ### mean temperature
  strg2 = ''.join(['t',strg0])
  exec(strg2+" = np.nanmean(newtemp27["+strg1+"])")
  if i==0: exec("tbinned27 = np.array(("+strg2+"))")
  elif i>0: exec("tbinned27 = np.append(tbinned27,"+strg2+")")
  ### std temperature
  strg3 = ''.join(['ts',strg0])
  exec(strg3+" = np.nanstd(newtemp27["+strg1+"])")
  if i==0: exec("tsbinned27 = np.array(("+strg3+"))")
  elif i>0: exec("tsbinned27 = np.append(tsbinned27,"+strg3+")")
  #### mean Qvap (dewpoint)
  strg4 = ''.join(['q',strg0,'_dp'])
  exec(strg4+" = np.nanmean(newq01_dp27["+strg1+"])")
  if i==0: exec("qdpbinned27 = np.array(("+strg4+"))")
  elif i>0: exec("qdpbinned27 = np.append(qdpbinned27,"+strg4+")")
  #### mean Qvap (frostpoint)
  strg5 = ''.join(['q',strg0,'_fp'])
  exec(strg5+" = np.nanmean(newq01_fp27["+strg1+"])")
  if i==0: exec("qfpbinned27 = np.array(("+strg5+"))")
  elif i>0: exec("qfpbinned27 = np.append(qfpbinned27,"+strg5+")")
  #### std Qvap (frostpoint)
  strg6 = ''.join(['qs',strg0,'_fp'])
  exec(strg6+" = np.nanstd(newq01_fp27["+strg1+"])")
  if i==0: exec("qsfpbinned27 = np.array(("+strg6+"))")
  elif i>0: exec("qsfpbinned27 = np.append(qsfpbinned27,"+strg6+")")  


# t100 = np.nanmean(newtemp27[ind100]); t200 = np.nanmean(newtemp27[ind200]); 
# t300 = np.nanmean(newtemp27[ind300]); t400 = np.nanmean(newtemp27[ind400]); 
# t500 = np.nanmean(newtemp27[ind500]); t600 = np.nanmean(newtemp27[ind600]); 
# t700 = np.nanmean(newtemp27[ind700]); t800 = np.nanmean(newtemp27[ind800]); 
# t900 = np.nanmean(newtemp27[ind900]); t1000 = np.nanmean(newtemp27[ind1000]); 
# tbinned27 = np.array((t100, t200, t300, t400, t500, t600, t700, t800, t900, t1000))

# ts100 = np.nanstd(newtemp27[ind100]); ts200 = np.nanstd(newtemp27[ind200]); 
# ts300 = np.nanstd(newtemp27[ind300]); ts400 = np.nanstd(newtemp27[ind400]); 
# ts500 = np.nanstd(newtemp27[ind500]); ts600 = np.nanstd(newtemp27[ind600]); 
# ts700 = np.nanstd(newtemp27[ind700]); ts800 = np.nanstd(newtemp27[ind800]); 
# ts900 = np.nanstd(newtemp27[ind900]); ts1000 = np.nanstd(newtemp27[ind1000]); 
# tsbinned27 = np.array((ts100, ts200, ts300, ts400, ts500, ts600, ts700, ts800, ts900, ts1000))

# tmax100 = np.nan; tmax200 = np.nanmax(newtemp27[ind200],axis=0); 
# tmax300 = np.nanmax(newtemp27[ind300],axis=0); tmax400 = np.nanmax(newtemp27[ind400],axis=0); 
# tmax500 = np.nanmax(newtemp27[ind500],axis=0); tmax600 = np.nanmax(newtemp27[ind600],axis=0); 
# tmax700 = np.nanmax(newtemp27[ind700],axis=0); tmax800 = np.nanmax(newtemp27[ind800],axis=0); 
# tmax900 = np.nanmax(newtemp27[ind900],axis=0); tmax1000 = np.nanmax(newtemp27[ind1000],axis=0); 
# tmaxbinned27 = np.array((tmax100,tmax200, tmax300, tmax400, tmax500, tmax600, tmax700, tmax800, tmax900, tmax1000))

# tmin100 = np.nan; tmin200 = np.nanmin(newtemp27[ind200]); 
# tmin300 = np.nanmin(newtemp27[ind300]); tmin400 = np.nanmin(newtemp27[ind400]); 
# tmin500 = np.nanmin(newtemp27[ind500]); tmin600 = np.nanmin(newtemp27[ind600]); 
# tmin700 = np.nanmin(newtemp27[ind700]); tmin800 = np.nanmin(newtemp27[ind800]); 
# tmin900 = np.nanmin(newtemp27[ind900]); tmin1000 = np.nanmin(newtemp27[ind1000]); 
# tminbinned27 = np.array((tmin100,tmin200, tmin300, tmin400, tmin500, tmin600, tmin700, tmin800, tmin900, tmin1000))

# q100_dp = np.nanmean(newq01_dp27[ind100]); q200_dp = np.nanmean(newq01_dp27[ind200]); 
# q300_dp = np.nanmean(newq01_dp27[ind300]); q400_dp = np.nanmean(newq01_dp27[ind400]); 
# q500_dp = np.nanmean(newq01_dp27[ind500]); q600_dp = np.nanmean(newq01_dp27[ind600]); 
# q700_dp = np.nanmean(newq01_dp27[ind700]); q800_dp = np.nanmean(newq01_dp27[ind800]); 
# q900_dp = np.nanmean(newq01_dp27[ind900]); q1000_dp = np.nanmean(newq01_dp27[ind1000]); 
# qdpbinned27 = np.array((q100_dp, q200_dp, q300_dp, q400_dp, q500_dp,
#                       q600_dp, q700_dp, q800_dp, q900_dp, q1000_dp))

# q100_fp = np.nanmean(newq01_fp27[ind100]); q200_fp = np.nanmean(newq01_fp27[ind200]); 
# q300_fp = np.nanmean(newq01_fp27[ind300]); q400_fp = np.nanmean(newq01_fp27[ind400]); 
# q500_fp = np.nanmean(newq01_fp27[ind500]); q600_fp = np.nanmean(newq01_fp27[ind600]); 
# q700_fp = np.nanmean(newq01_fp27[ind700]); q800_fp = np.nanmean(newq01_fp27[ind800]); 
# q900_fp = np.nanmean(newq01_fp27[ind900]); q1000_fp = np.nanmean(newq01_fp27[ind1000]); 
# qfpbinned27 = np.array((q100_fp, q200_fp, q300_fp, q400_fp, q500_fp,
#                       q600_fp, q700_fp, q800_fp, q900_fp, q1000_fp))

# qs100_fp = np.nanstd(newq01_fp27[ind100]); qs200_fp = np.nanstd(newq01_fp27[ind200]); 
# qs300_fp = np.nanstd(newq01_fp27[ind300]); qs400_fp = np.nanstd(newq01_fp27[ind400]); 
# qs500_fp = np.nanstd(newq01_fp27[ind500]); qs600_fp = np.nanstd(newq01_fp27[ind600]); 
# qs700_fp = np.nanstd(newq01_fp27[ind700]); qs800_fp = np.nanstd(newq01_fp27[ind800]); 
# qs900_fp = np.nanstd(newq01_fp27[ind900]); qs1000_fp = np.nanstd(newq01_fp27[ind1000]); 
# qsfpbinned27 = np.array((qs100_fp, qs200_fp, qs300_fp, qs400_fp, qs500_fp,
#                       qs600_fp, qs700_fp, qs800_fp, qs900_fp, qs1000_fp))

# altbinned = np.arange(100,1001,100)



###################################################
###################################################
##### 	28degW
###################################################
###################################################
#science219 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.0, data219['CORE']['Intp_time']<=22.5))
science28 = np.where(np.logical_and(data219['CORE']['Intp_time']>=20.45, data219['CORE']['Intp_time']<=22.5))
newalt28 = data219['CORE']['Intp_alt'][science28]
newnisg28 = data219['TWODS']['Intp_Nice'][science28]
newtemp28 = data219['CORE']['Intp_tempdi'][science28]
newq01_dp28 = data219['CORE']['Intp_q01_dp'][science28]*float(1e3)
newq01_fp28 = data219['CORE']['Intp_q01_fp'][science28]*float(1e3)

####  50m vertical resolution
setres = 50
for i in range(0,22):
  tempvar = (i+1)*setres
  if i==0: strg0 = "%2.f" % tempvar 
  elif i>20: strg0 = "%4.f" % tempvar
  elif i>0: strg0 = "%3.f" % tempvar
  strg1 = ''.join(['ind',strg0])
  exec(strg1+" = np.where(np.logical_and(newalt28>=tempvar-(setres/2), newalt28<tempvar+(setres/2)))")
  ### mean temperature
  strg2 = ''.join(['t',strg0])
  exec(strg2+" = np.nanmean(newtemp28["+strg1+"])")
  if i==0: exec("tbinned28 = np.array(("+strg2+"))")
  elif i>0: exec("tbinned28 = np.append(tbinned28,"+strg2+")")
  ### std temperature
  strg3 = ''.join(['ts',strg0])
  exec(strg3+" = np.nanstd(newtemp28["+strg1+"])")
  if i==0: exec("tsbinned28 = np.array(("+strg3+"))")
  elif i>0: exec("tsbinned28 = np.append(tsbinned28,"+strg3+")")
  #### mean Qvap (dewpoint)
  strg4 = ''.join(['q',strg0,'_dp'])
  exec(strg4+" = np.nanmean(newq01_dp28["+strg1+"])")
  if i==0: exec("qdpbinned28 = np.array(("+strg4+"))")
  elif i>0: exec("qdpbinned28 = np.append(qdpbinned28,"+strg4+")")
  #### mean Qvap (frostpoint)
  strg5 = ''.join(['q',strg0,'_fp'])
  exec(strg5+" = np.nanmean(newq01_fp28["+strg1+"])")
  if i==0: exec("qfpbinned28 = np.array(("+strg5+"))")
  elif i>0: exec("qfpbinned28 = np.append(qfpbinned28,"+strg5+")")
  #### std Qvap (frostpoint)
  strg6 = ''.join(['qs',strg0,'_fp'])
  exec(strg6+" = np.nanstd(newq01_fp28["+strg1+"])")
  if i==0: exec("qsfpbinned28 = np.array(("+strg6+"))")
  elif i>0: exec("qsfpbinned28 = np.append(qsfpbinned28,"+strg6+")")  


# ind100 = np.where(np.logical_and(newalt28>=50, newalt28<150)); ind200 = np.where(np.logical_and(newalt28>=150, newalt28<250))
# ind300 = np.where(np.logical_and(newalt28>=250, newalt28<350)); ind400 = np.where(np.logical_and(newalt28>=350, newalt28<450))
# ind500 = np.where(np.logical_and(newalt28>=450, newalt28<550)); ind600 = np.where(np.logical_and(newalt28>=550, newalt28<650))
# ind700 = np.where(np.logical_and(newalt28>=650, newalt28<750)); ind800 = np.where(np.logical_and(newalt28>=750, newalt28<850))
# ind900 = np.where(np.logical_and(newalt28>=850, newalt28<950)); ind1000 = np.where(np.logical_and(newalt28>=950, newalt28<1050))

# t100 = np.nanmean(newtemp28[ind100]); t200 = np.nanmean(newtemp28[ind200]); 
# t300 = np.nanmean(newtemp28[ind300]); t400 = np.nanmean(newtemp28[ind400]); 
# t500 = np.nanmean(newtemp28[ind500]); t600 = np.nanmean(newtemp28[ind600]); #t700 = np.nanmean(newtemp28[ind700]); t800 = np.nanmean(newtemp28[ind800]); 
# t900 = np.nanmean(newtemp28[ind900]); t1000 = np.nanmean(newtemp28[ind1000]);

# tbinned28 = np.array((t100, t200, t300, t400, t500, t600, t700, t800, t900, t1000))

# ts100 = np.nanstd(newtemp28[ind100]); ts200 = np.nanstd(newtemp28[ind200]); 
# ts300 = np.nanstd(newtemp28[ind300]); ts400 = np.nanstd(newtemp28[ind400]); 
# ts500 = np.nanstd(newtemp28[ind500]); ts600 = np.nanstd(newtemp28[ind600]); 
# ts700 = np.nanstd(newtemp28[ind700]); ts800 = np.nanstd(newtemp28[ind800]); 
# ts900 = np.nanstd(newtemp28[ind900]); ts1000 = np.nanstd(newtemp28[ind1000]); 
# tsbinned28 = np.array((ts100, ts200, ts300, ts400, ts500, ts600, ts700, ts800, ts900, ts1000))

# q100_dp = np.nanmean(newq01_dp28[ind100]); q200_dp = np.nanmean(newq01_dp28[ind200]); 
# q300_dp = np.nanmean(newq01_dp28[ind300]); q400_dp = np.nanmean(newq01_dp28[ind400]); 
# q500_dp = np.nanmean(newq01_dp28[ind500]); q600_dp = np.nanmean(newq01_dp28[ind600]); 
# q700_dp = np.nanmean(newq01_dp28[ind700]); q800_dp = np.nanmean(newq01_dp28[ind800]); 
# q900_dp = np.nanmean(newq01_dp28[ind900]); q1000_dp = np.nanmean(newq01_dp28[ind1000]); 
# qdpbinned28 = np.array((q100_dp, q200_dp, q300_dp, q400_dp, q500_dp,
#                       q600_dp, q700_dp, q800_dp, q900_dp, q1000_dp))

# q100_fp = np.nanmean(newq01_fp28[ind100]); q200_fp = np.nanmean(newq01_fp28[ind200]); 
# q300_fp = np.nanmean(newq01_fp28[ind300]); q400_fp = np.nanmean(newq01_fp28[ind400]); 
# q500_fp = np.nanmean(newq01_fp28[ind500]); q600_fp = np.nanmean(newq01_fp28[ind600]); 
# q700_fp = np.nanmean(newq01_fp28[ind700]); q800_fp = np.nanmean(newq01_fp28[ind800]); 
# q900_fp = np.nanmean(newq01_fp28[ind900]); q1000_fp = np.nanmean(newq01_fp28[ind1000]); 
# qfpbinned28 = np.array((q100_fp, q200_fp, q300_fp, q400_fp, q500_fp,
#                       q600_fp, q700_fp, q800_fp, q900_fp, q1000_fp))

# qs100_fp = np.nanstd(newq01_fp28[ind100]); qs200_fp = np.nanstd(newq01_fp28[ind200]); 
# qs300_fp = np.nanstd(newq01_fp28[ind300]); qs400_fp = np.nanstd(newq01_fp28[ind400]); 
# qs500_fp = np.nanstd(newq01_fp28[ind500]); qs600_fp = np.nanstd(newq01_fp28[ind600]); 
# qs700_fp = np.nanstd(newq01_fp28[ind700]); qs800_fp = np.nanstd(newq01_fp28[ind800]); 
# qs900_fp = np.nanstd(newq01_fp28[ind900]); qs1000_fp = np.nanstd(newq01_fp28[ind1000]); 
# qsfpbinned28 = np.array((qs100_fp, qs200_fp, qs300_fp, qs400_fp, qs500_fp,
#                       qs600_fp, qs700_fp, qs800_fp, qs900_fp, qs1000_fp))


###################################################
# ###################################################
# ##### 	29degW
# ###################################################
# ###################################################

# ## science period for M219 @ 29degW
# science29 = np.where(np.logical_and(data219['CORE']['Intp_time']>=21.3, data219['CORE']['Intp_time']<=22.5))
# newalt29 = data219['CORE']['Intp_alt'][science29]
# newnisg29 = data219['TWODS']['Intp_Nice'][science29]
# newtemp29 = data219['CORE']['Intp_tempdi'][science29]
# newq01_dp29 = data219['CORE']['Intp_q01_dp'][science29]*float(1e3)
# newq01_fp29 = data219['CORE']['Intp_q01_fp'][science29]*float(1e3)

# ####  50m vertical resolution
# setres = 50
# for i in range(0,22):
#   tempvar = (i+1)*setres
#   if i==0: strg0 = "%2.f" % tempvar 
#   elif i>20: strg0 = "%4.f" % tempvar
#   elif i>0: strg0 = "%3.f" % tempvar
#   strg1 = ''.join(['ind',strg0])
#   exec(strg1+" = np.where(np.logical_and(newalt29>=tempvar-(setres/2), newalt29<tempvar+(setres/2)))")
#   ### mean temperature
#   strg2 = ''.join(['t',strg0])
#   exec(strg2+" = np.nanmean(newtemp29["+strg1+"])")
#   if i==0: exec("tbinned29 = np.array(("+strg2+"))")
#   elif i>0: exec("tbinned29 = np.append(tbinned29,"+strg2+")")
#   ### std temperature
#   strg3 = ''.join(['ts',strg0])
#   exec(strg3+" = np.nanstd(newtemp29["+strg1+"])")
#   if i==0: exec("tsbinned29 = np.array(("+strg3+"))")
#   elif i>0: exec("tsbinned29 = np.append(tsbinned29,"+strg3+")")
#   #### mean Qvap (dewpoint)
#   strg4 = ''.join(['q',strg0,'_dp'])
#   exec(strg4+" = np.nanmean(newq01_dp29["+strg1+"])")
#   if i==0: exec("qdpbinned29 = np.array(("+strg4+"))")
#   elif i>0: exec("qdpbinned29 = np.append(qdpbinned29,"+strg4+")")
#   #### mean Qvap (frostpoint)
#   strg5 = ''.join(['q',strg0,'_fp'])
#   exec(strg5+" = np.nanmean(newq01_fp29["+strg1+"])")
#   if i==0: exec("qfpbinned29 = np.array(("+strg5+"))")
#   elif i>0: exec("qfpbinned29 = np.append(qfpbinned29,"+strg5+")")
#   #### std Qvap (frostpoint)
#   strg6 = ''.join(['qs',strg0,'_fp'])
#   exec(strg6+" = np.nanstd(newq01_fp29["+strg1+"])")
#   if i==0: exec("qsfpbinned29 = np.array(("+strg6+"))")
#   elif i>0: exec("qsfpbinned29 = np.append(qsfpbinned29,"+strg6+")")  


# ind100 = np.where(np.logical_and(newalt29>=50, newalt29<150)); ind200 = np.where(np.logical_and(newalt29>=150, newalt29<250))
# ind300 = np.where(np.logical_and(newalt29>=250, newalt29<350)); ind400 = np.where(np.logical_and(newalt29>=350, newalt29<450))
# ind500 = np.where(np.logical_and(newalt29>=450, newalt29<550)); ind600 = np.where(np.logical_and(newalt29>=550, newalt29<650))
# ind700 = np.where(np.logical_and(newalt29>=650, newalt29<750)); ind800 = np.where(np.logical_and(newalt29>=750, newalt29<850))
# ind900 = np.where(np.logical_and(newalt29>=850, newalt29<950)); ind1000 = np.where(np.logical_and(newalt29>=950, newalt29<1050))

# t100 = np.nanmean(newtemp29[ind100]); t200 = np.nanmean(newtemp29[ind200]); 
# t300 = np.nanmean(newtemp29[ind300]); t400 = np.nanmean(newtemp29[ind400]); 
# t500 = np.nanmean(newtemp29[ind500]); t600 = np.nanmean(newtemp29[ind600]); 
# t700 = np.nanmean(newtemp29[ind700]); t800 = np.nanmean(newtemp29[ind800]); 
# t900 = np.nanmean(newtemp29[ind900]); t1000 = np.nanmean(newtemp29[ind1000]); 
# tbinned29 = np.array((t100, t200, t300, t400, t500, t600, t700, t800, t900, t1000))

# ts100 = np.nanstd(newtemp29[ind100]); ts200 = np.nanstd(newtemp29[ind200]); 
# ts300 = np.nanstd(newtemp29[ind300]); ts400 = np.nanstd(newtemp29[ind400]); 
# ts500 = np.nanstd(newtemp29[ind500]); ts600 = np.nanstd(newtemp29[ind600]); 
# ts700 = np.nanstd(newtemp29[ind700]); ts800 = np.nanstd(newtemp29[ind800]); 
# ts900 = np.nanstd(newtemp29[ind900]); ts1000 = np.nanstd(newtemp29[ind1000]); 
# tsbinned29 = np.array((ts100, ts200, ts300, ts400, ts500, ts600, ts700, ts800, ts900, ts1000))

# q100_dp = np.nanmean(newq01_dp29[ind100]); q200_dp = np.nanmean(newq01_dp29[ind200]); 
# q300_dp = np.nanmean(newq01_dp29[ind300]); q400_dp = np.nanmean(newq01_dp29[ind400]); 
# q500_dp = np.nanmean(newq01_dp29[ind500]); q600_dp = np.nanmean(newq01_dp29[ind600]); 
# q700_dp = np.nanmean(newq01_dp29[ind700]); q800_dp = np.nanmean(newq01_dp29[ind800]); 
# q900_dp = np.nanmean(newq01_dp29[ind900]); q1000_dp = np.nanmean(newq01_dp29[ind1000]); 
# qdpbinned29 = np.array((q100_dp, q200_dp, q300_dp, q400_dp, q500_dp,
#                       q600_dp, q700_dp, q800_dp, q900_dp, q1000_dp))

# q100_fp = np.nanmean(newq01_fp29[ind100]); q200_fp = np.nanmean(newq01_fp29[ind200]); 
# q300_fp = np.nanmean(newq01_fp29[ind300]); q400_fp = np.nanmean(newq01_fp29[ind400]); 
# q500_fp = np.nanmean(newq01_fp29[ind500]); q600_fp = np.nanmean(newq01_fp29[ind600]); 
# q700_fp = np.nanmean(newq01_fp29[ind700]); q800_fp = np.nanmean(newq01_fp29[ind800]); 
# q900_fp = np.nanmean(newq01_fp29[ind900]); q1000_fp = np.nanmean(newq01_fp29[ind1000]); 
# qfpbinned29 = np.array((q100_fp, q200_fp, q300_fp, q400_fp, q500_fp,
#                       q600_fp, q700_fp, q800_fp, q900_fp, q1000_fp))

# qs100_fp = np.nanstd(newq01_fp29[ind100]); qs200_fp = np.nanstd(newq01_fp29[ind200]); 
# qs300_fp = np.nanstd(newq01_fp29[ind300]); qs400_fp = np.nanstd(newq01_fp29[ind400]); 
# qs500_fp = np.nanstd(newq01_fp29[ind500]); qs600_fp = np.nanstd(newq01_fp29[ind600]); 
# qs700_fp = np.nanstd(newq01_fp29[ind700]); qs800_fp = np.nanstd(newq01_fp29[ind800]); 
# qs900_fp = np.nanstd(newq01_fp29[ind900]); qs1000_fp = np.nanstd(newq01_fp29[ind1000]); 
# qsfpbinned29 = np.array((qs100_fp, qs200_fp, qs300_fp, qs400_fp, qs500_fp,
#                       qs600_fp, qs700_fp, qs800_fp, qs900_fp, qs1000_fp))


# define plot title
#title = 'TEMP_mean1530-1630'
# lonlabel = np.round(xlon[a2])
# strg2 = "%3.1f" % lonlabel[0]
# strg3 = '{:%H%M_%d-%B-%Y}'.format(wrf_time218[0])
# strg4 = '{:%H%M_%d-%B-%Y}'.format(wrf_time218[-1])
# # title = ''.join([strg1,'_Lon',strg2,'deg_',strg3,''])
# title = ''.join([strg1,'_XLON',strg2,'deg_',strg3,'-',strg4])

# ###################################################
# ###################################################
# ##### 	TEMPERATURE COMPARISON
# ###################################################
# ###################################################


# def range(x, axis=0): return np.max(x, axis=axis) - np.min(x, axis=axis)


SMALL_SIZE = 10
MED_SIZE = 12
LARGE_SIZE = 14

plt.rc('font',size=SMALL_SIZE)
plt.rc('axes',titlesize=MED_SIZE)
plt.rc('axes',labelsize=MED_SIZE)
plt.rc('xtick',labelsize=MED_SIZE)
plt.rc('ytick',labelsize=MED_SIZE)
plt.rc('legend',fontsize=MED_SIZE)
# plt.rc('figure',titlesize=LARGE_SIZE)

## create figure and axes instances
fig = plt.figure(figsize=(6,7))

###################################
## 	M218 - 27degW
###################################
ax  = fig.add_axes([0.13,0.55,0.4,0.35])	# left, bottom, width, height
# plt.errorbar(data['udom_t27']-273.16,data['udom_z27'],xerr=data['udom_ts27'],fmt='o-',color='k',ecolor='k',markersize=4,label='WRF 1km')
plt.errorbar(np.nanmean(data['udom_t'][0:3,:],0)-273.16,data['udom_z27'],xerr=np.nanstd(data['T'][0:3,0:16,:,:]),fmt='o-',color='k',ecolor='grey',markersize=4,label='WRF 1km')
# plt.errorbar(np.nanmean(data['udom_t'][0:3,:],0)-273.16,data['udom_z27'],xerr=np.nanstd(data['T27'][0:22,:,:]),fmt='o-',color='k',ecolor='grey',markersize=4,label='WRF 1km')
plt.plot(newtemp27-273.16,newalt27,'.',color='pink',markersize=2,label='De-iced T')
#plt.plot(tbinned27-273.16,altbinned,'rd-',markersize=5,label='Mean observed')
plt.errorbar(tbinned27-273.16,altbinned,xerr=tsbinned27,fmt='d-',color='r',ecolor='r',markersize=4,label='Mean observed')
plt.xlabel('Temperature [$^{\circ}C$]')
plt.ylabel('Z [m]')
plt.title('M218')
axes = plt.gca()
axes.set_ylim([0,1100])
axes.set_xlim([-13,2])
# plt.legend()
plt.grid('on')

###################################
## 	M219 - 28degW
###################################
ax  = fig.add_axes([0.58,0.55,0.4,0.35])	# left, bottom, width, height
# plt.errorbar(data['udom_t28']-273.16,data['udom_z28'],xerr=data['udom_ts28'],fmt='o-',color='k',ecolor='k',markersize=4,label='WRF 1km')
plt.errorbar(np.nanmean(data['udom_t'][3:9,:],0)-273.16,data['udom_z27'],xerr=np.nanstd(data['T'][3:9,0:16,:,:]),fmt='o-',color='k',ecolor='grey',markersize=4,label='WRF 1km')
plt.plot(newtemp28-273.16,newalt28,'.',color='pink',markersize=2,label='De-iced T')
#plt.plot(tbinned28-273.16,altbinned,'rd-',markersize=5,label='Mean observed')
plt.errorbar(tbinned28-273.16,altbinned,xerr=tsbinned28,fmt='d-',color='r',ecolor='r',markersize=4,label='Mean observed')
plt.xlabel('Temperature [$^{\circ}C$]')
# plt.ylabel('Z [m]')
plt.title('M219')
axes = plt.gca()
axes.set_ylim([0,1100])
axes.set_xlim([-13,2])
axes.axes.yaxis.set_ticklabels([])
#plt.legend()
plt.grid('on')

# ###################################
# ## 	M219 - 29degW
# ###################################
# ax  = fig.add_axes([0.66,0.55,0.25,0.35])	# left, bottom, width, height
# # plt.errorbar(data['udom_t29']-273.16,data['udom_z29'],xerr=data['udom_ts29'],fmt='o-',color='k',ecolor='k',markersize=4,label='WRF 1km')
# plt.errorbar(np.nanmean(data['udom_t'][6:9,:],0)-273.16,data['udom_z27'],xerr=np.nanstd(data['T'][6:9,0:16,:,:]),fmt='o-',color='k',ecolor='grey',markersize=4,label='u-Dom')
# plt.plot(newtemp29-273.16,newalt29,'.',color='pink',markersize=2,label='Obs')
# #plt.plot(tbinned29-273.16,altbinned,'rd-',markersize=5,label='Mean observed')
# plt.errorbar(tbinned29-273.16,altbinned,xerr=tsbinned29,fmt='d-',color='r',ecolor='r',markersize=4,label='Mean Obs')
# plt.xlabel('Temperature [$^{\circ}C$]')
# # plt.ylabel('Z [m]')
# plt.title('29$^{\circ}W$')
# axes = plt.gca()
# axes.set_ylim([0,1100])
# axes.set_xlim([-13,2])
# axes.axes.yaxis.set_ticklabels([])
# #plt.legend()
# plt.legend(bbox_to_anchor=(.55, 0.6, 1., .102), loc=3, ncol=1)
# plt.grid('on')

###################################
## 	M218 - 27degW
###################################
ax  = fig.add_axes([0.13,0.1,0.4,0.35])	# left, bottom, width, height
# plt.errorbar(data['udom_q0127']*float(1e3),data['udom_z27'],xerr=data['udom_q01s27']*float(1e3),fmt='o-',color='k',ecolor='k',markersize=4,label='WRF 1km')
plt.errorbar(np.nanmean(data['udom_q01'][0:3,:],0)*float(1e3),data['udom_z27'],xerr=np.nanstd(data['udom_qvap'][0:3,0:16,:,:])*float(1e3),
  fmt='o-',color='k',ecolor='grey',markersize=4,label='WRF 1km')
##plt.plot(newq01_dp27,newalt27,'.',color='palegreen',markersize=2,label='DPt Buck')
##plt.plot(qdpbinned27,altbinned,'gd-',markersize=5,label='Mean DPt Buck')
plt.plot(newq01_fp27,newalt27,'.',color='skyblue',markersize=2,label='Buck')
#plt.plot(qfpbinned27,altbinned,'bd-',markersize=5,label='Mean FPt Buck')
plt.errorbar(qfpbinned27,altbinned,xerr=qsfpbinned27,fmt='d-',color='b',ecolor='b',markersize=4,label='Mean Buck')
plt.xlabel('WVMR [g $kg^{-1}$]')
plt.ylabel('Z [m]')
#plt.title('27$^{\circ}W$')
axes = plt.gca()
axes.set_ylim([0,1100])
axes.set_xlim([1,3.5])
# plt.legend()
plt.grid('on')

###################################
## 	M219 - 28degW
###################################
ax  = fig.add_axes([0.58,0.1,0.4,0.35])	# left, bottom, width, height
# plt.errorbar(data['udom_q0128']*float(1e3),data['udom_z28'],xerr=data['udom_q01s28']*float(1e3),fmt='o-',color='k',ecolor='k',markersize=4,label='WRF 1km')
plt.errorbar(np.nanmean(data['udom_q01'][3:9,:],0)*float(1e3),data['udom_z27'],xerr=np.nanstd(data['udom_qvap'][3:9,0:16,:,:])*float(1e3),
  fmt='o-',color='k',ecolor='grey',markersize=4,label='WRF 1km')
##plt.plot(newq01_dp28,newalt28,'.',color='palegreen',markersize=2,label='DPt Buck')
##plt.plot(qdpbinned28,altbinned,'gd-',markersize=5,label='Mean DPt Buck')
plt.plot(newq01_fp28,newalt28,'.',color='skyblue',markersize=2,label='Buck')
#plt.plot(qfpbinned28,altbinned,'bd-',markersize=5,label='Mean FPt Buck')
plt.errorbar(qfpbinned28,altbinned,xerr=qsfpbinned28,fmt='d-',color='b',ecolor='b',markersize=4,label='Mean Buck')
plt.xlabel('WVMR [g $kg^{-1}$]')
# plt.ylabel('Z [m]')
#plt.title('28$^{\circ}W$')
axes = plt.gca()
axes.set_ylim([0,1100])
axes.set_xlim([1,3.5])
axes.axes.yaxis.set_ticklabels([])
#plt.legend()
plt.grid('on')

# ###################################
# ## 	M219 - 29degW
# ###################################
# ax  = fig.add_axes([0.66,0.1,0.25,0.35])	# left, bottom, width, height
# # plt.errorbar(data['udom_q0129']*float(1e3),data['udom_z29'],xerr=data['udom_q01s29']*float(1e3),fmt='o-',color='k',ecolor='k',markersize=4,label='WRF 1km')
# plt.errorbar(np.nanmean(data['udom_q01'][6:9,:],0)*float(1e3),data['udom_z27'],xerr=np.nanstd(data['udom_qvap'][6:9,0:16,:,:])*float(1e3),
#   fmt='o-',color='k',ecolor='grey',markersize=4,label='WRF 1km')
# ##plt.plot(newq01_dp29,newalt29,'.',color='palegreen',markersize=2,label='DPt Buck')
# ##plt.plot(qdpbinned29,altbinned,'gd-',markersize=5,label='Mean DPt Buck')
# plt.plot(newq01_fp29,newalt29,'.',color='skyblue',markersize=2,label='Obs')
# #plt.plot(qfpbinned29,altbinned,'bd-',markersize=5,label='Mean FPt Buck')
# plt.errorbar(qfpbinned29,altbinned,xerr=qsfpbinned29,fmt='d-',color='b',ecolor='b',markersize=4,label='Mean Obs')
# plt.xlabel('WVMR [g $kg^{-1}$]')
# # plt.ylabel('Z [m]')
# #plt.title('29$^{\circ}W$')
# axes = plt.gca()
# axes.set_ylim([0,1100])
# axes.set_xlim([1,3.5])
# axes.axes.yaxis.set_ticklabels([])
# #plt.legend()
# # plt.legend(bbox_to_anchor=(.7, 0.35, 1., .102), loc=3, ncol=1)
# plt.grid('on')


# show plot
# plt.savefig('../Figures/FigS3.svg')
plt.show()

