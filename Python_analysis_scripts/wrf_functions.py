### -------------------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------------------

def wrf_load(nc,dno):

	import constants
	import numpy as np
	from datetime import datetime

	data = {}

	###################################
	###################################
	### INDICES
	###################################
	###################################

	## set time index
	data['time218'] = np.array((31,32,33))
	data['time219'] = np.array((40,41,42,43,44,45))

	data['time27'] = np.array((31,32,33))
	data['time28'] = np.array((40,41,42))
	data['time29'] = np.array((43,44,45))

	data['time_sci'] = np.array((31,32,33,40,41,42,43,44,45))

	###################################
	###################################
	### STANDARD DOMAIN VARIABLES 
	###################################
	###################################

	# x_dim and y_dim are the x and y dimensions 
	# of the model domain in gridpoints
	data['x_dim'] = len(nc.dimensions['west_east'])
	data['y_dim'] = len(nc.dimensions['south_north'])

	# Get the grid spacing
	data['dx'] = float(nc.DX)
	data['dy'] = float(nc.DY)

	data['width_meters']  = data['dx'] * (data['x_dim'] - 1)
	data['height_meters'] = data['dy'] * (data['y_dim'] - 1)

	data['cen_lat']  = float(nc.CEN_LAT)
	data['cen_lon']  = float(nc.CEN_LON)
	data['truelat1'] = float(nc.TRUELAT1)
	data['truelat2'] = float(nc.TRUELAT2)
	data['standlon'] = float(nc.STAND_LON)
	# data['hgt'] = nc.variables['HGT_M']

	data['xlat'] = nc.variables['XLAT'][data['time218'][0]]
	data['xlon'] = nc.variables['XLONG'][data['time218'][0]]

	# data['lonindex'] = np.where(np.logical_and(xlon3>=-27.01, xlon3<=-26.99))

	# data['lonindex27'] = np.where(np.logical_and(data['xlon'][300,:]>=-27.075, data['xlon'][300,:]<=-26.925))
	# data['lonindex28'] = np.where(np.logical_and(data['xlon'][250,:]>=-28.20, data['xlon'][250,:]<=-28.05))
	# data['lonindex29'] = np.where(np.logical_and(data['xlon'][200,:]>=-29.075, data['xlon'][200,:]<=-28.925))

	# # data['lonindex27'] = np.where(np.logical_and(data['xlon']>=-27.075, data['xlon']<=-26.925))
	# # data['lonindex28'] = np.where(np.logical_and(data['xlon']>=-28.20, data['xlon']<=-28.05))
	# # data['lonindex29'] = np.where(np.logical_and(data['xlon']>=-29.075, data['xlon']<=-28.925))

	# # data['lonindex27'] = np.where(np.logical_and(data['xlon']>=-29.055, data['xlon']<=-26.95))
	# # data['lonindex28'] = np.where(np.logical_and(data['xlon']>=-29.055, data['xlon']<=-26.95))
	# # data['lonindex29'] = np.where(np.logical_and(data['xlon']>=-29.055, data['xlon']<=-26.95))

	# # data['lonindex27'] = np.where(np.logical_and(data['xlon']>=-29.5, data['xlon']<=-26.5))
	# # data['lonindex28'] = np.where(np.logical_and(data['xlon']>=-29.5, data['xlon']<=-26.5))
	# # data['lonindex29'] = np.where(np.logical_and(data['xlon']>=-29.5, data['xlon']<=-26.5))

	# data['lonindex_udom'] = np.where(np.logical_and(data['xlon']>=-29.5, data['xlon']<=-26.5))
	# # data['lonindex_udom2'] = np.where(np.logical_and(data['xlon'][250,:]>=-29.5, data['xlon']<=-26.5))

	###################################
	###################################

	data['times'] = nc.variables['Times']
	## get time index
	data['wrf_time218'] = datetime.strptime("".join(data['times'][data['time218'][0]]),'%Y-%m-%d_%H:%M:%S')
	data['wrf_time219'] = datetime.strptime("".join(data['times'][data['time219'][0]]),'%Y-%m-%d_%H:%M:%S')

	data['wrf_time27'] = datetime.strptime("".join(data['times'][data['time27'][0]]),'%Y-%m-%d_%H:%M:%S')
	data['wrf_time28'] = datetime.strptime("".join(data['times'][data['time28'][0]]),'%Y-%m-%d_%H:%M:%S')
	data['wrf_time29'] = datetime.strptime("".join(data['times'][data['time29'][0]]),'%Y-%m-%d_%H:%M:%S')

	###################################
	###################################
	### TD + MET VARIABLES 
	###################################
	###################################

	data['u10']  = nc.variables['U10']
	data['v10']  = nc.variables['V10']

	###################################
	# Z : 'Height' : 'm'
	###################################
	ph = nc.variables['PH'][data['time27']]
	phb = nc.variables['PHB'][data['time27']]
	tempvar1 = (ph + phb)/9.81
	data['Z27'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

	ph = nc.variables['PH'][data['time28']]
	phb = nc.variables['PHB'][data['time28']]
	tempvar1 = (ph+phb)/9.81
	data['Z28'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

	ph = nc.variables['PH'][data['time29']]
	phb = nc.variables['PHB'][data['time29']]
	tempvar1 = (ph+phb)/9.81
	data['Z29'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

	ph = nc.variables['PH'][data['time_sci']]
	phb = nc.variables['PHB'][data['time_sci']]
	tempvar1 = (ph+phb)/9.81
	data['Zsci'] = np.nanmean(0.5*(tempvar1[0:len(tempvar1)-2,:,:]+tempvar1[1:len(tempvar1)-1,:,:]),0)

	###################################
	# Tk : 'Temperature' : 'K'
	# need to convert to float BEFORE division
	###################################
	# data['theta218'] = nc.variables['T'][data['time218']]+300 # potential temperature in K
	# data['p218'] = (nc.variables['P'][data['time218'][0]]+nc.variables['PB'][data['time218'][0]]) # pressure in Pa
	# tempvar = constants.R/float(1005)
	# tempvar0 = (data['p218']/100000)**tempvar     
	# data['Tk218'] = tempvar0*np.nanmean(data['theta218'],0)

	# data['theta219'] = nc.variables['T'][data['time219']]+300 # potential temperature in K
	# data['p219'] = (nc.variables['P'][data['time219'][0]]+nc.variables['PB'][data['time219'][0]]) # pressure in Pa
	# tempvar = constants.R/float(1005)
	# tempvar0 = (data['p219']/100000)**tempvar     
	# data['Tk219'] = tempvar0*np.nanmean(data['theta219'],0)

	###################################
	data['theta27'] = nc.variables['T'][data['time27']]+300 # potential temperature in K
	data['p27'] = (nc.variables['P'][data['time27'][0]]+nc.variables['PB'][data['time27'][0]])   # pressure in Pa
	tempvar = constants.R/float(1005)
	tempvar0 = (data['p27']/100000)**tempvar       
	data['Tk27'] = tempvar0*np.nanmean(data['theta27'],0)

	data['theta28'] = nc.variables['T'][data['time28']]+300 # potential temperature in K
	data['p28'] = (nc.variables['P'][data['time28'][0]]+nc.variables['PB'][data['time28'][0]])   # pressure in Pa
	tempvar = constants.R/float(1005)
	tempvar0 = (data['p28']/100000)**tempvar       
	data['Tk28'] = tempvar0*np.nanmean(data['theta28'],0)

	data['theta29'] = nc.variables['T'][data['time29']]+300 # potential temperature in K
	data['p29'] = (nc.variables['P'][data['time29'][0]]+nc.variables['PB'][data['time29'][0]])   # pressure in Pa
	tempvar = constants.R/float(1005)
	tempvar0 = (data['p29']/100000)**tempvar       
	data['Tk29'] = tempvar0*np.nanmean(data['theta29'],0)

	data['theta'] = nc.variables['T'][data['time_sci']]+300 # potential temperature in K
	data['p'] = (nc.variables['P'][data['time_sci'][0]]+nc.variables['PB'][data['time_sci'][0]])   # pressure in Pa
	tempvar = constants.R/float(1005)
	tempvar0 = (data['p']/100000)**tempvar       
	data['Tk'] = tempvar0*data['theta']

	data['theta_all'] = nc.variables['T'][:]+300 # potential temperature in K
	data['p_all'] = (nc.variables['P'][:]+nc.variables['PB'][:])   # pressure in Pa
	tempvar = constants.R/float(1005)
	tempvar0 = (data['p_all']/100000)**tempvar       
	data['Tk_all'] = tempvar0*data['theta_all']

	tempvar1 = data['p_all']/9.81
	data['Z_all'] = np.nanmean(0.5*(tempvar1[:,0:np.size(tempvar1,1)-2,:,:]+tempvar1[:,1:np.size(tempvar1,1)-1,:,:]),0)


	###################################
	# rho : 'Air density' : 'kg/m3'
	###################################

	# data['rho218'] = data['p218']/(constants.R*data['Tk218'])
	# data['rho219'] = data['p219']/(constants.R*data['Tk219'])

	###################################

	data['rho27'] = data['p27']/(constants.R*data['Tk27'])
	data['rho28'] = data['p28']/(constants.R*data['Tk28'])
	data['rho29'] = data['p29']/(constants.R*data['Tk29'])
	data['rho'] = data['p']/(constants.R*data['Tk'])
	data['rho_all'] = data['p_all']/(constants.R*data['Tk_all'])

	###################################
	###################################
	### Q-FIELDS - Flight averaged
	###################################
	###################################

	# data['qcloud218'] = np.nanmean(nc.variables['QCLOUD'][data['time218'],:,:,:],0)	# Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
	# data['qcloud218'][data['qcloud218']<0]=0
	# data['qnisg218'] = np.nanmean((nc.variables['QNICE'][data['time218'],:,:,:]+
	#     nc.variables['QNSNOW'][data['time218'],:,:,:]+
	#     nc.variables['QNGRAUPEL'][data['time218'],:,:,:]),0)	# Qisg mean over M218 flight times @ lon=-27 [z,lat,lon]
	# data['qnisg218'][data['qnisg218']<0]=0
	# data['qvap218'] = np.nanmean(nc.variables['QVAPOR'][data['time218'],:,:,:],0)	# Qvapour mean over M218 flight times @ lon=-27 [z,lat,lon]
	# data['nisg80_218'] = np.nanmean(nc.variables['NISG80'][data['time218'],:,:,:],0)
	# data['qisg80_218'] = np.nanmean(nc.variables['QISG80'][data['time218'],:,:,:],0)
	# data['nig50_218'] = np.nanmean(nc.variables['NI50'][data['time218'],:,:,:],0) + np.nanmean(nc.variables['NG50'][data['time218'],:,:,:],0)
	# data['nsmic_218'] = data['qnisg218'] - data['nig50_218']

	# data['qcloud219'] = np.nanmean(nc.variables['QCLOUD'][data['time219'],:,:,:],0)	# Qcloud mean over M219 flight times @ lon=-27 [z,lat,lon]
	# data['qcloud219'][data['qcloud219']<0]=0
	# data['qnisg219'] = np.nanmean((nc.variables['QNICE'][data['time219'],:,:,:]+
	#     nc.variables['QNSNOW'][data['time219'],:,:,:]+
	#     nc.variables['QNGRAUPEL'][data['time219'],:,:,:]),0)	# Qisg mean over M219 flight times @ lon=-27 [z,lat,lon]
	# data['qnisg219'][data['qnisg219']<0]=0
	# data['qvap219'] = np.nanmean(nc.variables['QVAPOR'][data['time219'],:,:,:],0)	# Qvapour mean over M219 flight times @ lon=-27 [z,lat,lon]
	# data['nisg80_219'] = np.nanmean(nc.variables['NISG80'][data['time219'],:,:,:],0)
	# data['qisg80_219'] = np.nanmean(nc.variables['QISG80'][data['time219'],:,:,:],0)
	# data['nig50_219'] = np.nanmean(nc.variables['NI50'][data['time219'],:,:,:],0) + np.nanmean(nc.variables['NG50'][data['time219'],:,:,:],0)
	# data['nsmic_219'] = data['qnisg219'] - data['nig50_219']

	###################################
	###################################
	### Q-FIELDS - Lon averaged
	###################################
	###################################

	data['qcloud27'] = np.nanmean(nc.variables['QCLOUD'][data['time27'],:,:,:],0) # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
	data['qcloud27'][data['qcloud27']<0]=0
	data['qnisg27'] = np.nanmean((nc.variables['QNICE'][data['time27'],:,:,:]+
	nc.variables['QNSNOW'][data['time27'],:,:,:]+
	nc.variables['QNGRAUPEL'][data['time27'],:,:,:]),0)*data['rho']    # Qisg mean over M218 flight times @ lon=-27 [z,lat,lon]
	data['qnisg27'][data['qnisg27']<0]=0
	data['qvap27'] = np.nanmean(nc.variables['QVAPOR'][data['time27'],:,:,:],0)   # Qvapour mean over M218 flight times @ lon=-27 [z,lat,lon]
	data['nisg80_27'] = np.nanmean(nc.variables['NISG80'][data['time27'],:,:,:],0)*(data['rho'])*(data['rho'])
	data['qisg80_27'] = np.nanmean(nc.variables['QISG80'][data['time27'],:,:,:],0)
	data['nig50_27'] = np.nanmean(nc.variables['NI50'][data['time27'],:,:,:],0)*(data['rho'])*(data['rho']) + np.nanmean(nc.variables['NG50'][data['time27'],:,:,:],0)*(data['rho'])*(data['rho'])
	data['nsmic_27'] = data['qnisg27'] - data['nig50_27']

	data['qcloud28'] = np.nanmean(nc.variables['QCLOUD'][data['time28'],:,:,:],0) # Qcloud mean over M218 flight times @ lon=-28 [z,lat,lon]
	data['qcloud28'][data['qcloud28']<0]=0
	data['qnisg28'] = np.nanmean((nc.variables['QNICE'][data['time28'],:,:,:]+
	nc.variables['QNSNOW'][data['time28'],:,:,:]+
	nc.variables['QNGRAUPEL'][data['time28'],:,:,:]),0)*data['rho']    # Qisg mean over M218 flight times @ lon=-28 [z,lat,lon]
	data['qnisg28'][data['qnisg28']<0]=0
	data['qvap28'] = np.nanmean(nc.variables['QVAPOR'][data['time28'],:,:,:],0)   # Qvapour mean over M218 flight times @ lon=-28 [z,lat,lon]
	data['nisg80_28'] = np.nanmean(nc.variables['NISG80'][data['time28'],:,:,:],0)*(data['rho'])*(data['rho'])
	data['qisg80_28'] = np.nanmean(nc.variables['QISG80'][data['time28'],:,:,:],0)
	data['nig50_28'] = np.nanmean(nc.variables['NI50'][data['time28'],:,:,:],0)*(data['rho'])*(data['rho']) + np.nanmean(nc.variables['NG50'][data['time28'],:,:,:],0)*(data['rho'])*(data['rho'])
	data['nsmic_28'] = data['qnisg28'] - data['nig50_28']

	data['qcloud29'] = np.nanmean(nc.variables['QCLOUD'][data['time29'],:,:,:],0) # Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['qcloud29'][data['qcloud29']<0]=0
	data['qnisg29'] = np.nanmean((nc.variables['QNICE'][data['time29'],:,:,:]+
	nc.variables['QNSNOW'][data['time29'],:,:,:]+
	nc.variables['QNGRAUPEL'][data['time29'],:,:,:]),0)*data['rho']    # Qisg mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['qnisg29'][data['qnisg29']<0]=0
	data['qvap29'] = np.nanmean(nc.variables['QVAPOR'][data['time29'],:,:,:],0)   # Qvapour mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['nisg80_29'] = np.nanmean(nc.variables['NISG80'][data['time29'],:,:,:],0)*(data['rho'])*(data['rho'])
	data['qisg80_29'] = np.nanmean(nc.variables['QISG80'][data['time29'],:,:,:],0)
	data['nig50_29'] = np.nanmean(nc.variables['NI50'][data['time29'],:,:,:],0)*(data['rho'])*(data['rho']) + np.nanmean(nc.variables['NG50'][data['time29'],:,:,:],0)*(data['rho'])*(data['rho'])
	data['nsmic_29'] = data['qnisg29'] - data['nig50_29']

	data['qcloud'] = nc.variables['QCLOUD'][data['time_sci'],:,:,:]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['qcloud'][data['qcloud']<0]=0
	data['qliq'] = data['qcloud'] + nc.variables['QRAIN'][data['time_sci'],:,:,:]# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['qliq'][data['qliq']<0]=0
	data['qnisg'] = (nc.variables['QNICE'][data['time_sci'],:,:,:]+
	nc.variables['QNSNOW'][data['time_sci'],:,:,:]+
	nc.variables['QNGRAUPEL'][data['time_sci'],:,:,:])*data['rho']    # Qisg mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['qnisg'][data['qnisg']<0]=0
	data['qvap'] = nc.variables['QVAPOR'][data['time_sci'],:,:,:]   # Qvapour mean over M218 flight times @ lon=-29 [z,lat,lon]
	data['nisg80'] = nc.variables['NISG80'][data['time_sci'],:,:,:]*(data['rho'])*(data['rho'])
	data['qisg80'] = nc.variables['QISG80'][data['time_sci'],:,:,:]
	data['nig50'] = nc.variables['NI50'][data['time_sci'],:,:,:]*(data['rho'])*(data['rho']) + nc.variables['NG50'][data['time_sci'],:,:,:]*(data['rho'])*(data['rho'])
	data['nsmic'] = data['qnisg'] - data['nig50']
	data['qisg'] = (nc.variables['QICE'][data['time_sci'],:,:,:]+
	nc.variables['QSNOW'][data['time_sci'],:,:,:]+
	nc.variables['QGRAUP'][data['time_sci'],:,:,:])    # Qisg mean over M218 flight times @ lon=-29 [z,lat,lon]

	# data['qcloud27'][data['qcloud27']==0] = np.nan    # only include cloud data >=0.01 g/kg
	# data['qnisg27'][data['qnisg27']<=0.005*float(1000)] = np.nan
	# data['nisg80_27'][data['nisg80_27']<=0.005*float(1000)] = np.nan
	# data['qisg80_27'][data['qisg80_27']==0] = np.nan

	data['qnisg27'][data['qnisg27']<=0.0] = np.nan
	data['nisg80_27'][data['nisg80_27']<=0.0] = np.nan
	data['nsmic_27'][data['nsmic_27']==0] = np.nan

	# data['qcloud28'][data['qcloud28']==0] = np.nan
	# data['qnisg28'][data['qnisg28']<=0.005*float(1000)] = np.nan
	# data['nisg80_28'][data['nisg80_28']<=0.005*float(1000)] = np.nan
	# data['qisg80_28'][data['qisg80_28']==0] = np.nan

	data['qnisg28'][data['qnisg28']<=0.0] = np.nan
	data['nisg80_28'][data['nisg80_28']<=0.0] = np.nan
	data['nsmic_28'][data['nsmic_28']==0] = np.nan

	# data['qcloud29'][data['qcloud29']==0] = np.nan
	# data['qnisg29'][data['qnisg29']<=0.005*float(1000)] = np.nan
	# data['nisg80_29'][data['nisg80_29']<=0.005*float(1000)] = np.nan
	# data['qisg80_29'][data['qisg80_29']==0] = np.nan

	data['qnisg29'][data['qnisg29']<=0.0] = np.nan
	data['nisg80_29'][data['nisg80_29']<=0.0] = np.nan
	data['nsmic_29'][data['nsmic_29']==0] = np.nan

	# data['qcloud29'][data['qcloud29']==0] = np.nan
	# data['qnisg'][data['qnisg']<=0.005*float(1000)] = np.nan
	# data['nisg80'][data['nisg80']<=0.005*float(1000)] = np.nan

	data['qnisg'][data['qnisg']<=0] = np.nan
	data['nisg80'][data['nisg80']<=0] = np.nan
	data['qisg80'][data['qisg80']<=0] = np.nan


	###################################
	# rh : 'Relative Humidity' : '%'
	###################################

	# #es = 10*0.6112*np.exp(17.67*(Tk-273.15))/(Tk-29.65)
	# #qvs = (0.622*es)/(0.01*p-(1-0.622)*es)
	# #rh = (qvapor[time]/qvs)*100

	# #       dt = max(-80.,t-273.16)
	#    dT = data['Tk218'] - 273.16
	#    polysvp = constants.a0 + dT*(constants.a1+dT*(constants.a2+dT*(constants.a3+dT*
	# 	(constants.a4+dT*(constants.a5+dT*(constants.a6+dT*(constants.a7+constants.a8*dT)))))))
	#    data['evs218'] = polysvp*100. ## saturation vapour pressure over liquid (evs) in Pa
	#    data['qvs218'] = (0.622*data['evs218'])/(data['p218']-data['evs218']) ## saturation vapour mixing ratio over liquid (qvs)
	#    data['rh218'] = (data['qvap218']/data['qvs218'])*100

	# #       dt = max(-80.,t-273.16)
	#    dT = data['Tk219'] - 273.16
	#    polysvp = constants.a0 + dT*(constants.a1+dT*(constants.a2+dT*(constants.a3+dT*
	# 	(constants.a4+dT*(constants.a5+dT*(constants.a6+dT*(constants.a7+constants.a8*dT)))))))
	#    data['evs219'] = polysvp*100. ## saturation vapour pressure over liquid (evs) in Pa
	#    data['qvs219'] = (0.622*data['evs219'])/(data['p219']-data['evs219']) ## saturation vapour mixing ratio over liquid (qvs)
	#    data['rh219'] = (data['qvap219']/data['qvs219'])*100

	###################################
	dT = data['Tk27'] - 273.16
	polysvp = constants.a0 + dT*(constants.a1+dT*(constants.a2+dT*(constants.a3+dT*
	(constants.a4+dT*(constants.a5+dT*(constants.a6+dT*(constants.a7+constants.a8*dT)))))))
	data['evs27'] = polysvp*100. ## saturation vapour pressure over liquid (evs) in Pa
	data['qvs27'] = (0.622*data['evs27'])/(data['p27']-data['evs27']) ## saturation vapour mixing ratio over liquid (qvs)
	data['rh27'] = (data['qvap27']/data['qvs27'])*100

	dT = data['Tk28'] - 273.16
	polysvp = constants.a0 + dT*(constants.a1+dT*(constants.a2+dT*(constants.a3+dT*
	(constants.a4+dT*(constants.a5+dT*(constants.a6+dT*(constants.a7+constants.a8*dT)))))))
	data['evs28'] = polysvp*100. ## saturation vapour pressure over liquid (evs) in Pa
	data['qvs28'] = (0.622*data['evs28'])/(data['p28']-data['evs28']) ## saturation vapour mixing ratio over liquid (qvs)
	data['rh28'] = (data['qvap28']/data['qvs28'])*100

	dT = data['Tk29'] - 273.16
	polysvp = constants.a0 + dT*(constants.a1+dT*(constants.a2+dT*(constants.a3+dT*
	(constants.a4+dT*(constants.a5+dT*(constants.a6+dT*(constants.a7+constants.a8*dT)))))))
	data['evs29'] = polysvp*100. ## saturation vapour pressure over liquid (evs) in Pa
	data['qvs29'] = (0.622*data['evs29'])/(data['p29']-data['evs29']) ## saturation vapour mixing ratio over liquid (qvs)
	data['rh29'] = (data['qvap29']/data['qvs29'])*100

	dT = data['Tk'] - 273.16
	polysvp = constants.a0 + dT*(constants.a1+dT*(constants.a2+dT*(constants.a3+dT*
	(constants.a4+dT*(constants.a5+dT*(constants.a6+dT*(constants.a7+constants.a8*dT)))))))
	data['evs'] = polysvp*100. ## saturation vapour pressure over liquid (evs) in Pa
	data['qvs'] = (0.622*data['evs'])/(data['p']-data['evs']) ## saturation vapour mixing ratio over liquid (qvs)
	data['rh'] = (data['qvap']/data['qvs'])*100

	# ###################################
	# ###################################
	# ### HARD-CODED FOR D02 Make micro-(u)domain for flight comparison
	# ###################################
	# ###################################

	if dno == 'd02':

		# #!------ INDICES

		data['lonindex27'] = np.where(np.logical_and(data['xlon'][300,:]>=-27.075, data['xlon'][300,:]<=-26.925))
		data['lonindex28'] = np.where(np.logical_and(data['xlon'][250,:]>=-28.20, data['xlon'][250,:]<=-28.05))
		data['lonindex29'] = np.where(np.logical_and(data['xlon'][200,:]>=-29.075, data['xlon'][200,:]<=-28.925))
		
		data['lonindex_udom'] = np.where(np.logical_and(data['xlon']>=-29.5, data['xlon']<=-26.5))

		### min/max temps for 218
		#min218 = 260.0015351657724
		#max218 = 274.70261130222667
		min218 = 264
		max218 = 271
		T3D = np.arange(np.round(min218),np.round(max218),1)

		#######################################
		##!---- D02
		#######################################
		data['d02_z27'] = np.nanmean(np.nanmean(data['Z27'],1),1)   # mean alt over udomain height
		data['d02_z28'] = np.nanmean(np.nanmean(data['Z28'],1),1)   # mean alt over udomain height
		data['d02_z29'] = np.nanmean(np.nanmean(data['Z29'],1),1)   # mean alt over udomain height
		data['d02_t'] = np.nanmean(np.nanmean(data['Tk'],3),2)    # mean temperature over udomain height
		data['d02_q01'] = np.nanmean(np.nanmean(data['qvap'],3),2)    # mean temperature over udomain height    
		data['d02_ts'] = np.nanstd(np.nanstd(data['Tk'],3),2)    # stdev of temperature over udomain height
		data['d02_q01s'] = np.nanstd(np.nanmean(data['qvap'],3),2)    # stddev of q01 over udomain height

		##! ---- values for full udom
		data['d02_percind_t27'] = np.where(np.logical_and(data['d02_t'][0:3,:]>=T3D[0], data['d02_t'][0:3,:]<=T3D[-1]))
		data['d02_percind_t28'] = np.where(np.logical_and(data['d02_t'][3:6,:]>=T3D[0], data['d02_t'][3:6,:]<=T3D[-1]))
		data['d02_percind_t29'] = np.where(np.logical_and(data['d02_t'][6:9,:]>=T3D[0], data['d02_t'][6:9,:]<=T3D[-1]))

		# data['d02_nisg80'][data['d02_nisg80']==0] = np.nan    
		##! ---- Altitude(temperature)-binned
		data['d02_997_nisg'] = np.nanpercentile(np.nanpercentile(data['nisg80'],99.7,3),99.7,2)/float(1000) # nisg (total) outliers over d02ain
		data['d02_mean_nisg'] = np.nanmean(np.nanmean(data['nisg80'],3),2)/float(1000) # mean nisg (total) over d02ain
		data['d02_997_nsmic'] = np.nanpercentile(np.nanpercentile(data['nsmic'],99.7,3),99.7,2)/float(1000) # nisg (total) outliers over d02ain
		data['d02_mean_nsmic'] = np.nanmean(np.nanmean(data['nsmic'],3),2)/float(1000) # mean nisg (total) over d02ain

		# data['udom_qvap27'] = data['qvap218'][:,190:340,np.unique(data['lonindex27'][1])]	# Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qcloud27'] = data['qcloud218'][:,190:340,np.unique(data['lonindex27'][1])]	# Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qnisg27'] = data['qnisg218'][:,190:340,np.unique(data['lonindex27'][1])]
		# data['udom_nisg80_27'] = data['nisg80_218'][:,190:340,np.unique(data['lonindex27'][1])]
		# data['udom_nsmic27'] = data['nsmic_218'][:,190:340,np.unique(data['lonindex27'][1])]

		# data['udom_qvap28'] = data['qvap219'][:,190:340,np.unique(data['lonindex28'][1])]	# Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qcloud28'] = data['qcloud219'][:,190:340,np.unique(data['lonindex28'][1])]	# Qcloud mean over M218 flight times @ lon=-28 [z,lat,lon]
		# data['udom_qnisg28'] = data['qnisg219'][:,190:340,np.unique(data['lonindex28'][1])]
		# data['udom_nisg80_28'] = data['nisg80_219'][:,190:340,np.unique(data['lonindex28'][1])]
		# data['udom_nsmic28'] = data['nsmic_219'][:,190:340,np.unique(data['lonindex28'][1])]

		# data['udom_qvap29'] = data['qvap219'][:,190:340,np.unique(data['lonindex29'][1])]	# Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qcloud29'] = data['qcloud219'][:,190:340,np.unique(data['lonindex29'][1])]	# Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
		# data['udom_qnisg29'] = data['qnisg219'][:,190:340,np.unique(data['lonindex29'][1])]
		# data['udom_nisg80_29'] = data['nisg80_219'][:,190:340,np.unique(data['lonindex29'][1])]
		# data['udom_nsmic29'] = data['nsmic_219'][:,190:340,np.unique(data['lonindex29'][1])]   

		# data['T27'] = data['Tk218'][:,190:340,np.unique(data['lonindex27'][1])]	# -75 -> -73.5 degS, 27degW
		# data['T28'] = data['Tk219'][:,190:340,np.unique(data['lonindex28'][1])]	# -75 -> -73.5 degS, 27degW
		# data['T29'] = data['Tk219'][:,190:340,np.unique(data['lonindex29'][1])]	# -75 -> -73.5 degS, 27degW

		# data['zz27'] = data['Z218'][:-1,190:340,np.unique(data['lonindex27'][1])]	# [altitude, latitude, longitude]	# -75 -> -73.5 degS
		# data['zz28'] = data['Z219'][:-1,190:340,np.unique(data['lonindex28'][1])]	# [altitude, latitude, longitude]	# -75 -> -73.5 degS
		# data['zz29'] = data['Z219'][:-1,190:340,np.unique(data['lonindex29'][1])]	# [altitude, latitude, longitude]	# -75 -> -73.5 degS

		#!------ MICRO DOMAIN    

		# data['udom_qvap27'] = data['qvap27'][:,190:340,np.unique(data['lonindex27'][1])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qcloud27'] = data['qcloud27'][:,190:340,np.unique(data['lonindex27'][1])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qnisg27'] = data['qnisg27'][:,190:340,np.unique(data['lonindex27'][1])]
		# data['udom_nisg80_27'] = data['nisg80_27'][:,190:340,np.unique(data['lonindex27'][1])]
		# data['udom_nsmic27'] = data['nsmic_27'][:,190:340,np.unique(data['lonindex27'][1])]

		# data['udom_qvap28'] = data['qvap28'][:,190:340,np.unique(data['lonindex28'][1])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qcloud28'] = data['qcloud28'][:,190:340,np.unique(data['lonindex28'][1])]   # Qcloud mean over M218 flight times @ lon=-28 [z,lat,lon]
		# data['udom_qnisg28'] = data['qnisg28'][:,190:340,np.unique(data['lonindex28'][1])]
		# data['udom_nisg80_28'] = data['nisg80_28'][:,190:340,np.unique(data['lonindex28'][1])]
		# data['udom_nsmic28'] = data['nsmic_28'][:,190:340,np.unique(data['lonindex28'][1])]

		# data['udom_qvap29'] = data['qvap29'][:,190:340,np.unique(data['lonindex29'][1])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom_qcloud29'] = data['qcloud29'][:,190:340,np.unique(data['lonindex29'][1])]   # Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
		# data['udom_qnisg29'] = data['qnisg29'][:,190:340,np.unique(data['lonindex29'][1])]
		# data['udom_nisg80_29'] = data['nisg80_29'][:,190:340,np.unique(data['lonindex29'][1])]
		# data['udom_nsmic29'] = data['nsmic_29'][:,190:340,np.unique(data['lonindex29'][1])]   
		
		#######################################
		##!---- UDOM - BULK 
		#######################################
		data['udom_qvap27'] = data['qvap27'][:,190:340,np.squeeze(data['lonindex27'])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		data['udom_qcloud27'] = data['qcloud27'][:,190:340,np.squeeze(data['lonindex27'])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		data['udom_qnisg27'] = data['qnisg27'][:,190:340,np.squeeze(data['lonindex27'])]
		data['udom_nisg80_27'] = data['nisg80_27'][:,190:340,np.squeeze(data['lonindex27'])]
		data['udom_nsmic27'] = data['nsmic_27'][:,190:340,np.squeeze(data['lonindex27'])]

		data['udom_qvap28'] = data['qvap28'][:,190:340,np.squeeze(data['lonindex28'])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		data['udom_qcloud28'] = data['qcloud28'][:,190:340,np.squeeze(data['lonindex28'])]   # Qcloud mean over M218 flight times @ lon=-28 [z,lat,lon]
		data['udom_qnisg28'] = data['qnisg28'][:,190:340,np.squeeze(data['lonindex28'])]
		data['udom_nisg80_28'] = data['nisg80_28'][:,190:340,np.squeeze(data['lonindex28'])]
		data['udom_nsmic28'] = data['nsmic_28'][:,190:340,np.squeeze(data['lonindex28'])]

		data['udom_qvap29'] = data['qvap29'][:,190:340,np.squeeze(data['lonindex29'])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		data['udom_qcloud29'] = data['qcloud29'][:,190:340,np.squeeze(data['lonindex29'])]   # Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
		data['udom_qnisg29'] = data['qnisg29'][:,190:340,np.squeeze(data['lonindex29'])]
		data['udom_nisg80_29'] = data['nisg80_29'][:,190:340,np.squeeze(data['lonindex29'])]
		data['udom_nsmic29'] = data['nsmic_29'][:,190:340,np.squeeze(data['lonindex29'])]   

		data['udom_qvap'] = data['qvap'][:,:,190:340,np.unique(data['lonindex_udom'][1])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		data['udom_qcloud'] = data['qcloud'][:,:,190:340,np.unique(data['lonindex_udom'][1])]   # Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
		data['udom_qnisg'] = data['qnisg'][:,:,190:340,np.unique(data['lonindex_udom'][1])]
		data['udom_nisg80'] = data['nisg80'][:,:,190:340,np.unique(data['lonindex_udom'][1])]
		data['udom_qisg80'] = data['qisg80'][:,:,190:340,np.unique(data['lonindex_udom'][1])]
		data['udom_nsmic'] = data['nsmic'][:,:,190:340,np.unique(data['lonindex_udom'][1])]   
		data['udom_liqmass'] = data['qliq'][:,:,190:340,np.unique(data['lonindex_udom'][1])]
		data['udom_icemass'] = data['qisg80'][:,:,190:340,np.unique(data['lonindex_udom'][1])]

		# data['udom2_qvap'] = data['qvap'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]   # Qcloud mean over M218 flight times @ lon=-27 [z,lat,lon]
		# data['udom2_qcloud'] = data['qcloud'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]   # Qcloud mean over M218 flight times @ lon=-29 [z,lat,lon]
		# data['udom2_qnisg'] = data['qnisg'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]
		# data['udom2_nisg80'] = data['nisg80'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]
		# data['udom2_nsmic'] = data['nsmic'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]   
		# data['udom2_liqmass'] = data['qliq'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]
		# data['udom2_icemass'] = data['qisg80'][:,:,190:340,np.squeeze(data['lonindex_udom2'])]

		data['T27'] = data['Tk27'][:,190:340,np.squeeze(data['lonindex27'])] # -75 -> -73.5 degS, 27degW
		data['T28'] = data['Tk28'][:,190:340,np.squeeze(data['lonindex28'])] # -75 -> -73.5 degS, 27degW
		data['T29'] = data['Tk29'][:,190:340,np.squeeze(data['lonindex29'])] # -75 -> -73.5 degS, 27degW
		data['T'] = data['Tk'][:,:,190:340,np.unique(data['lonindex_udom'][1])] # -75 -> -73.5 degS, 27degW
		# data['T2'] = data['Tk'][:,:,190:340,np.squeeze(data['lonindex_udom2'])] # -75 -> -73.5 degS, 27degW

		data['zz27'] = data['Z27'][:-1,190:340,np.squeeze(data['lonindex27'])]   # [altitude, latitude, longitude]   # -75 -> -73.5 degS
		data['zz28'] = data['Z28'][:-1,190:340,np.squeeze(data['lonindex28'])]   # [altitude, latitude, longitude]   # -75 -> -73.5 degS
		data['zz29'] = data['Z29'][:-1,190:340,np.squeeze(data['lonindex29'])]   # [altitude, latitude, longitude]   # -75 -> -73.5 degS
		data['zz'] = data['Zsci'][:-1,190:340,np.unique(data['lonindex_udom'][1])]   # [altitude, latitude, longitude]   # -75 -> -73.5 degS   
		# data['zz2'] = data['Z219'][:-1,190:340,np.squeeze(data['lonindex_udom2'])]   # [altitude, latitude, longitude]   # -75 -> -73.5 degS   
		#         
		#######################################
		##!---- UDOM - STATISTICS
		#######################################
		data['udom_z'] = np.nanmean(np.nanmean(data['zz'],1),1)   # mean alt over udomain height
		data['udom_t'] = np.nanmean(np.nanmean(data['T'],3),2)    # mean temperature over udomain height
		data['udom_q01'] = np.nanmean(np.nanmean(data['udom_qvap'],3),2)    # mean temperature over udomain height    
		# data['udom_ts'] = np.nanstd(data['T'],3),2)    # stdev of temperature over udomain height
		# data['udom_q01s'] = np.nanstd(np.nanmean(data['udom_qvap'],3),2)    # stddev of q01 over udomain height

		##! ---- values for full udom
		data['udom_percind_t27'] = np.where(np.logical_and(data['udom_t'][0:3,:]>=T3D[0], data['udom_t'][0:3,:]<=T3D[-1]))
		data['udom_percind_t28'] = np.where(np.logical_and(data['udom_t'][3:6,:]>=T3D[0], data['udom_t'][3:6,:]<=T3D[-1]))
		data['udom_percind_t29'] = np.where(np.logical_and(data['udom_t'][6:9,:]>=T3D[0], data['udom_t'][6:9,:]<=T3D[-1]))

		# data['udom_nisg80'][data['udom_nisg80']==0] = np.nan    
		##! ---- Altitude(temperature)-binned
		data['udom_997_nisg'] = np.nanpercentile(np.nanpercentile(data['udom_nisg80'],99.7,3),99.7,2)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nisg'] = np.nanmean(np.nanmean(data['udom_nisg80'],3),2)/float(1000) # mean nisg (total) over udomain
		data['udom_997_nsmic'] = np.nanpercentile(np.nanpercentile(data['udom_nsmic'],99.7,3),99.7,2)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nsmic'] = np.nanmean(np.nanmean(data['udom_nsmic'],3),2)/float(1000) # mean nisg (total) over udomain

		# #######################################
		# ##!---- UDOM
		# #######################################
		# data['udom2_z'] = np.nanmean(np.nanmean(data['zz2'],1),1)   # mean alt over udomain height
		# data['udom2_t'] = np.nanmean(np.nanmean(data['T2'],3),2)    # mean temperature over udomain height
		# data['udom2_q01'] = np.nanmean(np.nanmean(data['udom2_qvap'],3),2)    # mean temperature over udomain height    
		# data['udom2_ts'] = np.nanstd(np.nanstd(data['T2'],3),2)    # stdev of temperature over udomain height
		# data['udom2_q01s'] = np.nanstd(np.nanmean(data['udom2_qvap'],3),2)    # stddev of q01 over udomain height

		# ##! ---- values for full udom
		# data['udom2_percind_t27'] = np.where(np.logical_and(data['udom2_t'][0:3,:]>=T3D[0], data['udom2_t'][0:3,:]<=T3D[-1]))
		# data['udom2_percind_t28'] = np.where(np.logical_and(data['udom2_t'][3:6,:]>=T3D[0], data['udom2_t'][3:6,:]<=T3D[-1]))
		# data['udom2_percind_t29'] = np.where(np.logical_and(data['udom2_t'][6:9,:]>=T3D[0], data['udom2_t'][6:9,:]<=T3D[-1]))

		# # data['udom_nisg80'][data['udom_nisg80']==0] = np.nan    
		# ##! ---- Altitude(temperature)-binned
		# data['udom2_997_nisg'] = np.nanpercentile(np.nanpercentile(data['udom2_nisg80'],99.7,3),99.7,2)/float(1000) # nisg (total) outliers over udomain
		# data['udom2_mean_nisg'] = np.nanmean(np.nanmean(data['udom2_nisg80'],3),2)/float(1000) # mean nisg (total) over udomain
		# data['udom2_997_nsmic'] = np.nanpercentile(np.nanpercentile(data['udom2_nsmic'],99.7,3),99.7,2)/float(1000) # nisg (total) outliers over udomain
		# data['udom2_mean_nsmic'] = np.nanmean(np.nanmean(data['udom2_nsmic'],3),2)/float(1000) # mean nisg (total) over udomain


		#######################################
		##!---- 27degW
		data['udom_z27'] = np.nanmean(np.nanmean(data['zz27'],1),1)   # mean alt over udomain height
		data['udom_t27'] = np.nanmean(np.nanmean(data['T27'],1),1)    # mean temperature over udomain height
		data['udom_q0127'] = np.nanmean(np.nanmean(data['udom_qvap27'],1),1)    # mean temperature over udomain height    
		# data['udom_ts27'] = np.nanstd(np.nanstd(data['T27'],1),1)    # stdev of temperature over udomain height
		# data['udom_q01s27'] = np.nanstd(np.nanmean(data['udom_qvap27'],1),1)    # stddev of q01 over udomain height    

		##! ---- values for full udom
		data['udom_percind27'] = np.where(np.logical_and(data['udom_t27']>=T3D[0], data['udom_t27']<=T3D[-1]))
		# data['iceindex27'] = np.where(data['udom_nisg80_27']>0.005)
		# data['udom_nisg80_27'][data['udom_nisg80_27']==0] = np.nan

		##! ---- Altitude(temperature)-binned
		data['udom_997_nisg27'] = np.nanpercentile(np.nanpercentile(data['udom_nisg80_27'],99.7,1),99.7,1)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nisg27'] = np.nanmean(np.nanmean(data['udom_nisg80_27'],1),1)/float(1000) # mean nisg (total) over udomain
		data['udom_997_nsmic27'] = np.nanpercentile(np.nanpercentile(data['udom_nsmic27'],99.7,1),99.7,1)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nsmic27'] = np.nanmean(np.nanmean(data['udom_nsmic27'],1),1)/float(1000) # mean nisg (total) over udomain

		#######################################
		##!---- 28degW
		data['udom_z28'] = np.nanmean(np.nanmean(data['zz28'],1),1)   # mean alt over udomain height
		data['udom_t28'] = np.nanmean(np.nanmean(data['T28'],1),1)    # mean temperature over udomain height
		data['udom_q0128'] = np.nanmean(np.nanmean(data['udom_qvap28'],1),1)    # mean temperature over udomain height    
		# data['udom_ts28'] = np.nanstd(np.nanstd(data['T28'],1),1)    # stdev of temperature over udomain height
		# data['udom_q01s28'] = np.nanstd(np.nanmean(data['udom_qvap28'],1),1)    # stddev of q01 over udomain height

		##! ---- values for full udom
		data['udom_percind28'] = np.where(np.logical_and(data['udom_t28']>=T3D[0], data['udom_t28']<=T3D[-1]))
		data['udom_nisg80_28'][data['udom_nisg80_28']==0] = np.nan
		##! ---- Altitude(temperature)-binned
		data['udom_997_nisg28'] = np.nanpercentile(np.nanpercentile(data['udom_nisg80_28'],99.7,1),99.7,1)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nisg28'] = np.nanmean(np.nanmean(data['udom_nisg80_28'],1),1)/float(1000) # mean nisg (total) over udomain
		data['udom_997_nsmic28'] = np.nanpercentile(np.nanpercentile(data['udom_nsmic28'],99.7,1),99.7,1)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nsmic28'] = np.nanmean(np.nanmean(data['udom_nsmic28'],1),1)/float(1000) # mean nisg (total) over udomain

		#######################################
		##!---- 29degW
		data['udom_z29'] = np.nanmean(np.nanmean(data['zz29'],1),1)   # mean alt over udomain height
		data['udom_t29'] = np.nanmean(np.nanmean(data['T29'],1),1)    # mean temperature over udomain height
		data['udom_q0129'] = np.nanmean(np.nanmean(data['udom_qvap29'],1),1)    # mean temperature over udomain height    
		# data['udom_ts29'] = np.nanstd(np.nanstd(data['T29'],1),1)    # stdev of temperature over udomain height
		# data['udom_q01s29'] = np.nanstd(np.nanmean(data['udom_qvap29'],1),1)    # stddev of q01 over udomain height

		##! ---- values for full udom
		data['udom_percind29'] = np.where(np.logical_and(data['udom_t29']>=T3D[0], data['udom_t29']<=T3D[-1]))
		data['udom_nisg80_29'][data['udom_nisg80_29']==0] = np.nan    
		##! ---- Altitude(temperature)-binned
		data['udom_997_nisg29'] = np.nanpercentile(np.nanpercentile(data['udom_nisg80_29'],99.7,1),99.7,1)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nisg29'] = np.nanmean(np.nanmean(data['udom_nisg80_29'],1),1)/float(1000) # mean nisg (total) over udomain
		data['udom_997_nsmic29'] = np.nanpercentile(np.nanpercentile(data['udom_nsmic29'],99.7,1),99.7,1)/float(1000) # nisg (total) outliers over udomain
		data['udom_mean_nsmic29'] = np.nanmean(np.nanmean(data['udom_nsmic29'],1),1)/float(1000) # mean nisg (total) over udomain


	return data

### -------------------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------------------
### -------------------------------------------------------------------------------------------------

def params(T3D):

    import constants
    import numpy as np
    
    prams = {}

    ###################################
    # M218 - Average Grimm out of cloud 0.5<d<1.6, CAS>/3micron, 2DS<=0.01/L
    ###################################
    
    ## np.nanmean(np.sum(np.sum(data218['GRIMM']['Intp_conc_binned'][outofcloud218,6:14],0),1))/1000
    M218_GRIMM = 0.56429092673704939	# scm-3

    ## np.nanmean(np.sum(np.sum(data218['GRIMM']['Intp_conc_binned'][outofcloud218,6:14],0),1))/1000
    M219_GRIMM = 0.4081844790333829	# scm-3

    ###################################
    # Inter-flight average
    ###################################
    
    prams['GRIMM_MEAN'] = np.mean((M218_GRIMM,M219_GRIMM))
    
    prams['D10_GRIMM'] = constants.demott_a*(np.power(273.15-
        T3D,constants.demott_b))*(np.power(prams['GRIMM_MEAN'],(constants.demott_c*(273.15-
                               T3D)+constants.demott_d)))*1000 # Num of ice crystals M-3
  
    ###################################
    # C86
    ###################################
    
    prams['C86'] = 0.005*np.exp(0.304*(273.15-T3D))*1000 	# convert from L-1 to m-3 
    prams['T'] = T3D
    return prams
