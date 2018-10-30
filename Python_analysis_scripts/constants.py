###################################
###################################
### CONSTANTS
###################################
###################################

## specific gas constant for dry air (287.05 J/kg/K)
R = float(287.05) 

## gravitational acceleration [m/s2]
g = float(9.81)

## constants for saturation vapour pressure calculation
a0 = 6.11239921
a1 = 0.443987641
a2 = 0.142986287e-1
a3 = 0.264847430e-3
a4 = 0.302950461e-5
a5 = 0.206739458e-7
a6 = 0.640689451e-10
a7 = -0.952447341e-13
a8 = -0.976195544e-15


## constants for DeMott et al., 2010 calculation
demott_a = 5.94e-5
demott_b = 3.33
demott_c = 0.0264
demott_d = 0.0033

## constants for Tobo et al., 2013 calculation
tobo_alpha = -0.074
tobo_beta = 3.8
tobo_gamma = 0.414
tobo_delta = -9.671

