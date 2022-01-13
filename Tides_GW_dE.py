from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
import sys
from itertools import combinations
import itertools
from scipy import linalg
from sympy import mpmath
from subprocess import call
from mpl_toolkits.mplot3d import Axes3D
import subprocess
from scipy.optimize import fsolve
from scipy import interpolate
import matplotlib as mpl


#----------------------------------------------------------
#Units and conversions:
#----------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI       = 299792458.0        #m/s
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
#----------------------------------------------------------


#----------------------------------------------------------
#Settings:
#----------------------------------------------------------
nres                = 100
rp_unitRstar_arr    = np.logspace(np.log10(1.0), np.log10(100.0), num=nres, base=10)

#tidal obj:
b1_gas_n    = 1.5
MTO = 1.4 #0.6                                                  #SET THIS!
RTO = 1.7246e-5 #0.013*((1.43/MTO)**(1./3.))*((1.-MTO/1.43)**(0.447))  #SET THIS!
#compact obj:
MCO = 1.4
RCO = 1.7246e-5

MTOsi   = MTO*M_sun_SI
RTOsi   = RTO*R_sun_SI

MCOsi   = MCO*M_sun_SI
RCOsi   = RCO*R_sun_SI

#----------------------------------------------------------


#----------------------------------------------------------
#Energy loss: Tides
#----------------------------------------------------------
eta_arr         = ((MTOsi/(MTOsi + MCOsi))**(1./2.))*(rp_unitRstar_arr)**(3./2.)
log10eta_arr    = np.log10(eta_arr) 
#For n=1.5:
if (b1_gas_n == 1.5):
    print 'PT: n=', b1_gas_n
    fA = -0.397
    fB = 1.678
    fC = 1.277
    fD = -12.42
    fE = 9.446
    fF = -5.550
#For n=3.0:
if (b1_gas_n == 3.):
    print 'PT: n=', b1_gas_n
    fA = -1.124
    fB = 0.877
    fC = -13.37
    fD = 21.55
    fE = -16.48
    fF = 4.124
log10T2         = fA + fB*(log10eta_arr**(1.)) + fC*(log10eta_arr**(2.)) + fD*(log10eta_arr**(3.)) + fE*(log10eta_arr**(4.)) + fF*(log10eta_arr**(5))
T2_PT           = 10.**(log10T2)
dE_arr_tides    = ((G_new_SI*MCOsi**(2.))/RTOsi)*((rp_unitRstar_arr)**(-6.))*T2_PT
#----------------------------------------------------------
#Energy loss: GWs
#----------------------------------------------------------
#assuming high ecc limit:
dE_arr_GWs      = (85.*np.pi/(12.*np.sqrt(2.)))*((G_new_SI**(7./2.))/(c_SI**(5.)))*((MTOsi**2)*(MCOsi**2)*((MTOsi+MCOsi)**(1./2.)))/((RTOsi*rp_unitRstar_arr)**(7./2.))
#----------------------------------------------------------
#scale dE:
#----------------------------------------------------------
EGm2R   = ((G_new_SI*MTOsi**(2.))/RTOsi)
dE_UGm2R_arr_tides  = dE_arr_tides/EGm2R
dE_UGm2R_arr_GWs    = dE_arr_GWs/EGm2R
#take ratio:
dEtides_over_dEGW   = dE_UGm2R_arr_tides/dE_UGm2R_arr_GWs
#adjust scale:
pos_ok  = np.where(eta_arr > 1.0)[0]
rp_unitRstar_arr    = rp_unitRstar_arr[pos_ok]
dE_UGm2R_arr_tides  = dE_UGm2R_arr_tides[pos_ok]
dE_UGm2R_arr_GWs    = dE_UGm2R_arr_GWs[pos_ok]
dEtides_over_dEGW   = dEtides_over_dEGW[pos_ok]
#----------------------------------------------------------


#----------------------------------------------------------
#PLOT:
#----------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

#ax1.plot(rp_unitRstar_arr, dE_UGm2R_arr_tides*(rp_unitRstar_arr**6.), linewidth=1.0, linestyle='-')
#ax1.plot(rp_unitRstar_arr, dE_UGm2R_arr_GWs, linewidth=1.0, linestyle='-')

#test:
x = 1./((np.sqrt((2./5.)*(rp_unitRstar_arr)**3.)))
ax1.plot(x, np.log10(2.*dEtides_over_dEGW), linewidth=1.0, linestyle=':')
ax1.set_xlabel(r'$\nu_{p}/\nu(f-mode)$')
ax1.set_ylabel(r'$log_{10}(dE_t/dE_GW)$')

print x
print np.log10(eta_arr[pos_ok])

#plt.xscale('log')
#plt.yscale('log')

plt.show()





























