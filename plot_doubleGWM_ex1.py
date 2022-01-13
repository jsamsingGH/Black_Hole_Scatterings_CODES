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
from matplotlib import colors as ccolor
from matplotlib import cm

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
yr_sec     = 31536000.0    
#----------------------------------------------------------


#Setting general font:
font = {'family' : 'serif',
        'size'   : 25}
mpl.rc('font', **font)










#----------------------------------------------------------
#plot:
#----------------------------------------------------------
fig = plt.figure(figsize=(12,9))
ax  = plt.subplot(111)

vkick_arr   = np.array([10.**(0.0), 10.**(0.5), 10.**(1.0), 10.**(1.5), 10.**(2.0), 10.**(2.5), 10.**(3.0)])

P_set_cop       = np.array([0.228966319274, 0.228653218231, 0.219931117771, 0.161560137764, 0.0856554993962, 0.0208435836651, 0.00241535089681])
Perr_set_cop    = np.array([0.00319646094211, 0.00319583498285, 0.003128773, 0.0026985916, 0.00197466, 0.000988089, 0.00031627])

P_set_iso       = np.array([0.0222222222222, 0.0233333333333, 0.02, 0.0244444444444, 0.0233333333333, 0.00666666666667, 0.00111111111111])
Perr_set_iso    = np.array([0.005, 0.005, 0.005, 0.005, 0.003, 0.003, 0.001])


#N-body res:
ax.errorbar(vkick_arr, P_set_cop, yerr=Perr_set_cop,  marker='',  linestyle='', linewidth=2.0, color='black')
ax.plot(vkick_arr, P_set_cop,        marker='',  linestyle='-', linewidth=2.0, color='black', label=r'$\gamma = 0$')
ax.plot(vkick_arr, P_set_cop,        marker='o', markersize=10.0, linestyle='', color='black')

ax.errorbar(vkick_arr, P_set_iso, yerr=Perr_set_iso,  marker='',  linestyle='', linewidth=2.0, color='red')
ax.plot(vkick_arr, P_set_iso,        marker='',  linestyle='-', linewidth=2.0, color='red', label=r'$\gamma = iso$')
ax.plot(vkick_arr, P_set_iso,        marker='o', markersize=10.0, linestyle='', color='red')

#labels etc:
ax.set_xlabel(r'kick velocity $v_{\rm k}$ [kms$^{-1}$]')
ax.set_ylabel(r'$P(t_{12}<10^{8}$ yrs.)')
ax.set_title(r'Prompt 2G Mergers')

#ax.set_title(r'Double GW Merger Probability')
ax.legend(loc='upper right', numpoints = 1, ncol = 1, fontsize = 40.0, frameon = False)

#axis:
ax.set_xlim([10.**(-0.25), 10.**(3.25)])
ax.set_ylim([10**(-3.0), 1e0])
ax.set_xscale('log')
ax.set_yscale('log')


#save and show fig:
save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/doubleGWmergers/'
plt.savefig(save_data_folder + 'tlife_prob.eps', bbox_inches='tight')
plt.show()
exit()
#----------------------------------------------------------








#----------------------------------------------------------
#input, calc, define:
#----------------------------------------------------------
Mmass   = 20.
#Example 1:
#vkick_kmsec = 0.0
##GW freq threshold:
#fGW_limit   = 50.0  #in Hz
##ecc limit:
#ecc_limit   = 0.1
##time-limit for second GW merger:
#t_limit     = 10.0  #in years
Ex1_prob_arr        = np.array([0.180245301497,    0.121821656464,    0.0683186007147,      0.0375850220882,    0.0183301546938])
Ex1_proberr_arr     = np.array([0.0020500960351,   0.00198642802532,   0.00179228122555,    0.00162342569284,   0.00139361585484])
Ex1_ISO_prob_arr    = np.array([0.0158394931362,   0.00548081714001,   0.00271370420624,    0.0,   0.0])
Ex1_ISO_proberr_arr = np.array([0.00166962917644,  0.00165252854527,   0.00191887864637,    0.0,   0.0])
#Example 2:
#vkick_kmsec = 10.0
##GW freq threshold:
#fGW_limit   = 50.0  #in Hz
##ecc limit:
#ecc_limit   = 0.1
##time-limit for second GW merger:
#t_limit     = 10.0  #in years
Ex2_prob_arr    = np.array([0.179802266474,     0.119295177016,      0.0545420349821,    0.0149358390015,    0.00275482093664])
Ex2_proberr_arr = np.array([0.00204757496204,    0.00196572167534,   0.00160140928775,   0.00102338682556,   0.000540264835091])
Ex2_ISO_prob_arr    = np.array([0.0154875043999,   0.00647732934728,   0.00135685210312,    0.0,   0.0])
Ex2_ISO_proberr_arr = np.array([0.00165097351631,  0.00179648792998,   0.00135685210312,    0.0,   0.0])

#Example 3:
#vkick_kmsec = 100.0
##GW freq threshold:
#fGW_limit   = 50.0  #in Hz
##ecc limit:
#ecc_limit   = 0.1
##time-limit for second GW merger:
#t_limit     = 10.0  #in years
Ex3_prob_arr    = np.array([0.15590169286,     0.0575907751109,    0.0126481098364,    0.00210363929598,   0.000317863954228])
Ex3_proberr_arr = np.array([0.00190663536788,   0.00136579980192,   0.000771168867165,      0.000384070231755,   0.000183518839539])
Ex3_ISO_prob_arr    = np.array([0.0158394931362,   0.00647732934728,   0.00135685210312,    0.0,   0.0])
Ex3_ISO_proberr_arr = np.array([0.00166962917644,  0.00179648792998,   0.00135685210312,    0.0,   0.0])





SMAa0_arr       = np.array([10**(-2.0), 10**(-1.5), 10**(-1.0), 10**(-0.5), 10**(0.0)])





##calc v_c:
#M_msun      = 20.0   #equal mass case.
#fHB         = 1./0.1
#v_c_fHB_kms = np.sqrt((3./2.)*(G_new_SI*(M_msun*M_sun_SI)/(fHB*SMAa0_arr*AU_SI)))/(1000.)
#print v_c_fHB_kms

##test print:
#t_limit_yr  = 1.0
#M_msun      = 20.0   #equal mass case.
#a0_AU       = 0.01
#dbp = (2.*(85./(4.*np.sqrt(2.)))**(1./7.))*(G_new_SI**(3./7.)/(c_SI**(5./7.)))*((t_limit_yr*yr_sec)**(1./7.))*(1.**(-1./14.))*((a0_AU*AU_SI)**(-4./7.))*((M_msun*M_sun_SI)**(3./7.))
#print dbp, (np.sqrt(3.)/2.), (np.sqrt(3.)/2.) + dbp, (np.sqrt(3.)/2.) - dbp
#exit()
#----------------------------------------------------------



#----------------------------------------------------------
#plot:
#----------------------------------------------------------
fig = plt.figure(figsize=(12,9))
ax  = plt.subplot(111)

#N-body res:
ax.errorbar(SMAa0_arr, Ex1_prob_arr, yerr=Ex1_proberr_arr,  marker='',  linestyle='', linewidth=2.0, color='black')
ax.plot(SMAa0_arr, Ex1_prob_arr,        marker='',  linestyle='-', linewidth=2.0, color='black', label=r'$v_{\rm K} = 0, \gamma = 0$')
ax.plot(SMAa0_arr, Ex1_prob_arr,        marker='o', markersize=10.0, linestyle='', color='black')
ax.errorbar(SMAa0_arr, Ex1_ISO_prob_arr, yerr=Ex1_ISO_proberr_arr,  marker='',  linestyle='', linewidth=2.0, color='red')
ax.plot(SMAa0_arr, Ex1_ISO_prob_arr,        marker='',  linestyle='-', linewidth=2.0, color='red', label=r'$v_{\rm K} = 0, \gamma = iso$')
ax.plot(SMAa0_arr, Ex1_ISO_prob_arr,        marker='o', markersize=10.0, linestyle='', color='red')

#N-body res:
ax.errorbar(SMAa0_arr, Ex2_prob_arr, yerr=Ex2_proberr_arr,  marker='',  linestyle='', linewidth=2.0, color='black')
ax.plot(SMAa0_arr, Ex2_prob_arr, marker='',  linestyle='--', linewidth=2.0, color='black', label=r'$v_{\rm K} = 10^{1}$ kms$^{-1}$, $\gamma = 0$')
ax.plot(SMAa0_arr, Ex2_prob_arr, marker='o', markersize=10.0, linestyle='', color='black')
ax.errorbar(SMAa0_arr, Ex2_ISO_prob_arr, yerr=Ex2_ISO_proberr_arr,  marker='',  linestyle='', linewidth=2.0, color='red')
ax.plot(SMAa0_arr, Ex2_ISO_prob_arr,        marker='',  linestyle='--', linewidth=2.0, color='red', label=r'$v_{\rm K} = 10^{1}$ kms$^{-1}$, $\gamma = iso$')
ax.plot(SMAa0_arr, Ex2_ISO_prob_arr,        marker='o', markersize=10.0, linestyle='', color='red')

#N-body res:
ax.errorbar(SMAa0_arr, Ex3_prob_arr, yerr=Ex3_proberr_arr,  marker='',  linestyle='', linewidth=2.0, color='black')
ax.plot(SMAa0_arr, Ex3_prob_arr, marker='',  linestyle='-.', linewidth=2.0, color='black', label=r'$v_{\rm K} = 10^{2}$ kms$^{-1}$, $\gamma = 0$')
ax.plot(SMAa0_arr, Ex3_prob_arr, marker='o', markersize=10.0, linestyle='', color='black')
ax.errorbar(SMAa0_arr, Ex3_ISO_prob_arr, yerr=Ex3_ISO_proberr_arr,  marker='',  linestyle='', linewidth=2.0, color='red')
ax.plot(SMAa0_arr, Ex3_ISO_prob_arr,        marker='',  linestyle='-.', linewidth=2.0, color='red', label=r'$v_{\rm K} = 10^{2}$ kms$^{-1}$, $\gamma = iso$')
ax.plot(SMAa0_arr, Ex3_ISO_prob_arr,        marker='o', markersize=10.0, linestyle='', color='red')


#analytical calc:
a0_pas      = 10.**(np.arange(-2.5,1,0.01))
prob_pas    = 0.26*(a0_pas/(10.**(-2.)))**(-4./7.) 
ax.plot(a0_pas, prob_pas, linestyle=':', color='black', label=r'$\propto a_{0}^{-4/7}, v_{\rm K} = 0$')

a0_pas      = 10.**(np.arange(-2.5,1,0.01))
prob_pas    = 0.02*(a0_pas/(10.**(-2.)))**(-8./7.) 
ax.plot(a0_pas, prob_pas, linestyle=':', color='red', label=r'$\propto a_{0}^{-8/7}, v_{\rm K} = 0$')

#labels etc:
ax.set_xlabel(r'semi-major axis $a_{0}$ [AU]')
ax.set_ylabel(r'$P_{\rm 12}(t_{12}<10$ yrs$|e_{1}>0.1(50 Hz))$')
#ax.set_title(r'Double GW Merger Probability')
ax.legend(loc='upper right', numpoints = 1, ncol = 2, fontsize = 20.0, frameon = False)

#axis:
x_lim_SMA_AU    = np.array([10.**(-2.15), 10.**(0.25)])
ax.set_xlim(x_lim_SMA_AU)
ax.set_ylim([10**(-3.0), 1e0])
ax.set_xscale('log')
ax.set_yscale('log')

#add upper x-axis showing log t_life:
Mmass_SI        = Mmass*M_sun_SI
t_beta          = (64./5.)*(G_new_SI**3./c_SI**5.)*Mmass_SI*Mmass_SI*(Mmass_SI+Mmass_SI)
x_lim_t_insp    = (((x_lim_SMA_AU*AU_SI)**4./(4.*t_beta))/yr_sec)
ax3 = ax.twiny()
ax3.set_xlim(x_lim_t_insp)
ax3.set_xscale('log')
ax3.set_xlabel(r'initial binary life time ${t_{0}}$ [year]', labelpad = 15)


#save and show fig:
save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/doubleGWmergers/'
plt.savefig(save_data_folder + 'doubleGWM_prob_ex.eps', bbox_inches='tight')
plt.show()
#----------------------------------------------------------












