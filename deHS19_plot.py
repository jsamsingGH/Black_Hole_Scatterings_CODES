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
import re    
import csv





#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI       = 299792458.0        #m/s
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
yr_sec     = 31536000.0
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI 
#------------------------------------------------------------------------------------------


#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


rp_min  = 3.0
rp_max  = 20.0
rp_arr  = 10.**(np.arange(np.log10(rp_min), np.log10(rp_max), 0.01))


ecc_ini = 0.999
m1      = 20.0
m2      = 20.0
m3      = 20.0
m12     = m1+m2

#PLOT 1:
fig1 = plt.figure(figsize=(5, 8))
f1ax1 = fig1.add_subplot(211)
#angles:
w_rot                   = 0.0
i_rot                   = np.pi/2.
Omega_rot_arr_fracpi    = [-(1./4.), -(1./8.), -(1./16.), (1./16.), (1./8.), (1./4.)]
Omega_rot_arr           = [i*np.pi for i in Omega_rot_arr_fracpi]
#colors/names:
color_arr               = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
legend_arr              = ['-1/4', '-1/8', '-1/16', '1/16', '1/8', '1/4']
#plot:
nras            = len(Omega_rot_arr)
for ac in range (0,nras):
    #calc de (HS19):
    Omega_rot   = Omega_rot_arr[ac]
    epsSA       = (((m3**2.)/(m12*(m12 + m3)))*((1./rp_arr)**3.)*((1.+1.)**(-3.)))**(1./2.) 
    de_HS19     = epsSA*(15.*np.pi/4.)*ecc_ini*np.sqrt(1.-ecc_ini**2.)*np.sin(2.*(-Omega_rot))*(np.sin(i_rot)**2.) + (epsSA**2.)*(3.*np.pi*ecc_ini/512.)*(100.*(1.-ecc_ini**2.)*np.sin(2.*(-Omega_rot))*((5.*np.cos(i_rot)+3.*np.cos(3.*i_rot))*np.cos(2.*(-w_rot)) + 6.*np.sin(i_rot)*np.sin(2.*i_rot)) + 4.*np.cos(2.*i_rot)*(3.*np.pi*(81.*(ecc_ini**2.) - 56.) + 200.*(1.-ecc_ini**2.)*np.cos(2.*(-Omega_rot))*np.sin(2.*(-w_rot))) + 3.*np.pi*(200.*(ecc_ini**2.)*(np.sin(i_rot)**4.)*np.cos(4.*(-Omega_rot)) + 8.*(16.*(ecc_ini**2.) + 9.)*(np.sin(2.*i_rot)**2.)*np.cos(2.*(-Omega_rot)) + (39.*(ecc_ini**2.) + 36.)*np.cos(4.*i_rot) - 299.*(ecc_ini**2.) + 124.))
    #plot:
    pos = np.where(de_HS19[:] > 0.0)[0]
    f1ax1.plot(rp_arr[pos], abs(de_HS19[pos]), linestyle='--', linewidth=2.0, alpha = 0.75, color=color_arr[ac])
    pos = np.where(de_HS19[:] < 0.0)[0]
    f1ax1.plot(rp_arr[pos], abs(de_HS19[pos]), linestyle=':',  linewidth=2.0, alpha = 0.75, color=color_arr[ac])
    #dummy plots for legends:
    f1ax1.plot(-1,-1, linestyle='-', linewidth=2.0, alpha = 0.75, color=color_arr[ac], label=r'$\Omega/\pi = $'+legend_arr[ac])
#dummy plots for legends:
f1ax1.plot(-1,-1, linestyle='--', linewidth=2.0, alpha = 0.75, color='black', label=r'$\Delta{e} > 0$')
f1ax1.plot(-1,-1, linestyle=':',  linewidth=2.0, alpha = 0.75, color='black', label=r'$\Delta{e} < 0$')
#settings:
f1ax1.set_xlim(rp_min, rp_max)
f1ax1.set_ylim(1e-6, 1e-1)
f1ax1.set_xscale("log")
f1ax1.set_yscale("log")
f1ax1.set_xlabel(r'$r_p/a_0$')
f1ax1.set_ylabel(r'$|\Delta{e}|$')
plt.text(6.0,0.05, r'[$e_0 = 0.999$, $m_1 = m_2 = m_3$, $i = \pi/2$, $\omega = 0.0$]', fontsize=8)
f1ax1.legend(loc='lower left', numpoints = 1, fontsize = 9.0, ncol = 3, frameon = False)


#PLOT 2:
f1ax2 = fig1.add_subplot(212)
#angles:
w_rot                   = 0.0
Omega_rot               = -np.pi/4.
i_rot_arr_fracpi        = [(1./8.), (2./8.), (3./8.), (4./8.)]
i_rot_arr               = [i*np.pi for i in i_rot_arr_fracpi]
#colors/names:
color_arr               = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
legend_arr              = ['1/8', '2/8', '3/8', '4/8']
#plot:
nras            = len(i_rot_arr)
for ac in range (0,nras):
    #calc de (HS19):
    i_rot       = i_rot_arr[ac]
    epsSA       = (((m3**2.)/(m12*(m12 + m3)))*((1./rp_arr)**3.)*((1.+1.)**(-3.)))**(1./2.) 
    de_HS19     = epsSA*(15.*np.pi/4.)*ecc_ini*np.sqrt(1.-ecc_ini**2.)*np.sin(2.*(-Omega_rot))*(np.sin(i_rot)**2.) + (epsSA**2.)*(3.*np.pi*ecc_ini/512.)*(100.*(1.-ecc_ini**2.)*np.sin(2.*(-Omega_rot))*((5.*np.cos(i_rot)+3.*np.cos(3.*i_rot))*np.cos(2.*(-w_rot)) + 6.*np.sin(i_rot)*np.sin(2.*i_rot)) + 4.*np.cos(2.*i_rot)*(3.*np.pi*(81.*(ecc_ini**2.) - 56.) + 200.*(1.-ecc_ini**2.)*np.cos(2.*(-Omega_rot))*np.sin(2.*(-w_rot))) + 3.*np.pi*(200.*(ecc_ini**2.)*(np.sin(i_rot)**4.)*np.cos(4.*(-Omega_rot)) + 8.*(16.*(ecc_ini**2.) + 9.)*(np.sin(2.*i_rot)**2.)*np.cos(2.*(-Omega_rot)) + (39.*(ecc_ini**2.) + 36.)*np.cos(4.*i_rot) - 299.*(ecc_ini**2.) + 124.))
    #plot:
    pos = np.where(de_HS19[:] > 0.0)[0]
    f1ax2.plot(rp_arr[pos], abs(de_HS19[pos]), linestyle='--', linewidth=2.0, alpha = 0.75, color=color_arr[ac])
    pos = np.where(de_HS19[:] < 0.0)[0]
    f1ax2.plot(rp_arr[pos], abs(de_HS19[pos]), linestyle=':',  linewidth=2.0, alpha = 0.75, color=color_arr[ac])
    #dummy plots for legends:
    f1ax2.plot(-1,-1, linestyle='-', linewidth=2.0, alpha = 0.75, color=color_arr[ac], label=r'$i/\pi = $'+legend_arr[ac])
#dummy plots for legends:
f1ax2.plot(-1,-1, linestyle='--', linewidth=2.0, alpha = 0.75, color='black', label=r'$\Delta{e} > 0$')
f1ax2.plot(-1,-1, linestyle=':',  linewidth=2.0, alpha = 0.75, color='black', label=r'$\Delta{e} < 0$')
#settings:
f1ax2.set_xlim(rp_min, rp_max)
f1ax2.set_ylim(1e-6, 1e-1)
f1ax2.set_xscale("log")
f1ax2.set_yscale("log")
f1ax2.set_xlabel(r'$r_p/a_0$')
f1ax2.set_ylabel(r'$|\Delta{e}|$')
plt.text(6.0,0.05, r'[$e_0 = 0.999$, $m_1 = m_2 = m_3$, $\Omega = -\pi/4$, $\omega = 0.0$]', fontsize=8)
f1ax2.legend(loc='lower left', numpoints = 1, fontsize = 9.0, ncol = 2, frameon = False)


#save fig:
plt.savefig('deHS19_ill.eps', bbox_inches='tight')        


plt.show()
exit()





















