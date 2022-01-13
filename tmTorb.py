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
import matplotlib as mpl
from matplotlib import cm
from matplotlib.patches import Rectangle
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)



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
tU_tSI     = np.sqrt(R_sun_SI**3./(G_new_SI*M_sun_SI))    
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
m_parsec    = 3.086*(10**16.)   #m
#----------------------------------------------------------


mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'
mpl.rcParams['figure.autolayout'] = 'true'



#----------------------------------------------------------
#Input:
#----------------------------------------------------------
fGW_HZ      = 1             #in Hz (linear between 1 and 10 Hz: 10**(0, 0.5, 1.0) = 1, 3, 10)
v_disp_C1   = 10.0*1000.    #in m/s
v_disp_C2   = 100.0*1000.   #in m/s
f_disp_esc  = 5.
n_BH_m3     = (10.**5)/(m_parsec**3.)
#----------------------------------------------------------




nr_bins         = 1000

log10_R_AU_arr  = np.linspace(-2, 2, nr_bins)
M_msun_arr      = np.linspace(1.0, 100., nr_bins)

#define:
R_AU_arr        = 10.**(log10_R_AU_arr)
tm_Torb_arr     = np.zeros((nr_bins,nr_bins), dtype='d')
Ftid_Fbin_arr   = np.zeros((nr_bins,nr_bins), dtype='d')

for Rc in range(0,nr_bins):
    for Mc in range(0,nr_bins):
        #m and R:
        R_si    = (R_AU_arr[Rc])*AU_SI
        M_si    = (M_msun_arr[Mc])*M_sun_SI
        #tm/Torb: 
        const   = ((768./425.)*(5./256.)*(2.**(11./4.)/np.sqrt(85.))*(np.sqrt(6.)/(np.pi**(5.)))*np.sqrt(6.))*(c_SI**(15./2.)/(G_new_SI**2.))
        tm_Torb = const*(M_si**(-2.))*(fGW_HZ**(-7./2.))*(R_si**(-3./2.))
        tm_Torb_arr[Rc,Mc]  = tm_Torb
        #Ftid/Fbin:
        const_a0    = (6.*(2.**(7./6.)))/(85.*(np.pi**(10./3.)))*((c_SI**(5.))/(G_new_SI**(4./3.))) 
        a0_si       = const_a0*(M_si**(-4./3.))*(fGW_HZ**(-7./3.))
        Ftid_Fbin   = (a0_si/R_si)**(3.)
        Ftid_Fbin_arr[Rc,Mc]  = Ftid_Fbin

#----------------------------------------------------------
#Calc:
#----------------------------------------------------------
d3b         = 7./9.
Deltafac    = 1.0 - d3b
Rmax_AU_C1  = np.zeros(nr_bins, dtype='d')
Rmin_AU_C1  = np.zeros(nr_bins, dtype='d')
Rmax_AU_C2  = np.zeros(nr_bins, dtype='d')
Rmin_AU_C2  = np.zeros(nr_bins, dtype='d')

R_M_mask_C1 = np.zeros((nr_bins,nr_bins), dtype='i')
R_M_mask_C2 = np.zeros((nr_bins,nr_bins), dtype='i')

for nc in range(0,nr_bins):
    
    Ac      = (7.**(7./2.))*85.*(G_new_SI**2.)/(((10.*Deltafac)**(7./2.))*9.*np.pi*(c_SI**5.))
    M_si    = M_msun_arr[nc]*M_sun_SI
    
    aHB_AU_C1   = ((3./2.)*(G_new_SI*M_si)/((v_disp_C1)**(2.)))/AU_SI
    aEJ_AU_C1   = ((1./6.)*((1./d3b) - 1.)*(G_new_SI*M_si)/((v_disp_C1*f_disp_esc)**(2.)))/AU_SI
    aGW_AU_C1   = ((Ac**(1./5.))*(M_si**(2./5.))*(v_disp_C1**(1./5.))/(n_BH_m3**(1./5.)))/AU_SI
    
    aHB_AU_C2   = ((3./2.)*(G_new_SI*(M_si))/((v_disp_C2)**(2.)))/AU_SI
    aEJ_AU_C2   = ((1./6.)*((1./d3b) - 1.)*(G_new_SI*M_si)/((v_disp_C2*f_disp_esc)**(2.)))/AU_SI
    aGW_AU_C2   = ((Ac**(1./5.))*(M_si**(2./5.))*(v_disp_C2**(1./5.))/(n_BH_m3**(1./5.)))/AU_SI
    
    Rmax_AU_C1[nc]  = aHB_AU_C1
    Rmin_AU_C1[nc]  = max([aEJ_AU_C1, aGW_AU_C1])
    
    Rmax_AU_C2[nc]  = aHB_AU_C2
    Rmin_AU_C2[nc]  = max([aEJ_AU_C2, aGW_AU_C2])
    
    #mask array 1:
    pos_GTRmin      = np.where(R_AU_arr[:] > Rmin_AU_C1[nc])[0]
    pos_LTRmax      = np.where(R_AU_arr[:] < Rmax_AU_C1[nc])[0]
    pos_RminRmax    = list(set(pos_GTRmin).intersection(pos_LTRmax))    
    R_M_mask_C1[pos_RminRmax,nc]    = 1

    #mask array 2:
    pos_GTRmin      = np.where(R_AU_arr[:] > Rmin_AU_C2[nc])[0]
    pos_LTRmax      = np.where(R_AU_arr[:] < Rmax_AU_C2[nc])[0]
    pos_RminRmax    = list(set(pos_GTRmin).intersection(pos_LTRmax))    
    R_M_mask_C2[pos_RminRmax,nc]    = 1

#----------------------------------------------------------




#----------------------------------------------------------
#PLOTs:
#----------------------------------------------------------

fig = plt.figure(figsize=(6.0, 3.5))
ax  = fig.add_subplot(111)

#oplot contours:
X, Y        = np.meshgrid(log10_R_AU_arr, M_msun_arr)

#HB/EJ limits:
#cluster model C1:
ax.plot(np.log10(Rmax_AU_C1[:]), M_msun_arr[:], linestyle = '--', linewidth = 1, color='red', alpha=1)
ax.plot(np.log10(Rmin_AU_C1[:]), M_msun_arr[:], linestyle = '--', linewidth = 1, color='red', alpha=1)
ZC  = np.transpose(R_M_mask_C1)
CS  = ax.contourf(X, Y, ZC, 3, colors=['white', 'red'], alpha=0.5)
#cluster model C2:
ax.plot(np.log10(Rmax_AU_C2[:]), M_msun_arr[:], linestyle = '--', linewidth = 1, color='orange', alpha=1)
ax.plot(np.log10(Rmin_AU_C2[:]), M_msun_arr[:], linestyle = '--', linewidth = 1, color='orange', alpha=1)
ZC  = np.transpose(R_M_mask_C2)
CS  = ax.contourf(X, Y, ZC, 3, colors=['white', 'orange'], alpha=0.5)


#Tidal threshold (mark where Ftid/Fbin > 1):
Z           = np.transpose(Ftid_Fbin_arr)
[pi,pj]     = np.where(Z > 1.0)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '\\\\\\'], alpha=0)
CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '//////'], alpha=0)


#tm vs Torb:
Z           = np.transpose(np.log10(tm_Torb_arr))
#quad        = plt.pcolormesh(X, Y, Z, cmap='cubehelix', linewidth=0, rasterized=True, vmin=-2, vmax=0)
CS          = plt.contour(X, Y, Z, list([np.linspace(-2, 0.0, num=11)]), colors='teal', linestyles=':', linewidths = 2)
CS          = plt.contour(X, Y, Z, list([-2.0, -1.0, 0.0]), colors=['teal', 'teal', 'teal'], linestyles=['-', '-', '-'], linewidths = 2)
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.2f')


#plot settings etc.:
if (fGW_HZ == 1.):
    ax.text(0.75, 0.9, r'$f_{GW} = 1\ Hz$', ha='left', fontsize = 12, bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'grey', 'pad': 2}, transform=ax.transAxes)
    ax.set_xlim([min(log10_R_AU_arr), max(log10_R_AU_arr)])
    ax.set_ylim([0.,100.])

if (fGW_HZ == 3.):
    ax.text(0.75, 0.9, r'$f_{GW} = 3\ Hz$', ha='left', fontsize = 12, bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'grey', 'pad': 2}, transform=ax.transAxes)
    ax.set_xlim([min(log10_R_AU_arr), max(log10_R_AU_arr)])
    ax.set_ylim([0.,100.])

if (fGW_HZ == 10.):
    ax.text(0.75, 0.9, r'$f_{GW} = 10\ Hz$', ha='left', fontsize = 12, bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'grey', 'pad': 2}, transform=ax.transAxes)
    ax.set_xlim([min(log10_R_AU_arr), max(log10_R_AU_arr)])
    ax.set_ylim([0.,100.])
        
    
    
ax.set_xlabel(r'Binary Seperation $R$ $[\log\ AU]$', fontsize = 10)
ax.set_ylabel(r'Black Hole Mass $M$ $[M_{\odot}]$ ', fontsize = 10)



#save and show:
plt.savefig('FIG_tm_Torb_fHZ_' + str(int(fGW_HZ)) + '.pdf', bbox_inches='tight')     
plt.show()



exit()
#----------------------------------------------------------
































