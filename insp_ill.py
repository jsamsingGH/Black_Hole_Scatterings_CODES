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
from scipy.integrate import quad
import matplotlib as mpl


#-----------------------------------------------------------------
#integral function:
#-----------------------------------------------------------------
def func_ap_int(a, b):
    return (a**(1./b - 1.))*((a-1.)**(-3./(2.*b)))
#-----------------------------------------------------------------

    
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


c_SI       = 299792458.0        #m/s
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)


#input:
m1      = 0.4
m2      = 1.4
m3      = 1.4 

mi      = 1.4
mj      = 1.4
mk      = 0.4

a0_Rsun = (1e-4)*AU_U

#define:
m_bs    = m1+m2+m3
m_ij    = mi+mj
mu_ij   = mi*mj/m_ij

#for GW inspirals: 
beta    = 7./2.
Ms      = mu_ij
Rs      = (2.*G_new_SI*(m_ij*M_sun_SI)/(c_SI**2.))/R_sun_SI
Es      = np.pi*(85./96.)
Mfac    = ((m1*m2)/(mi*mj))*(((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
#apu     = (mk/mu_ij)*((m_ij/m_bs)**(1./3.)) + 1.
apu     = ((mk/mu_ij)**(2./3.))*0.75 + 1.

print apu, 'apu'

#plot results:
ap_arr  = np.linspace(1.0, apu, num=10000)
ecc     = 1. - Mfac*(Es**(1./beta))*((a0_Rsun/Rs)**(1./beta - 1.))*((ap_arr**(1./beta-1.))*(ap_arr-1.)**(-3./(2.*beta)))

fig = plt.figure(figsize=(5, 5))
fig.add_subplot(111).plot(ap_arr,ecc)
fig.add_subplot(111).plot(ap_arr,((ap_arr**(1./beta-1.))*(ap_arr-1.)**(-3./(2.*beta))) )
fig.add_subplot(111).plot([apu,apu],[-10,10])
fig.add_subplot(111).set_ylim(0.0,1.0)


#plt.show()


#-----------------------------------------------------------------
#inspiral integral:
#-----------------------------------------------------------------
fig = plt.figure(figsize=(5, 3))
ax  = fig.add_subplot(111)

#settings:
nres        = 1000
apu_arr     = np.linspace(1.0001, 8.0, num=nres)
beta_arr    = np.array([7./2., 6., 10.0]) 
beta_txt    = ['7/2 (GWs)', '6.0 (tides)', '10.0 (tides)'] 
nbeta       = len(beta_arr)
int_arr     = np.zeros(nres, dtype=np.float64)
plot_beta_colors = ['dodgerblue','brown', 'orange']

#calc int per beta:
for bc in range(0,nbeta):
    for ic in range(0,nres):
        b           = beta_arr[bc]
        apl         = 1.0
        apu         = apu_arr[ic]
        ap_int      = quad(func_ap_int, apl, apu, args=(b))[0]
        int_arr[ic] = ap_int          
    #functional app:
    int_func_app    = (1.05/(1.-1.7/b))*np.log(apu_arr)*(apu_arr-1.)**(-(1./(b+1.)))
    #plot:
    ax.plot(apu_arr, int_arr, linestyle='-',        color=plot_beta_colors[bc])#, label=r'numerical solution ($\beta = 7/2$)')
    ax.plot(apu_arr, int_func_app, linestyle=':',   color=plot_beta_colors[bc])#, label=r'functional approximation')
    #text for different beta:
    ax.text(0.7, 0.35-0.1*bc, r'$\beta$ = ' + beta_txt[bc],
            verticalalignment='bottom', horizontalalignment='left',
            transform=ax.transAxes,
            color=plot_beta_colors[bc], fontsize=10)

ax.plot(apu_arr, np.log(apu_arr), linestyle='--',   color='black')#, label=r'ln$(a^{\prime}_{u})$')
ax.text(0.7, 0.35-0.1*3, r'$\beta = \infty$ (collision)',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10)

#fake plot for making line legends:           
#ax.plot([0,0], [0,0], linestyle='-',    color='black', label=r'numerical solution')
#ax.plot([0,0], [0,0], linestyle=':',    color='black', label=r'approximation')

#labels, limits etc:
ax.set_xlabel(r'$a^{\prime}_{u}$')
ax.set_ylabel(r'$\mathscr{I}(a^{\prime}_{u}, \beta)$')
ax.set_ylim(0.0,3.0)
ax.set_xlim(1.0,8.)
ax.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)

#save and show:
plt.savefig('insp_integral.eps', bbox_inches='tight')        
plt.show()
exit()
#-----------------------------------------------------------------










nres    = 1000
b_arr = np.linspace(3.5, 10., num=nres)
int_arr = np.zeros(nres, dtype=np.float64)
#set:
#calc int:
for ic in range(0,nres):
    apl         = 1.0
    apu         = 10.0
    ap_int      = quad(func_ap_int, apl, apu, args=(b_arr[ic]))[0]/(np.log(apu)*(apu-1.)**(-(1./(b_arr[ic]+1.))))
    int_arr[ic] = ap_int
#functional app:
#int_func_app    = 2.0*np.log(apu_arr)*(apu_arr-1)**(-1./4.7)
#int_func_app    = 1.275*np.log(apu_arr)*(apu_arr-1)**(-(1./(b+1.)))
int_func_app    = 1.08/(1.-1.7/b_arr)

fig = plt.figure(figsize=(5, 3))
fig.add_subplot(111).plot(b_arr, int_arr, linestyle='-', color='black')
fig.add_subplot(111).plot(b_arr, int_func_app, linestyle=':', color='black')

plt.show()
exit()

 













