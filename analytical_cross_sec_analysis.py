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



#mass ratio dependence:

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

fig = plt.figure(figsize=(5.0, 4.0))

qmin    = 1./5.0
qmax    = 5.0/1.
q_arr = np.linspace(qmin, qmax, num=10000)

termA_qarr      = q_arr*(3.*q_arr/(2.+q_arr))**(1./7.)
termB_qarr      = ((2.+q_arr)/3.)
termC_qarr      = ((1.+q_arr)/(2.*q_arr))**(1./3.)

beta            = 7./2.
f_tid           = 0.5
#EM case:
apu_EM          = ((f_tid/2.)**(1./3.))*((2.)**(2./3.)) + 1.
insp_I_EM       = (1.05/(1.-1.7/beta))*np.log(apu_EM)*(apu_EM-1.)**(-(1./(beta+1.)))
#as a function of q:
apu_qarr        = ((f_tid/2.)**(1./3.))*((2.*q_arr)**(2./3.)) + 1.
insp_I_qarr     = (1.05/(1.-1.7/beta))*np.log(apu_qarr)*(apu_qarr-1.)**(-(1./(beta+1.)))
#scaled:
insp_I_qarr_UEM = insp_I_qarr/insp_I_EM


termTOT_qarr    = insp_I_qarr_UEM*termA_qarr*termB_qarr*termC_qarr


#oplot guideline:
fig.add_subplot(111).plot([1e-10,1e10], [1.0,1.0],  marker='', linestyle=':', alpha=0.75, linewidth=1.0, color='black')

fig.add_subplot(111).plot(q_arr, termA_qarr,        linestyle='-', alpha=0.75, linewidth=2.0, color='black',        label=r'kinematic term')
fig.add_subplot(111).plot(q_arr, termB_qarr,        linestyle='-', alpha=0.75, linewidth=2.0, color='blue',         label=r'grav. foc. term')
fig.add_subplot(111).plot(q_arr, termC_qarr,        linestyle='-', alpha=0.75, linewidth=2.0, color='red',          label=r'tidal radius term')
fig.add_subplot(111).plot(q_arr, insp_I_qarr_UEM,   linestyle='-', alpha=0.75, linewidth=2.0, color='green',        label=r'insp. integral ${\mathscr{I}}/\tilde{{\mathscr{I}}}$')
fig.add_subplot(111).plot(q_arr, termTOT_qarr,      linestyle='-', alpha=0.75, linewidth=2.0, color='purple',       label=r'combined')


fig.add_subplot(111).set_title(r'GW cross section weight terms')
fig.add_subplot(111).set_xlabel(r'$q \equiv m_{\rm A}/m_{\rm B}$')
fig.add_subplot(111).set_ylabel(r'$value$')


fig.add_subplot(111).legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)

fig.add_subplot(111).set_xlim([qmin,qmax])    
fig.add_subplot(111).set_ylim([0.1,10.])    

fig.add_subplot(111).set_xscale('log')
fig.add_subplot(111).set_yscale('log')

fig.savefig('unequal_mass_case_ex_1_qterms.eps', bbox_inches='tight')

plt.show()
exit()











