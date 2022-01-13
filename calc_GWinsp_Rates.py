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



#EQUAL MASS CASE:

betaGW      = 7./2.
a_amin_amax = np.array([1e-2, 1e0])     #AU amin,amax
vel_cs      = 10.0                      #km/s
mEM         = 30.0                      #Msun

#define:
m1      = M_sun_SI*mEM
m2      = M_sun_SI*mEM
m3      = M_sun_SI*mEM
mi      = M_sun_SI*mEM
mj      = M_sun_SI*mEM  
mk      = M_sun_SI*mEM      
mu_ij   = mi*mj/(mi+mj)
mu_12   = m1*m2/(m1+m2)
m_ij    = mi+mj
m_bs    = m1 + m2 + m3
eps     = 85.*np.pi/96.
Rsch    = 2.*G_new_SI*m_ij/(c_SI**2.)

Mpfac   = (((mu_ij/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./betaGW)
f_tid   = 0.5
apu     = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
insp_Ip = ((1.05/(1.-1.7/betaGW))*np.log(apu)*(apu-1.)**(-(1./(betaGW+1.))))/np.log(apu)
Nfac    = 4.0  #for EM N \approx 5.
cs_ijk_R_AU2    = (Nfac*((((m3/mu_12)**(1./3.))*(2.*np.pi*G_new_SI*m_bs*Rsch)/((1000*vel_cs)**2.))*((m1*m2)/(mi*mj)))/(AU_SI**2.))*np.log(apu)
cs_ijk_insp_AU2 = 3.*cs_ijk_R_AU2*(eps**(1./betaGW))*insp_Ip*Mpfac*((a_amin_amax*AU_SI/Rsch)**(1./betaGW))
cs_ijk_insp_m2  = cs_ijk_insp_AU2*(AU_SI**2.)

#calc rates:
Ntot        = 1000.0
fbin        = 0.5
pc_m        = 3.086*(10.**(16.))            #pc in meters
V_GC        = (0.1*pc_m)**(3.)              #typical vol of a GC core      
nrGC_Gpc3   = 1.*(1000.*1000.*1000)         #nr GC per Gpc^-3
sec_year    = 31536000.


print cs_ijk_insp_AU2

print (fbin*(1.-fbin)/2.)*(betaGW*(vel_cs*1000.)*(1./V_GC)*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]/a_amin_amax[0])))*(sec_year*(10.**(6)))
print np.sqrt((fbin*(1.-fbin)/2.)*(betaGW*(vel_cs*1000.)*(1./V_GC)*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]/a_amin_amax[0])))*(sec_year*(10.**(6))))**(-1.)

print (fbin*(1.-fbin)/2.)*(betaGW*(vel_cs*1000.)*(1./V_GC)*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]/a_amin_amax[0])))*(nrGC_Gpc3*sec_year)
print (np.sqrt((fbin*(1.-fbin)/2.)*(betaGW*(vel_cs*1000.)*(1./V_GC)*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]/a_amin_amax[0])))*(nrGC_Gpc3*sec_year)))**(-1.)

print (Ntot**2.)*(fbin*(1.-fbin)/2.)*(betaGW*(vel_cs*1000.)*(1./V_GC)*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]/a_amin_amax[0])))*(nrGC_Gpc3*sec_year)

print 10.**(3./2.)

print (fbin*(1.-fbin)/2.)

print ((Ntot/80.)**2.)*(1.4/30.)**(12./7.)













