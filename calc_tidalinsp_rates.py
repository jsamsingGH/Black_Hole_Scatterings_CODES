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



c_SI       = 299792458.0        #m/s
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)

#INSPIRALS + COLL:
a_amin_amax = np.array([1e-3, 1e0])     #AU amin,amax
vel_cs      = 10.0                       #km/s

mWD         = 0.5
RWD_Rsun    = 0.013*((1.43/mWD)**(1./3.))*((1.-mWD/1.43)**(0.447))

m1      = mWD*M_sun_SI
m2      = 1.4*M_sun_SI
m3      = 1.4*M_sun_SI

mi      = m1
mj      = m2
mk      = m3     
mu_ij   = mi*mj/(mi+mj)
mu_12   = m1*m2/(m1+m2)
m_ij    = mi+mj
m_bs    = m1+m2+m3

eps     = 0.5
x_tide  = 0.0   #T_2(eta) \approx Es*\eta**(-x_tide)
beta    = 6. + 3.*x_tide/2.
Ms      = mj*((mj/mu_ij)**(x_tide/4.))
Rs      = RWD_Rsun*R_sun_SI

Mpfac   = (((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
f_tid   = 0.5
apu     = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
insp_Ip = ((1.05/(1.-1.7/beta))*np.log(apu)*(apu-1.)**(-(1./(beta+1.))))/np.log(apu)
Nfac    = 2.0  #for WD=0.5, NS=1.4 N \approx 2.
cs_ijk_R_AU2    = (Nfac*((((m3/mu_12)**(1./3.))*(2.*np.pi*G_new_SI*m_bs*Rs)/((1000*vel_cs)**2.))*((m1*m2)/(mi*mj)))/(AU_SI**2.))*np.log(apu)
cs_ijk_insp_AU2 = 2.*cs_ijk_R_AU2*(eps**(1./beta))*insp_Ip*Mpfac*((a_amin_amax*AU_SI/Rs)**(1./beta))
cs_ijk_insp_m2  = cs_ijk_insp_AU2*(AU_SI**2.)

#ONLY COLL:
cs_ijk_coll_AU2 = 2.*(Nfac*((((m3/mu_12)**(1./3.))*(2.*np.pi*G_new_SI*m_bs*Rs)/((1000*vel_cs)**2.))*((m1*m2)/(mi*mj)))/(AU_SI**2.))*np.log(apu)
cs_ijk_coll_m2  = cs_ijk_coll_AU2*(AU_SI**2.)

print cs_ijk_insp_AU2, cs_ijk_coll_AU2

#calc rates:
Ntot        = 1000.0
pc_m        = 3.086*(10.**(16.))            #pc in meters
V_GC        = (0.1*pc_m)**(3.)              #typical vol of a GC core      
nrGC_Gpc3   = 1.*(1000.*1000.*1000)         #nr GC per Gpc^-3
sec_year    = 31536000.

fbin        = 0.3
Nbin        = fbin*Ntot
Nsin        = (1.-fbin)*Ntot

print 'INSPIRAL rate: (per Gpc^3 per year)'
print (nrGC_Gpc3*sec_year)*(((vel_cs*1000.)*Nbin*Nsin/(V_GC))*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]**(1./beta)) - np.log(a_amin_amax[0]**(1./beta))))
print 'INSPIRAL rate: (per galaxy per year)'
print (250.*sec_year)*(((vel_cs*1000.)*Nbin*Nsin/(V_GC))*(cs_ijk_insp_m2[1]-cs_ijk_insp_m2[0])/(np.log(a_amin_amax[1]**(1./beta)) - np.log(a_amin_amax[0]**(1./beta))))


print 'COLL rate: (per Gpc^3 per year)'
print (nrGC_Gpc3*sec_year)*((vel_cs*1000.)*Nbin*Nsin/(V_GC))*cs_ijk_coll_m2
print 'COLL rate: (per galaxy per year)'
print (250.*sec_year)*((vel_cs*1000.)*Nbin*Nsin/(V_GC))*cs_ijk_coll_m2








