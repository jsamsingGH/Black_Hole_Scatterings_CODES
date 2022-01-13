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
c_SI        = 299792458.0        #m/s
M_sun_SI    = 1.989*(10.**30.)   #kg
R_sun_SI    = 695800000.         #m
AU_SI       = 149597871000.      #m 
G_new_SI    = 6.67*(10.**(-11.))
AU_U        = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U     = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
sec_year    = 31536000.
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
#----------------------------------------------------------


q       = np.arange(0.1, 1.0, 0.001)
eta     = q/((1.+q)**2.)
Afac    = 1.2*(10.**4)
Bfac    = -0.93 
vm      = Afac*(eta**2.)*((1.-q)/(1.+q))*(1.+Bfac*eta)

fig = plt.figure(figsize=(6, 5))
fig.add_subplot(111).plot(q, vm)
plt.show()
exit()



Pcap    = 2.5*(2.*20.)*(((10.)*Rsch_1Msun_unitRsun*R_sun_SI)/((10.**(-2.))*AU_SI))**(5./7.)
Pcoll   = (4.*20.)*(((10.)*Rsch_1Msun_unitRsun*R_sun_SI)/((10.**(-2.))*AU_SI))**(1./1.)
print Pcap, Pcoll, Pcap/Pcoll
exit() 


Nims    = 20.0
m_BH    = 20.0
a_BBH   = 0.1
f_ob    = 10.0
e_ob    = 0.1

Fe_ob   = (e_ob**(12./19.)/(1.+e_ob))*(1. + (121./304.)*(e_ob**2.))**(870./2299.)

rf      = ((2.*G_new_SI*(m_BH*M_sun_SI)/((f_ob**2.)*(np.pi**2.)))**(1./3.))/AU_SI
rEC     = rf*(1./(2.*Fe_ob))*((425./304.)**(870./2299.))

print rf
exit()  

PEM     = Nims*(2.*(rEC/a_BBH))
print PEM, PEM*(9./2.)
print (1./(2.*Fe_ob))*((425./304.)**(870./2299.))
print (9./2.)

print np.log10(((1./6.)*(G_new_SI*(m_BH*M_sun_SI)/((50.*1000)**2.))*(1./(7./9.) - 1.0))/AU_SI)

delta   = 7./9.
vesc    = 50.*1000.     #m/s
print Nims*(delta/((1.-delta)**2.))*(12.*(2.**(1./3.))/(np.pi**(2./3.)))*((1./(2.*Fe_ob))*((425./304.)**(870./2299.)))*(1./(G_new_SI**(2./3.)))*(1./(f_ob**(2./3.)))*(1./((m_BH*M_sun_SI)**(2./3.)))*(vesc**(2.))
print ((1./6.)*((1./delta) - 1.)*(G_new_SI*(m_BH*M_sun_SI)/(vesc**2.)))/AU_SI

exit()

print (((6.*2**(5./3.))/(85.*np.pi**(10./3.)))*((c_SI**5.)/(G_new_SI**(4./3.)))*((m_BH*M_sun_SI)**(-4./3.))*((f_ob)**(-7./3.)))/AU_SI


print ((85.*np.pi/(6.*np.sqrt(64.)))**(2./7.))*(1.0*AU_SI/(20.*Rsch_1Msun_unitRsun*R_sun_SI))**(2./7.)
print ((85.*np.pi/(6.*np.sqrt(64.)))**(2./7.))
print rEC/(20.*Rsch_1Msun_unitRsun/AU_U)
exit()


#calc a0t:
m_SI    = 20.0*M_sun_SI
tau_SI  = (10.**10.)*sec_year 
aot_AU     = ((tau_SI*4.*((64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.)))**(1./4.))/AU_SI
print aot_AU



exit()






print ((3./2.)*(G_new_SI*(20.*M_sun_SI))/(20000.0**2.))/AU_SI
print 9./(1./(0.8) - 1.0)
print (1./(1.-0.9))*(np.log(50.))
exit()

T0  = 2.*np.pi*np.sqrt(((0.1*AU_SI)**3.)/(2.*G_new_SI*(30.*M_sun_SI)))
NI  = 10.
fGW = 10.0

print (32.**(1./3.))*NI*(fGW*T0)**(-2./3.)
print T0

eG = 0.1
F_eG = (eG**(12./19.)/(1.+eG))*(1. + (121./304.)*(eG**2.))**(870./2299.) 

print NI*((2.*G_new_SI/(np.pi**2.))**(1./3.))*((425./304.)**(870./2299.))*((30.*M_sun_SI)**(1./3.))*(1./fGW**(2./3.))*(1./F_eG)*(1./(0.1*AU_SI))

print (1./(2.*F_eG))*(1. + (121./304.))**(870./2299.) 

print 0.00268*2.669

exit()


m_SI    = 1.0*M_sun_SI
a0_SI   = 1.0*AU_SI
v_SI    = 10000.

C_GW    = (85.*np.pi/(3.*np.sqrt(3.0)))**(2./7.)
F_fac   = 6.0
betaGW  = 7./2.
apu     = 2.0
insp_I = ((1.05/(1.-1.7/betaGW))*np.log(apu)*(apu-1.)**(-(1./(betaGW+1.))))

cs_GWinsp_all_AU2   = 3.0*(C_GW*6.*np.pi*insp_I*F_fac)*((G_new_SI**(12./7.))*(m_SI**(12./7.))*(a0_SI**(2./7.))/((c_SI**(10./7.))*(v_SI**2.)))/(AU_SI**2.)
print cs_GWinsp_all_AU2

#calc rates:
pc_m        = 3.086*(10.**(16.))            #pc in meters
ac_AU       = 1.0
V_GC        = (0.1*pc_m)**(3.)              #typical vol of a GC core      
nrGC_Gpc3   = 1.*(1000.*1000.*1000)         #nr GC per Gpc^-3
sec_year    = 31536000.
vel_SI      = 10000.
s_sdev      = 0.5
f_bin       = 0.5
cs_I_norm   = ((0.025*(20.0)**(12./7.)*(ac_AU)**(2./7.)))*(AU_SI**2.) #for 10 km/s



print ((f_bin*(1.-f_bin)/2.)*(vel_SI/V_GC)*cs_I_norm*np.exp(((np.log(10.)**2.)/2.)*(s_sdev/betaGW)**2.)*(nrGC_Gpc3*sec_year))**(-1./2.)



#calc a0t:
m_SI    = 1.0*M_sun_SI
tau_SI  = (10.**10.)*sec_year 

aot_AU     = ((tau_SI*4.*((64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.)))**(1./4.))/AU_SI

print aot_AU






m_SI    = 20.0*M_sun_SI
v_SI    = 10000.
aHB     = ((3./2.)*G_new_SI*m_SI/(v_SI**2.))/AU_SI

print  aHB

Rsch = (2.*G_new_SI*m_SI)/(c_SI**2.)

print ((10.**10.)*((3.3*10**(-2.))**(7./2.)))*(20.**(-1./2.))*(25.)**(3./2.)

print (1.-(1.-(30.*Rsch)/(0.1*AU_SI))**2.)*100.


print 2./(2.*np.pi*np.sqrt(((30.*Rsch)**3.)/(G_new_SI*2.*m_SI)))
print (1./np.pi)*np.sqrt(G_new_SI*(2.*m_SI)/((30.*Rsch)**3.))*(2.0**(1.1954-1.5))


print (((2.9*(10**(-4.)))/0.01)**(7./5.))*(20.)

print (2.*10**(6.))*(0.005/25.)**(3./2.)

print 0.15*((0.1/0.01)**(-7./5.))

fGW = 10.
print (((1./(fGW*np.pi))*np.sqrt(G_new_SI*(2.*m_SI))*(2.0)**(1.1954-1.5))**(2./3.))/AU_SI


RC  = 1e-5
print 3.0*((RC*AU_SI)*(6.*(2.*np.pi*G_new_SI*(3.*m_SI)/(v_SI**2.))*np.log(2.0)))/(AU_SI**2.)

print 0.4/(0.025*(20.)**(12./7.))

print (0.1)**(-2./7.)




















