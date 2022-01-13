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
from mpl_toolkits.axes_grid1 import make_axes_locatable
import matplotlib as mpl
import numpy.ma as ma



def func_insp_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, Rtid_i_Rsun, Es, f_tid, beta, Ms, ap):
    #1,2 is in ini binary - 3 is incoming
    #i,j is in IMS binary, and i is the tidal object (star).
    m_bs        = m1+m2+m3
    m_ij        = mi+mj
    mu_ij       = mi*mj/m_ij
    mu_12       = m1*m2/(m1+m2)
    apu         = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
    
    Mfac        = ((m1*m2)/(mi*mj))*(((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
    insp_I      = (1.05/(1.-1.7/beta))*(np.log(apu)/((apu-1.)**((1./(beta+1.)))))
    Einsp       = (Es**(1./beta))*Mfac*((a0_Rsun/Rtid_i_Rsun)**(1./beta - 1.))*(ap**(1./beta -1.))*((ap-1.)**(-3./(2.*beta)))
    Ainsp       = (Es**(1./beta))*insp_I*Mfac*((a0_Rsun/Rtid_i_Rsun)**(1./beta - 1.))
    return [Einsp, Ainsp]

    

def func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, Rcoll_i_Rsun, f_tid, ap):
    #1,2 is in ini binary - 3 is incoming
    #i,j is in IMS binary, and i is the tidal object (star).
    m_bs        = m1+m2+m3
    m_ij        = mi+mj
    mu_ij       = mi*mj/m_ij
    mu_12       = m1*m2/(m1+m2)
    apu         = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
    
    coll_I      = np.log(apu)
    Ecoll       = (Rcoll_i_Rsun/a0_Rsun)*((m1*m2)/(mi*mj))*(1./ap)
    Acoll       = (Rcoll_i_Rsun/a0_Rsun)*((m1*m2)/(mi*mj))*coll_I
    return [Ecoll, Acoll]





#-----------------------------------------------------------------
#Units and conversions:
#-----------------------------------------------------------------
c_SI       = 299792458.0        #m/s
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
#-----------------------------------------------------------------





#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)







#-----------------------------------------------------------------
#Example 4: (BC,NS,NS)
#-----------------------------------------------------------------
nres            = 100
loga0_AU_arr    = np.linspace(-3.0, 1.0, nres)
a0_Rsun_arr     = (10.**(loga0_AU_arr))*AU_U
#mass:
mA_Msun_arr     = np.linspace(0.5, 3.0, nres)   #m1
mB_Msun         = 1.4                           #m2
mC_Msun         = 1.4                           #m3
#radius:
rA_Rsun         = 1e-5
rB_Rsun         = 12000./R_sun_SI
rC_Rsun         = 12000./R_sun_SI

insp_over_coll_arr  = np.zeros((nres,nres), dtype=np.float64)

for ac in range(0,nres):        #sma a
    for mc in range(0,nres):    #mass
        #read:
        a0_Rsun     = a0_Rsun_arr[ac]       
        mA_Msun     = mA_Msun_arr[mc]
        #initial config:
        m1      = mA_Msun                   #bin 1
        m2      = mB_Msun                   #bin 2
        m3      = mC_Msun                   #sin 3
        #[i,j]-k endstate config:
        mi      = mB_Msun                   #insp i
        mj      = mC_Msun                   #insp j
        mk      = mA_Msun                   #bound single k

        R_collij    = rB_Rsun + rC_Rsun 

        #GW Inspirals:
        f_tid       = 0.5
        m_bs        = m1+m2+m3
        m_ij        = mi+mj
        mu_ij       = mi*mj/m_ij
        mu_12       = m1*m2/(m1+m2)
        apu         = 1. + ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.))
        Es          = 85.*np.pi/96.
        Ms          = mu_ij
        Rs          = (2.*G_new_SI*(M_sun_SI*m_ij)/(c_SI**2.))/R_sun_SI 
        beta        = 7./2.
        p_Mfac      = (((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
        p_inspI     = ((1.05/(1.-1.7/beta))*(np.log(apu)/((apu-1.)**((1./(beta+1.))))))/np.log(apu)
        insp_over_coll_arr[ac,mc]   = np.log10((Rs/R_collij)*((Es**(1./beta))*p_inspI*p_Mfac*((a0_Rsun/Rs)**(1./beta))))

print insp_over_coll_arr[:,99]

fig     = plt.figure(figsize=(5.5, 4))
#The formating here is quite a mess: It seems like we have to flip the two axis around to get it to work.
#Prepare x,y,z for plotting:
X, Y    = np.meshgrid(loga0_AU_arr, mA_Msun_arr)
Z       = np.transpose(insp_over_coll_arr)
#plot colors:
ax = plt.subplot(111)
ax.set_xlabel(r'log $a_{\rm 0}$ [AU]')
ax.set_ylabel(r'$m_{\rm 1}$ $[M_{\odot}]$')
ax.set_title(r'[BC$_{1}$,NS$_{2}$] $\leftarrow$ NS$_{3}$, log($\sigma_{I_{23}}/\sigma_{C_{23}}$)')
quad = plt.pcolormesh(X, Y, Z, cmap='bone', linewidth=0, rasterized=True, vmin=0.0, vmax=3.0)
plt.colorbar()
#oplot contours:
CS = plt.contour(X, Y, Z, list(np.linspace(1.0, 3.0, num=10+1)), colors='yellow')
plt.clabel(CS, fontsize=9, inline=1)
#save and show:
plt.savefig('GWinsp_coll_contour_BC1NS2NS3.eps', bbox_inches='tight')     
plt.show()

exit()
#-----------------------------------------------------------------








#-----------------------------------------------------------------
#Example 5: (BH,BH,BH)
#-----------------------------------------------------------------
nres            = 100
loga0_AU_arr    = np.linspace(-3.0, 1.0, nres)
a0_Rsun_arr     = (10.**(loga0_AU_arr))*AU_U
#mass:
mA_Msun         = 30.                           #m1
mB_Msun         = 20.                           #m2
mC_Msun_arr     = np.linspace(10.0, 40.0, nres) #m3
#radius:
Rsch_1Msun_Rsun = (2.*G_new_SI*M_sun_SI/(c_SI**2.))/R_sun_SI
rA_Rsun         = Rsch_1Msun_Rsun*mA_Msun
rB_Rsun         = Rsch_1Msun_Rsun*mB_Msun
rC_Rsun_arr     = Rsch_1Msun_Rsun*mC_Msun_arr

insp_over_coll_arr  = np.zeros((nres,nres), dtype=np.float64)

for ac in range(0,nres):        #sma a
    for mc in range(0,nres):    #mass
        #read:
        a0_Rsun     = a0_Rsun_arr[ac]       
        mC_Msun     = mC_Msun_arr[mc]
        rC_Rsun     = rC_Rsun_arr[mc]
        #initial config:
        m1      = mA_Msun                   #bin 1
        m2      = mB_Msun                   #bin 2
        m3      = mC_Msun                   #sin 3
        #[i,j]-k endstate config:
        mi      = mA_Msun                   #insp i
        mj      = mC_Msun                   #insp j
        mk      = mB_Msun                   #bound single k

        R_collij    = rA_Rsun + rC_Rsun

        #GW Inspirals:
        f_tid       = 0.5
        m_bs        = m1+m2+m3
        m_ij        = mi+mj
        mu_ij       = mi*mj/m_ij
        mu_12       = m1*m2/(m1+m2)
        apu         = 1. + ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.))
        Es          = 85.*np.pi/96.
        Ms          = mu_ij
        Rs          = (2.*G_new_SI*(M_sun_SI*m_ij)/(c_SI**2.))/R_sun_SI 
        beta        = 7./2.
        p_Mfac      = (((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
        p_inspI     = ((1.05/(1.-1.7/beta))*(np.log(apu)/((apu-1.)**((1./(beta+1.))))))/np.log(apu)
        insp_over_coll_arr[ac,mc]   = np.log10((Rs/R_collij)*((Es**(1./beta))*p_inspI*p_Mfac*((a0_Rsun/Rs)**(1./beta))))

print insp_over_coll_arr[:,99]

fig     = plt.figure(figsize=(5.5, 4))
#The formating here is quite a mess: It seems like we have to flip the two axis around to get it to work.
#Prepare x,y,z for plotting:
X, Y    = np.meshgrid(loga0_AU_arr, mC_Msun_arr)
Z       = np.transpose(insp_over_coll_arr)
#plot colors:
ax = plt.subplot(111)
ax.set_xlabel(r'log $a_{\rm 0}$ [AU]')
ax.set_ylabel(r'$m_{\rm 3}$ $[M_{\odot}]$')
ax.set_title(r'[BH$_{1}$,BH$_{2}$] $\leftarrow$ BH$_{3}$, log($\sigma_{I_{13}}/\sigma_{C_{13}}$)')
quad = plt.pcolormesh(X, Y, Z, cmap='bone', linewidth=0, rasterized=True, vmin=0.0, vmax=3.0)
plt.colorbar()
#oplot contours:
CS = plt.contour(X, Y, Z, list(np.linspace(1.0, 3.0, num=10+1)), colors='yellow')
plt.clabel(CS, fontsize=9, inline=1)
#save and show:
plt.savefig('GWinsp_coll_contour_BH1BH2BH3.eps', bbox_inches='tight')     
plt.show()

exit()
#-----------------------------------------------------------------





#-----------------------------------------------------------------
#Example 6: (MS,NS,NS)
#-----------------------------------------------------------------
Rcoll_yesno = 1
Rtde_yesno  = 0

nres            = 100
loga0_Rsun_arr  = np.linspace(-5.0, 2.0, nres)
a0_Rsun_arr     = (10.**(loga0_Rsun_arr))*AU_U
#mass:
mMS_Msun_arr    = np.linspace(0.5, 3.0, nres)   #m1 (i)
mP1_Msun        = 1.4                           #m2 (j)
mP2_Msun        = 1.4                           #m3 (k)
#MS radius:
RMS_Rsun_arr    = mMS_Msun_arr**(0.8)

insp_over_coll_arr  = np.zeros((nres,nres), dtype=np.float64)
insp_over_tde_arr   = np.zeros((nres,nres), dtype=np.float64)

for ac in range(0,nres):        #sma a
    for mc in range(0,nres):    #wd mass, radius
        #read:
        a0_Rsun     = a0_Rsun_arr[ac]       
        mMS         = mMS_Msun_arr[mc]      #bin 1
        R_i_Rsun    = RMS_Rsun_arr[mc]
        #initial config:
        m1      = mMS                       #bin 1
        m2      = mP1_Msun                  #bin 2
        m3      = mP2_Msun                  #sin 3
        #[i,j]-k endstate config:
        mi      = mMS                       #IMS i  STAR
        mj      = mP1_Msun                  #IMS j  COMPACT OBJ
        mk      = mP2_Msun                  #bound single k

        ap          = 1.0+0.5   #plays no role in this example
        f_tid       = 0.5
        #Tidal Inspirals:
        Es          = 0.15
        x_tide      = 3.5   #T_2(eta) \approx Es*\eta**(-x_tide)
        beta        = 6. + 3.*x_tide/2.
        Ms          = mj*(((mi+mj)/mi)**(x_tide/4.))
        A_insp          = func_insp_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_i_Rsun, Es, f_tid, beta, Ms, ap)[1]
        #Solid Sphere (SS) Collisions:
        R_SS_coll_i_Rsun    = R_i_Rsun  
        A_SScoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_SS_coll_i_Rsun, f_tid, ap)[1]
        #Tidal Disruption (TD) Collisions:
        R_TD_coll_i_Rsun    = max([R_i_Rsun, R_i_Rsun*(mj/mi)**(1./3.)])          
        A_TDcoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_TD_coll_i_Rsun, f_tid, ap)[1]  
        #Relative cs:
        insp_over_coll_arr[ac,mc]   = A_insp/A_SScoll
        insp_over_tde_arr[ac,mc]    = A_insp/A_TDcoll
        #MASK:
        if (a0_Rsun < 1.0*R_i_Rsun):    insp_over_coll_arr[ac,mc]   = np.ma.masked
        

#PLOT results:
if (Rcoll_yesno == 1):   
   
    #mask and chose colors:
    insp_over_coll_arr = ma.masked_invalid(insp_over_coll_arr)
    cmap = plt.cm.magma
    cmap.set_bad('grey',1) 

    fig     = plt.figure(figsize=(5.5, 4))
    #The formating here is quite a mess: It seems like we have to flip the two axis around to get it to work.
    #Prepare x,y,z for plotting:
    X, Y    = np.meshgrid(loga0_Rsun_arr, mMS_Msun_arr)
    Z       = np.transpose(insp_over_coll_arr)
    #plot colors:
    ax = plt.subplot(111)
    ax.set_xlabel(r'log $a_{\rm 0}$ [AU]')
    ax.set_ylabel(r'$m_{\rm MS}$ $[M_{\odot}]$')
    ax.set_title(r'[MS,NS]-NS')
    quad = plt.pcolormesh(X, Y, Z, cmap=cmap, linewidth=0, rasterized=True, vmin=1.0, vmax=5.0)
    plt.colorbar()
    #oplot black contours:
    CS = plt.contour(X, Y, Z, [1.0,1.5,2.0,2.5,3.0],
                     colors='lime',  # negative contours will be dashed by default
                     )
    plt.clabel(CS, fontsize=12, inline=1)
    #save and show:
    plt.savefig('Pinsp_Pcoll_contour_MSNSNS.eps', bbox_inches='tight')     
    plt.show()
    exit()


exit()
#-----------------------------------------------------------------


#-----------------------------------------------------------------
#Example 3: (WD,NS,NS)
#-----------------------------------------------------------------
Rcoll_yesno = 1
Rtde_yesno  = 0

nres            = 100
loga0_Rsun_arr  = np.linspace(-5.0, 2.0, nres)
a0_Rsun_arr     = (10.**(loga0_Rsun_arr))*AU_U
#mass:
mWD_Msun_arr    = np.linspace(0.4, 1.0, nres)   #m1 (i)
mP1_Msun        = 1.4                           #m2 (j)
mP2_Msun        = 1.4                           #m3 (k)
#WD radius:
RWD_Rsun_arr    = 0.013*((1.43/mWD_Msun_arr)**(1./3.))*((1.-mWD_Msun_arr/1.43)**(0.447))

insp_over_coll_arr  = np.zeros((nres,nres), dtype=np.float64)
insp_over_tde_arr   = np.zeros((nres,nres), dtype=np.float64)

for ac in range(0,nres):        #sma a
    for mc in range(0,nres):    #wd mass, radius
        #read:
        a0_Rsun     = a0_Rsun_arr[ac]       
        mWD         = mWD_Msun_arr[mc]      #bin 1
        R_i_Rsun    = RWD_Rsun_arr[mc]
        #initial config:
        m1      = mWD                       #bin 1
        m2      = mP1_Msun                  #bin 2
        m3      = mP2_Msun                  #sin 3
        #[i,j]-k endstate config:
        mi      = mWD                       #IMS i  STAR
        mj      = mP1_Msun                  #IMS j  COMPACT OBJ
        mk      = mP2_Msun                  #bound single k

        ap          = 1.0+0.5   #plays no role in this example
        f_tid       = 0.5
        #Tidal Inspirals:
        Es          = 0.5
        x_tide      = 0.0   #T_2(eta) \approx Es*\eta**(-x_tide)
        beta        = 6. + 3.*x_tide/2.
        Ms          = mj*(((mi+mj)/mi)**(x_tide/4.))
        A_insp          = func_insp_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_i_Rsun, Es, f_tid, beta, Ms, ap)[1]
        #Solid Sphere (SS) Collisions:
        R_SS_coll_i_Rsun    = R_i_Rsun  
        A_SScoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_SS_coll_i_Rsun, f_tid, ap)[1]
        #Tidal Disruption (TD) Collisions:
        R_TD_coll_i_Rsun    = max([R_i_Rsun, R_i_Rsun*(mj/mi)**(1./3.)])          
        A_TDcoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_TD_coll_i_Rsun, f_tid, ap)[1]  
        #Relative cs:
        insp_over_coll_arr[ac,mc]   = A_insp/A_SScoll
        insp_over_tde_arr[ac,mc]    = A_insp/A_TDcoll
        #MASK:
        if (a0_Rsun < 1.0*R_i_Rsun):    insp_over_coll_arr[ac,mc]   = np.ma.masked
        

#PLOT results:
if (Rcoll_yesno == 1):
   
    #mask and chose colors:
    insp_over_coll_arr = ma.masked_invalid(insp_over_coll_arr)
    cmap = plt.cm.magma
    cmap.set_bad('grey',1)
    
    fig     = plt.figure(figsize=(5.5, 4))
    #The formating here is quite a mess: It seems like we have to flip the two axis around to get it to work.
    #Prepare x,y,z for plotting:
    X, Y    = np.meshgrid(loga0_Rsun_arr, mWD_Msun_arr)
    Z       = np.transpose(insp_over_coll_arr)
    #plot colors:
    ax = plt.subplot(111)
    ax.set_xlabel(r'log $a_{\rm 0}$ [AU]')
    ax.set_ylabel(r'$m_{\rm WD}$ $[M_{\odot}]$')
    ax.set_title(r'[WD,NS]-NS')
        
    quad = plt.pcolormesh(X, Y, Z, cmap=cmap, linewidth=0, rasterized=True, vmin=1.0, vmax=5.0)
    plt.colorbar()
    #oplot black contours:
    CS = plt.contour(X, Y, Z, [1,2,3,4,5],
                     colors='lime',  # negative contours will be dashed by default
                     )
    plt.clabel(CS, fontsize=12, inline=1)
    #save and show:
    plt.savefig('Pinsp_Pcoll_contour_WDNSNS.eps', bbox_inches='tight')     
    plt.show()
    exit()

exit()
#-----------------------------------------------------------------

























#-----------------------------------------------------------------
#Example 1:
#-----------------------------------------------------------------
#1,2 is in ini binary - 3 is incoming
#i,j is in IMS binary, and i is the tidal object (star).
m1  = 0.6    #bin 1
m2  = 1.4    #bin 2
m3  = 1.4    #incoming 3
mi  = m1    #IMS i  STAR
mj  = m2    #IMS j  COMPACT OBJ
mk  = m3    #bound single k
R_i_Rsun    = 0.013*((1.43/mi)**(1./3.))*((1.-mi/1.43)**(0.447))       #radius of STAR (WD) (obj i)
f_tid       = 0.5
a0_AU       = 0.01           #AU
a0_Rsun     = a0_AU*AU_U    #Rsun

#Tidal Inspirals:
Es          = 0.1
x_tide      = 0.0   #3.5   #T_2(eta) \approx Es*\eta**(-x_tide)
beta        = 6. + 3.*x_tide/2.
Ms          = mj*((mi+mj)/mi)**(x_tide/4.)

#Solid Sphere (SS) Collisions:
R_SS_coll_i_Rsun    = R_i_Rsun  

#Tidal Disruption (TD) Collisions:
R_TD_coll_i_Rsun    = R_i_Rsun*(mj/mi)**(1./3.)  #we have removed th fac of 4 - people do that and refer to this as the r_t


#PLOTS:
ap_arr  = np.linspace(1.00001, 2.0, 1000)

#inspirals:
Einsp_arr       = func_insp_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_i_Rsun, Es, f_tid, beta, Ms, ap_arr)[0]
A_insp          = func_insp_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_i_Rsun, Es, f_tid, beta, Ms, ap_arr)[1]
eccinsp_arr     = 1.-Einsp_arr

#SS-collisions:
ESScoll_arr     = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_SS_coll_i_Rsun, f_tid, ap_arr)[0]
A_SScoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_SS_coll_i_Rsun, f_tid, ap_arr)[1]
eccSScoll_arr   = 1.-ESScoll_arr

#TD-collisions:
ETDcoll_arr     = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_TD_coll_i_Rsun, f_tid, ap_arr)[0]
A_TDcoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_TD_coll_i_Rsun, f_tid, ap_arr)[1]
eccTDcoll_arr   = 1.-ETDcoll_arr

print A_insp, A_SScoll, A_TDcoll
print A_insp/A_SScoll, A_insp/A_TDcoll

#PLOT insp boundary:
fig = plt.figure(figsize=(5, 4))
ax1 = fig.add_subplot(1,1,1)

arr_1   = 1.0+0.0*ap_arr
ax1.fill_between(ap_arr, eccinsp_arr,   eccTDcoll_arr, where=(eccinsp_arr   < eccTDcoll_arr),   facecolor='lightblue',  interpolate=True)
ax1.fill_between(ap_arr, eccTDcoll_arr, eccSScoll_arr, where=(eccTDcoll_arr < eccSScoll_arr),   facecolor='purple',     interpolate=True)
ax1.fill_between(ap_arr, eccSScoll_arr, arr_1, where=(eccSScoll_arr < arr_1),                   facecolor='orange',     interpolate=True)

#legend:
fig.add_subplot(111).plot(-1, -1, linestyle='-', linewidth=5.0, color = 'lightblue',    label='Tidal Inspiral')  
fig.add_subplot(111).plot(-1, -1, linestyle='-', linewidth=5.0, color = 'purple',       label='Tidal Disruption')  
fig.add_subplot(111).plot(-1, -1, linestyle='-', linewidth=5.0, color = 'orange',       label='Solid Sphere Collision')  
plt.legend(loc='lower right', numpoints = 1, fontsize = 10.0, frameon = False)


fig.add_subplot(111).set_ylim(0.925,1.0)
fig.add_subplot(111).set_xlim(1.0,2.0)
            
plt.xlabel(r'$a_{\rm IMS}/a_{\rm c}$')
plt.ylabel(r'$e_{\rm IMS}$')
plt.title(r'Inspirals, Disruptions, Collisions')

plt.savefig('boundary_insp_coll_td.eps', bbox_inches='tight')     
     
plt.show()
exit()
#-----------------------------------------------------------------
















#-----------------------------------------------------------------
#Example 2:
#-----------------------------------------------------------------
Rcoll_yesno = 0
Rtde_yesno  = 1

v_inf_kmsec = 10.

resm    = 100
m1_arr  = np.linspace(1.0, 10.0, resm)
m2_arr  = np.linspace(1.0, 10.0, resm)

insp_over_coll = np.zeros((resm,resm), dtype=np.float64)
insp_over_tde = np.zeros((resm,resm), dtype=np.float64)

m3          = 1.0            #incoming 3
#R_i_Rsun    = 0.013*((1.43/m3)**(1./3.))*((1.-m3/1.43)**(0.447))    #radius of STAR (obj i)
R_i_Rsun    = m3**(0.8)                                             #radius of STAR (obj i)

for mc_1 in range(0,resm):
    for mc_2 in range(0,resm):
        m1      = m1_arr[mc_1]  #bin 1
        m2      = m2_arr[mc_2]  #bin 2
        m_bs    = m1+m2+m3
        
        mi  = m3    #IMS i  STAR
        mj  = m2    #IMS j  COMPACT OBJ
        mk  = m1    #bound single k
        
        f_tid       = 0.5
        
        a_HB_Rsun   = G_new_SI*M_sun_SI*(m1*m2*m_bs/(m3*(m1+m2)))*(1./(v_inf_kmsec*1000.)**(2.))/R_sun_SI                
        a0_Rsun     = a_HB_Rsun    #Rsun
        
        ap          = 1.0+0.5   #plays no role here
        #Tidal Inspirals:
        Es          = 0.1
        x_tide      = 3.5   #T_2(eta) \approx Es*\eta**(-x_tide)
        beta        = 6. + 3.*x_tide/2.
        Ms          = mj*((mi+mj)/mi)**(x_tide/4.)
        A_insp          = func_insp_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_i_Rsun, Es, f_tid, beta, Ms, ap)[1]
        #Solid Sphere (SS) Collisions:
        R_SS_coll_i_Rsun    = R_i_Rsun  
        A_SScoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_SS_coll_i_Rsun, f_tid, ap)[1]
        #Tidal Disruption (TD) Collisions:
        R_TD_coll_i_Rsun    = R_i_Rsun*(mj/mi)**(1./3.)          
        A_TDcoll        = func_coll_A_E(m1,m2,m3, mi,mj,mk, a0_Rsun, R_TD_coll_i_Rsun, f_tid, ap)[1]
                
        
        #save:
        insp_over_coll[mc_1,mc_2] = A_insp/A_SScoll 
        insp_over_tde[mc_1,mc_2] = A_insp/A_TDcoll 
        

Pinsp_over_PSScoll = (1./2.)*insp_over_coll
Pinsp_over_PTDcoll = (1./2.)*insp_over_tde



if (Rcoll_yesno == 1):
    matrix = np.matrix(Pinsp_over_PSScoll)  
    fig     = plt.figure(figsize=(4, 4))
    a       = fig.add_subplot(1,1,1)
    im      = plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.cubehelix)#, extent=(0.5,10.5,0.5,10.5))
    a.set_xlabel(r'$m_{1}$ $[M_{\odot}]$')
    a.set_ylabel(r'$m_{2}$ $[M_{\odot}]$')
    a.set_title(r'$P_{insp}/P_{coll}$')
    #set color bar:
    ax      = plt.gca()
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    #plt.colorbar()#ticks=[0.1,0.5,1.5,2.0])
    plt.savefig('Pinsp_Pcoll_contour_1.eps', bbox_inches='tight')     
    plt.show()
        
if (Rtde_yesno == 1):
    matrix = np.matrix(Pinsp_over_PTDcoll)  
    fig     = plt.figure(figsize=(4, 4))
    a       = fig.add_subplot(1,1,1)
    im      = plt.imshow(matrix, interpolation='nearest', cmap=plt.cm.cubehelix)#, extent=(0.5,10.5,0.5,10.5))
    a.set_xlabel(r'$m_{1}$ $[M_{\odot}]$')
    a.set_ylabel(r'$m_{2}$ $[M_{\odot}]$')
    a.set_title(r'$P_{insp}/P_{td}$')
    #set color bar:
    ax      = plt.gca()
    divider = make_axes_locatable(ax)
    cax     = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    #plt.colorbar()#ticks=[0.1,0.5,1.5,2.0])
    plt.savefig('Pinsp_Ptde_contour_1.eps', bbox_inches='tight')     
    plt.show()


exit()
#-----------------------------------------------------------------























