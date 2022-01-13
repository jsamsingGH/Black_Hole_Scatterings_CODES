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
from scipy.integrate import solve_ivp
from matplotlib import rcParams


#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI        = 299792458.0       #m/s
M_sun_SI    = 1.989*(10.**30.)  #kg
R_sun_SI    = 695800000.        #m
AU_SI       = 149597871000.     #m 
G_new_SI    = 6.67*(10.**(-11.))
AU_U        = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U     = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
sec_year    = 31536000.
m_parsec    = 3.086*(10**16.)   #m
#------------------------------------------------------------------------------------------

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'



Bc      = (9.*np.pi/((21.**5)*85.))*((20./63.)**(7./2.))
fed     = 2.*np.sqrt(3.)
delta   = 7./9.
tH_sec  = (13.7*(10.**9.))*sec_year  
print tH_sec*((((4.*np.pi/63.)**(7./2.))*((G_new_SI**4.)*(fed**3.)/(Bc*(c_SI**5.))))**(2./7.))*(((10.**5)/(m_parsec**3.))**(5./7.))*((1.4*M_sun_SI)**(8./7.))*(10000.**(1./7.))
print tH_sec/((delta*(fed**2.)/(np.pi*(G_new_SI**2.)*((1.-delta)**2.)))*((10000.**3.)/(((10.**5)/(m_parsec**3.))*((1.4*M_sun_SI)**2.))))
exit()



#--------------------------------------------------------------
#N2/N1, m fig:
#--------------------------------------------------------------

v_kms   = 10.0                          #DISPERSION vel
n_arr   = np.array([10**4., 10**5])        #per 1/pc3
f_arr   = np.array([0.01, 0.1])           #BBH fraction
m_arr   = np.linspace(1.0, 50.0, 1000)  #mass
nrn     = len(n_arr)
nrf     = len(f_arr)
nrm     = len(m_arr)

N2oN1_nfm   = np.zeros((nrn,nrf,nrm), dtype='d')
info_nfm    = np.zeros((nrn,nrf,nrm), dtype='d')

for nc in range(0,nrn):
    for fc in range(0,nrf):
        for mc in range(0,nrm):
        
            
            #--------------------------------------
            #INPUT:
            #--------------------------------------            
            n_nr_pc3    = n_arr[nc]     #per 1/pc3
            fbin        = f_arr[fc]     #0.01              #BBH fraction
            m_Msun      = m_arr[mc]     #30.
            #define/calc:
            fed         = 2.*np.sqrt(3.)    # = v_e/v_d
            vesc_kms    = v_kms*fed         #ESCAPE     vel
            tH_SI       = (10.**10.)*sec_year
            Nims        = 20
            delta_bs    = 7./9.
            v_SI        = v_kms*1000.
            vesc_SI     = vesc_kms*1000.
            m_SI        = m_Msun*M_sun_SI
            n_nr_SI     = n_nr_pc3/(m_parsec**3.)
            #calc a_HB:
            a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(v_SI**2.))
            #calc a_ej:
            a_ej_SI     = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(vesc_SI**2.))
            #--------------------------------------
            #--------------------------------------
            #CALC:
            #--------------------------------------
            #ASSUME a_m = a_ej:
            a_max       = a_HB_SI
            a_min       = a_ej_SI
            #cycle time Tc:
            Tc_SI   = ((((6.*np.pi*G_new_SI)**(-1.))/(1.-delta_bs))*(v_SI/n_nr_SI)*(1./m_SI)*((1./a_min) - (1./a_max)))    
            #Prob merger:
            tint    = (1./(6.*np.pi*G_new_SI))*(v_SI/(n_nr_SI*m_SI*a_min))
            tcinsp  = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((a_min**4.)/(m_SI**3.))
            Pm      = (tint/tcinsp)**(2./7.)
            Pmtot   = Pm*((7./10.)/(1.-delta_bs))
            #input nr single ejections per bin:
            Ns_ej   = 4.0
            Nt_ej   = 2.0 + Ns_ej
            #pi_2 = 0 (obj "2" dont interact):
            Nc_tH       = tH_SI/Tc_SI
            alpha_as    = (Nt_ej - Pmtot*Ns_ej)
            beta_as     = Pmtot
            N2oN1_as    = (beta_as/alpha_as)*(np.exp(alpha_as*fbin*Nc_tH) - 1.0)
            N2oN1_nfm[nc, fc, mc]   = N2oN1_as
            #info array:
            if (Nc_tH > 1.0):   inf1 =  1
            if (Nc_tH < 1.0):   inf1 = -1
            info_nfm[nc, fc, mc]    = inf1
            #--------------------------------------


#--------------------------------------
#PLOT/analyze:
#--------------------------------------
#define:
fig = plt.figure(figsize=(4.0, 7.0))
ax1  = fig.add_subplot(211)
ax2  = fig.add_subplot(212)

#v_kms   = 10.0                          #DISPERSION vel
#n_arr   = np.array([10**4., 10**6])        #per 1/pc3
#f_arr   = np.array([0.01, 0.1])           #BBH fraction
#m_arr   = np.linspace(1.0, 50.0, 1000)  #mass

#PLOT for n1:
#fb1:
nc  = 0
fc  = 0
info_arr    = info_nfm[nc,fc,:]
N2oN1_arr   = N2oN1_nfm[nc,fc,:]
pos    = np.where(info_arr[:] ==  1)[0] #Nc > 1
ax1.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}>1, f_{b} = $' + str(f_arr[fc]),       color = 'steelblue',     linestyle = '-', linewidth = 3.0, zorder=99)
pos    = np.where(info_arr[:] == -1)[0] #Nc < 1
ax1.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}<1, f_{b} = $' + str(f_arr[fc]),       color = 'steelblue',     linestyle = ':', linewidth = 3.0, zorder=99)
#fb2:
nc  = 0
fc  = 1
info_arr    = info_nfm[nc,fc,:]
N2oN1_arr   = N2oN1_nfm[nc,fc,:]
pos    = np.where(info_arr[:] ==  1)[0] #Nc > 1
ax1.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}>1, f_{b} = $' + str(f_arr[fc]),       color = 'seagreen',     linestyle = '-', linewidth = 3.0, zorder=99)
pos    = np.where(info_arr[:] == -1)[0] #Nc < 1
ax1.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}<1, f_{b} = $' + str(f_arr[fc]),       color = 'seagreen',     linestyle = ':',  linewidth = 3.0, zorder=99)
ax1.set_xscale("log")
ax1.set_yscale("log")
ax1.set_xlim(min(m_arr), max(m_arr))
ax1.set_ylim(1e-3,1.0)
#guide lines:
#to LMG:
ax1.plot([1.4, 1.4],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax1.plot([2.5, 2.5],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax1.axvspan(1.4, 2.5, alpha=0.5, color='grey')
#the LMG
ax1.plot([3.0, 3.0],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax1.plot([5.0, 5.0],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax1.axvspan(3.0, 5.0, alpha=0.5, color='red')
ax1.text(3.1, 0.08, 'LMG', ha='left')
#arrow:
af   = 0.003
ax1.arrow(2.0, 0.06, 1.5, 0, width = af, head_width = af*3, color = 'black', head_length = af*65, zorder=100)
#to UMG:
ax1.plot([25.0, 25.0],  [0.0, 1.0],          color = 'black',     linestyle = '--')
ax1.axvspan(25.0, 50., alpha=0.5, color='grey')
#labels etc:
ax1.legend(loc='lower right', numpoints = 1, fontsize = 7.5, frameon = True, facecolor = 'white', ncol=1, framealpha = 1.0)
ax1.text(1.15, 0.55, r'$v_{d} = 10\ kms^{-1}, log(n) = $' + str(np.log10(n_arr[nc])), ha='left', fontsize = 9, bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'grey', 'pad': 2})
ax1.set_xlabel(r'$m$ [$M_\odot$]')
ax1.set_ylabel(r'$max(N_2/N_1),\ t = t_{H}$')



#PLOT for n2:
#fb1:
nc  = 1
fc  = 0
info_arr    = info_nfm[nc,fc,:]
N2oN1_arr   = N2oN1_nfm[nc,fc,:]
pos    = np.where(info_arr[:] ==  1)[0] #Nc > 1
ax2.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}>1, f_{b} = $' + str(f_arr[fc]),       color = 'steelblue',     linestyle = '-', linewidth = 3.0, zorder=99)
pos    = np.where(info_arr[:] == -1)[0] #Nc < 1
ax2.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}<1, f_{b} = $' + str(f_arr[fc]),       color = 'steelblue',     linestyle = ':', linewidth = 3.0, zorder=99)
#fb2:
nc  = 1
fc  = 1
info_arr    = info_nfm[nc,fc,:]
N2oN1_arr   = N2oN1_nfm[nc,fc,:]
pos    = np.where(info_arr[:] ==  1)[0] #Nc > 1
ax2.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}>1, f_{b} = $' + str(f_arr[fc]),       color = 'seagreen',     linestyle = '-', linewidth = 3.0, zorder=99)
pos    = np.where(info_arr[:] == -1)[0] #Nc < 1
ax2.plot(m_arr[pos], N2oN1_arr[pos],           label = r'$N_{c}<1, f_{b} = $' + str(f_arr[fc]),       color = 'seagreen',     linestyle = ':',  linewidth = 3.0, zorder=99)
ax2.set_xscale("log")
ax2.set_yscale("log")
ax2.set_xlim(min(m_arr), max(m_arr))
ax2.set_ylim(1e-3,1.0)
#guide lines:
#to LMG:
ax2.plot([1.4, 1.4],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax2.plot([2.5, 2.5],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax2.axvspan(1.4, 2.5, alpha=0.5, color='grey')
#the LMG
ax2.plot([3.0, 3.0],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax2.plot([5.0, 5.0],    [0.0, 1.0],          color = 'black',     linestyle = '--')
ax2.axvspan(3.0, 5.0, alpha=0.5, color='red')
ax2.text(3.1, 0.08, 'LMG', ha='left')
#arrow:
af   = 0.003
ax2.arrow(2.0, 0.06, 1.5, 0, width = af, head_width = af*3, color = 'black', head_length = af*65, zorder=100)
#to UMG:
ax2.plot([25.0, 25.0],  [0.0, 1.0],          color = 'black',     linestyle = '--')
ax2.axvspan(25.0, 50., alpha=0.5, color='grey')
#labels etc:
ax2.legend(loc='lower right', numpoints = 1, fontsize = 7.5, frameon = True, facecolor = 'white', ncol=1, framealpha = 1.0)
ax2.text(1.15, 0.55, r'$v_{d} = 10\ kms^{-1}, log(n) = $' + str(np.log10(n_arr[nc])), ha='left', fontsize = 9, bbox={'facecolor': 'white', 'alpha': 1, 'edgecolor': 'grey', 'pad': 2}, zorder = 100)
ax2.set_xlabel(r'$m$ [$M_\odot$]')
ax2.set_ylabel(r'$max(N_2/N_1),\ t = t_{H}$')





#PLOT for n2:




#save and plot:
plt.savefig('N2N1marr' + str(int(10*m_Msun)) + '.pdf', bbox_inches='tight')
plt.show()
#--------------------------------------
exit()
#--------------------------------------------------------------
#--------------------------------------------------------------









#--------------------------------------------------------------
#Overview figure:
#--------------------------------------------------------------

#--------------------------------------------------------------
#SETTINGS:
#--------------------------------------------------------------
#System:
m_Msun      = 30.               #MASS in M_sun

#model params/values:
fed         = 2.*np.sqrt(3.)    # = v_e/v_d
Nims        = 20
delta_bs    = 7./9.
tH_SI       = (10.**10.)*sec_year

#binnings etc:
nrb_vnfig       = 100       #nr bins.
vdis_mm         = [2,50]    #vdis min, max
logn_mm         = [2,8]     #logn min, max

#define:
vdis_kms_arr    = np.linspace(vdis_mm[0], vdis_mm[1], nrb_vnfig)
nden_pc3_arr    = 10.**(np.linspace(logn_mm[0],logn_mm[1], nrb_vnfig))
vn_info_arr     = np.zeros((nrb_vnfig,nrb_vnfig,10), dtype='d')


for vc in range(0,nrb_vnfig):
    for dc in range(0,nrb_vnfig):

        
        #--------------------------------------
        #INPUT:
        #--------------------------------------
        #input:
        v_kms       = vdis_kms_arr[vc]    #DISPERSION vel
        n_nr_pc3    = nden_pc3_arr[dc]    #n per 1/pc3
        #calc/convert:
        vesc_kms    = fed*v_kms             #v_e
        v_SI        = v_kms*1000.
        vesc_SI     = vesc_kms*1000.
        m_SI        = m_Msun*M_sun_SI
        n_nr_SI     = n_nr_pc3/(m_parsec**3.)
        #--------------------------------------
        #CALC/DEFINE:
        #--------------------------------------
        #a_HB:
        a_HB_SI         = (3./2.)*(G_new_SI*m_SI/(v_SI**2.))
        tint_HB_SI      = (1./(6.*np.pi*G_new_SI))*(v_SI/(n_nr_SI*m_SI*a_HB_SI))
        
        #a_ej:
        a_ej_SI         = (1./6.)*((1./delta_bs) - 1.)*(G_new_SI*m_SI/(vesc_SI**2.))
        
        #a_GW:
        tcinsp_HB_SI    = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((a_HB_SI**4.)/(m_SI**3.))
        Pm_ai           = (tint_HB_SI/tcinsp_HB_SI)**(2./7.)    #P for merging at ai = a_HB
        a_GW_SI         = a_HB_SI*((10./7.)*((1.-delta_bs)/Pm_ai) + 1.0)**(-7./10.) 
        
        #a_tH:
        a_tH_SI         = a_HB_SI*((1.-delta_bs)*(tH_SI/tint_HB_SI) + 1.0)**(-1.)
        
        #a limit:
        alim_SI         = max([a_ej_SI,a_GW_SI,a_tH_SI]) 
        
        #Nr burning cycles in tH:
        tlim_SI         = tint_HB_SI*(((a_HB_SI/alim_SI) - 1.0)/(1.-delta_bs))
        N_bc            = tH_SI/tlim_SI
        
        #P2 merging before alim (in between 3-body int):
        Pfac        = (85./(9*np.pi))*((G_new_SI**2.)/(c_SI**5.))*((m_SI**2.)*v_SI/n_nr_SI)
        P2_alim     = (1./(1.-delta_bs))*(7./10.)*((Pfac**(2./7.))*(alim_SI**(-10./7.)))

        #P3 merging before alim (during 3-body int):
        Rm_SI       = (2.*G_new_SI*m_SI/(c_SI**2.))
        P3_alim     = ((1./(1.-delta_bs))*(7./5.))*(2.*Nims*((Rm_SI/alim_SI)**(5./7.)))
        
        #P23 - total prob for merging from 3-body int:
        P23_alim    = P2_alim + P3_alim
        
        #bs and ss merger rates:        
        #bs mergers:
        Rate_23_SI  = P23_alim/tlim_SI
        #ss mergers:
        Rc_ss       = (((85.*np.pi)/(24.*np.sqrt(2.)))**(2./7.))*Rm_SI*((c_SI/v_SI)**(4./7.))
        sigma_ss    = ((2.*np.pi*G_new_SI)*(2.*m_SI/(v_SI**2.)))*Rc_ss
        Rate_ss_SI  = n_nr_SI*sigma_ss*v_SI
        
        #Other params:
        Nbc_P23     = N_bc*P23_alim
        
        #time evolving cluster model:
        Ns_ej       = 4.0
        Nt_ej       = 2.0 + Ns_ej
        fb_evm      = 0.01
        alpha_evm   = (Nt_ej - P2_alim*Ns_ej)
        beta_evm    = P2_alim
        N1oN10_tH   = np.exp(-alpha_evm*fb_evm*N_bc)
        N2oN1_tH    = (beta_evm/alpha_evm)*(np.exp(alpha_evm*fb_evm*N_bc) - 1.0)                 
        #--------------------------------------
        #save:
        #--------------------------------------        
        vn_info_arr[vc,dc,0]                            = N_bc                      #0
        if (alim_SI == a_ej_SI): vn_info_arr[vc,dc,1]   = 0                         #1
        if (alim_SI == a_GW_SI): vn_info_arr[vc,dc,1]   = 1                         #1
        if (alim_SI == a_tH_SI): vn_info_arr[vc,dc,1]   = 2                         #1
        vn_info_arr[vc,dc,2]                            = P23_alim                  #2
        vn_info_arr[vc,dc,3]                            = Nbc_P23                   #3
        vn_info_arr[vc,dc,4]                            = P3_alim/P23_alim          #4
        vn_info_arr[vc,dc,5]                            = Rate_ss_SI/Rate_23_SI     #5
        vn_info_arr[vc,dc,6]                            = P3_alim/P2_alim           #6
        vn_info_arr[vc,dc,7]                            = N2oN1_tH                  #7
        vn_info_arr[vc,dc,8]                            = N1oN10_tH                 #8
        #--------------------------------------



#--------------------------------------
#PLOT:
#--------------------------------------
mpl.rcParams['hatch.linewidth'] = 1

#--------------------------------------
#FIG 1:
#--------------------------------------
fig = plt.figure(figsize=(5.0, 4.0))
ax  = fig.add_subplot(111)

X, Y    = np.meshgrid(vdis_kms_arr, np.log10(nden_pc3_arr))

#oplot regions:
Z           = np.transpose(vn_info_arr[:,:,1])
cMap        = ccolor.ListedColormap(['dodgerblue', 'crimson', 'grey', 'green'])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
quad        = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=3.0)

#oplot contours:
Z           = np.transpose(vn_info_arr[:,:,5])
Fbs         = 0.05  #set binary fraction cut.
CS          = plt.contour(X, Y, Z, list([2.*Fbs]), colors='black', alpha=1.0, linestyles='dotted')
[pi,pj]     = np.where(Z > 2.*Fbs)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '..'],  alpha=0)

Z           = np.transpose(vn_info_arr[:,:,6])
CS          = plt.contour(X, Y, Z, list([1]), colors='black', alpha=1.0, linestyles='dotted')
[pi,pj]     = np.where(Z > 1)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '\\\\'], alpha=0)

Z           = np.transpose(vn_info_arr[:,:,0])
CS          = plt.contour(X, Y, Z, list([2,5,10,25, 50, 100, 250, 500, 1000]), colors='yellow', alpha=0.75, linestyles='dashed')
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.1f')

Z           = np.transpose(vn_info_arr[:,:,2])
CS          = plt.contour(X, Y, Z, list([0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9]), colors='crimson', alpha=1.0, linestyles='solid')
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.2f')

#oplot lines:
vdis_parr_kms   = np.linspace(1, 100, 1000)
vdis_parr_SI    = vdis_parr_kms*1000.
Dfac            = (9./(21**5.))*(np.pi/85.)*((20./63.)**(7./2.))
Efac            = ((20./63.)**(-7./2.))*(85./(9.*np.pi))*(((4./3.)*np.pi)**5.)
Kfac            = (63./5.)*(42.**(5./7.))
n_pc3_ejGW      = ((1./Dfac)*((fed**10.)/((G_new_SI**3.)*(c_SI**5.)))*((vdis_parr_SI**(11.))/(m_SI**3.)))*(m_parsec**3.)
n_pc3_ejtH      = ((63./(4.*np.pi))*((fed**2)/((G_new_SI**2.)*tH_SI))*((vdis_parr_SI**3.)/(m_SI**2.)))*(m_parsec**3.)
n_pc3_GWtH      = (((Efac*(G_new_SI**7.)*(tH_SI**5.)/(c_SI**5.))**(-1./4.))*(m_SI**(-7./4.))*vdis_parr_SI)*(m_parsec**3.)
n_pc3_23m       = ((vdis_parr_SI**(6))*(m_SI**(-3.))*((Kfac**(7./2.))*(G_new_SI**3.)*Dfac*(Nims**(7./2.))*(fed**(-5.)))**(-1.))*(m_parsec**3.)
vkms_BP         = (((63./(4.*np.pi))*((Dfac*G_new_SI*(c_SI**5.))/((fed**8)*(tH_SI)))*m_SI)**(1./8.))/1000.
npc3_BP         = (((((63./(4.*np.pi))**(11.))*(((Dfac**3.)*(c_SI**15.))/((G_new_SI**(13.))*(tH_SI**(11.))*(fed**8.))))**(1./8.))*(m_SI**(-13./8.)))*(m_parsec**3.)
pos_LTBP        = np.where(vdis_parr_kms[:] <= vkms_BP)[0]
pos_GTBP        = np.where(vdis_parr_kms[:] >= vkms_BP)[0]
ax.plot(vdis_parr_kms[pos_GTBP], np.log10(n_pc3_ejGW[pos_GTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vdis_parr_kms[pos_LTBP], np.log10(n_pc3_ejtH[pos_LTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vdis_parr_kms[pos_GTBP], np.log10(n_pc3_GWtH[pos_GTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vkms_BP, np.log10(npc3_BP), marker='o', markersize=6.0, linestyle='', color='lime', markeredgewidth=1, markeredgecolor = 'black', alpha = 0.9)

#axis settings etc:
ax.set_xlim([min(vdis_kms_arr),max(vdis_kms_arr)])
ax.set_ylim([min(np.log10(nden_pc3_arr)),max(np.log10(nden_pc3_arr))])
ax.set_xlabel(r'$v_{\rm d}$ $[kms^{-1}]$', fontsize = 12)
ax.set_ylabel(r'$log(n)$ $[pc^{-3}]$', fontsize = 12)

plt.text(37, 2.5, r'$m = $'+ str(m_Msun) + r'$M_{\odot}$', size=12, ha="left", va="center", bbox=dict(boxstyle="square", ec='lightgrey', fc='lightgrey',), zorder=100)

#save and show:
plt.savefig('alim_fig_' + str(int(10*m_Msun)) + '.pdf', bbox_inches='tight')     
plt.show()
#--------------------------------------

#--------------------------------------
#FIG 2:
#--------------------------------------
fig = plt.figure(figsize=(5.0, 4.0))
ax  = fig.add_subplot(111)

#oplot contours:
X, Y    = np.meshgrid(vdis_kms_arr, np.log10(nden_pc3_arr))

Z           = np.transpose(vn_info_arr[:,:,3])
quad        = plt.pcolormesh(X, Y, np.log10(Z), cmap='cubehelix', linewidth=0, rasterized=True, vmin=-2, vmax=4)

Z           = np.transpose(vn_info_arr[:,:,5])
Fbs         = 0.05  #set binary fraction cut.
CS          = plt.contour(X, Y, Z, list([2.*Fbs]), colors='black', alpha=1.0, linestyles='dotted')
[pi,pj]     = np.where(Z > 2.*Fbs)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '..'],  alpha=0)

Z           = np.transpose(vn_info_arr[:,:,6])
CS          = plt.contour(X, Y, Z, list([1]), colors='black', alpha=1.0, linestyles='dotted')
[pi,pj]     = np.where(Z > 1)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '\\\\'], alpha=0)

Z           = np.transpose(vn_info_arr[:,:,3])
CS          = plt.contour(X, Y, Z, list(10.**np.linspace(-2, -0.1, num=10)), colors='white', linestyles='dotted')
CS          = plt.contour(X, Y, Z, list([1, 2,5,10, 25, 50, 100, 250, 500, 1000]), colors='black', linestyles='solid')
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.0f')

#oplot lines:
ax.plot(vdis_parr_kms[pos_GTBP], np.log10(n_pc3_ejGW[pos_GTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vdis_parr_kms[pos_LTBP], np.log10(n_pc3_ejtH[pos_LTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vdis_parr_kms[pos_GTBP], np.log10(n_pc3_GWtH[pos_GTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vkms_BP, np.log10(npc3_BP), marker='o', markersize=6.0, linestyle='', color='lime', markeredgewidth=1, markeredgecolor = 'black', alpha = 0.9)

#axis settings etc:
ax.set_xlim([min(vdis_kms_arr),max(vdis_kms_arr)])
ax.set_ylim([min(np.log10(nden_pc3_arr)),max(np.log10(nden_pc3_arr))])
ax.set_xlabel(r'$v_{\rm d}$ $[kms^{-1}]$', fontsize = 12)
ax.set_ylabel(r'$log(n)$ $[pc^{-3}]$', fontsize = 12)

plt.text(37, 2.5, r'$m = $'+ str(m_Msun) + r'$M_{\odot}$', size=12, ha="left", va="center", bbox=dict(boxstyle="square", ec='lightgrey', fc='lightgrey',), zorder=100)

#save and show:
plt.savefig('Nmbin_' + str(int(10*m_Msun)) + '.pdf', bbox_inches='tight')     
plt.show()
#--------------------------------------




#--------------------------------------
#FIG 3:
#--------------------------------------
fig = plt.figure(figsize=(5.0, 4.0))
ax  = fig.add_subplot(111)

#oplot contours:
X, Y        = np.meshgrid(vdis_kms_arr, np.log10(nden_pc3_arr))

#highlight Nc<1 area in grey:
Z           = np.transpose(vn_info_arr[:,:,1])
cMap        = ccolor.ListedColormap(['white', 'white', 'grey', 'white'])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
quad        = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=3.0)

#plot: N2/N1
#clean for inf/NAN:
Z_N2oN1         = np.transpose(vn_info_arr[:,:,7])
[pi,pj]         = np.where(np.isinf(Z_N2oN1))
Z_N2oN1[pi,pj]  = 1e10
[pi,pj]         = np.where(np.isnan(Z_N2oN1))
Z_N2oN1[pi,pj]  = 1e10
#mask Nc < 1 area:
Z_Nc        = np.transpose(vn_info_arr[:,:,0])
[pi,pj]     = np.where(Z_Nc < 1.001)
Z_N2oN1_M           = Z_N2oN1
Z_N2oN1_M[pi,pj]    = float('NaN')
#define Z for easy plotting and plot:
Z           = Z_N2oN1_M
quad        = plt.pcolormesh(X, Y, np.log10(Z), cmap='viridis', linewidth=0, rasterized=True, vmin=-2, vmax=0)
CS          = plt.contour(X, Y, Z, list([0.01, 0.025, 0.05, 0.1, 0.25, 1.0]), colors='red', linestyles='solid')
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.2f')
[pi,pj]     = np.where(Z > 1.0)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '\\\\'], alpha=0)
[pi,pj]     = np.where(Z < 0.01)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', '+'], alpha=0)

#plot: N1/N1(0)
Z           = np.transpose(vn_info_arr[:,:,8])
CS          = plt.contour(X, Y, Z, list([1e-4]), colors='black', linestyles='solid')
[pi,pj]     = np.where(Z < 1e-4)
ZC          = np.zeros(Z.shape)
ZC[pi,pj]   = 1
if (len(pi) > 0):   CS          = ax.contourf(X, Y, ZC, 3, hatches=['', 'XX'], alpha=0)

#oplot lines:
ax.plot(vdis_parr_kms[pos_GTBP], np.log10(n_pc3_ejGW[pos_GTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vdis_parr_kms[pos_LTBP], np.log10(n_pc3_ejtH[pos_LTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vdis_parr_kms[pos_GTBP], np.log10(n_pc3_GWtH[pos_GTBP]), linestyle = '--', linewidth = 1.5, color='lime')
ax.plot(vkms_BP, np.log10(npc3_BP), marker='o', markersize=6.0, linestyle='', color='lime', markeredgewidth=1, markeredgecolor = 'black', alpha = 0.9)

#axis settings etc:
ax.set_xlim([min(vdis_kms_arr),max(vdis_kms_arr)])
ax.set_ylim([min(np.log10(nden_pc3_arr)),max(np.log10(nden_pc3_arr))])
ax.set_xlabel(r'$v_{\rm d}$ $[kms^{-1}]$', fontsize = 12)
ax.set_ylabel(r'$log(n)$ $[pc^{-3}]$', fontsize = 12)

plt.text(37, 2.5, r'$m = $'+ str(m_Msun) + r'$M_{\odot}$', size=12, ha="left", va="center", bbox=dict(boxstyle="square", ec='lightgrey', fc='lightgrey',), zorder=100)

#save and show:
plt.savefig('N2oN1_' + str(int(10*m_Msun)) + '.pdf', bbox_inches='tight')     
plt.show()
#--------------------------------------
exit()
#--------------------------------------------------------------
#--------------------------------------------------------------









#--------------------------------------------------------------
#N_C >> 1 model:
#--------------------------------------------------------------

#--------------------------------------
#INPUT:
#--------------------------------------
#input/edit:
N1_t0       = 1.0               #nr m1 at t0
N2_t0       = 0.0               #nr m2 at t0
incl_N2_int = 1                 #incl N2 interactions (pi_2 > 0). (#yes=1, no=0)

#fix bin frac:
fix_binfrac = 1 #yes=1, no=0
fbin        = 0.01              #BBH fraction
#fix bin nr:
fix_binnr   = 0 #yes=1, no=0
Nb_t0       = 5.0               #nr bin at t0

#cluster + m:
m_Msun      = 30.
v_kms       = 10.0              #DISPERSION vel
n_nr_pc3    = (10.**(4.))       #per 1/pc3

#define/calc:
fed         = 2.*np.sqrt(3.)    # = v_e/v_d
vesc_kms    = v_kms*fed         #ESCAPE     vel
tH_SI       = (10.**10.)*sec_year
Nims        = 20
delta_bs    = 7./9.
v_SI        = v_kms*1000.
vesc_SI     = vesc_kms*1000.
m_SI        = m_Msun*M_sun_SI
n_nr_SI     = n_nr_pc3/(m_parsec**3.)

#calc a_HB:
a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(v_SI**2.))
#calc a_ej:
a_ej_SI     = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(vesc_SI**2.))
#--------------------------------------

#--------------------------------------
#CALC/DEFINE:
#--------------------------------------
#ASSUME a_m = a_ej:
a_max       = a_HB_SI
a_min       = a_ej_SI
#cycle time Tc:
Tc_SI   = ((((6.*np.pi*G_new_SI)**(-1.))/(1.-delta_bs))*(v_SI/n_nr_SI)*(1./m_SI)*((1./a_min) - (1./a_max)))    
#Prob merger:
tint    = (1./(6.*np.pi*G_new_SI))*(v_SI/(n_nr_SI*m_SI*a_min))
tcinsp  = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((a_min**4.)/(m_SI**3.))
Pm      = (tint/tcinsp)**(2./7.)
Pmtot   = Pm*((7./10.)/(1.-delta_bs))
#input nr single ejections per bin:
Ns_ej   = 4.0
Nt_ej   = 2.0 + Ns_ej
#define:
f21     = 1.0
w21     = 1.0
Rm_11   = 1.0
pes_211 = (2./3.)*w21
pes_112 = 1.0-pes_211
B_fac   = (2.*f21*w21)/(3.-2.*w21) 
A_fac   = (1.+B_fac)*(pes_112*Ns_ej + pes_211)
C_fac   = A_fac - Pmtot*(A_fac - B_fac)
#--------------------------------------

#--------------------------------------
#EVOLVE SYSTEM:
#--------------------------------------
Nc_min  = 0.0
Nc_max  = 1000
Nc_eval = np.linspace(Nc_min, Nc_max, 1000)

def rhs(t, v):
    
    #nr N1, N2:
    N1      = v[0]
    N2      = v[1]
    N12     = N1+N2
    
    #N2 interaction probability:
    if (incl_N2_int == 1):  pi_2 = (N2/N12)*f21
    if (incl_N2_int == 0):  pi_2 = 0.0
    
    #define:
    pb_21   = pi_2*B_fac
    pb_11   = 1.0 - pb_21
    
    #nr binaries:
    if (fix_binfrac == 1):  Nb  = fbin*N1       #keep tot frac of bin constant
    if (fix_binnr   == 1):  Nb  = Nb_t0         #keep tot nr of bin constant
    
    #evol eq:
    dN1 = Nb*(-(Nt_ej - Pmtot*Ns_ej) + pi_2*C_fac)
    dN2 = Nb*((pb_11*Pmtot*Rm_11)    - pi_2*C_fac)
    return [dN1, dN2]

#evolve:
res = solve_ivp(rhs, (Nc_min, Nc_max), [N1_t0, N2_t0], t_eval = Nc_eval)
#--------------------------------------

#--------------------------------------
#PLOT/analyze:
#--------------------------------------
#define:
fig = plt.figure(figsize=(5.0, 3.0))
ax  = fig.add_subplot(111)

def Nc2ttH(x):
    return x*(Tc_SI/tH_SI)
def ttH2Nc(x):
    return x*(tH_SI/Tc_SI)

t_tH_eval   = Nc2ttH(Nc_eval) 

#from evolv code:
N1_ev       = res.y[0, :]  
N2_ev       = res.y[1, :]
N2oN1_ev    = N2_ev/N1_ev
ax.plot(t_tH_eval, N1_ev,              label = r'(A): $N_1/N_1(0)$',         color = 'black',   linestyle = '--', zorder=100)
ax.plot(t_tH_eval, N2_ev,              label = r'(A): $N_2/N_1(0)$',         color = 'red',     linestyle = '--')
ax.plot(t_tH_eval, N2oN1_ev,           label = r'(A): $N_2/N_1$',            color = 'blue',    linestyle = '--')

#pi_2 = 0 (obj "2" dont interact):
alpha_as    = (Nt_ej - Pmtot*Ns_ej)
beta_as     = Pmtot
N1_as       = N1_t0*np.exp(-alpha_as*fbin*Nc_eval)
N2_as       = N1_t0*(beta_as/alpha_as)*(1.0 - np.exp(-alpha_as*fbin*Nc_eval))
N2_inf      = N1_t0*(beta_as/alpha_as)
N2oN1_as    = (beta_as/alpha_as)*(np.exp(alpha_as*fbin*Nc_eval) - 1.0)
ax.plot(t_tH_eval, N1_as,              label = r'(B): $N_1/N_1(0)$',    color = 'black',    linestyle = '-', zorder=100)
ax.plot(t_tH_eval, N2_as,              label = r'(B): $N_2/N_1(0)$',    color = 'red',      linestyle = '-')
ax.plot(t_tH_eval, N2oN1_as,           label = r'(B): $N_2/N_1$',       color = 'blue',     linestyle = '-')
#ax.plot([1e-3, 1e1], [N2_inf,N2_inf], linestyle = ':', color = 'black')

#fill inbetween lines:
ax.fill_between(t_tH_eval, N1_ev,       N1_as,      color = 'black',    alpha=0.25)
ax.fill_between(t_tH_eval, N2_ev,       N2_as,      color = 'red',      alpha=0.25)
ax.fill_between(t_tH_eval, N2oN1_ev,    N2oN1_as,   color = 'blue',     alpha=0.25)

#oplot guidelines:
ax.plot([1,1], [1e-10, 1e10], linestyle = ':', color = 'black', alpha = 0.25)

#legend/labels:
plt.legend(loc='center left', numpoints = 1, fontsize = 8, frameon = False, ncol=1)
plt.xlabel(r'time $(t/t_{\rm H})$')
plt.ylabel(r'$N_1, N_2$ evolution')
plt.text(0.1, 10**(-3.65), r'$m = $'+ str(m_Msun) + r'$M_{\odot}$' + ',' + r'$\ log(n/pc^{3}) =$' + str(np.log10(n_nr_pc3)), size=10, ha="left", va="center", bbox=dict(boxstyle="square", ec='black', fc='lightgrey',), zorder=100)

#axes settings and limites:
ax.set_xscale('log')
ax.set_yscale('log')
ax_xmm = [5e-3,3] 
ax_ymm = [1e-4,2] 
ax.set_xlim(ax_xmm)
ax.set_ylim(ax_ymm)

#create top x-axes:
ax2 = ax.twiny()
ax2.set_xlabel(r'nr. interaction cycles ($N_c$)')
ax2.plot(Nc_eval, N1_ev, color = '', linestyle = '')
ax2.set_xscale('log')
ax2.set_yscale('log')
ax2_xmm = ttH2Nc(np.array(ax_xmm))
ax2_ymm = ax_ymm
ax2.set_xlim(ax2_xmm)
ax2.set_ylim(ax2_ymm)

#save and plot:
plt.savefig('Nevol_' + str(int(10*m_Msun)) + '.pdf', bbox_inches='tight')
plt.show()
#--------------------------------------
exit()
#--------------------------------------------------------------
#--------------------------------------------------------------































