import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import integrate as integrate
import matplotlib as mpl

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
#------------------------------------------------------------------------------------------

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

tsec_year   = 31536000.

#[1,2]-3 initial config:
m1      = M_sun_SI*30.
m2      = M_sun_SI*30.
m3      = M_sun_SI*30.
#[i,j]-k endstate config:
mi      = m1
mj      = m2
mk      = m3
#GW Inspirals:
f_tid       = 0.5
m_bs        = m1+m2+m3
m_ij        = mi+mj
m_12        = m1+m2
mu_ij       = mi*mj/m_ij
mu_12       = m1*m2/m_12
apu         = 1. + ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.))
Es          = 85.*np.pi/96.
Ms          = mu_ij
beta        = 7./2.
p_Mfac      = (((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
inspI       = ((1.05/(1.-1.7/beta))*(np.log(apu)/((apu-1.)**((1./(beta+1.))))))
#approx values:
PBSij       = 0.3
Nfac        = 5.   #for EM: N \approx 5 (from fit to our first paper)

tau_arr_year    = (10.**(np.linspace(0, 10, num=1000)))
tau_arr_secs    = tsec_year*tau_arr_year
V_disp_1_SI     = 1.*1000.
V_disp_2_SI     = 10.*1000.
V_disp_3_SI     = 100.*1000.

sigmaI_sigmaM_aHB_V1   = (Nfac*inspI*p_Mfac/PBSij)*((3.*(np.pi**2.)/32.)*((G_new_SI*m_bs*mu_12/((V_disp_1_SI**3.)*tau_arr_secs*mu_ij))**2.)*(m_ij/m_12)*(mi*mj/(m3*m3))*(m_bs/(3.*m3)))**(1./7.)
sigmaI_sigmaM_aHB_V2   = (Nfac*inspI*p_Mfac/PBSij)*((3.*(np.pi**2.)/32.)*((G_new_SI*m_bs*mu_12/((V_disp_2_SI**3.)*tau_arr_secs*mu_ij))**2.)*(m_ij/m_12)*(mi*mj/(m3*m3))*(m_bs/(3.*m3)))**(1./7.)
sigmaI_sigmaM_aHB_V3   = (Nfac*inspI*p_Mfac/PBSij)*((3.*(np.pi**2.)/32.)*((G_new_SI*m_bs*mu_12/((V_disp_3_SI**3.)*tau_arr_secs*mu_ij))**2.)*(m_ij/m_12)*(mi*mj/(m3*m3))*(m_bs/(3.*m3)))**(1./7.)


fig = plt.figure(figsize=(5, 4.0))
    
#plot ratio:
fig.add_subplot(111).plot(tau_arr_year, sigmaI_sigmaM_aHB_V1, linestyle='-', linewidth=1.5, label=r'$v_{\infty} = 1$ kms$^{-1}$', color = 'grey')    
fig.add_subplot(111).plot(tau_arr_year, sigmaI_sigmaM_aHB_V2, linestyle='-', linewidth=1.5, label=r'$v_{\infty} = 10$ kms$^{-1}$', color = 'purple')    
fig.add_subplot(111).plot(tau_arr_year, sigmaI_sigmaM_aHB_V3, linestyle='-', linewidth=1.5, label=r'$v_{\infty} = 100$ kms$^{-1}$', color = 'orange')    

#oplot guidelines:
fig.add_subplot(111).plot([1e-10, 1e10], [1.0,1.0], linestyle=':', linewidth=1.0, color = 'black')    

fig.add_subplot(111).set_xscale('log')
fig.add_subplot(111).set_yscale('log')

fig.add_subplot(111).legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = False, ncol=1)

fig.add_subplot(111).set_xlabel(r'$\tau$ [years]')
fig.add_subplot(111).set_ylabel(r'${\sigma_{\rm  I,ij}}/{\sigma_{\rm  M,ij}^{<\tau}}$ $(a_{0} = a_{\rm HB})$')

ax = fig.add_subplot(111)

ax.text(0.1, 0.125, r'Interaction:',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=12)

ax.text(0.1, 0.05, r'[BH(30$M_{\odot}$),BH(30$M_{\odot}$)]-BH(30$M_{\odot}$)',
        verticalalignment='bottom', horizontalalignment='left',
        transform=ax.transAxes,
        color='black', fontsize=10)

fig.add_subplot(111).set_xlim(1e3, 1e10)
fig.add_subplot(111).set_ylim(1e-3, 1e3)

#save and show fig:
plt.savefig('sigmaIsigmaM_aHB.eps', bbox_inches='tight')
    
plt.show()


exit()













