import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import integrate as integrate
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
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


nres                = 200
m_Msun_arr          = (10.**(np.linspace(0.0,  3.0, nres)))
a0_AU_arr           = (10.**(np.linspace(-3.0, 3.0, nres)))

GWinsp_over_PImerg  = np.zeros((nres,nres), dtype=np.float64)

tlimit_yr           = 10.**(10.)    #in years
CIfac               = 0.5

#calc t_insptime:
#tinsp_yr            = (a0_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))/sec_year

#fig = plt.figure(figsize=(5, 4))

for mc in range(0,nres):        #mass m
    for ac in range(0,nres):    #SMA  a
        
        #read m,a:
        m_Msun  = m_Msun_arr[mc]
        a0_AU   = a0_AU_arr[ac]
        #different units:
        m_SI    = m_Msun*M_sun_SI
        a0_SI   = a0_AU*AU_SI
                
        #calc GWinsp cs (AU2):
        cs_AU2_10kms_GWinsp_allp        = 0.025*(m_Msun**(12./7.))*(a0_AU**(2./7.))
                
        #calc CI cs: 
        cs_AU2_10kms_CI_allp            = (CIfac*2.*np.pi*G_new_SI*(3.0*m_SI)*a0_SI/(10000.**2.))/(AU_SI**2.)
         
        #calc info:
        a0_AU_tlimit                    = (((sec_year*tlimit_yr)*(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.)))**(1./4.))/AU_SI
        cs_AU2_10kms_CI_a0_AU_tlimit    = (CIfac*2.*np.pi*G_new_SI*(3.0*m_SI)*(a0_AU_tlimit*AU_SI)/(10000.**2.))/(AU_SI**2.)
        
        if (a0_AU < a0_AU_tlimit):
            cs_AU2_10kms_PImerg_allp    = cs_AU2_10kms_CI_allp
        if (a0_AU > a0_AU_tlimit):
            cs_AU2_10kms_PImerg_allp    = 10.**((-1./7.)*(np.log10(a0_AU) - np.log10(a0_AU_tlimit)) + np.log10(cs_AU2_10kms_CI_a0_AU_tlimit))
        
        #TEST:
        #fig.add_subplot(111).plot(np.log10(a0_AU), np.log10(cs_AU2_10kms_PImerg_allp), marker = '+', color='black')
        #fig.add_subplot(111).plot(np.log10(a0_AU), np.log10(cs_AU2_10kms_GWinsp_allp), marker = '+', color='red')
        #fig.add_subplot(111).plot(np.log10(a0_AU), np.log10(cs_AU2_10kms_GWinsp_allp/cs_AU2_10kms_PImerg_allp), marker = '+', color='green')
        #print cs_AU2_10kms_CI_allp
        #print a0_AU_tlimit, a0_AU
        #print 10.**((-1./7.)*(np.log10(a0_AU) - np.log10(a0_AU_tlimit)) + np.log10(cs_AU2_10kms_CI_a0_AU_tlimit))
        #print 0.75*m_Msun**(13./7.)*a0_AU**(-1./7.)*(tlimit_yr/(10.**10.))**(2./7.)
        #print cs_AU2_10kms_GWinsp_allp/cs_AU2_10kms_CI_allp
        #exit()
            
        GWinsp_over_PImerg[ac,mc]   = cs_AU2_10kms_GWinsp_allp/cs_AU2_10kms_PImerg_allp
        
        
fig     = plt.figure(figsize=(5.5, 4))
#The formating here is quite a mess: It seems like we have to flip the two axis around to get it to work.
#Prepare x,y,z for plotting:

log10_a0_AU_arr             = np.log10(a0_AU_arr) 
log10_m_Msun_arr            = np.log10(m_Msun_arr)
log10_GWinsp_over_PImerg    = np.log10(GWinsp_over_PImerg)

X, Y    = np.meshgrid(log10_a0_AU_arr, log10_m_Msun_arr)
Z       = np.transpose(log10_GWinsp_over_PImerg)
#plot colors:
ax = plt.subplot(111)
ax.set_xlabel(r'log $a_{\rm 0}$ [AU]')
ax.set_ylabel(r'log $m_{\rm}$ $[M_{\odot}]$')
ax.set_title(r'log($\sigma_{\rm I}/\sigma_{\rm M}^{<\tau}$), $\tau$ = 10 Gyrs')
quad = plt.pcolormesh(X, Y, Z, cmap='bone', linewidth=0, rasterized=True, vmin=-3, vmax=0.0)
plt.colorbar()
#oplot contours:
matplotlib.rcParams['contour.negative_linestyle'] = 'solid'
CS = plt.contour(X, Y, Z, list(np.linspace(-2.0, 0.0, num=10+1)), colors='red')
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.1f')

#oplot line:
fbeta           = 2.*64.*(G_new_SI**3.)/(5.*c_SI**5.)
log10_m_SI      = (1./3.)*(4.*np.log10(a0_AU_arr*AU_SI) - np.log10(4.*fbeta*tlimit_yr*sec_year)) 
log10_m_Msun    = np.log10((10.**(log10_m_SI))/M_sun_SI)
ax.plot(log10_a0_AU_arr, log10_m_Msun, linestyle = '--', linewidth = 2, color='yellow')

ax.set_xlim([-3,3])
ax.set_ylim([0,3])

#save and show:
plt.savefig('GWinsp_over_GWmerg.eps', bbox_inches='tight')     
plt.show()

exit()








