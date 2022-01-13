import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import integrate as integrate
import matplotlib as mpl
from scipy.optimize import fsolve


#----------------------------------------------------------
#Units and conversions:
#----------------------------------------------------------
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
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#figure settings:
fignr = 2
fig, ax1 = plt.subplots(figsize=(5, 4))    

#input settings:
m_Msun_arr  = [30.0, 20.0, 10.0, 5.0]
m_color_arr = ['#1f78b4', '#33a02c', '#e31a1c', '#ff7f00']
nr_m        = len(m_Msun_arr)

for mc in range(0,nr_m):

    #----------------------------------------------------------
    #input/define:
    #----------------------------------------------------------
    #input:
    vdis_kms    = 12.0      #DISPERSION vel
    m_Msun      = m_Msun_arr[mc]
    n_nr_SI     = (10.**(4.))/(m_parsec**3.)
    #define:
    vdis_SI     = vdis_kms*1000.
    m_SI        = m_Msun*M_sun_SI
    delta_bs    = 7./9.
    Eps         = delta_bs**(-10./7.)
    #----------------------------------------------------------
    #----------------------------------------------------------
    #calc:
    #----------------------------------------------------------
    #a_HB:
    a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(vdis_SI**2.))
    #prob for merger at HB:
    tbs_HB_SI   = 1./(n_nr_SI*(2.*np.pi*G_new_SI*(3.*m_SI)*a_HB_SI/(vdis_SI**2.))*vdis_SI)
    tcl_HB_SI   = (a_HB_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
    Pm_HB       = (tbs_HB_SI/tcl_HB_SI)**(2./7.)
    #T orb at HB:
    THB_yrs     = (2.*np.pi*np.sqrt(a_HB_SI**3./(G_new_SI*(m_SI + m_SI))))/sec_year
    THB_dys     = 365.*THB_yrs
    #----------------------------------------------------------
    
    
    
    #----------------------------------------------------------
    #TEST:
    #----------------------------------------------------------
    #print np.log10(THB_yrs*(((30./21.)*Pm_HB*Eps/(Eps-1.))**(21./20.)))
    #print (THB_yrs*(((30./21.)*Pm_HB*Eps/(Eps-1.))**(21./20.)))*365.
    #print (THB_yrs*(2./63.)**(3./2.)*(5.0**(-3.)))*365.
    print np.log10((1./Pm_HB))/np.log10(Eps)
    nmax    = -np.log10(Pm_HB)/np.log10(Eps)
    amin    = (delta_bs**(nmax))*a_HB_SI
    Tmin    = (2.*np.pi*np.sqrt(amin**3./(G_new_SI*(m_SI + m_SI))))
    print Tmin
    
    print ((((512.*(32.**(1./3.)))*(np.pi**(7./3.))*(G_new_SI**(1./3.))/(30.*c_SI**5.))*(m_SI**(1./3.)*vdis_SI/n_nr_SI))**(3./10.))
    
    Ccont   = 512.*G_new_SI**2./(30.*np.pi*c_SI**5.)
    print ((Ccont*m_SI**2.*vdis_SI/n_nr_SI)**(3./5.)*(2.*np.pi**2./(G_new_SI*m_SI)))**(1./2.)
    
    print (2.*np.pi*np.sqrt(a_HB_SI**3./(G_new_SI*(m_SI + m_SI))))*(Pm_HB)**(21./20.)
    
    #print time limits:
    print '1', THB_dys*(Pm_HB)**(21./20.) 
    print '2', THB_dys*(5.0)**(-3.0)*((1.-delta_bs)/(9.*delta_bs))**(3./2.)    
    print '3', THB_dys*(tbs_HB_SI/(sec_year*(10.**10)))**(3./2.)*(delta_bs/(1.0 - delta_bs))**(3./2.)
    print '4', THB_dys*(Pm_HB)**(21./20.)*((10./7.)*(Eps)/(Eps - 1.0))**(21./20.)
    print '5', np.log10(THB_dys*(Pm_HB)**(21./20.)*((10./7.)*(Eps)/(Eps - 1.0))**(21./20.))

    #exit()
    
    
    print ((10./7.)*Eps/(Eps-1.0))**(21./20.)
    print nmax
    
    Ms      = M_sun_SI
    Rs      = R_sun_SI
    Rt      = Rs*(m_SI/Ms)**(1./3.)
    n_s     = (10.**(4.))/(m_parsec**3.)
    rate_TD = (n_nr_SI*(2.*np.pi*G_new_SI*(2.*m_SI)*(2.*Rt)/(vdis_SI**2.))*vdis_SI)*sec_year
    print rate_TD*5.*200.
    
    #exit()
    
    vesc    = 50000.0 #m/s
    aej     = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(vesc**2.))
    tbs_ej  = 1./(n_nr_SI*(2.*np.pi*G_new_SI*(3.*m_SI)*aej/(vdis_SI**2.))*vdis_SI)
    print (n_s*(2.*np.pi*G_new_SI*(2.*m_SI)*(2.*Rt)/(vdis_SI**2.))*vdis_SI)*tbs_ej 
    print aej/R_sun_SI
    print aej/AU_SI
    print a_HB_SI/AU_SI
    print 365.*((2.*np.pi*np.sqrt(aej**3./(G_new_SI*(m_SI + m_SI))))/sec_year)
    print (2./((2.*np.pi*np.sqrt(aej**3./(G_new_SI*(m_SI + m_SI))))))*(1.-0.99)**(-3./2.)
    print 2.*Rt/(0.1*AU_SI)
    print ((2.*Rt)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))/sec_year

    adyn    = 0.5   #AU
    t3b     = 0.1   #year
    t2b     = 1e7   #year
    tH      = 1e10  #year
    
    print (7./5.)*(1./(1. - delta_bs))*20.*((t3b*sec_year)/(((adyn*AU_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))))**(2./7.)
    print (7./10.)*(1./(1. - delta_bs))*1.*((t2b*sec_year)/(((adyn*AU_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))))**(2./7.)
    print 1.*((tH*sec_year)/(((adyn*AU_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))))**(2./7.)
    
    rcap    = ((85*np.pi/(6.*np.sqrt(2.)))**(2./7.))*(G_new_SI*(m_SI**(2./7.))*(m_SI**(2./7.))*((m_SI + m_SI)**(3./7.)))/((c_SI**(10./7.))*(vdis_SI**(4./7.)))
    fcap    = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rcap**3.)) 
    print fcap
    
    print ((2.*np.pi*np.sqrt((0.2*AU_SI)**(3.)/(G_new_SI*(20.*M_sun_SI + 20.*M_sun_SI))))/sec_year)*365.
    
    #Aconst  = 1./(6.*np.pi*G_new_SI)
    #Bconst  = (768./425.)*(5.*(c_SI**5.))/(512.*(G_new_SI**3.))
    #Cconst  = Bconst*(2.**(7./2.))
    #Econst  = (1./6.)*(1./delta_bs - 1.)*G_new_SI
    #print (((Aconst/Bconst)*(vdis_SI*(m_SI**2.))/n_nr_SI)**(1./5.))/R_sun_SI
    #print (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(vesc**2.))/R_sun_SI
    #print (((Aconst/Cconst)*(vdis_SI/n_nr_SI)*((m_SI**(1./2.))/(Econst**(3./2.)))*(vesc**(3.)))**(2./7.))/R_sun_SI
    #print 1.0 - (1.-2./50.)**2.
    
    
    #exit()
    #----------------------------------------------------------
    
    
    #----------------------------------------------------------
    #PLOT:
    #----------------------------------------------------------
    #----------------------------------------------------------
    if (fignr == 1):
        #------------------------------------------------------
        #FIG 1:
        #------------------------------------------------------
        arr_TTHB    = 10.**(np.arange(-10.0,0.0,0.01))
        Pnm_TTHB    = np.exp(-Pm_HB*((1.0 - Eps*(arr_TTHB**(-20./21.)))/(1.0 - Eps)))
    
        PTD_GRno    = arr_TTHB**(-2./3.)
        PTD_GRyes   = PTD_GRno*Pnm_TTHB
    
        Pnorm       = 10.0  #rnd norm
        PTD_GRno    = PTD_GRno/Pnorm
        PTD_GRyes   = PTD_GRyes/Pnorm
        
        #lower limits:
        if (mc == 0):
            #lim 1:
            ve  = 50.0      #esc vel in km/s
            fed         = ve/vdis_kms #2.*np.sqrt(3.0)
            Tej_THB     = ((2./63.)**(3./2.))*(fed**(-3.))
            ax1.plot([np.log10(Tej_THB), np.log10(Tej_THB)], [0.0,100.0], linestyle = ':', linewidth = 1, color = 'grey', label = r'lower lim. $v_{\rm e} = 50$ kms$^{-1}$')
            #lim 2:
            ve  = 100.0     #esc vel in km/s
            fed         = ve/vdis_kms
            Tej_THB     = ((2./63.)**(3./2.))*(fed**(-3.))
            ax1.plot([np.log10(Tej_THB), np.log10(Tej_THB)], [0.0,100.0], linestyle = ':', linewidth = 1, color = 'grey', label = r'lower lim. $v_{\rm e} = 100$ kms$^{-1}$')
    
        #no GR:
        if (mc == 0): ax1.plot(np.log10(arr_TTHB), PTD_GRno,    linestyle = '--',   linewidth = 2, color = 'black',         label   = r'Newtonian prediction')
        #with GR:
        if (mc >= 0): ax1.plot(np.log10(arr_TTHB), PTD_GRyes,   linestyle = '-',    linewidth = 2, color = m_color_arr[mc], label   = r'+GR, m = ' + str(m_Msun) + r'$M_{\odot}$' + r', $T_{\rm HB} \approx$' + str(int(THB_yrs)) + ' yrs')
        #ax1.fill_between(np.log10(arr_TTHB), 0, PTD_GRyes, color = m_color_arr[mc], alpha = 0.1)
        #------------------------------------------------------
    #---------------------------------------------------------- 
    if (fignr == 2):
    #----------------------------------------------------------
        #------------------------------------------------------
        #FIG 2:
        #------------------------------------------------------        
        #FULL T scale:
        arr_Torb_dys    = THB_dys*(10.**(np.arange(0.0, -10.0, -0.01)))
        Pnm_Torb_dys    = np.exp(-Pm_HB*((1.0 - Eps*((arr_Torb_dys/THB_dys)**(-20./21.)))/(1.0 - Eps)))
        PTD_GRyes       = (arr_Torb_dys**(-2./3.))*Pnm_Torb_dys
        PTD_GRyes       = PTD_GRyes/PTD_GRyes[0]    #normalize to HB limit value
        PTD_GRno        = (arr_Torb_dys**(-2./3.))*1.0
        PTD_GRno        = PTD_GRno/PTD_GRno[0]      #normalize to HB limit value
        #define norm factor:
        if (mc == 0): normfac = max(PTD_GRyes)
        #plot:
        ax1.plot(np.log10(arr_Torb_dys), PTD_GRyes/normfac, linestyle = '-', linewidth = 2, color = m_color_arr[mc], label = r'$m_{\rm BH} = $' + str(m_Msun) + r'$M_{\odot}$')
        ax1.plot(np.log10(arr_Torb_dys), PTD_GRno/normfac,  linestyle = ':', linewidth = 1, color = m_color_arr[mc])#, label = r'$m_{\rm BH} = $' + str(m_Msun) + r'$M_{\odot}$ (Newt.)')
        #plot limits:
        fed = 5.0
        Tmin    = THB_dys*(fed**(-3.0))*((1.0 - delta_bs)/(9.0*delta_bs))**(3./2.)
        pTmin_x = [np.log10(Tmin), np.log10(Tmin)]
        pTmin_y = [1.05, 1.15]
        ax1.plot(pTmin_x, pTmin_y, linestyle = '-', linewidth = 1, color = m_color_arr[mc])
        if (mc == 0):
            ax1.text(-1.45, pTmin_y[0], r'$T_{\rm min,ej} \rightarrow$',
            horizontalalignment='left',
            verticalalignment='bottom',
            rotation='horizontal')
        #plot limits:
        Tmin    = THB_dys*(Pm_HB**(21./20))
        pTmin_x = [np.log10(Tmin), np.log10(Tmin)]
        pTmin_y = [1.15, 1.25]
        ax1.plot(pTmin_x, pTmin_y, linestyle = '-', linewidth = 1, color = m_color_arr[mc])
        if (mc == 0):
            ax1.text(-1.45, pTmin_y[0], r'$T_{\rm min,m} \rightarrow$',
            horizontalalignment='left',
            verticalalignment='bottom',
            rotation='horizontal')
        #plot limits:
        Tmin    = THB_dys*((tbs_HB_SI/(sec_year*(10.**10)))**(3./2.))*(delta_bs/(1.0 - delta_bs))**(3./2.)
        pTmin_x = [np.log10(Tmin), np.log10(Tmin)]
        pTmin_y = [1.25, 1.35]
        ax1.plot(pTmin_x, pTmin_y, linestyle = '-', linewidth = 1, color = m_color_arr[mc])
        if (mc == 0):
            ax1.text(-1.45, pTmin_y[0], r'$T_{\rm min,\tau}\ \rightarrow$',
            horizontalalignment='left',
            verticalalignment='bottom',
            rotation='horizontal')        
        #------------------------------------------------------
    #----------------------------------------------------------
    #---------------------------------------------------------- 
    if (fignr == 3 and mc == 0):
    #----------------------------------------------------------
        #------------------------------------------------------
        #FIG 3:
        #------------------------------------------------------        
        #define:
        arr_state_n     = np.arange(0.0, 40.0, 1)
        Pnm_state_n     = np.exp(-Pm_HB*((1.0 - Eps**(arr_state_n + 1.0))/(1.0 - Eps)))
        Pobs_state_n    = 0.0001*delta_bs**(-arr_state_n)
        Pfin_state_n    = Pnm_state_n*Pobs_state_n
        #plot data:
        ax1.plot(arr_state_n, (Pnm_state_n),    linestyle = '--',   linewidth = 1, color = 'black', label   = r'$\mathscr{P}(n)$')
        ax1.plot(arr_state_n, (Pobs_state_n),   linestyle = ':',    linewidth = 1, color = 'black', label   = r'$\mathscr{W}_{n}$')
        ax1.plot(arr_state_n, (Pfin_state_n),   linestyle = '-',    linewidth = 2, color = 'black', label   = r'$\mathscr{W}_{n} \times \mathscr{P}(n)$')
        #plot settings:
        ax1.set_xlim([0.0, 40.0])
        ax1.set_ylim([1e-5, 1e1])
        ax1.set_xlabel(r'BBH hardening step $n$')
        ax1.set_ylabel(r'[rnd norm.]')
        plt.legend(loc='upper right', numpoints = 1, fontsize = 12.0, ncol = 3, frameon = False)
        ax1.set_yscale('log')
        #save/show fig:
        fig.savefig('PnmPobs_ex.eps', bbox_inches='tight')
        plt.show()
        exit()
        #------------------------------------------------------
    #----------------------------------------------------------
    #----------------------------------------------------------
    
#----------------------------------------------------------
#plot settings: FIG. 2
#----------------------------------------------------------
if (fignr == 2):
    #dummy plot for adding label:
    ax1.plot(-1000, -1000,  linestyle = '-', linewidth = 2, color = 'black', label = r'$+$ GW mergers')
    ax1.plot(-1000, -1000,  linestyle = ':', linewidth = 1, color = 'black', label = r'$-$ GW mergers')
    #plot settings:
    ax1.set_xlim([-1.6, 3.5])
    ax1.set_ylim([0.0, 1.4])
    ax1.set_xlabel(r'log $T$ [days]')
    ax1.set_ylabel(r'$d\Gamma_{\rm TD}/d\log T$ [rnd. norm.]')
    #ax1.set_title(r'TDE Rates from Binary Black Holes')
    plt.legend(loc='center left', numpoints = 1, fontsize = 9.0, ncol = 1, frameon = False)

    #save/show fig:
    fig.savefig('Torb_dist.eps', bbox_inches='tight')
    plt.show()
    #exit()
#----------------------------------------------------------




#exit()




#----------------------------------------------------------
#plot settings: FIG. 4
#----------------------------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 4))    
#----------------------------------------------------------
#PLOT data:
#----------------------------------------------------------
cmc_data_file   = '/Users/jsamsing/Desktop/TIDES_PROJ/BBH_TDE/CMC_data/'
tf = open(cmc_data_file + 'BBH_2e6.dat',  "r")
cmc_data        = np.loadtxt(tf, dtype=float)
cmc_Torb_dys    = cmc_data[:,4]
cmc_m1_Msun     = cmc_data[:,0]
cmc_m2_Msun     = cmc_data[:,1]
tf.close()
#BH mass cut:
cmc_ml  = 10.0
cmc_mu  = 20.0
pos_m1l  = np.where(cmc_m1_Msun[:] > cmc_ml)[0]
pos_m1u  = np.where(cmc_m1_Msun[:] < cmc_mu)[0]
pos_m1lu = list(set(pos_m1l).intersection(pos_m1u))
pos_m2l  = np.where(cmc_m2_Msun[:] > cmc_ml)[0]
pos_m2u  = np.where(cmc_m2_Msun[:] < cmc_mu)[0]
pos_m2lu = list(set(pos_m2l).intersection(pos_m2u))
pos_m12lu   = list(set(pos_m1lu).intersection(pos_m2lu))
#print info:
print cmc_Torb_dys
print cmc_m1_Msun[pos_m12lu]
print cmc_m2_Msun[pos_m12lu]
#make Torb hist:
plot_data   = cmc_Torb_dys[pos_m12lu]
cmc_weight  = np.zeros(len(plot_data)) + 0.000079
ax1.hist(np.log10(plot_data), weights = cmc_weight, range=(-2,5), bins=28, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'grey', label=r'CMC data', hatch = '\\\\\\\\')
#----------------------------------------------------------
#PLOT models:
#----------------------------------------------------------
#PARAM 1:
vdis_kms    = 12.0
m_Msun      = 15.0
n_nr_SI     = (10.**(4.))/(m_parsec**3.)
vdis_SI     = vdis_kms*1000.
m_SI        = m_Msun*M_sun_SI
a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(vdis_SI**2.))
tbs_HB_SI   = 1./(n_nr_SI*(2.*np.pi*G_new_SI*(3.*m_SI)*a_HB_SI/(vdis_SI**2.))*vdis_SI)
tcl_HB_SI   = (a_HB_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
Pm_HB       = (tbs_HB_SI/tcl_HB_SI)**(2./7.)
THB_yrs     = (2.*np.pi*np.sqrt(a_HB_SI**3./(G_new_SI*(m_SI + m_SI))))/sec_year
THB_dys     = 365.*THB_yrs
fed = 5.0
Tmin_dys        = THB_dys*(fed**(-3.0))*((1.0 - delta_bs)/(9.0*delta_bs))**(3./2.)  #Tej
arr_Torb_dys    = THB_dys*(10.**(np.arange(0.0, np.log10(Tmin_dys/THB_dys), -0.01)))    #added lower cutoff
Pnm_Torb_dys    = np.exp(-Pm_HB*((1.0 - Eps*((arr_Torb_dys/THB_dys)**(-20./21.)))/(1.0 - Eps)))
PTD_GRyes       = (arr_Torb_dys**(-2./3.))*Pnm_Torb_dys
PTD_GRyes       = PTD_GRyes/max(PTD_GRyes)
#adjust arrays to include Tmin:
arr_Torb_dys    = np.append(arr_Torb_dys, Tmin_dys)
PTD_GRyes       = np.append(PTD_GRyes, 0.0)
#plot:
normfac = 1.0
ax1.plot(np.log10(arr_Torb_dys), PTD_GRyes*normfac, linestyle = ':', linewidth = 3, color = 'black', label = r'analytical model (params. 1)')
#PARAM 2:
vdis_kms    = 10.0
m_Msun      = 15.0
n_nr_SI     = (10.**(3.))/(m_parsec**3.)
vdis_SI     = vdis_kms*1000.
m_SI        = m_Msun*M_sun_SI
a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(vdis_SI**2.))
tbs_HB_SI   = 1./(n_nr_SI*(2.*np.pi*G_new_SI*(3.*m_SI)*a_HB_SI/(vdis_SI**2.))*vdis_SI)
tcl_HB_SI   = (a_HB_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
Pm_HB       = (tbs_HB_SI/tcl_HB_SI)**(2./7.)
THB_yrs     = (2.*np.pi*np.sqrt(a_HB_SI**3./(G_new_SI*(m_SI + m_SI))))/sec_year
THB_dys     = 365.*THB_yrs
fed = 5.0
Tmin_dys        = THB_dys*(fed**(-3.0))*((1.0 - delta_bs)/(9.0*delta_bs))**(3./2.)  #Tej
arr_Torb_dys    = THB_dys*(10.**(np.arange(0.0, np.log10(Tmin_dys/THB_dys), -0.01)))    #added lower cutoff
Pnm_Torb_dys    = np.exp(-Pm_HB*((1.0 - Eps*((arr_Torb_dys/THB_dys)**(-20./21.)))/(1.0 - Eps)))
PTD_GRyes       = (arr_Torb_dys**(-2./3.))*Pnm_Torb_dys
PTD_GRyes       = PTD_GRyes/max(PTD_GRyes)
#adjust arrays to include Tmin:
arr_Torb_dys    = np.append(arr_Torb_dys, Tmin_dys)
PTD_GRyes       = np.append(PTD_GRyes, 0.0)
#plot:
normfac = 0.7
ax1.plot(np.log10(arr_Torb_dys), PTD_GRyes*normfac, linestyle = '-.', linewidth = 3, color = 'black', label = r'analytical model (params. 2)')

#plot settings:
ax1.set_xlim([0.0, 5.0])
ax1.set_ylim([-0.01, 1.01])
ax1.set_xlabel(r'log $T$ [days]')
ax1.set_ylabel(r'$dN_{\rm BBH}/d\log T$ [rnd. norm.]')
plt.legend(loc='upper right', numpoints = 1, fontsize = 9.0, ncol = 1, frameon = False)

#save/show fig:
fig.savefig('Torb_dist_CMCdata.eps', bbox_inches='tight')
plt.show()
exit()
#----------------------------------------------------------





















exit()
#----------------------------------------------------------
