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


#----------------------------------------------------------
#calc t_insptime:
#m_SI        = 10.0*M_sun_SI
#a_SI        = (1./6.)*((1./(7./9.) - 1.0)*G_new_SI*m_SI/((500.*1000.)**2.)) #0.005*AU_SI
#tinsp_yr    = (a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))/sec_year
#print (10.0/tinsp_yr)**(2./7.)
#exit()

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

nres                = 100
v_kms_arr           = (10.**(np.linspace(np.log10(6.0), np.log10(110.0), nres)))
a_AU_arr            = (10.**(np.linspace(-1.4, 1.3, nres)))
m_Msun_arr          = (10.**(np.linspace(np.log10(5.0),  np.log10(100.0), nres)))

Pecc_over_Pcir_am   = np.zeros((nres,nres), dtype=np.float64)
Pecc_over_Pcir_vm   = np.zeros((nres,nres), dtype=np.float64)
P_IsolatedMerger_vm = np.zeros((nres,nres), dtype=np.float64)
Pcap_over_Pcir_vm   = np.zeros((nres,nres), dtype=np.float64)

test1   = np.zeros((nres,nres), dtype=np.float64)
test1[:,:]  = -1.0

#input:
tlimit_yr   = 10.**(10.)    #in years
f_ob        = 10.0
e_ob        = 0.1
Nims        = 20.0
delta_bs    = 7./9.
#calc/define:
Fe_ob   = (e_ob**(12./19.)/(1.+e_ob))*(1. + (121./304.)*(e_ob**2.))**(870./2299.)
#----------------------------------------------------------

#----------------------------------------------------------
for yc in range(0,nres):        #m
    for xc in range(0,nres):    #a or v
        
        
        #calc Pecc_over_Pcir:
        #read a,m:
        a_AU    = a_AU_arr[xc]
        m_Msun  = m_Msun_arr[yc]
        #different units:
        a_SI    = a_AU*AU_SI
        m_SI    = m_Msun*M_sun_SI  
        #calc P(ecc M):
        rf_SI   = (2.*G_new_SI*m_SI/((f_ob**2.)*(np.pi**2.)))**(1./3.)
        rEM_SI  = rf_SI*(1./(2.*Fe_ob))*((425./304.)**(870./2299.))        
        PeccM           = (1./(1.-delta_bs))*(Nims*(2.*(rEM_SI/a_SI)))
        #calc P(cir M):
        tinsp_e0_yr = (a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))/sec_year
        if (tinsp_e0_yr < tlimit_yr):
            PcirM       = 1.0
        if (tinsp_e0_yr > tlimit_yr):
            PcirM       = ((tlimit_yr/tinsp_e0_yr))**(2./7.)    #SKIPPED HERE (425./768) OTHERWISE IT DOES NOT FIT WITH THE OTHER LIMIT.
        #save:    
        Pecc_over_Pcir_am[xc,yc]   = PeccM/PcirM


        #calc Pecc_over_Pcir:
        #read v,m:
        v_kms   = v_kms_arr[xc]     #v_esc
        m_Msun  = m_Msun_arr[yc]    #mass
        #convert to corresponding a_ej:
        a_AU    = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*(m_Msun*M_sun_SI)/((v_kms*1000.)**(2.)))/AU_SI  #a_ej
        #different units:
        a_SI    = a_AU*AU_SI
        m_SI    = m_Msun*M_sun_SI  
        #calc P(ecc M):
        rf_SI   = (2.*G_new_SI*m_SI/((f_ob**2.)*(np.pi**2.)))**(1./3.)
        rEM_SI  = rf_SI*(1./(2.*Fe_ob))*((425./304.)**(870./2299.))        
        PeccM           = (1./(1.-delta_bs))*(Nims*(2.*(rEM_SI/a_SI)))
        #calc P(cir M):
        tinsp_e0_yr = (a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))/sec_year
        if (tinsp_e0_yr < tlimit_yr):
            PcirM       = 1.0
        if (tinsp_e0_yr > tlimit_yr):
            PcirM       = ((tlimit_yr/tinsp_e0_yr))**(2./7.)    #SKIPPED HERE (425./768) OTHERWISE IT DOES NOT FIT WITH THE OTHER LIMIT.
        #save:
        Pecc_over_Pcir_vm[xc,yc]   = PeccM/PcirM
        
        
        #calc Pcap_over_Pcir_vm:
        #read v,m:
        v_SI    = v_kms_arr[xc]*1000.       #v_esc
        m_SI    = m_Msun_arr[yc]*M_sun_SI   #mass        
        #convert to corresponding a_ej:
        a_SI    = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(v_SI**2.))  #a_ej        
        #calc P(cap M):
        Rsch_SI     = (2.*G_new_SI*m_SI)/(c_SI**2.)
        rcap_aej    = Rsch_SI*(a_SI/Rsch_SI)**(2./7.)
        Pcap_aej    = (2.*rcap_aej/a_SI)*Nims
        Pcap_ainaej = (7./5.)*(1./(1.-delta_bs))*Pcap_aej        
        #calc P(cir M):
        tinsp_e0_yr = (a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))/sec_year
        if (tinsp_e0_yr < tlimit_yr):
            PcirM       = 1.0
        if (tinsp_e0_yr > tlimit_yr):
            PcirM       = ((tlimit_yr/tinsp_e0_yr))**(2./7.)    #SKIPPED HERE (425./768) OTHERWISE IT DOES NOT FIT WITH THE OTHER LIMIT.
        #save:
        Pcap_over_Pcir_vm[xc,yc]   = Pcap_ainaej/PcirM
        
        
        #calc P_IM (int P for isolated merger between interactions):
        #assume v_esc and v_disp are similar.
        #read v,m:
        v_SI    = v_kms_arr[xc]*1000.       #v_esc
        m_SI    = m_Msun_arr[yc]*M_sun_SI   #mass        
        #convert to corresponding a_ej:
        a_SI    = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(v_SI**2.))  #a_ej
        #calc binary-single (bs) interaction time:
        n_nr_SI     = (10.**(6.))/(m_parsec**3.) 
        sigma_bs_SI = 2.*np.pi*G_new_SI*(3.*m_SI)*a_SI/(v_SI**2.)
        Gamma_bs_SI = n_nr_SI*sigma_bs_SI*v_SI
        time_bs_SI  = 1./Gamma_bs_SI
        #calc binary GW life time e=0:
        tinsp_e0_SI = (a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
        #calc P_IM:
        P_IM_aej    = (time_bs_SI/tinsp_e0_SI)**(2./7.)
        P_IM        = (7./10.)*(1./(1.-delta_bs))*P_IM_aej
        #save:
        P_IsolatedMerger_vm[xc,yc] = P_IM
     
     
     
     
     
        
     
        #calc P_IM (int P for isolated merger between interactions):
        #assume v_esc and v_disp are similar.
        #read v,m:
        v_SI    = v_kms_arr[xc]*1000.       #v_esc
        m_SI    = m_Msun_arr[yc]*M_sun_SI   #mass        
        #input:
        n_nr_SI = (10.**(6.))/(m_parsec**3.)
        tH_SI   = (10.**10.)*sec_year
        #define:
        eps_cap = delta_bs**(-5./7.)
        eps_EM  = delta_bs**(-1./1.)
        eps_IM  = delta_bs**(-10./7.)
        eps_t   = delta_bs**(-1./1.)
        #calc a_HB:
        a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(v_SI**2.))
        #calc a_ej:
        a_ej_SI     = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(v_SI**2.))  #a_ej
        #calc a_Ht:
        k_bs        = (2.*np.pi*G_new_SI*(3.*m_SI)/(v_SI**2.))
        C_bs        = (1./(n_nr_SI*k_bs*v_SI))
        a_Ht_SI     = (C_bs/(tH_SI*(eps_t - 1.0)))
        #calc a_IM:
        a_IM_SI     = (C_bs/(1./(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))))**(1./5.)
        #define a_max, a_min (the system transitions from a_max to a_min):
        a_max_SI    = a_HB_SI       
        a_min_SI    = max([a_ej_SI, a_Ht_SI, a_IM_SI])    
        #calc number of bin-sin int from a_max_SI to a_min_SI:
        Nbs         = int((np.log(a_min_SI/a_max_SI)/np.log(delta_bs)))
        #redefine a_min:
        a_min_SI    = a_max_SI*(delta_bs**(Nbs))
        if (Nbs > 1):
            #time to reach a_min (time inside Cluster):
            T_inC_SI    = (C_bs/(a_min_SI*(eps_t - 1.0)))
            #if ejected, the BBH has this time outside (time outside Cluster):
            T_outC_SI   = tH_SI - T_inC_SI
            if (T_outC_SI <= 0.0): T_outC_SI = 0.0
            #calc P_cap at a_max:
            Rsch_SI     = (2.*G_new_SI*m_SI)/(c_SI**2.)
            rcap_amax   = Rsch_SI*(a_max_SI/Rsch_SI)**(2./7.)
            Pcap_amax   = (2.*rcap_amax/a_max_SI)*Nims              #Pcap_amax
            #calc P_EM at a_max:
            rf_SI       = (2.*G_new_SI*m_SI/((f_ob**2.)*(np.pi**2.)))**(1./3.)
            rEM_SI      = rf_SI*(1./(2.*Fe_ob))*((425./304.)**(870./2299.))        
            PEM_amax    = (2.*rEM_SI/a_max_SI)*Nims                 #PEM_amax
            #calc P_IM at a_max:
            sigma_bs_amax_SI    = 2.*np.pi*G_new_SI*(3.*m_SI)*a_max_SI/(v_SI**2.)
            Gamma_bs_amax_SI    = n_nr_SI*sigma_bs_amax_SI*v_SI
            tbs_amax_SI         = 1./Gamma_bs_amax_SI                                                           #binary-single (bs) interaction time at a_ej.
            tinsp_e0_amax_SI    = (a_max_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))  #binary GW life time e=0 at a_ej.
            PIM_amax    = (tbs_amax_SI/tinsp_e0_amax_SI)**(2./7.)   #PIM_amax
            #calc probabilities:
            lnPNm_sum_n = 0.0
            P_narr_EM   = np.zeros(Nbs+1, dtype=np.float64)
            #loop over interactions from a_max(n=0) to a_min(n=Nbs):
            for n in range(0,Nbs+1):
                #calc probabilities at this n(=a_n):
                Pm_n_cap    = Pcap_amax*(eps_cap**n)    #prob at n(=a_n) for cap merger.
                Pm_n_EM     = PEM_amax*(eps_EM**n)      #prob at n(=a_n) for ecc Merger.
                Pm_n_IM     = PIM_amax*(eps_IM**n)      #prob at n(=a_n) for isolated Merger.
                Pm_n_tot    = Pm_n_cap + Pm_n_IM        #tot prob at n(=a_n) for merger.
                PNm_n_tot   = 1. - Pm_n_tot             #tot prob at n(=a_n) for NOT merger.
                #correct:
                if (PNm_n_tot < 0.0): PNm_n_tot = 1e-10
                #calc prob at this n(=a_n) INCLUDING previous n:
                lnPNm_sum_n     = lnPNm_sum_n + np.log(PNm_n_tot)
                PNm_sum_n       = np.exp(lnPNm_sum_n)
                P_narr_EM[n]    = Pm_n_EM*PNm_sum_n
            #define:
            PNm_amax_amin   = np.exp(lnPNm_sum_n)             

            #The following cuts are not be necessary, its just to make it look better:
            if (a_IM_SI > max([a_ej_SI,a_Ht_SI])): PNm_amax_amin = 0.0
                        
            Pm_amax_amin    = 1.0 - PNm_amax_amin            
            #calc total prob for merger <tH INSIDE cluster:            
            Pm_tot_inC          = Pm_amax_amin
            #calc total prob for merger <tH OUTSIDE cluster:
            tinsp_e0_amin_SI    = (a_min_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))  #binary GW life time e=0 at a_ej.
            Pm_tot_outC         = (1.0 - Pm_tot_inC)*((T_outC_SI/tinsp_e0_amin_SI)**(2./7.))
            #total probability for that BBH merge within a hubble time.
            Pm_tot              = Pm_tot_inC + Pm_tot_outC    
            if (Pm_tot > 1.0): Pm_tot = 1.0
            #calc total EM prob:
            PEM_tot = sum(P_narr_EM[:])

            #calc final info:
            # ...
            #save:
            test1[xc,yc] = PEM_tot/Pm_tot
            #PEM_tot/Pm_tot#Pm_tot_inC/Pm_tot #max([a_Ht_SI, a_IM_SI, a_ej_SI])/a_IM_SI #Pm_tot#PEM_tot/Pm_tot#Pm_tot_inC/Pm_tot #max([a_Ht_SI, a_IM_SI, a_ej_SI])/a_ej_SI #Pm_tot_inC/Pm_tot #PEM_tot/Pm_tot

            #ALSO: incl. varying density n_s. 
            #chekc for prob that exceeds one.
            #problem that N is not always >> 1?  
            #CHECK everything again.     
            #is Pm_tot_outC correct?
            #are all the N factors OK? are we always considering the same amin limit?
            #CHECK ALL N things and limits again!! Make an assumption and follow that. Its not perfect, but ok.
#----------------------------------------------------------                
       
       
       




  



       
       
#----------------------------------------------------------      
#PLOT:  
#----------------------------------------------------------      
fig     = plt.figure(figsize=(6, 5))

#make plot window:
ax = plt.subplot(111)

#define x,y:
x_input = v_kms_arr
y_input = m_Msun_arr
X, Y    = np.meshgrid(x_input, y_input)


#contour plot P_IM:
Z       = np.transpose(test1*100.)    #P in percent
#Z       = np.transpose(test1)

CS = plt.contourf(X, Y, Z, list([0,1]), colors='blue',      alpha = 0.75, hatches=[None], extend='lower')
CS = plt.contourf(X, Y, Z, list([1,10]), colors='green',      alpha = 0.75, hatches=[None], extend='lower')
CS = plt.contourf(X, Y, Z, list([10,99]), colors='grey',      alpha = 0.75, hatches=[None], extend='lower')

CS = plt.contourf(X, Y, Z, list([99,1000000]), colors='red', alpha = 0.75, hatches=[''], extend='lower')
CS_IM = plt.contour(X, Y, Z, list([1,10,50,100]), linestyles = ':', linewidths=2, colors='blue')
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
# Recast levels to new class
CS_IM.levels = [nf(val) for val in CS_IM.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r %%'
plt.clabel(CS_IM, CS_IM.levels, inline=True, fmt=fmt, fontsize=12)


#contour plot P_EM:
Z       = np.transpose(test1*100)       #P in percent
#Z       = np.transpose(test1)

CS_EM = plt.contour(X, Y, Z, list([1,2,3,4,5,6,7,8,9,10,20,30,40,50,60, 70, 75, 80, 90, 100]), linewidths=2, colors='black')
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
# Recast levels to new class
CS_EM.levels = [nf(val) for val in CS_EM.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r %%'
plt.clabel(CS_EM, CS_EM.levels, inline=True, fmt=fmt, fontsize=12)


#axis:
ax.set_xlabel(r'escape velocity $v_{\rm esc}$ [kms$^{-1}$]')
ax.set_ylabel(r'black hole mass $m_{\rm}$ $[M_{\odot}]$')
ax.set_title(r'GW Merger Probability')

ax.set_xlim([min(x_input),  max(x_input)])
ax.set_ylim([min(y_input),  max(y_input)])

plt.xscale('log')
plt.yscale('log')

#legend:
#lines = [CS_IM.collections[0], CS_EM.collections[0]]
#labels = [r'$P_{\rm IM}$ ($n_{\rm s} = 10^{6}$ pc$^{-3}$)',r'$\Gamma_{\rm EM}/\Gamma_{\rm CM}$ ($P_{\rm IM} \ll 1$)']
#plt.legend(lines, labels, loc='upper right', numpoints = 1, fontsize = 10.0, frameon = True)





#save and show:
plt.savefig('test.pdf', bbox_inches='tight')     
plt.show()

exit()
#----------------------------------------------------------




       
       














#----------------------------------------------------------      
#PLOT:  
#----------------------------------------------------------      
fig     = plt.figure(figsize=(6, 5))

#make plot window:
ax = plt.subplot(111)

#define x,y:
x_input = v_kms_arr
y_input = m_Msun_arr
X, Y    = np.meshgrid(x_input, y_input)




#contour plot P_CM:
Z       = np.transpose(Pcap_over_Pcir_vm*100)       #P in percent
CS_CM   = plt.contour(X, Y, Z, list([1,2,3,4,5,6,8,10,12,14,16,18,20,25,50]), linewidths=2, colors='black')
CSTST   = plt.pcolormesh(X, Y, Z, cmap='summer', linewidth=0, rasterized=True, vmin=0, vmax=20.0)
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
# Recast levels to new class
CS_CM.levels = [nf(val) for val in CS_CM.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r %%'
plt.clabel(CS_CM, CS_CM.levels, inline=True, fmt=fmt, fontsize=12)



#contour plot P_IM:
Z       = np.transpose(P_IsolatedMerger_vm*100.)    #P in percent
#CS = plt.contourf(X, Y, Z, list([0,10]), colors='green',      alpha = 0.75, hatches=[None], extend='lower')
#CS = plt.contourf(X, Y, Z, list([10,100]), colors='grey',      alpha = 0.75, hatches=[None], extend='lower')
CS = plt.contourf(X, Y, Z, list([10,1000000]), colors='black', alpha = 0.75, hatches=['/'], extend='lower')
CS_IM = plt.contour(X, Y, Z, list([1,10,50,100]), linestyles = ':', linewidths=2, colors='blue')
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
# Recast levels to new class
CS_IM.levels = [nf(val) for val in CS_IM.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r %%'
plt.clabel(CS_IM, CS_IM.levels, inline=True, fmt=fmt, fontsize=12)





#annotations:
ax.text(35, 6, r'$P_{\rm IM} \approx 1$ ($n_{\rm s} = 10^{6}$ pc$^{-3}$):', color='yellow', fontsize = 10)
ax.text(35, 5.25, 'BBH merge before ejection', color='yellow', fontsize = 11)

#axis:
ax.set_xlabel(r'escape velocity $v_{\rm esc}$ [kms$^{-1}$]')
ax.set_ylabel(r'black hole mass $m_{\rm}$ $[M_{\odot}]$')
#ax.set_title(r'$N_{\rm cap}/N_{\rm ej}$')

ax.set_xlim([min(x_input),  max(x_input)])
ax.set_ylim([min(y_input),  max(y_input)])

plt.xscale('log')
plt.yscale('log')

#legend:
lines = [CS_IM.collections[0], CS_CM.collections[0]]
labels = [r'$P_{\rm IM}$ ($n_{\rm s} = 10^{6}$ pc$^{-3}$)',r'$P_{\rm bs}/P_{\rm ej}$ ($P_{\rm IM} \ll 1$)']
plt.legend(lines, labels, loc='upper right', numpoints = 1, fontsize = 10.0, frameon = True)



#save and show:
plt.savefig('Pcap_over_PHM.pdf', bbox_inches='tight')     
plt.show()

exit()
#----------------------------------------------------------











#----------------------------------------------------------      
#PLOT:  
#----------------------------------------------------------      
fig     = plt.figure(figsize=(6, 5))

#make plot window:
ax = plt.subplot(111)

#define x,y:
x_input = v_kms_arr
y_input = m_Msun_arr
X, Y    = np.meshgrid(x_input, y_input)


#contour plot P_IM:
Z       = np.transpose(P_IsolatedMerger_vm*100.)    #P in percent
CS = plt.contourf(X, Y, Z, list([0,10]), colors='forestgreen',      alpha = 0.75, hatches=[None], extend='lower')
CS = plt.contourf(X, Y, Z, list([10,100]), colors='grey',      alpha = 0.75, hatches=[None], extend='lower')

CS = plt.contourf(X, Y, Z, list([100,1000000]), colors='indianred', alpha = 0.75, hatches=[None], extend='lower')
CS_IM = plt.contour(X, Y, Z, list([1,10,50,100]), linestyles = ':', linewidths=2, colors='blue')
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
# Recast levels to new class
CS_IM.levels = [nf(val) for val in CS_IM.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r %%'
plt.clabel(CS_IM, CS_IM.levels, inline=True, fmt=fmt, fontsize=12)


#contour plot P_EM:
Z       = np.transpose(Pecc_over_Pcir_vm*100)       #P in percent
CS_EM = plt.contour(X, Y, Z, list([1,2,3,4,5,6,7,8,9,10,25,50,75,100]), linewidths=2, colors='black')
class nf(float):
    def __repr__(self):
        str = '%.1f' % (self.__float__(),)
        if str[-1] == '0':
            return '%.0f' % self.__float__()
        else:
            return '%.1f' % self.__float__()
# Recast levels to new class
CS_EM.levels = [nf(val) for val in CS_EM.levels]
# Label levels with specially formatted floats
if plt.rcParams["text.usetex"]:
    fmt = r'%r \%%'
else:
    fmt = '%r %%'
plt.clabel(CS_EM, CS_EM.levels, inline=True, fmt=fmt, fontsize=12)


#annotations:
#ax.text(70, 9, r'$P_{\rm IM}$ = 1 ($n_{\rm s} = 10^{6}$ pc$^{-3}$):', color='yellow', fontsize = 10)
#ax.text(70, 8, 'BBH merge before ejection', color='yellow', fontsize = 10)

#axis:
ax.set_xlabel(r'escape velocity $v_{\rm esc}$ [kms$^{-1}$]')
ax.set_ylabel(r'black hole mass $m_{\rm}$ $[M_{\odot}]$')
ax.set_title(r'GW Merger Probability')

ax.set_xlim([min(x_input),  max(x_input)])
ax.set_ylim([min(y_input),  max(y_input)])

plt.xscale('log')
plt.yscale('log')

#legend:
lines = [CS_IM.collections[0], CS_EM.collections[0]]
labels = [r'$P_{\rm IM}$ ($n_{\rm s} = 10^{6}$ pc$^{-3}$)',r'$\Gamma_{\rm EM}/\Gamma_{\rm CM}$ ($P_{\rm IM} \ll 1$)']
plt.legend(lines, labels, loc='upper right', numpoints = 1, fontsize = 10.0, frameon = True)





#save and show:
plt.savefig('data_xy_input.eps', bbox_inches='tight')     
plt.show()

exit()
#----------------------------------------------------------







