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
Rsch_1Msun_SI   = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))
#----------------------------------------------------------



#FOR KAZAs proj: make sure that the HB cross sections are ok for high velocity disp!!
#exit()

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


#----------------------------------------------------------
#settings:
#----------------------------------------------------------
filename = raw_input('filename: ')
#input/edit:
Nsamp       = 1000

v_kms       = 10.0      #DISPERSION vel
vesc_kms    = 50.0      #ESCAPE     vel
m_Msun      = 10.0
n_nr_SI     = (10.**(5.))/(m_parsec**3.)

fGW_limit   = 1e-10

#define/calc:
tH_SI       = (10.**10.)*sec_year
Nims        = 20
delta_bs    = 7./9.
v_SI        = v_kms*1000.
vesc_SI     = vesc_kms*1000.
m_SI        = m_Msun*M_sun_SI

#calc a_HB:
a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(v_SI**2.))
#calc a_ej:
a_ej_SI     = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(vesc_SI**2.))
#define:
a_max_SI = a_HB_SI
a_min_SI = a_ej_SI
#calc nr. int:
Nbs         = int((np.log(a_min_SI/a_max_SI)/np.log(delta_bs)))
#redefine a_min:
a_min_SI    = a_max_SI*(delta_bs**(Nbs))

#arrays:
output_INT_arr      = np.zeros((Nsamp,10), dtype=int)
output_REAL_arr     = np.zeros((Nsamp,10), dtype=np.float64)

#for background calculations (BGcalc):
tot_nr_bins = (Nbs+1)*Nsamp 
bin_info_arr_REAL   = np.zeros((tot_nr_bins,5), dtype=np.float64)
bin_info_arr_INT    = np.zeros((tot_nr_bins,5), dtype=int)
tttc        = 0
#----------------------------------------------------------
#----------------------------------------------------------






#TEST PRINT INFO:
#----------------

print Nbs
print a_max_SI
print a_min_SI

kfac    = (2.*np.pi*G_new_SI*(3.*m_SI)/(v_SI**2.))
Cfac    = (1./(n_nr_SI*kfac*v_SI))
Tej     = ((Cfac/(delta_bs**(-1.0) -1.0))*(1./a_ej_SI))/tH_SI

sigma_bs_SI = 2.*np.pi*G_new_SI*(3.*m_SI)*a_ej_SI/(v_SI**2.)
Gamma_bs_SI = n_nr_SI*sigma_bs_SI*v_SI
tbs_SI      = 1./Gamma_bs_SI
tc_SI       = (a_ej_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
rm          = (1./2.)*a_ej_SI*(tbs_SI/tc_SI)**(2./7.)
fm          = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rm**3.))
print np.log10(fm)
print np.log10(((v_kms/50.)**(-12./7.))*((n_nr_SI/((10.**(5.))/(m_parsec**3.)))**(3./7.))*(10**(-3.837248)))
print 'Tej', Tej
print (7./10.)*(1./(1.-delta_bs))*(tbs_SI/tH_SI)**(2./7.)
print 0.79827431101*(10.**(-2./7.))

print (tbs_SI/tc_SI)**(2./7.), 'ES1'
print (((512.*G_new_SI**2.)/(30.*np.pi*c_SI**5.))*((m_SI**2.)*v_SI/((a_ej_SI**5.)*n_nr_SI)))**(2./7.), 'ES2'
print ((8192.*v_SI**11.)/(3645.*np.pi*G_new_SI**3.*c_SI**5.*m_SI**3.*n_nr_SI))**(2./7.), 'ESHB at HB limit'


a   = 1.0*AU_SI
sigma_bs_SI = 2.*np.pi*G_new_SI*(3.*m_SI)*a/(v_SI**2.)
Gamma_bs_SI = n_nr_SI*sigma_bs_SI*v_SI
tbs_SI      = 1./Gamma_bs_SI
tc_SI       = (a**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
print np.sqrt(1.-(tbs_SI/tc_SI)**(2./7.)), 'ES3'


print ((768./425.)*((((0.5*AU_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.)))))*(1-(0.8**2.))**(7./2.)

e_0     = 0.99
Omega   = -np.pi/4.
iang    = np.pi/2.   

print ((e_0*np.sqrt(1-e_0**2.)/(e_0 - 1.))*np.sin(2.*Omega)*np.sin(iang)**2.*15.*np.pi/(16.*np.sqrt(3.)))**(2./3.)

mBH     = 20.*M_sun_SI
aBBH    = 0.5*AU_SI
eBBH    = 0.99
print  (((c_SI**4.)/(G_new_SI**2.))*(aBBH**2./(mBH**2.))*(1.-eBBH**2.)**(2.))**(1./3.)



#HAVE TO CHEKC WITH SIM!!!

t0_int  = (10**5.)*sec_year
a0      = 10.0*AU_SI 
ak      = 0.1*AU_SI
tk_GW   = (ak**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
t0_GW   = (a0**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))

#co-planar:
delta_cpl   = 2./3.
t0_int  = t0_int*np.sqrt(a0/a0)
tk_int  = t0_int*np.sqrt(a0/ak)
print 'planar:'
p2_0    = (t0_int/t0_GW)**(1./7.)
p2_k    = (tk_int/tk_GW)**(1./7.)
print p2_k, 'p2_k'
P2_0k   = (1./(1.-delta_cpl))*(7./5.)*(p2_k - p2_0)
Rcm     = 2.*G_new_SI*m_SI/(c_SI**2.)
#p3_0    = (np.sqrt(2.)*20.)*(Rcm/a0)**(5./14.)
#p3_k    = (np.sqrt(2.)*20.)*(Rcm/ak)**(5./14.)
p3_0    = (np.sqrt(2.)*20.)*(Rcm/a0)**(5./28.)
p3_k    = (np.sqrt(2.)*20.)*(Rcm/ak)**(5./28.)
print p3_k, 'p3_k'
P3_0k   = (1./(1.-delta_cpl))*(14./5.)*(p3_k-p3_0)
print P2_0k, P3_0k, 'P2_0k, P3_0k'
print a0/m_parsec, ak/m_parsec

#isotropic:
delta_iso   = 7./9.
t0_int  = t0_int*(a0/a0)
tk_int  = t0_int*(a0/ak)
print 'isotropic:'
p2_0    = (t0_int/t0_GW)**(2./7.)
p2_k    = (tk_int/tk_GW)**(2./7.)
print p2_k, 'p2_k'
P2_0k   = (1./(1.-delta_iso))*(7./10.)*(p2_k - p2_0)
Rcm     = 2.*G_new_SI*m_SI/(c_SI**2.)
p3_0    = (2.*20.)*(Rcm/a0)**(5./7.)
p3_k    = (2.*20.)*(Rcm/ak)**(5./7.)
print p3_k, 'p3_k'
P3_0k   = (1./(1.-delta_iso))*(7./5.)*(p3_k - p3_0)
print P2_0k, P3_0k, 'P2_0k, P3_0k'
print a0/m_parsec, ak/m_parsec

fig, ax1 = plt.subplots(figsize=(5, 4))
aplot   = 10.**np.array([-2.0, -1.5, -1.0, -0.5, 0.0])
pid5p   = np.array([0.256096, 0.219344, 0.18312, 0.14888, 0.118424])
pan5p   = 0.12*aplot**(-5./28.)
ax1.plot(np.log10(aplot), np.log10(pid5p), linestyle = '--', linewidth = 2, color = 'red')
ax1.plot(np.log10(aplot), np.log10(pan5p), linestyle = '--', linewidth = 2, color = 'black')

plt.show()

exit()





a   = 0.5*AU_SI
tau = tH_SI#0.1*sec_year
tc  = (a**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
rp  = (a/2.)*(tau/tc)**(2./7.)
fp  = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp**3.))

print ((1. - 0.8**2.)**(7./2.))*tc/(tH_SI)

print tc
print (768./425.)*(5.*(c_SI**5.)/(512.*G_new_SI**3.))*(a**4.)/(m_SI**3.)
#exit()

print (2.*np.pi*np.sqrt(a**(3.0)/(2.0*G_new_SI*m_SI)))/sec_year
print tbs_SI/sec_year, 'here'
print a_ej_SI/AU_SI
print fp
print np.log10(fp)

print 0.2*AU_SI/(R_sun_SI*(30)**(1./3.))
#print 

Eps = delta_bs**(-10./7.)
print (Eps**(-7./10.) - 1.0)/(Eps**(-17./10.) - 1.0)
print Eps*(Eps**(-7./10.) - 1.0)/(Eps**(-17./10.) - 1.0)




nstar   = (10.**(5.))/(m_parsec**3.)
vdis    = 20.0*(1000.0) 
SMAbbh  = 1.0*AU_SI
Mstar   = 1.0*M_sun_SI
Mbh     = 30.0*M_sun_SI
Rstar   = 1.0*R_sun_SI
Rtde    = Rstar*(Mbh/Mstar)**(1./3.)
Ptde    = 1.0*(2.*Rtde/SMAbbh)
sigbs   = 2.*np.pi*G_new_SI*(2.*Mbh + Mstar)*SMAbbh/(vdis**2.)
sigtde  = Ptde*sigbs
Gtde    = nstar*sigtde*vdis

nrbbhs  = 5.0
nrGCs   = 200.0
Gtde_year_MWgal   = (Gtde*sec_year)*nrbbhs*nrGCs

print Gtde*sec_year
print Gtde_year_MWgal
print Ptde

rcap        = ((85*np.pi/(6.*np.sqrt(2.)))**(2./7.))*(G_new_SI*(m_SI**(2./7.))*(m_SI**(2./7.))*((m_SI + m_SI)**(3./7.)))/((c_SI**(10./7.))*(v_SI**(4./7.)))
print rcap/AU_SI

#exit()
#----------------








#----------------------------------------------------------
#binary-single:
#----------------------------------------------------------
merger_c            = 0
#--------------------------------------
for b in range(0,Nsamp):
#--------------------------------------
    
    #----------------------------------
    #initialize BBH:
    #----------------------------------
    merge_yesno = 0
    id_merger   = -1
    #----------------------------------
    
    #----------------------------------
    #loop over nr hardening steps 'n':
    #----------------------------------
    for n in range(0,Nbs+1):
       
        
        #----------------------
        #General info for this 'n' state:
        #----------------------
        a_n_SI  = a_max_SI*(delta_bs**n)
        #----------------------
        
        
        #----------------------
        #(2-body merger) Isolated state 'n':
        #----------------------
                
        #now define bin a,e:
        abin_SI = a_n_SI                                    #SI
        ebin    = np.sqrt(np.random.random_sample())        #draw from P(e)=2e
                
        #calc binary info:
        #time before next encounter:
        sigma_bs_SI = 2.*np.pi*G_new_SI*(3.*m_SI)*abin_SI/(v_SI**2.)
        Gamma_bs_SI = n_nr_SI*sigma_bs_SI*v_SI
        tbs_SI      = 1./Gamma_bs_SI
        #other params:
        rp      = abin_SI*(1.-ebin)
        fGWbin  = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp**3.))
        #save info:
        bin_info_arr_REAL[tttc,0]  = abin_SI        #SMA a
        bin_info_arr_REAL[tttc,1]  = ebin           #ecc e     
        bin_info_arr_REAL[tttc,2]  = fGWbin         #GW peak  
        bin_info_arr_REAL[tttc,3]  = tbs_SI         #time before next encounter
        
        #check for merger 'before next encounter':
        if (merge_yesno == 0):
            tinsp_e0_SI = (abin_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
            e_IM        = (1.0 - (tbs_SI/tinsp_e0_SI)**(2./7.))**(1./2.)
            #print n, e_IM
            #not merge:
            if (ebin < e_IM):
                bin_info_arr_INT[tttc,0] = 1
                texist_SI   = tbs_SI
            #will merge:
            if (ebin >= e_IM):
                bin_info_arr_INT[tttc,0] = -1
                texist_SI   = tinsp_e0_SI*(1.0-ebin**(2.))**(7./2.)
            #save existance time:
            bin_info_arr_REAL[tttc,4]  = texist_SI      #existence time
            #----------------
            #merger:
            #----------------
            if (ebin >= e_IM):
                #merger, save:
                am  = abin_SI
                em  = ebin
                merge_yesno = 1
                id_merger   = 1                     #MERGER ID 1 (isolated merger)
           #---------------- 
        
        #update counter:
        tttc = tttc+1
        #----------------------
        #----------------------
        
        #----------------------
        #(3-body merger) check for merger 'during encounter':
        #----------------------
        if (merge_yesno == 0):
            Rsch_SI = (2.*G_new_SI*m_SI)/(c_SI**2.)
            rcap_SI = 1.8*Rsch_SI*(abin_SI/Rsch_SI)**(2./7.) #should we add here front factor of 1.8?
            e_cap   = 1.0 - (rcap_SI/abin_SI)
            for rc in range(0,Nims):
                if (merge_yesno == 0):
                    #draw from P(e)=2e
                    e_binsin = np.sqrt(np.random.random_sample())
                    #check for merger:
                    if (e_binsin >= e_cap):
                        #merger, save:
                        am    = abin_SI
                        em    = e_binsin
                        merge_yesno = 1   
                        id_merger   = 2             #MERGER ID 2 (capture merger)
        #----------------------
        #----------------------
                                        
    #loop over nr hardening steps 'n'    
    #----------------------------------
    
    #----------------------------------
    #if the BBH escapes:
    #----------------------------------
    if (merge_yesno == 0):
        tinsp_e0_SI = (a_min_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
        if (tinsp_e0_SI > tH_SI):   e_Ht    = (1.0 - (tH_SI/tinsp_e0_SI)**(2./7.))**(1./2.)
        if (tinsp_e0_SI < tH_SI):   e_Ht    = 0.0
        #draw from P(e)=2e MC:
        e_ej = np.sqrt(np.random.random_sample())
        #check for merger:
        if (e_ej >= e_Ht):
             #merger, save:
             am    = a_min_SI
             em    = e_ej
             merge_yesno    = 1   
             id_merger      = 3                     #MERGER ID 3 (ejected merger)    
    #----------------------------------
    
    #----------------------------------
    #if BBH merged:
    #----------------------------------
    if (merge_yesno == 1):
        
        #calc initial fGW values:
        rp      = am*(1.-em)                                    #rp at formation of merger
        fGW0    = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp**3.)) #fGW at formation of merger
        
        #CASE A: propagate ecc to fGW:
        setA_yesno  = 0
        ecc_A       = 0.0
        a_A         = 0.0
        if (fGW0 < fGW_limit):
            setA_yesno    = 1
            #propagating ini a,e to a,e(fGW_limit)
            c0      = rp*(1.+em)*(em**(-12./19.))*((1.+(121./304.)*em**2.)**(-870./2299.)) #am/((em**(12./19.)/(1.-em**2.))*((1.+(121./304.)*em**2.)**(870./2299.)))
            func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(m_SI+m_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 1e-8
            ecc_A               = fsolve(func, ecc_initial_guess)   #eccentricity at fGW_limit                        
            a_A                 = c0*((ecc_A**(12./19.))/(1.-ecc_A**2.))*((1.+(121./304.)*(ecc_A**2.))**(870./2299.))
        #save info:
        output_INT_arr[merger_c,:]     = [id_merger,setA_yesno,0,0,0,0,0,0,0,0]
        output_REAL_arr[merger_c,:]    = [am,em,fGW0,ecc_A,a_A,0,0,0,0,0]
        
        #update merger counter:
        merger_c = merger_c + 1
    #----------------------------------   
#trim output arrays:
output_INT_arr  = output_INT_arr[0:merger_c,:]
output_REAL_arr = output_REAL_arr[0:merger_c,:]
#--------------------------------------         
#--------------------------------------
#save data:
#--------------------------------------
tf = open(filename + '_BBH_merger_data_INT.txt', "w")
np.savetxt(tf, output_INT_arr,   fmt='%5i')
tf.close()

tf = open(filename + '_BBH_merger_data_REAL.txt', "w")
np.savetxt(tf, output_REAL_arr,   fmt='%5f')
tf.close()
#----------------------------------------------------------
#----------------------------------------------------------



#----------------------------------------------------------
#single-single:
#----------------------------------------------------------
#----------------------------------
#set:
#----------------------------------
Fbinsin     = 0.01
nr_sins_tot = 1./Fbinsin
nr_Rbins    = 1000   #1000
nr_sampl    = 100    #100

minlogR     = -2.
maxlogR     =  2.
#----------------------------------
#calc/define:
#----------------------------------
delogR                      = (maxlogR-minlogR)/(1.*nr_Rbins)
logR_arr                    = np.arange(minlogR, maxlogR, delogR)+(0.5*delogR)
save_fGW_dist_sinsin        = np.zeros((nr_Rbins*nr_sampl,10), dtype=np.float64)
save_a0rp0_PLUMMER_sinsin   = np.zeros((nr_Rbins*nr_sampl,2), dtype=np.float64)
save_eccinfo_temp           = np.zeros((nr_sampl,2), dtype=np.float64)
rsch                        = (2.*G_new_SI*m_SI)/(c_SI**2.)
#----------------------------------
#central properties:
#----------------------------------
n_R     = n_nr_SI   #SI
v_R     = v_SI      #SI
rcap    = (((85.*np.pi)/(24.*np.sqrt(2.)))**(2./7.))*rsch*((c_SI**2./v_R**2.)**(2./7.))
#solve for rpmax:
Ac  = (85.*np.pi/12.)*(G_new_SI**(7./2.))/(c_SI**5.)
Bc  = G_new_SI/2.
Cc  = (768./425.)*(5.*(c_SI**5.)/(512.*(G_new_SI**3.)))*(2.**(7./2.))
Dc  = (1./(6.*np.pi*G_new_SI))*(v_R/n_R) 
func        = lambda rm : Ac*(m_SI**(9./2.))*(rm**(-7./2.)) - (0.5*(m_SI/2.)*(v_R**2.)) - Bc*(m_SI**2.)*(((Cc/Dc)*(1./(m_SI**2.)))**(2./3.))*(rm**(7./3.))
rm_IG       = rcap
rm_solve    = fsolve(func, rm_IG)                        
rpmax       = rm_solve[0]
#calc bmax and sample fGW:
bmax    = np.sqrt((2.*G_new_SI*(2.*m_SI)*rpmax)/(v_R**2.))
#define:
rcap_0  = rcap
rpmax_0 = rpmax
bmax_0  = bmax
#----------------------------------


#----------------------------------
#loop over radial bins:
#----------------------------------
for nR in range(0,nr_Rbins):
    
    print nR, nr_Rbins
    #bin val:
    logR    = logR_arr[nR]
    valR    = 10.**(logR)
    
    #----------------------
    #PLUMMER SPHERE:
    #----------------------
    #n,v:
    n_R     = n_nr_SI*(1.+valR**2.)**(-5./2.)   #SI
    v_R     = v_SI*(1.+valR**2.)**(-1./4.)      #SI
    #calc rcap,bcap,...:
    rcap    = (((85.*np.pi)/(24.*np.sqrt(2.)))**(2./7.))*rsch*((c_SI**2./v_R**2.)**(2./7.))
    #solve for rpmax:
    Ac  = (85.*np.pi/12.)*(G_new_SI**(7./2.))/(c_SI**5.)
    Bc  = G_new_SI/2.
    Cc  = (768./425.)*(5.*(c_SI**5.)/(512.*(G_new_SI**3.)))*(2.**(7./2.))
    Dc  = (1./(6.*np.pi*G_new_SI))*(v_R/n_R) 
    func        = lambda rm : Ac*(m_SI**(9./2.))*(rm**(-7./2.)) - (0.5*(m_SI/2.)*(v_R**2.)) - Bc*(m_SI**2.)*(((Cc/Dc)*(1./(m_SI**2.)))**(2./3.))*(rm**(7./3.))
    rm_IG       = rcap
    rm_solve    = fsolve(func, rm_IG)                        
    rpmax       = rm_solve[0]
    #calc bmax and sample fGW:
    bmax    = np.sqrt((2.*G_new_SI*(2.*m_SI)*rpmax)/(v_R**2.))
    brnd    = bmax*np.sqrt(np.random.random(nr_sampl))
    rp_rnd  = ((brnd**2.)*(v_R**2.))/(2.*G_new_SI*(2.*m_SI))
    fGW_rnd = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp_rnd**3.))
    #CASE A: propagate ecc to fGW:
    for nS in range(0,nr_sampl):
        #encounter properties:
        rp0         = rp_rnd[nS]
        fGW0        = fGW_rnd[nS]
        #a0,e0 right after capture:
        Eorb0       = (1./2.)*(m_SI*m_SI/(m_SI+m_SI))*(v_R**2.) - (85.*np.pi/12.)*(G_new_SI**(7./2.))*(c_SI**(-5.))*(m_SI**(9./2.))*(rp0**(-7./2.))
        a0          = -G_new_SI*m_SI*m_SI/(2.*Eorb0)    #a0 right after capture
        e0          = 1. - rp0/a0                       #e0 right after capture
        save_a0rp0_PLUMMER_sinsin[nR*nr_sampl+nS,0:2]  = [a0, rp0]
        setA_yesno  = 0
        ecc_A       = 0.0        
        if (fGW0 < fGW_limit):
            setA_yesno    = 1   
            #propagating ini a,e to a,e(fGW_limit)
            c0      = rp0*(1.+e0)*(e0**(-12./19.))*((1.+(121./304.)*e0**2.)**(-870./2299.))#    a0/((e0**(12./19.)/(1.-e0**2.))*((1.+(121./304.)*e0**2.)**(870./2299.)))
            func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(m_SI+m_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 1e-8
            ecc_A               = fsolve(func, ecc_initial_guess)   #eccentricity at fGW_limit                        
        #save:
        save_eccinfo_temp[nS,0] = setA_yesno
        save_eccinfo_temp[nS,1] = ecc_A
    #calc rate,weight fac,...:
    sigma_R     = (2.*np.pi*G_new_SI)*(2.*m_SI*rpmax/(v_R**2.))     #cross section
    rate_ps     = n_R*sigma_R*v_R                                   #rate per single
    nr_sins_R   = 3.*nr_sins_tot*(valR**2.)*((1.+valR**2.)**(-5./2.))*(valR*delogR/np.log10(np.e))
    rate_tot    = (rate_ps*nr_sins_R)/2.
    weight_fac  = (1./(1.*nr_sampl))*rate_tot
    #save:
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 0]   =  fGW_rnd
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 1]   =  rate_tot
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 2]   =  weight_fac
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 3]   =  save_eccinfo_temp[:,0]    #setA_yesno
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 4]   =  save_eccinfo_temp[:,1]    #ecc_A
    #----------------------
    
    #----------------------
    #UNIFORM SPHERE:
    #----------------------
    #n,v:
    n_R     = n_nr_SI   #SI
    v_R     = v_SI      #SI
    #calc rcap,bcap,...:
    #define:
    rcap    = rcap_0
    rpmax   = rpmax_0
    bmax    = bmax_0
    brnd    = bmax*np.sqrt(np.random.random(nr_sampl))
    rp_rnd  = ((brnd**2.)*(v_R**2.))/(2.*G_new_SI*(2.*m_SI))
    fGW_rnd = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp_rnd**3.))
    #CASE A: propagate ecc to fGW:
    for nS in range(0,nr_sampl):
        #encounter properties:
        rp0         = rp_rnd[nS]
        fGW0        = fGW_rnd[nS]
        #a0,e0 right after capture:
        Eorb0   = (1./2.)*(m_SI*m_SI/(m_SI+m_SI))*(v_R**2.) - (85.*np.pi/12.)*(G_new_SI**(7./2.))*(c_SI**(-5.))*(m_SI**(9./2.))*(rp0**(-7./2.))
        a0      = -G_new_SI*m_SI*m_SI/(2.*Eorb0)    #a0 right after capture
        e0      = 1. - rp0/a0                       #e0 right after capture  
        setA_yesno  = 0
        ecc_A       = 0.0
        if (fGW0 < fGW_limit):
            setA_yesno    = 1
            #propagating ini a,e to a,e(fGW_limit)
            c0      = rp0*(1.+e0)*(e0**(-12./19.))*((1.+(121./304.)*e0**2.)**(-870./2299.))#    a0/((e0**(12./19.)/(1.-e0**2.))*((1.+(121./304.)*e0**2.)**(870./2299.)))
            func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(m_SI+m_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 1e-8
            ecc_A               = fsolve(func, ecc_initial_guess)   #eccentricity at fGW_limit                      
        #save:
        save_eccinfo_temp[nS,0] = setA_yesno
        save_eccinfo_temp[nS,1] = ecc_A
    #calc rate,weight fac,...:
    sigma_R     = (2.*np.pi*G_new_SI)*(2.*m_SI*rpmax/(v_R**2.))     #cross section
    rate_ps     = n_R*sigma_R*v_R                                   #rate per single
    nr_sins_R   = 3.*nr_sins_tot*(valR**2.)*(valR*delogR/np.log10(np.e))
    rate_tot    = (rate_ps*nr_sins_R)/2.
    if (valR <  1.0):    weight_fac  = (1./(1.*nr_sampl))*rate_tot
    if (valR >= 1.0):    weight_fac  = 0.0
    #save:
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 5]   =  fGW_rnd
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 6]   =  rate_tot
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 7]   =  weight_fac
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 8]   =  save_eccinfo_temp[:,0]    #setA_yesno
    save_fGW_dist_sinsin[(nR+0)*nr_sampl:(nR+1)*nr_sampl, 9]   =  save_eccinfo_temp[:,1]    #ecc_A
    #----------------------
    
#----------------------------------
#----------------------------------------------------------
#----------------------------------------------------------











#TEST:
#pos_mN      = np.where(bin_info_arr_INT[:,0] == 1)[0]
#a_arr       = bin_info_arr_REAL[pos_mN,0]
#wfac_arr    = bin_info_arr_REAL[pos_mN,3]
#T_arr       = 365.*(2.*np.pi*np.sqrt(a_arr**3./(G_new_SI*(2.*m_SI)))/sec_year) 
#logT_arr    = np.log10(T_arr)
#fig, ax1 = plt.subplots(figsize=(5, 4))
#ax1.hist(logT_arr, weights = wfac_arr, range=(-2.0,3.0), bins=50, histtype='step', fill=False, color = 'red',   linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')
#plt.show()
#exit()
 



#exit()

#--------------------------------------
#analyze data:
#--------------------------------------
pos_id1 = np.where(output_INT_arr[:,0] == 1)[0] # = IM      (1 - isolated merger)
pos_id2 = np.where(output_INT_arr[:,0] == 2)[0] # = capM    (2 - capture merger)
pos_id3 = np.where(output_INT_arr[:,0] == 3)[0] # = ejM     (3 - ejected merger)

pos_A   = np.where(output_INT_arr[:,1] == 1)[0]

pos_id1_A   = list(set(pos_id1).intersection(pos_A))
pos_id2_A   = list(set(pos_id2).intersection(pos_A))
pos_id3_A   = list(set(pos_id3).intersection(pos_A))

fGW0_all    = output_REAL_arr[:,2]
fGW0_id1    = fGW0_all[pos_id1] 
fGW0_id2    = fGW0_all[pos_id2] 
fGW0_id3    = fGW0_all[pos_id3] 

ecc_all_A   = output_REAL_arr[pos_A,3]
ecc_id1_A   = output_REAL_arr[pos_id1_A,3]
ecc_id2_A   = output_REAL_arr[pos_id2_A,3]
ecc_id3_A   = output_REAL_arr[pos_id3_A,3]

a_all_A   = output_REAL_arr[pos_A,4]
a_id1_A   = output_REAL_arr[pos_id1_A,4]
a_id2_A   = output_REAL_arr[pos_id2_A,4]
a_id3_A   = output_REAL_arr[pos_id3_A,4]

#weight at steady state:
T_aHB_aej       = (1./n_nr_SI)*(v_SI/(1.-delta_bs))*(1./(6.*np.pi*G_new_SI*m_SI))*(1./a_ej_SI)
rate_binform    = 1./T_aHB_aej
prob_2merg      = (1.*len(pos_id1))/(1.*Nsamp)
rate_2merg      = rate_binform*prob_2merg
prob_3merg      = (1.*len(pos_id2))/(1.*Nsamp)
rate_3merg      = rate_binform*prob_3merg

#create pos arrays:
pos_mN              = np.where(bin_info_arr_INT[:,0] ==  1)[0]
pos_mY              = np.where(bin_info_arr_INT[:,0] == -1)[0]
pos_mNY             = np.array(list(pos_mN)+list(pos_mY))
pos_e99cut          = np.where(bin_info_arr_REAL[:,1] < 0.99)[0]
pos_e9cut           = np.where(bin_info_arr_REAL[:,1] < 0.9)[0]

#no GR effets:
arr_all_a_SI        = bin_info_arr_REAL[:,0]
arr_all_e_SI        = bin_info_arr_REAL[:,1]
arr_all_fGW_SI      = bin_info_arr_REAL[:,2]
arr_all_texist_SI   = bin_info_arr_REAL[:,3]
arr_all_weightfac   = (arr_all_texist_SI/tH_SI)

#calc:
arr_all_Torb_SI     = (2.*np.pi)*np.sqrt(arr_all_a_SI**3./(G_new_SI*(m_Msun+m_Msun)*M_sun_SI))
Torb_SI_HB          = (2.*np.pi)*np.sqrt(a_HB_SI**3./(G_new_SI*(m_Msun+m_Msun)*M_sun_SI))
arr_all_Torb_dys    = arr_all_Torb_SI/(60.*60.*24.)
arr_all_Torb_THB    = arr_all_Torb_SI/Torb_SI_HB

#with GR effects:
arr_mN_a_SI         = bin_info_arr_REAL[pos_mN,0]
arr_mN_e_SI         = bin_info_arr_REAL[pos_mN,1]
arr_mN_fGW_SI       = bin_info_arr_REAL[pos_mN,2]
arr_mN_texist_SI    = bin_info_arr_REAL[pos_mN,4]
arr_mN_weightfac    = (arr_mN_texist_SI/tH_SI)

arr_mY_a_SI         = bin_info_arr_REAL[pos_mY,0]
arr_mY_e_SI         = bin_info_arr_REAL[pos_mY,1]
arr_mY_fGW_SI       = bin_info_arr_REAL[pos_mY,2]
arr_mY_texist_SI    = bin_info_arr_REAL[pos_mY,4]
arr_mY_weightfac    = (arr_mY_texist_SI/tH_SI)

arr_mNY_a_SI        = bin_info_arr_REAL[pos_mNY,0]
arr_mNY_e_SI        = bin_info_arr_REAL[pos_mNY,1]
arr_mNY_fGW_SI      = bin_info_arr_REAL[pos_mNY,2]
arr_mNY_texist_SI   = bin_info_arr_REAL[pos_mNY,4]
arr_mNY_weightfac   = (arr_mNY_texist_SI/tH_SI)
#--------------------------------------
#--------------------------------------



#--------------------------------------
#TEST FOR KAZE PROJ:
#--------------------------------------
#input:
nr_x_bins   = 100   #t binning
logt_min    = -10   #log yr
logt_max    = 10    #log yr
#set:
pos_inC     = np.array(list(pos_id1) + list(pos_id2))
nr_inC      = len(pos_inC)
#calc:
dlogt       = (logt_max-logt_min)/(1.*nr_x_bins)
binarr_logt = np.arange(logt_min,logt_max,dlogt) + dlogt/2.
binarr_obsW = (10.**(binarr_logt + dlogt/2.)) - (10.**(binarr_logt - dlogt/2.))
#define:
arr_obsINFO = np.zeros((nr_inC, nr_x_bins,5), dtype=np.float64)
tC_yr = ((768./425.)*(5.*(c_SI**5.)/(512.*G_new_SI**3.)))/sec_year
#loop:
for nIC in range(0,nr_inC):
    i_inC   = pos_inC[nIC]
    #print info:
    print nIC, nr_inC
    #initial val:
    a0      = output_REAL_arr[i_inC,0]
    e0      = output_REAL_arr[i_inC,1]
    c0      = a0*(1-(e0**2.))*(e0**(-12./19.))*(1. + (121./304.)*(e0**2.))**(-870./2299.)
    t0_yr   = (tC_yr*(a0**4.)/(m_SI**3.))*((1.-e0**2.)**(7./2.))
    #calc ecc at each tbin:
    posmax      = max(np.where(binarr_logt[:] < np.log10(t0_yr))[0])
    eguess      = e0
    for nb in range(posmax,-1,-1):  #counts backwards from posmax to 0.
        logtnb_yr           = binarr_logt[nb]
        func                = lambda ecc : logtnb_yr - np.log10(tC_yr*((((c0*(ecc**(12./19.))/(1-(ecc**2.)))*(1. + (121./304.)*(ecc**2.))**(870./2299.))**(4.))/(m_SI**3.))*((1.-ecc**2.)**(7./2.)))
        ecc_initial_guess   = eguess
        ecc_nb              = fsolve(func, ecc_initial_guess)
        #update:
        eguess  = ecc_nb
        #save:
        arr_obsINFO[nIC, nb, 0] = logtnb_yr
        arr_obsINFO[nIC, nb, 1] = ecc_nb
        arr_obsINFO[nIC, nb, 2] = binarr_obsW[nb]

fig, ax1 = plt.subplots(figsize=(5,4))
for nIC in range(0,nr_inC):
    posmax  = max(np.where(arr_obsINFO[nIC,:,2] > 0.0)[0])
    ax1.plot(arr_obsINFO[nIC,0:posmax+1,0], arr_obsINFO[nIC,0:posmax+1,1], color='black', linestyle = "-")

plt.show()
#--------------------------------------












#--------------------------------------
#PLOT: fGW dist:
#--------------------------------------
#-------------------
#fGW at formation:
#-------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

Hmin    = -4
Hmax    = 2   
Hnrbins = 80

#Binary-Single:
#2-body mergers:
data_arr    = np.log10(fGW0_id1)
weight_fac  = (1./(1.*len(data_arr)))*rate_2merg
weight_arr  = np.full(len(data_arr),weight_fac)
ax1.hist(data_arr[:],   weights = weight_arr[:], range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=True, facecolor = 'forestgreen', ec = 'black', lw=1, alpha=0.5, label='2-body')
#3-body mergers:
data_arr    = np.log10(fGW0_id2)
weight_fac  = (1./(1.*len(data_arr)))*rate_3merg
weight_arr  = np.full(len(data_arr),weight_fac)
ax1.hist(data_arr[:],   weights = weight_arr[:], range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=True, facecolor = 'red', ec = 'black', lw=1, alpha=0.5, label='3-body')

#Single-Single:
#PLUMMER:
data_arr    = np.log10(save_fGW_dist_sinsin[:,0])
weight_arr  = save_fGW_dist_sinsin[:,2]
ax1.hist(data_arr[:],   weights = weight_arr[:], range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=True, facecolor = 'grey',       ec = 'black', lw=1, linestyle = '-', alpha=0.50,    label='sin-sin (Plummer)')
#UNIFORM:
data_arr    = np.log10(save_fGW_dist_sinsin[:,5])
weight_arr  = save_fGW_dist_sinsin[:,7]
ax1.hist(data_arr[:],   weights = weight_arr[:], range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, facecolor = 'lightgrey',  ec = 'black', lw=1, linestyle = ':', alpha=0.50,   label='sin-sin (uniform)')
    

ax1.set_xlabel(r'log $f_{\rm GW}$ [Hz]')
ax1.set_ylabel(r'Rate $\Gamma$ [rnd. norm.]')
ax1.set_title(r'GW frequency $f_{\rm GW}$ at formation')
ax1.set_yticklabels([])
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
fig.savefig(filename + '_hist_fGW0_BBH_bs_sinsin.pdf', bbox_inches='tight')
plt.show()
#--------------------------------------
#--------------------------------------


#--------------------------------------
#PLOT: ecc at X Hz:
#--------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 3))

Hmin    = -3.0
Hmax    =  0.0   
Hnrbins =  40

#Binary-Single:
#2-body mergers:
data_id1        = np.array([])
weight_id1      = np.array([])
data_arr        = np.log10(ecc_id1_A)
if (len(data_arr) > 0):
    weight_fac  = (1./(1.*len(data_arr)))*(1.*len(data_arr)/(1.*len(pos_id1)))*rate_2merg
    weight_arr  = np.full(len(data_arr),weight_fac)
    ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, color = 'forestgreen', linewidth = 3, alpha=0.95, label='2-body', zorder=40)
    data_id1    = data_arr
    weight_id1  = weight_arr   
#3-body mergers:
data_id2        = np.array([])
weight_id2      = np.array([])
data_arr        = np.log10(ecc_id2_A)
if (len(data_arr) > 0):
    weight_fac  = (1./(1.*len(data_arr)))*(1.*len(data_arr)/(1.*len(pos_id2)))*rate_3merg
    weight_arr  = np.full(len(data_arr),weight_fac)
    ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, color = 'red', linewidth = 3, alpha=0.95, label='3-body', zorder=30)
    data_id2    = data_arr
    weight_id2  = weight_arr 

#Single-Single:
#PLUMMER:
data_PLUM       = np.array([])
weight_PLUM     = np.array([])
setA_yesno_arr  = save_fGW_dist_sinsin[:,3]
setA_yesno_arr  = np.array([int(i) for i in setA_yesno_arr])
pos             = np.where(setA_yesno_arr == 1)[0]
setA_ecc_arr    = save_fGW_dist_sinsin[pos,4]
data_arr        = np.log10(setA_ecc_arr)
weight_arr      = save_fGW_dist_sinsin[pos,2]
ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, color = 'grey', linewidth = 3, alpha=0.95, label='sin-sin (Plummer)', zorder=20)
data_PLUM       = data_arr
weight_PLUM     = weight_arr 
#UNIFORM:
#setA_yesno_arr  = save_fGW_dist_sinsin[:,8]
#setA_yesno_arr  = np.array([int(i) for i in setA_yesno_arr])
#pos             = np.where(setA_yesno_arr == 1)[0]
#setA_ecc_arr    = save_fGW_dist_sinsin[pos,9]
#data_arr        = np.log10(setA_ecc_arr)
#weight_arr      = save_fGW_dist_sinsin[pos,7]
#ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, facecolor = 'lightgrey',  ec = 'black', lw=1, linestyle = ':', alpha=0.50,   label='sin-sin (uniform)')

#2-body + 3-body:
data_arr        = np.array(list(data_id1)+list(data_id2))
weight_arr      = np.array(list(weight_id1)+list(weight_id2))
ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=True, facecolor = 'steelblue', ec = 'black', lw=0, linestyle = ':', alpha=0.75, label='2-body + 3-body', zorder=0)

#All:
data_arr        = np.array(list(data_id1)+list(data_id2)+list(data_PLUM))
weight_arr      = np.array(list(weight_id1)+list(weight_id2)+list(weight_PLUM))
ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=True, facecolor = 'dodgerblue', ec = 'black', lw=2, linestyle = ':', alpha=0.25, label='All', zorder=10)
ax1.hist(data_arr[:], weights = weight_arr[:],  range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False,facecolor = 'dodgerblue', ec = 'black', lw=2, linestyle = ':', alpha=1.00, zorder=50)


ax1.set_xlim([Hmin, Hmax])
ax1.set_xlabel(r'log $e$')
ax1.set_ylabel(r'Rate $\Gamma$ [rnd. norm.]')
#ax1.set_title(r'ecc. dist. at $f_{\rm GW} = 1$ Hz')
ax1.set_title(r'ecc. dist. at $f_{\rm GW} = 10$ Hz')
ax1.set_yticklabels([])
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 2, frameon = False)

#save/show fig:
fig.savefig(filename + '_hist_eccfGW_BBHdm_sinsin.pdf', bbox_inches='tight')
plt.show()
#--------------------------------------
#--------------------------------------





#--------------------------------------
#SAVE DATA FOR DAN's PLOT:
#--------------------------------------
#This routine IS A MESS, but it should work.
#you have to choose the right number of initical samplings
#to make sure there is enough numbers to take from ... 
savedata = 0
if (savedata == 1):
    #weights and numbers:
    norm_weightarr_PLUMMER      = save_fGW_dist_sinsin[:,2]/(max(save_fGW_dist_sinsin[:,2]))
    rnd_nr_arr                  = np.random.random(nr_Rbins*nr_sampl)
    PLUMMER_weighted_dist_arp   = np.zeros((nr_Rbins*nr_sampl,2), dtype=np.float64) 
    sc = 0
    for i in range(nr_Rbins*nr_sampl):
        weight_i    = norm_weightarr_PLUMMER[i]
        rnd_i       = rnd_nr_arr[i]
        if (rnd_i < weight_i):
            PLUMMER_weighted_dist_arp[sc,0] = save_a0rp0_PLUMMER_sinsin[i,0]
            PLUMMER_weighted_dist_arp[sc,1] = save_a0rp0_PLUMMER_sinsin[i,1]
            sc = sc+1
    #trim array:
    PLUMMER_weighted_dist_arp = PLUMMER_weighted_dist_arp[0:sc,:]
    #relative weights and numbers:
    rate_P      = sum(save_fGW_dist_sinsin[:,2])
    ratew2      = rate_2merg/rate_P
    ratew3      = rate_3merg/rate_P
    ratewP      = rate_P/rate_P
    nrP         = len(PLUMMER_weighted_dist_arp[:,0])
    nrsave_2    = int(nrP*ratew2)
    nrsave_3    = int(nrP*ratew3)
    nrsave_P    = int(nrP*ratewP)
    totnrsave   = nrsave_2 + nrsave_3 + nrsave_P
    print nrsave_2, len(pos_id1)
    print nrsave_3, len(pos_id2)
    print nrsave_P
    #data:
    SMA_a_id2   = output_REAL_arr[pos_id1[0:nrsave_2],0]
    ECC_e_id2   = output_REAL_arr[pos_id1[0:nrsave_2],1]
    SMA_a_id3   = output_REAL_arr[pos_id2[0:nrsave_3],0]
    ECC_e_id3   = output_REAL_arr[pos_id2[0:nrsave_3],1]
    SMA_a_PLU   = PLUMMER_weighted_dist_arp[:,0]
    rPC_PLU     = PLUMMER_weighted_dist_arp[:,1]
    #define:
    rPC_id2     = SMA_a_id2*(1.-ECC_e_id2)
    rPC_id3     = SMA_a_id3*(1.-ECC_e_id3)
    #save:
    saveinfoarr = np.zeros((totnrsave,3), dtype=np.float64)
    saveinfoarr[0:nrsave_2, 0]                      = 2
    saveinfoarr[0:nrsave_2, 1]                      = SMA_a_id2[:]
    saveinfoarr[0:nrsave_2, 2]                      = rPC_id2[:]
    saveinfoarr[nrsave_2:(nrsave_2 + nrsave_3), 0]  = 3
    saveinfoarr[nrsave_2:(nrsave_2 + nrsave_3), 1]  = SMA_a_id3[:]
    saveinfoarr[nrsave_2:(nrsave_2 + nrsave_3), 2]  = rPC_id3[:]
    saveinfoarr[(nrsave_2 + nrsave_3):totnrsave, 0] = 1
    saveinfoarr[(nrsave_2 + nrsave_3):totnrsave, 1] = SMA_a_PLU[:]
    saveinfoarr[(nrsave_2 + nrsave_3):totnrsave, 2] = rPC_PLU[:]
    #save:
    tf = open('dansdata_id_a_rp.txt', "w")
    np.savetxt(tf, saveinfoarr,   fmt='%8f')
    tf.close()
    #test fig:
    fig, ax1 = plt.subplots(figsize=(5, 4))
    ax1.hist(np.log10(rPC_id2**(-3./2.)), bins=50, alpha = 0.25)
    ax1.hist(np.log10(rPC_id3**(-3./2.)), bins=50, alpha = 0.25)
    ax1.hist(np.log10(rPC_PLU**(-3./2.)), bins=50, alpha = 0.25)
    plt.show()
    exit()
#--------------------------------------


#CHECK ALL THINGS AGAIN!!!!! REPLOT!!!!!


exit()






#--------------------------------------
#PLOTs: PAPER 1.
#--------------------------------------
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#-------------------
#fGW at formation:
#-------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

Hmin    = -5
Hmax    = 2   
Hnrbins = 100


#binary-single:
ax1.hist(np.log10(fGW0_id3),      range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, color = 'royalblue',        label='ejected mergers (outside cluster)', hatch = '\\\\\\')
ax1.hist(np.log10(fGW0_id1),      range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, color = 'forestgreen',      label='2-body mergers (inside cluster)', hatch = '/////')
ax1.hist(np.log10(fGW0_id2),      range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, color = 'red',              label='3-body mergers (inside cluster)', hatch = '|||||')


#single-single:
#rmin: cap
rcap        = ((85*np.pi/(6.*np.sqrt(2.)))**(2./7.))*(G_new_SI*(m_SI**(2./7.))*(m_SI**(2./7.))*((m_SI + m_SI)**(3./7.)))/((c_SI**(10./7.))*(v_SI**(4./7.)))
fcap        = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rcap**3.))
print np.log10(fcap)
#rmin: merger
Ac  = (85.*np.pi/12.)*(G_new_SI**(7./2.))/(c_SI**5.)
Bc  = G_new_SI/2.
Cc  = (768./425.)*(5.*(c_SI**5.)/(512.*(G_new_SI**3.)))*(2.**(7./2.))
Dc  = (1./(6.*np.pi*G_new_SI))*(v_SI/n_nr_SI) 
func        = lambda rm : Ac*(m_SI**(9./2.))*(rm**(-7./2.)) - (0.5*(m_SI/2.)*(v_SI**2.)) - Bc*(m_SI**2.)*(((Cc/Dc)*(1./(m_SI**2.)))**(2./3.))*(rm**(7./3.))
rm_IG       = rcap
rm_solve    = fsolve(func, rm_IG)                        
fm          = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rm_solve**3.))
#Eorb        = (1./2.)*(m_SI/2.)*(v_SI**2.) - (85.*np.pi/12.)*(G_new_SI**(7./2.))*(c_SI**(-5.))*(m_SI**(9./2.))*(rm_solve**(-7./2.))
#am          = -G_new_SI*m_SI*m_SI/(2.*Eorb)
#print am/AU_SI
#print rcap/AU_SI, rm_solve/AU_SI
#exit()


print np.log10(fm)
#choose fmin:
fmin = fm
#make data and plot single-single dist:
logf_arr    = np.arange(np.log10(fmin), Hmax, ((Hmax - Hmin)/(1.*Hnrbins)))
fbs_1       = 0.02#1./(6.*Nims) #0.02
normfac     = 0.3*45.*(Nsamp/10000.)
plt.fill_between(logf_arr, normfac*(1./(fbs_1*(6.*Nims)))*(10.**((-2./3.)*logf_arr)), color="none", hatch="xxxx", edgecolor="grey", linewidth=1.0, step="post", label = r'sin-sin ($f_{\rm bs} = 2\%$)')


#print info:
print len(fGW0_id1)
print len(fGW0_id2)
print len(fGW0_id3)


ax1.set_xlabel(r'log $f_{\rm GW}$ [Hz]')
ax1.set_ylabel(r'Formation Rate $\Gamma$ [rnd. norm.]')
#ax1.set_title(r'GW frequency $f_{\rm GW}$ at formation')
ax1.set_yticklabels([])
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
#fig.savefig('hist_fGW0_BBHdm.pdf', bbox_inches='tight')

plt.show()
exit()
#-------------------



#-------------------
#ecc at X Hz:
#-------------------
#histogram:
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(ecc_id3_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'royalblue',    label='ejected mergers (outside cluster)', hatch = '\\\\\\\\')
ax1.hist(np.log10(ecc_id1_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'forestgreen',   label='2-body mergers (inside cluster)', hatch = '//////')

ax1.set_xlim([-3.5, 0.0])
ax1.set_xlabel(r'log $e$')
ax1.set_ylabel(r'nr events')
ax1.set_title(r'Eccentricity $e$ at $f_{\rm GW} = 10^{-2}$ Hz')
ax1.set_yticklabels([])
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
fig.savefig('hist_eccfGW_BBHdm.eps', bbox_inches='tight')
plt.show()



#cumulative histogram:
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(ecc_id3_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'royalblue',    cumulative=True, normed=1, label='ejected mergers', hatch = '\\\\\\\\')
ax1.hist(np.log10(ecc_all_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'black',   cumulative=True, normed=1, label='ejected + 2-body mergers', hatch = '//////', linewidth=2.0)

ax1.set_xlim([-3.5, 0.0])
ax1.set_ylim([0.0, 1.1])
ax1.set_xlabel(r'log $e$')
ax1.set_ylabel(r'cumulative nr events')
ax1.set_title(r'Eccentricity $e$ at $f_{\rm GW} = 10^{-2}$ Hz')
plt.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
fig.savefig('cumuhist_eccfGW_BBHdm.eps', bbox_inches='tight')
plt.show()
#-------------------
#--------------------------------------


exit()






#--------------------------------------
#Figure 2:
#--------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(arr_all_Torb_THB[:]),         weights = arr_all_weightfac[:],         range=(-5.0,0.0), bins=50, histtype='step', fill=False, color = 'red',   linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')
ax1.hist(np.log10(arr_all_Torb_THB[pos_mN]),    weights = arr_all_weightfac[pos_mN],    range=(-5.0,0.0), bins=50, histtype='step', fill=False, color = 'black', linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')
#ax1.hist(np.log10(arr_all_Torb_THB[pos_mY]),    weights = arr_all_weightfac[pos_mY],  range=(-5.0,0.0), bins=50, histtype='step', fill=False, color = 'purple', linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')

tbs_HB_SI   = 1./(n_nr_SI*(2.*np.pi*G_new_SI*(3.*m_SI)*a_HB_SI/(v_SI**2.))*v_SI)
tcl_HB_SI   = (a_HB_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
Pm_HB       = (tbs_HB_SI/tcl_HB_SI)**(2./7.)
Eps         = delta_bs**(-10./7.)

arr_TTHB    = 10.**(np.arange(-5.0,0.0,0.01))
Pnm_TTHB    = np.exp(-Pm_HB*((1.0 - Eps*arr_TTHB**(-20./21.))/(1.0 - Eps)))
PTD_GRno    = 0.0375*(arr_TTHB**(-2./3.))
PTD_GRyes   = PTD_GRno*Pnm_TTHB

ax1.plot(np.log10(arr_TTHB), PTD_GRno, linestyle = '--', linewidth = 2, color = 'red', label = r'analytical')
ax1.plot(np.log10(arr_TTHB), PTD_GRyes, linestyle = '--', linewidth = 2, color = 'blue', label = r'analytical')

fed         = vesc_kms/v_kms
Tej_THB     = ((2./63.)**(3./2.))*(fed**(-3.))
ax1.plot([np.log10(Tej_THB), np.log10(Tej_THB)], [0.0,100.0], linestyle = ':', linewidth = 1, color = 'grey', label = r'analytical')

fed         = 50.0/20.0
Tej_THB     = ((2./63.)**(3./2.))*(fed**(-3.))
ax1.plot([np.log10(Tej_THB), np.log10(Tej_THB)], [0.0,100.0], linestyle = ':', linewidth = 1, color = 'grey', label = r'analytical')


ax1.set_xlim([-5.0, 0.0])
ax1.set_ylim([0.0, 100.0])
#ax1.set_yticklabels([])

#ax1.set_yscale('log')

ax1.set_xlabel(r'log $f$ [Hz]')
ax1.set_ylabel(r'$d\Gamma_{\rm TD}/d\log(T/T_{\rm HB})$ [rnd. norm.]')
ax1.set_title(r'GW Frequency Distributions')
plt.legend(loc='upper right', numpoints = 1, fontsize = 9.0, ncol = 2, frameon = True)

#save/show fig:
fig.savefig(filename+'_TDE_Torb.eps', bbox_inches='tight')

plt.show()

exit()
#--------------------------------------










#--------------------------------------
#Figure 1:
#--------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

tbs_HB_SI   = 1./(n_nr_SI*(2.*np.pi*G_new_SI*(3.*m_SI)*a_HB_SI/(v_SI**2.))*v_SI)
tcl_HB_SI   = (a_HB_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
p_HB        = (tbs_HB_SI/tcl_HB_SI)**(2./7.)
Eps         = delta_bs**(-10./7.)
i_arr       = 1.0*np.arange(0,Nbs+1,0.001)
PnM_i       = np.exp(-p_HB*((1.-Eps**(i_arr + 1))/(1.-Eps)))
#PnM_i       = 1.0 + (-p_HB*((1.-Eps**(i_arr + 1))/(1.-Eps)))   #provides ok simple fit. just taylor of exp(x)

y_hist  = np.log10(arr_mN_a_SI/a_max_SI)
x_bins  = np.log10(delta_bs**(np.arange(0.0, Nbs+2,1) - 0.5))[::-1]
w_bins  = (1./(1.*Nsamp)) + np.zeros(len(y_hist), dtype=np.float64)

ax1.plot([max(x_bins), min(x_bins)], [1.0, 1.0], linestyle = ':', linewidth = 2, color = 'black', label = r'- GR evol.')
ax1.hist(y_hist, bins=x_bins, weights = w_bins, range=(-10,10), histtype='stepfilled', color='grey', label = r'+ GR evol.')
ax1.plot(np.log10(delta_bs**i_arr), PnM_i, linestyle = '--', linewidth = 2, color = 'red', label = r'analytical')

ax1.set_xlim([-3.0, 0.25])
ax1.set_ylim([0.0, 1.5])
ax1.set_xlabel(r'log $a/a_{\rm HB}$')
ax1.set_ylabel(r'Histogram')
plt.legend(loc='upper right', numpoints = 1, fontsize = 9.0, ncol = 3, frameon = False)

#save/show fig:
fig.savefig(filename+'_dist_SMA_mN.eps', bbox_inches='tight')

plt.show()

exit()
#--------------------------------------








#--------------------------------------
#TEST!
#--------------------------------------
print np.log10((1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(a_min_SI**3.)))
print np.log10((1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(a_max_SI**3.)))


fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(arr_all_fGW_SI),    range=(-11,1), bins=200, histtype='step', fill=False, color = 'red', linestyle = '-', linewidth = 2, alpha = 0.5)
ax1.hist(np.log10(arr_mN_fGW_SI),     range=(-11,1), bins=200, histtype='step', fill=False, color = 'black',  linestyle = '--',   linewidth = 2, alpha = 1.0)
ax1.hist(np.log10(arr_mY_fGW_SI),     range=(-11,1), bins=200, histtype='step', fill=False, color = 'green',  linestyle = '--',   linewidth = 2, alpha = 1.0)
ax1.hist(np.log10(arr_mNY_fGW_SI),    range=(-11,1), bins=200, histtype='step', fill=False, color = 'blue',  linestyle = '--',   linewidth = 2, alpha = 1.0)

xvals   = np.arange(-11,1,0.01)
yvals   = 100.*0.1*(10.**(xvals*(-2./3.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 2, color = 'red', label = r'analytical')
xvals   = np.arange(-11,1,0.01)
yvals   = 100.*0.0000003*(10.**(xvals*(-20./9.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 2, color = 'red', label = r'analytical')


ax1.set_yscale('log')

ax1.set_xlim([-6.0, 0.0])
ax1.set_ylim([1e-3, 1e5])

#save/show fig:
fig.savefig(filename+'_dist_fGW.eps', bbox_inches='tight')

plt.show()

#exit()
#--------------------------------------









#--------------------------------------
#Figure 2:
#--------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(arr_all_fGW_SI),  weights = arr_all_weightfac,    range=(-11,1), bins=200, histtype='step', fill=False, color = 'orangered',  linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')
ax1.hist(np.log10(bin_info_arr_REAL[pos_e99cut,2]), weights = bin_info_arr_REAL[pos_e99cut,3]/tH_SI, range=(-11,1), bins=200, histtype='step', fill=False, color = 'orangered', linestyle = ':', linewidth = 1, alpha = 1.0, label = r'Newt., e<0.99')
#ax1.hist(np.log10(bin_info_arr_REAL[pos_e9cut,2]), weights = bin_info_arr_REAL[pos_e9cut,3]/tH_SI, range=(-11,1), bins=100, histtype='step', fill=False, color = 'orange', linestyle = '-', linewidth = 1, alpha = 0.5)

ax1.hist(np.log10(arr_mY_fGW_SI),   weights = arr_mY_weightfac,     range=(-11,1), bins=200, histtype='step', fill=False, color = 'darkgray',  linestyle = '-',    linewidth = 2, alpha = 1.0, hatch = '////', label = r'GR, merg. pop.')
ax1.hist(np.log10(arr_mN_fGW_SI),   weights = arr_mN_weightfac,     range=(-11,1), bins=200, histtype='step', fill=False, color = 'royalblue',  linestyle = '-',    linewidth = 2, alpha = 1.0, hatch = '\\\\\\\\', label = r'GR, non-merg. pop.')
#ax1.hist(np.log10(arr_mNY_fGW_SI),  weights = arr_mNY_weightfac,    range=(-11,1), bins=200, histtype='step', fill=False, color = 'black',      linestyle = '-',    linewidth = 2, alpha = 1.0, label = r'+GR, all.')

xvals   = np.arange(-11,1,0.01)
yvals   = 0.1*0.2*(10.**(xvals*(-2./3.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 1, color = 'orangered', label = r'$\propto f^{-2/3}$')
xvals   = np.arange(-11,1,0.01)
yvals   = 0.1*0.000000000002*(10.**(xvals*(-34./9.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 1, color = 'royalblue', label = r'$\propto f^{-34/9}$')
xvals   = np.arange(-11,1,0.01)
yvals   = 0.1*0.000000003*(10.**(xvals*(-3.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 1, color = 'darkgray', label = r'$\propto f^{-3}$')

ax1.set_yscale('log')

ax1.set_xlim([-6.0, -1.0])
ax1.set_ylim([1e-2, 1e3])
#ax1.set_yticklabels([])

ax1.set_xlabel(r'log $f$ [Hz]')
ax1.set_ylabel(r'number events [rnd. norm.]')
ax1.set_title(r'GW Frequency Distributions')
plt.legend(loc='upper right', numpoints = 1, fontsize = 9.0, ncol = 2, frameon = True)

#save/show fig:
fig.savefig(filename+'_obsdist_fGW.eps', bbox_inches='tight')

plt.show()

exit()
#--------------------------------------









fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.hist(arr_mN_e_SI, range=(0,1), bins=100, histtype='step', fill=False)
plt.show()

fig, ax = plt.subplots(figsize=(5, 4))
plt.plot(np.log10(arr_mN_a_SI/AU_SI), np.log10(arr_mN_fGW_SI), marker='o', linestyle = '', fillstyle='full', markeredgewidth=0.0, markersize=2, alpha=0.75, c='black')
plt.show()

fig, ax = plt.subplots(figsize=(5, 4))
plt.plot(np.log10(arr_mN_a_SI/AU_SI), arr_mN_e_SI, marker='o', linestyle = '', fillstyle='full', markeredgewidth=0.0, markersize=2, alpha=0.75, c='black')
plt.show()

fig, ax1 = plt.subplots(figsize=(5, 4))
ax1.hist(np.log10(arr_mN_fGW_SI), weights = arr_mN_weightfac, range=(-11,1), bins=100, histtype='step', fill=False)
ax1.hist(np.log10(arr_mY_fGW_SI), weights = arr_mY_weightfac, range=(-11,1), bins=100, histtype='step', fill=False)
ax1.hist(np.log10(bin_info_arr_REAL[:,2]), weights = bin_info_arr_REAL[:,3]/tH_SI, range=(-11,1), bins=100, histtype='step', fill=False)

ax1.set_yscale('log')
plt.show()


exit()


 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 
 

#--------------------------------------
#PLOTs: PAPER 1.
#--------------------------------------
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#-------------------
#fGW at formation:
#-------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(fGW0_id3),      range=(-5,2), bins=50, histtype='step', stacked=True, fill=False, color = 'royalblue',     label='ejected mergers (outside cluster)', hatch = '\\\\\\\\')
ax1.hist(np.log10(fGW0_id1),      range=(-5,2), bins=50, histtype='step', stacked=True, fill=False, color = 'forestgreen',    label='2-body mergers (inside cluster)', hatch = '//////')
ax1.hist(np.log10(fGW0_id2),      range=(-5,2), bins=50, histtype='step', stacked=True, fill=False, color = 'red',      label='3-body mergers (inside cluster)', hatch = '||||||')

ax1.set_xlabel(r'log $f_{\rm GW}$ [Hz]')
ax1.set_ylabel(r'nr events')
ax1.set_title(r'GW frequency $f_{\rm GW}$ at formation')
ax1.set_yticklabels([])
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
fig.savefig('hist_fGW0_BBHdm.eps', bbox_inches='tight')

plt.show()
#-------------------



#-------------------
#ecc at X Hz:
#-------------------
#histogram:
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(ecc_id3_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'royalblue',    label='ejected mergers (outside cluster)', hatch = '\\\\\\\\')
ax1.hist(np.log10(ecc_id1_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'forestgreen',   label='2-body mergers (inside cluster)', hatch = '//////')

ax1.set_xlim([-3.5, 0.0])
ax1.set_xlabel(r'log $e$')
ax1.set_ylabel(r'nr events')
ax1.set_title(r'Eccentricity $e$ at $f_{\rm GW} = 10^{-2}$ Hz')
ax1.set_yticklabels([])
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
fig.savefig('hist_eccfGW_BBHdm.eps', bbox_inches='tight')
plt.show()



#cumulative histogram:
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(ecc_id3_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'royalblue',    cumulative=True, normed=1, label='ejected mergers', hatch = '\\\\\\\\')
ax1.hist(np.log10(ecc_all_A),   range=(-5,0), bins=50, alpha=1.0, histtype='step', stacked=True, fill=False, color = 'black',   cumulative=True, normed=1, label='ejected + 2-body mergers', hatch = '//////', linewidth=2.0)

ax1.set_xlim([-3.5, 0.0])
ax1.set_ylim([0.0, 1.1])
ax1.set_xlabel(r'log $e$')
ax1.set_ylabel(r'cumulative nr events')
ax1.set_title(r'Eccentricity $e$ at $f_{\rm GW} = 10^{-2}$ Hz')
plt.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)

#save/show fig:
fig.savefig('cumuhist_eccfGW_BBHdm.eps', bbox_inches='tight')
plt.show()
#-------------------
#--------------------------------------






#--------------------------------------
#Figure 2:
#--------------------------------------
fig, ax1 = plt.subplots(figsize=(5, 4))

ax1.hist(np.log10(arr_all_fGW_SI),  weights = arr_all_weightfac,    range=(-11,1), bins=200, histtype='step', fill=False, color = 'orangered',  linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')
ax1.hist(np.log10(bin_info_arr_REAL[pos_e99cut,2]), weights = bin_info_arr_REAL[pos_e99cut,3]/tH_SI, range=(-11,1), bins=200, histtype='step', fill=False, color = 'orangered', linestyle = ':', linewidth = 1, alpha = 1.0, label = r'Newt., e<0.99')
#ax1.hist(np.log10(bin_info_arr_REAL[pos_e9cut,2]), weights = bin_info_arr_REAL[pos_e9cut,3]/tH_SI, range=(-11,1), bins=100, histtype='step', fill=False, color = 'orange', linestyle = '-', linewidth = 1, alpha = 0.5)

ax1.hist(np.log10(arr_mY_fGW_SI),   weights = arr_mY_weightfac,     range=(-11,1), bins=200, histtype='step', fill=False, color = 'darkgray',  linestyle = '-',    linewidth = 2, alpha = 1.0, hatch = '////', label = r'2.5PN, merg. pop.')
ax1.hist(np.log10(arr_mN_fGW_SI),   weights = arr_mN_weightfac,     range=(-11,1), bins=200, histtype='step', fill=False, color = 'royalblue',  linestyle = '-',    linewidth = 2, alpha = 1.0, hatch = '\\\\\\\\', label = r'2.5PN, non-merg. pop.')
#ax1.hist(np.log10(arr_mNY_fGW_SI),  weights = arr_mNY_weightfac,    range=(-11,1), bins=200, histtype='step', fill=False, color = 'black',      linestyle = '-',    linewidth = 2, alpha = 1.0, label = r'+GR, all.')

xvals   = np.arange(-11,1,0.01)
yvals   = 0.1*0.2*(10.**(xvals*(-2./3.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 1, color = 'orangered', label = r'$\propto f^{-2/3}$')
xvals   = np.arange(-11,1,0.01)
yvals   = 0.1*0.000000000002*(10.**(xvals*(-34./9.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 1, color = 'royalblue', label = r'$\propto f^{-34/9}$')
xvals   = np.arange(-11,1,0.01)
yvals   = 0.1*0.000000003*(10.**(xvals*(-3.)))
ax1.plot(xvals, yvals, linestyle = '--', linewidth = 1, color = 'darkgray', label = r'$\propto f^{-3}$')

ax1.text(-2.35, 9, 'LISA', fontsize=15)
ax1.text(-3.0, 5, '////////////////////////////////', fontsize=10)

ax1.set_yscale('log')

ax1.set_xlim([-6.0, -1.0])
ax1.set_ylim([1e-2, 2e3])
#ax1.set_yticklabels([])

ax1.set_xlabel(r'log $f$ [Hz]')
ax1.set_ylabel(r'number events [rnd. norm.]')
ax1.set_title(r'GW Frequency Distributions')
plt.legend(loc='upper right', numpoints = 1, fontsize = 9.0, ncol = 2, frameon = True)

#save/show fig:
fig.savefig(filename+'_obsdist_fGW.eps', bbox_inches='tight')

plt.show()

exit()
#--------------------------------------





exit()
        






