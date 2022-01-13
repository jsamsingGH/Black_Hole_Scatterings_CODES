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

#----------------------------------------------------------
#Input:
#----------------------------------------------------------
filename = raw_input('filename: ')          #filename for output files
Nsamp       = 10000                          #number of samplings
v_kms       = 15.0                          #cluster DISPERSION vel
vesc_kms    = 5.*v_kms                      #cluster ESCAPE     vel
m_Msun      = 30.0                          #BH mass in units of Msun
n_nr_SI     = (10.**(4.))/(m_parsec**3.)    #nr density of single BHs in the cluster
#----------------------------------------------------------

#----------------------------------------------------------
#settings:
#----------------------------------------------------------
tH_SI       = (10.**10.)*sec_year
Nims        = 20
delta_bs    = 7./9.

v_SI    = v_kms*1000.
vesc_SI = vesc_kms*1000.
m_SI    = m_Msun*M_sun_SI

#calc a_HB (hard binary limit):
a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(v_SI**2.))
#calc a_ej (dynamical ejection limit):
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

tot_nr_bins = (Nbs+1)*Nsamp 
bin_info_arr_REAL   = np.zeros((tot_nr_bins,5), dtype=np.float64)
bin_info_arr_INT    = np.zeros((tot_nr_bins,5), dtype=int)

#initialize counters:
tttc        = 0
merger_c    = 0
#----------------------------------------------------------






#----------------------------------------------------------
#generate BBH merger distributions:
#----------------------------------------------------------
#--------------------------------------
#loop over number of samplings:
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
        #BBH orbital params at state 'n':
        #----------------------
        #semi-major axis:
        a_n_SI  = a_max_SI*(delta_bs**n)
        #---------------------
        #draw from P(e)=2e:
        #---------------------
        hit_yn = 0
        while (hit_yn == 0):
            #pick rnd x,y:
            xr = np.random.uniform(0.0, 1.0)
            yr = np.random.uniform(0.0, 2.0)
            #check for hit:
            if (yr <= 2.0*xr):
                #we have a hit:
                hit_yn = 1
        #---------------------
        #now define bin a,e:
        abin_SI = a_n_SI    #SI
        ebin    = xr        #SI
        #----------------------
        
        #----------------------
        #2-body state:
        #----------------------                
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
            #calc minimum ecc for merger (e_IM):
            tinsp_e0_SI = (abin_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
            e_IM        = (1.0 - (tbs_SI/tinsp_e0_SI)**(2./7.))**(1./2.)
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
                id_merger   = 1                     #MERGER ID 1 (2-body merger)
            #---------------- 
        
        #update counter:
        tttc = tttc+1
        #----------------------
        
        #----------------------
        #3-body state:
        #----------------------
        if (merge_yesno == 0):
            #calc minimum ecc for merger (e_cap):
            Rsch_SI = (2.*G_new_SI*m_SI)/(c_SI**2.)
            rcap_SI = 1.8*Rsch_SI*(abin_SI/Rsch_SI)**(2./7.)
            e_cap   = 1.0 - (rcap_SI/abin_SI)
            for rc in range(0,Nims):
                if (merge_yesno == 0):
                    #---------------------
                    #draw from P(e)=2e:
                    #---------------------
                    hit_yn = 0
                    while (hit_yn == 0):
                        #pick rnd x,y:
                        xr = np.random.uniform(0.0, 1.0)
                        yr = np.random.uniform(0.0, 2.0)
                        #check for hit:
                        if (yr <= 2.0*xr):
                            #we have a hit:
                            hit_yn = 1
                            e_bs = xr
                    #---------------------
                    #check if e_bs > e_c:
                    #---------------------
                    if (e_bs >= e_cap):
                        #merger, save:
                        am    = abin_SI
                        em    = e_bs
                        merge_yesno = 1   
                        id_merger   = 2             #MERGER ID 2 (3-body merger)
                    #---------------------
        #----------------------                                           
    #loop over nr hardening steps 'n'    
    #----------------------------------
    
    #----------------------------------
    #if the BBH escapes:
    #----------------------------------
    if (merge_yesno == 0):
        #calc minimum ecc for merger (e_Ht):
        tinsp_e0_SI = (a_min_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.))
        if (tinsp_e0_SI > tH_SI):   e_Ht    = (1.0 - (tH_SI/tinsp_e0_SI)**(2./7.))**(1./2.)
        if (tinsp_e0_SI < tH_SI):   e_Ht    = 0.0
        #---------------------
        #draw from P(e)=2e:
        #---------------------
        hit_yn = 0
        while (hit_yn == 0):
            #pick rnd x,y:
            xr = np.random.uniform(0.0, 1.0)
            yr = np.random.uniform(0.0, 2.0)
            #check for hit:
            if (yr <= 2.0*xr):
                #we have a hit:
                hit_yn = 1
                e_ej = xr
        #---------------------
        #check if e_ej > e_c:
        #---------------------
        if (e_ej >= e_Ht):
             #merger, save:
             am    = a_min_SI
             em    = e_ej
             merge_yesno    = 1   
             id_merger      = 3                     #MERGER ID 3 (ejected merger)    
        #---------------------
    #----------------------------------
    
    #----------------------------------
    #if BBH merged:
    #----------------------------------
    if (merge_yesno == 1):
        
        #calc initial fGW values:
        rp      = am*(1.-em)                                    #rp at formation of merger
        fGW0    = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp**3.)) #fGW at formation of merger        
        #save info:
        output_INT_arr[merger_c,:]     = [id_merger,0,0,0,0,0,0,0,0,0]
        output_REAL_arr[merger_c,:]    = [am,em,fGW0,0,0,0,0,0,0,0]
        
        #update merger counter:
        merger_c = merger_c + 1
    #----------------------------------   
#--------------------------------------
#trim output arrays:
#--------------------------------------
output_INT_arr  = output_INT_arr[0:merger_c,:]
output_REAL_arr = output_REAL_arr[0:merger_c,:]
#--------------------------------------
#----------------------------------------------------------         


#----------------------------------------------------------
#save data:
#----------------------------------------------------------
tf = open(filename + '_BBH_merger_data_INT.txt', "w")
np.savetxt(tf, output_INT_arr,   fmt='%5i')
tf.close()

tf = open(filename + '_BBH_merger_data_REAL.txt', "w")
np.savetxt(tf, output_REAL_arr,   fmt='%5f')
tf.close()
#----------------------------------------------------------


#----------------------------------------------------------
#EXAMPLE: analyze data
#----------------------------------------------------------
nr_merg = len(output_INT_arr[:,0])
pos_id1 = np.where(output_INT_arr[:,0] == 1)[0] # = IM      (1 - 2-body merger)
pos_id2 = np.where(output_INT_arr[:,0] == 2)[0] # = capM    (2 - 3-body merger)
pos_id3 = np.where(output_INT_arr[:,0] == 3)[0] # = ejM     (3 - ejected merger)

fGW0_all    = output_REAL_arr[:,2]
fGW0_id1    = fGW0_all[pos_id1] 
fGW0_id2    = fGW0_all[pos_id2] 
fGW0_id3    = fGW0_all[pos_id3] 

print len(pos_id1)/(1.*nr_merg), len(pos_id2)/(1.*nr_merg), len(pos_id3)/(1.*nr_merg) 

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)
fig, ax1 = plt.subplots(figsize=(5, 5))
Hmin    = -5
Hmax    = 2   
Hnrbins = 100
ax1.hist(np.log10(fGW0_id1), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, color = 'green', hatch = '\\\\\\')
ax1.hist(np.log10(fGW0_id2), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, color = 'red', hatch = '\\\\\\')
ax1.hist(np.log10(fGW0_id3), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, color = 'blue', hatch = '\\\\\\')
ax1.set_xlabel(r'$\log f_{\rm GW}$')
ax1.set_ylabel(r'nr counts')

plt.savefig('fGW_fig.eps', bbox_inches='tight')        

plt.show()
#----------------------------------------------------------


exit()


 





