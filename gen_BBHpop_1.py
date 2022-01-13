import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import integrate as integrate
import matplotlib as mpl
from scipy.optimize import fsolve
from numpy.random import choice


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


sim_YN      = 0 #simulate
sam_YN      = 0 #sample
analyze_YN  = 1 #analyze


#--------------------------------------------------------------
#Simulate:
#--------------------------------------------------------------
if (sim_YN == 1):
#--------------------------------------------------------------
    #----------------------------------------------------------
    #settings:
    #----------------------------------------------------------
    filename = raw_input('save filename: ')
    #input/edit:
    Nsamp       = 10000

    v_kms       = 2.75*10.0     #DISPERSION vel
    vesc_kms    = 2.75*50.0     #ESCAPE     vel
    m_Msun      = 20.0
    n_nr_pc3    = (10.**(5.))   #per 1/pc3

    fGW_limit   = 1e-10

    #define/calc:
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
    #define:
    a_max_SI = a_HB_SI
    a_min_SI = a_ej_SI
    #calc nr. int:
    Nbs         = int((np.log(a_min_SI/a_max_SI)/np.log(delta_bs)))
    #redefine a_min:
    a_min_SI    = a_max_SI*(delta_bs**(Nbs))

    #test info:
    #-------------------
    #amin_SI         = a_min_SI
    #T_BBHprod    = ((1./n_nr_SI)*(v_SI/(1.-delta_bs))*(1./(6.*np.pi*G_new_SI*m_SI))*(1./amin_SI))/(sec_year*10**10)
    #print T_BBHprod
    #exit()
    #abin_SI     = a_HB_SI
    #tinsp_e0_SI = ((abin_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_SI*m_SI*(m_SI+m_SI)/(c_SI**5.)))*(1.-0.0**2)**(7./2.)
    #print tinsp_e0_SI/sec_year
    #exit()
    #-------------------
    
    
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
    #----------------------------------------------------------
    #binary-single:
    #----------------------------------------------------------
    merger_c            = 0
    #--------------------------------------
    for b in range(0,Nsamp):
    #--------------------------------------
        
        print b, Nsamp
        
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
                if (tinsp_e0_SI > tbs_SI):  e_IM        = (1.0 - (tbs_SI/tinsp_e0_SI)**(2./7.))**(1./2.)
                if (tinsp_e0_SI < tbs_SI):  e_IM        = 0.0
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
    np.savetxt(tf, output_INT_arr,      fmt='%5i')
    tf.close()

    tf = open(filename + '_BBH_merger_data_REAL.txt', "w")
    np.savetxt(tf, output_REAL_arr,     fmt='%10.10f')
    tf.close()

    info_sim_arr = np.array([Nsamp, v_kms, vesc_kms, m_Msun, n_nr_pc3, delta_bs, a_HB_SI, a_ej_SI, a_max_SI, a_min_SI])
    tf = open(filename + '_info_sim_arr.txt', "w")
    np.savetxt(tf, info_sim_arr,        fmt='%10.10f')
    tf.close()
    #----------------------------------------------------------
    #----------------------------------------------------------
    exit()
#--------------------------------------------------------------
#--------------------------------------------------------------


#--------------------------------------------------------------
#Create arrays for Sampling:
#--------------------------------------------------------------
if (sam_YN == 1):
#--------------------------------------------------------------
    #----------------------------------------------------------
    #Create observable data set:
    #----------------------------------------------------------
    #--------------------------------------
    #open data:
    #--------------------------------------
    filename = raw_input('open filename: ')

    #open data:
    tf = open(filename + '_BBH_merger_data_INT.txt', "r")
    output_INT_arr  = np.loadtxt(tf, dtype=int)
    tf.close()

    tf = open(filename + '_BBH_merger_data_REAL.txt', "r")
    output_REAL_arr = np.loadtxt(tf, dtype=float)
    tf.close()

    tf = open(filename + '_info_sim_arr.txt', "r")
    info_sim_arr    = np.loadtxt(tf, dtype=float)
    tf.close()
    [Nsamp, v_kms, vesc_kms, m_Msun, n_nr_pc3, delta_bs, a_HB_SI, a_ej_SI, a_max_SI, a_min_SI] = info_sim_arr[:]
    v_SI    = v_kms*1000.
    vesc_SI = vesc_kms*1000.
    m_SI    = m_Msun*M_sun_SI
    n_nr_SI = n_nr_pc3/(m_parsec**3.)
    #--------------------------------------
    #--------------------------------------
    #Define:
    #--------------------------------------
    pos_id1 = np.where(output_INT_arr[:,0] == 1)[0] # = IM      (1 - isolated merger)
    pos_id2 = np.where(output_INT_arr[:,0] == 2)[0] # = capM    (2 - capture merger)
    pos_id3 = np.where(output_INT_arr[:,0] == 3)[0] # = ejM     (3 - ejected merger)
    #--------------------------------------

    #--------------------------------------
    #create obs dist:
    #--------------------------------------
    amin_SI         = a_ej_SI
    T_BBHprod_yr    = ((1./n_nr_SI)*(v_SI/(1.-delta_bs))*(1./(6.*np.pi*G_new_SI*m_SI))*(1./amin_SI))/sec_year
    #input:
    nr_t_bins   = 100                       #t binning
    logt_min    = -20                       #log yr
    logt_max    = np.log10(T_BBHprod_yr)    #log yr
    #set:
    pos_inC     = np.array(list(pos_id1) + list(pos_id2))
    nr_inC      = len(pos_inC)
    #calc:
    dlogt       = (logt_max-logt_min)/(1.*nr_t_bins)
    binarr_logt = np.arange(logt_min,logt_max,dlogt) + dlogt/2.
    binarr_obsW = (10.**(binarr_logt + dlogt/2.)) - (10.**(binarr_logt - dlogt/2.)) #in yr
    posarr_tbin = np.arange(0,nr_t_bins) 
    #define/initialize:
    arr_obsINFO         = np.zeros((nr_inC, nr_t_bins,6), dtype=np.float64)
    arr_obsINFO[:,:,2]  = -1
    tC_yr               = ((768./425.)*(5.*(c_SI**5.)/(512.*G_new_SI**3.)))/sec_year
    #loop:
    for nIC in range(0,nr_inC): #loop over nr in-cluster (id1 + id2) mergers.
        i_inC       = pos_inC[nIC]
        merg_id     = output_INT_arr[i_inC,0] 
        #print info:
        print nIC, nr_inC
        #initial val:
        a0      = output_REAL_arr[i_inC,0]
        e0      = output_REAL_arr[i_inC,1]
        c0      = a0*(1-(e0**2.))*(e0**(-12./19.))*(1. + (121./304.)*(e0**2.))**(-870./2299.)
        t0_yr   = (tC_yr*(a0**4.)/(m_SI**3.))*((1.-e0**2.)**(7./2.))
        rp0     = ((t0_yr/(tC_yr/(m_SI**3.)))/(((1.+e0)**4.)*((1.-e0**2.)**(-1./2.))))**(1./4.) 
        fGW0    = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp0**3.))
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
            #calc:
            tnb_yr  = 10.**(logtnb_yr)
            rp      = ((tnb_yr/(tC_yr/(m_SI**3.)))/(((1.+ecc_nb)**4.)*((1.-ecc_nb**2.)**(-1./2.))))**(1./4.) 
            fGW     = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp**3.))
            #save:
            if (nb == posmax): sflag    = 10
            if (nb <  posmax): sflag    = 1
            arr_obsINFO[nIC, nb, 0] = tnb_yr
            arr_obsINFO[nIC, nb, 1] = ecc_nb
            arr_obsINFO[nIC, nb, 2] = sflag
            arr_obsINFO[nIC, nb, 3] = fGW
            #DONT USE arr_obsINFO[nIC, nb, 4]: it will be used/overwritten later.
            arr_obsINFO[nIC, nb, 5] = merg_id
    #--------------------------------------

    #--------------------------------------
    #Make trimmed arrays:
    #--------------------------------------
    fGW_min         = 10**(-2.0)      #min peak fGW
    tinsp_arr       = np.reshape(arr_obsINFO[:,:,0], nr_inC*nr_t_bins)
    fGW_arr         = np.reshape(arr_obsINFO[:,:,3], nr_inC*nr_t_bins)
    pos             = np.where(fGW_arr > fGW_min)[0]
    tmax_fGW_min    = max(tinsp_arr[pos])
    maxpos_tbin     = max(np.where(binarr_logt[:] < np.log10(tmax_fGW_min))[0])
    max_tbin        = 10.**(binarr_logt[maxpos_tbin])
    nr_tbins_trim   = maxpos_tbin + 1
    #make arr obsINFO:
    arr_obsINFO_fGWtrim = arr_obsINFO[:,0:nr_tbins_trim,:]   #np.zeros((nr_inC, nr_tbins_trim, 5), dtype=np.float64)
    #make arr obsW:
    arr_obsW_fGWtrim    = np.zeros((nr_inC, nr_tbins_trim), dtype=np.float64)
    for nIC in range(0,nr_inC):
        arr_obsW_fGWtrim[nIC,:] = binarr_obsW[0:nr_tbins_trim]
    #--------------------------------------

    #--------------------------------------
    #Make sampling files:
    #--------------------------------------
    #if (nb == posmax): sflag    = 10
    #if (nb <  posmax): sflag    = 1
    #arr_obsINFO[nIC, nb, 0] = tnb_yr
    #arr_obsINFO[nIC, nb, 1] = ecc_nb
    #arr_obsINFO[nIC, nb, 2] = sflag
    #arr_obsINFO[nIC, nb, 3] = fGW
    #arr_obsINFO[nIC, nb, 5] = merg_id
    
    nrbins_1dtrim           = nr_inC*nr_tbins_trim
    nr_obsINFO              = len(arr_obsINFO_fGWtrim[0,0,:])
    obsINFO_fGWtrim_1darr   = np.zeros((nrbins_1dtrim, nr_obsINFO), dtype=np.float64)
    obsW_fGWtrim_1darr      = np.zeros(nrbins_1dtrim, dtype=np.float64)
    
    #obsINFO_fGWtrim_1darr:
    for nobs in range(0,nr_obsINFO):
        obsINFO_fGWtrim_1darr[:,nobs]   = np.reshape(arr_obsINFO_fGWtrim[:,:,nobs], nrbins_1dtrim)     
    
    #obsW_fGWtrim_1darr incl. obs time:
    obsW_fGWtrim_1darr[:]               = np.reshape(arr_obsW_fGWtrim[:,:], nrbins_1dtrim) 
    
    #make final sampling arrays:
    pos_obsBBH          = np.where(obsINFO_fGWtrim_1darr[:,2] > 0)[0]
    nr_obsBBH           = len(pos_obsBBH)
    obsINFO_1darr_samp  = obsINFO_fGWtrim_1darr[pos_obsBBH,:]    
    obsW_1darr_samp     = obsW_fGWtrim_1darr[pos_obsBBH]
    #merge into one array:
    obs_INFO_W_samp_arr = np.zeros((nr_obsBBH, nr_obsINFO), dtype=np.float64)        
    obs_INFO_W_samp_arr[:,:]    = obsINFO_1darr_samp[:,:]
    obs_INFO_W_samp_arr[:,4]    = obsW_1darr_samp[:]
    
    #write to file:
    saveoutput_arr  = obs_INFO_W_samp_arr  
    np.savez(filename + '_saveoutput_arr', saveoutput_arr)
    #--------------------------------------
    exit()
#--------------------------------------------------------------
#--------------------------------------------------------------    


#--------------------------------------------------------------
#Sample and Analyze:
#--------------------------------------------------------------
if (analyze_YN == 1):
#--------------------------------------------------------------
    #--------------------------------------
    #INPUT:
    #--------------------------------------
    #input numbers:
    obs_time_yrs    = 5.0       #input in YEARS
    Ns_tot          = 10000000  #nr of samplings
    #input data file:
    filename = raw_input('open filename: ')
    saveoutput_arr      = np.load(filename + '_saveoutput_arr' + '.npz')['arr_0']
    obs_INFO_W_samp_arr = saveoutput_arr
    tf = open(filename + '_info_sim_arr.txt', "r")
    info_sim_arr        = np.loadtxt(tf, dtype=float)
    tf.close()
    #define:
    [Nsamp, v_kms, vesc_kms, m_Msun, n_nr_pc3, delta_bs, a_HB_SI, a_ej_SI, a_max_SI, a_min_SI] = info_sim_arr[:]
    nr_obsbins          = len(obs_INFO_W_samp_arr[:,0])
    #--------------------------------------
    #--------------------------------------
    #incl. obs. time in weight array:
    #--------------------------------------
    pos_BBHform     = np.where(obs_INFO_W_samp_arr[:,2] == 10)[0]
    obs_INFO_W_samp_arr[pos_BBHform,4]  = obs_INFO_W_samp_arr[pos_BBHform,4] + obs_time_yrs
    obs_INFO_W_samp_arr[:,4]            = obs_INFO_W_samp_arr[:,4]/sum(obs_INFO_W_samp_arr[:,4])    #normalize 
    #--------------------------------------
    #SAMPLE:
    #--------------------------------------    
    pos_obsbins_arr = np.arange(nr_obsbins)
    sampweight_arr  = obs_INFO_W_samp_arr[:,4]
    pos_weightsampl = np.random.choice(pos_obsbins_arr, Ns_tot, p=sampweight_arr) 
    #output arrays:
    tinsp_WS    = obs_INFO_W_samp_arr[pos_weightsampl,0]    #in yrs
    ecc_WS      = obs_INFO_W_samp_arr[pos_weightsampl,1]
    fGW_WS      = obs_INFO_W_samp_arr[pos_weightsampl,3]    #in Hz
    mID_WS      = obs_INFO_W_samp_arr[pos_weightsampl,5]    #1 = (2-body merger), 2 = (3-body merger).
    #--------------------------------------
    #test PLOTs:
    #--------------------------------------
    fig, ax1 = plt.subplots(figsize=(5,4))
    
    pos_mID1    = np.where(mID_WS == 1)[0]
    pos_mID2    = np.where(mID_WS == 2)[0]
    
    #2-body mergers:
    dataplot    = fGW_WS[pos_mID1]
    ax1.hist(np.log10(dataplot), range=(-3,3), bins=50, histtype='step', stacked=True, fill=False, color = 'black')
    #3-body mergers:
    dataplot    = fGW_WS[pos_mID2]
    ax1.hist(np.log10(dataplot), range=(-3,3), bins=50, histtype='step', stacked=True, fill=False, color = 'red')

    ax1.set_yscale('log')
    plt.show()
    #--------------------------------------
    exit()
#--------------------------------------------------------------
#--------------------------------------------------------------


#NOTES:
#T_BBHprod_yr
#incl: ejected mergers
#source formed at high fGW are rare when sampling, but contributes at order unity to the LIGO rate. How do we sample them as good as possible?
#do such that the initial data is saved in a file so we don't have to run it again.
#sometimes this happens: posmax      = max(np.where(binarr_logt[:] < np.log10(t0_yr))[0])
#think again about if the new window binning is ok (?)
#are all the weights correct: do I get the correct powerlaws?
#ejected sources?

#NOTE: to sample more efficiently we sample
#only up to `max_tbin' (= `maxpos_tbin' in pos):
#ONLY CONSIDER fGW down to fGW_min. In the above routine we are sure to sample
#all sources with fGW > fGW_min, however, it will find sources with fGW < fGW_min but
#the sample is incomplete, so don't use it.


exit()

        






