import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import integrate as integrate
import matplotlib as mpl





#----------------------------------------------------------
#functions:
#----------------------------------------------------------
def Pdays_in_aAU_mbin(a_arr, mbin):
    Pdays = 2.*np.pi*np.sqrt((a_arr*AU_SI)**3./(G_new_SI*mbin*M_sun_SI))/(24.*3600.)
    return Pdays

def func_GWinsp_cs_fit(a, mWD, mNS):    #AU, Msun
    corrfac_D       = (1.+(mWD/fp_mb)**fp_r1)**(fp_r2/fp_r1)
    corrfac_A       = (1.+(a/fp_ab)**fp_g1)**(fp_g2/fp_g1)                        
    cs              = fp_csnorm*(corrfac_D*corrfac_A)*((mWD*(mWD+mNS+mNS))/(mWD+mNS))*(a**fp_beta)    
    GWinsp_cs_fit   = cs
    return GWinsp_cs_fit   
def func_aHB_AU(m1, m2, m3, vinf):      #Msun, km/sec
    dimfac_HB   = ((G_new_SI*M_sun_SI)/(1000.**2.))/AU_SI
    aHB_AU      = dimfac_HB*(m1*m2*(m1+m2+m3))/(m3*(m1+m2))/(vinf**2.)
    return aHB_AU
def func_aRL_AU(mWD, mNS, RWD):      #Msun, Rsun
    qfac        = mWD/mNS 
    aRL_AU      = RWD*((0.6*qfac**(2./3.) + np.log(1.+qfac**(1./3.)))/(0.49*qfac**(2./3.)))/AU_U
    return aRL_AU
def func_RWD(mWD):      #Msun
    RWD_Rsun    = 0.0041176640*(mWD/1.0)**(-1./3.)   #in units of Rsun
    return RWD_Rsun
def func_GWinsp_rate_fit(a, mWD, mNS):    #AU, Msun
    dPda            = 1./a                      
    rate            = dPda*func_GWinsp_cs_fit(a, mWD, mNS)    
    GWinsp_rate_fit = rate
    return GWinsp_rate_fit   
def func_GWinsp_rate_simpleana_fit(a, mWD, mNS):    #AU, Msun
    corrfac_D   = (1.+(mWD/fp_mb)**fp_r1)**(fp_r2/fp_r1)
    rate        = corrfac_D*((mWD*(mWD+mNS+mNS))/(mWD+mNS))**(1. + 2./7.) 
    return rate   
#----------------------------------------------------------


#compare cross sections for two datasets (data1, data2):

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

#--------------------------------------------------------------------------------------------------
#open data files: 
#--------------------------------------------------------------------------------------------------
#WD(1.2)-NS(1.2)-NS(1.2):
#DATASET     = 1
#name_data1  = 'SIM3_WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a_at6'
#name_data2  = 'SIM3_WD12pp12pp12_NT_0001to01AU_1kmsec_50000_rp1a_at6'
#name_dataset = [name_data1, name_data2]

##WD-NS-NS:
#DATASET     = 1
#name_data1  = 'WD14NSNS_WT_0001to01AU_1kmsec'
#name_data2  = 'WD14NSNS_NT_0001to01AU_1kmsec_50000'#'WD14NSNS_NT_0001to01AU_1kmsec'
#name_dataset = [name_data1, name_data2]

##MS-MS-MS:
#DATASET     = 2
#name_data1  = '1MS1MS1MS_WTWGR_05to5AU_1kmsec_50000'
#name_data2  = '1MS1MS1MS_NTNGR_05to5AU_1kmsec_50000'
#name_dataset = [name_data1, name_data2]

##MS-NS-NS:
#DATASET     = 3
#name_data1  = 'average_MS10pp10pp10_WT_05to5AU_1kmsec_50000_rp1a_at6'#'MS14NS14NS14_WTWGR_05to5AU_1kmsec_50000'#'MS10NS10NS10_WTWGR_05to5AU_1kmsec_50000'
#name_data2  = 'SIM1_MS10pp10pp10_NT_05to5AU_1kmsec_50000_rp1a_at6'
#name_dataset = [name_data1, name_data2]

#WD-WD-WD:
#DATASET     = 4
#name_data1  = 'SIM4_WDWDWD12_WT_NGR_0001_01_1kmsec_rp1a_at6'
#name_data2  = 'SIM3_WDWDWD12_NT_NGR_0001_01_1kmsec_rp1a_at6'
#name_dataset = [name_data1, name_data2]

#Grindlay analysis:
#DATASET     = 5
#name_data1  = 'JG_02WD_WGR_NT_tde3_insp2'
#name_data2  = 'GR1_WD02_NSNS_WT_WGR_001to100d_10kmsec_10000'
#name_dataset = [name_data1, name_data2]

##Grindlay analysis:
#DATASET     = 6
#name_data1  = 'JG_005WD_WGR_NT_tde3_insp2'
#name_data2  = 'JG_005WD_WGR_WT_tde3_insp2'
#name_dataset = [name_data1, name_data2]

##Grindlay analysis:
#DATASET     = 7
#name_data1  = 'JG_04MS_WGR_NT_tde3_insp2'
#name_data2  = 'JG_04MS_WGR_WT_tde3_insp2'
#name_dataset = [name_data1, name_data2]

##Grindlay analysis:
#DATASET     = 8
#name_data1  = 'JG_08MS_WGR_NT_tde3_insp2'
#name_data2  = 'GR1_MS08_NSNS_WT_WGR_001to100d_10kmsec_10000'
#name_dataset = [name_data1, name_data2]


##WD-NS-NS many dataset:
#DATASET         = 10

##name_dataset    = ['WD02_NSNS_WT_WGR_00003_0001_01_1kmsec', 'WD06_NSNS_WT_WGR_00003_0001_01_1kmsec', 'WD10_NSNS_WT_WGR_00003_0001_01_1kmsec', 'WD14_NSNS_WT_WGR_00003_0001_01_1kmsec']
##name_dataset    = ['TEST_WD02_NSNS_NT_WGR_00003_01_1kmsec', 'TEST_WD04_NSNS_NT_WGR_00003_01_1kmsec', 'TEST2_WD06_NSNS_NT_WGR_00003_01_1kmsec', 'WD10_NSNS_NT_WGR_00003_01_1kmsec', 'WD14_NSNS_NT_WGR_00003_01_1kmsec']
#name_dataset    = ['final_WD02_NSNS_NT_WGR_0001_01_1kmsec', 'final_WD04_NSNS_NT_WGR_0001_01_1kmsec', 'final_WD06_NSNS_NT_WGR_0001_01_1kmsec', 'final_WD10_NSNS_NT_WGR_0001_01_1kmsec', 'final_WD14_NSNS_NT_WGR_0001_01_1kmsec']
#name_dataset    = ['SIM1_WD04_NSNS_NT_WGR_0001_01_1kmsec', 'SIM1_WD06_NSNS_NT_WGR_0001_01_1kmsec', 'SIM1_WD08_NSNS_NT_WGR_0001_01_1kmsec', 'SIM1_WD10_NSNS_NT_WGR_0001_01_1kmsec', 'SIM1_WD12_NSNS_NT_WGR_0001_01_1kmsec', 'final_WD14_NSNS_NT_WGR_0001_01_1kmsec']



#calc using new rptid = 1 threshold:

#name_dataset    = ['average_pp04_NSNS_NT_WGR_0001_01_1kmsec', 'average_pp06_NSNS_NT_WGR_0001_01_1kmsec', 'average_pp10_NSNS_NT_WGR_0001_01_1kmsec', 'average_pp14_NSNS_NT_WGR_0001_01_1kmsec']
#WDmass_arr      = np.array([0.4, 0.6, 1.0, 1.4])       #in units of Msun
#plot_type       = 1

#name_dataset    = ['average_WD04_NSNS_NT_WGR_0001_01_1kmsec', 'average_WD06_NSNS_NT_WGR_0001_01_1kmsec', 'average_WD10_NSNS_NT_WGR_0001_01_1kmsec']
#WDmass_arr      = np.array([0.4, 0.6, 1.0])       #in units of Msun
#plot_type       = 2

#name_dataset    = ['SIM2_WD04_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'SIM2_WD06_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'SIM1_WD10_NSNS_WT_WGR_0001_01_1kmsec']
#WDmass_arr      = np.array([0.4, 0.6, 1.0])       #in units of Msun
#plot_type       = 3

#name_dataset        = ['average_WD04_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'average_WD04_NSNS_NT_WGR_0001_01_1kmsec', 'average_WD06_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'average_WD06_NSNS_NT_WGR_0001_01_1kmsec', 'average_WD10_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'average_WD10_NSNS_NT_WGR_0001_01_1kmsec']
#tidesyesno_dataset  = [1,0, 1,0, 1,0]
#WDmass_arr          = np.array([0.4,0.4, 0.6,0.6, 1.0,1.0])       #in units of Msun
#plot_type           = 31

#name_dataset        = ['average_WD06_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'average_WD06_NSNS_NT_WGR_0001_01_1kmsec']
#tidesyesno_dataset  = [1,0]
#WDmass_arr          = np.array([0.6,0.6])       #in units of Msun
#plot_type           = 31

#name_dataset    = ['SIM1_WD04_NSNS_WT_WGR_001_1_1kmsec_T1000', 'SIM1_WD06_NSNS_WT_WGR_001_1_1kmsec_T1000', 'SIM1_WD10_NSNS_WT_WGR_001_1_1kmsec_T1000']
#WDmass_arr      = np.array([0.4, 0.6, 1.0])       #in units of Msun
#plot_type       = 1



#BH analysis:
#DATASET     = 50
#name_data1  = 'SIM1_BH102020_NT_WGR_0001_10_1kmsec_T2500'
#name_data1  = 'average_BH102020_NT_WGR_0001_10_1kmsec_T2500'
#name_dataset = [name_data1]

DATASET     = 111
name_data1  = 'NIC_MS1CO1CO1_A1'#'average_NIC_WD06NS14NS14_A12345'#'average_Test_NIC_NS14NS14NS14_C123'#'Test_NIC_NS14NS14NS14_PN1'#'average_Test_NIC_NS14NS14NS14_A1234'#'average_Test_NIC_NS14NS14NS14_C123'#'Test_NIC_NS14NS14NS14_B2'
#'average_Test_NIC_NS14NS14NS14'#'average_Test_NIC_NS14NS14NS14_A1234'#'average_Test_NIC_NS14NS14NS14'#'Test_NIC_NS14NS14NS14_2'
name_dataset = [name_data1]









#define:
nr_dataset      = len(name_dataset) 
cs_arr_dataset       = []
cs_err_arr_dataset   = []

for nd in range(0,nr_dataset):
    #open data:
    tf = open('cs_data_'        + name_dataset[nd], "r")
    cs_arr_data        = np.loadtxt(tf, dtype=float)
    tf.close()
    tf = open('cs_err_data_'    + name_dataset[nd], "r")
    cs_err_arr_data    = np.loadtxt(tf, dtype=float)
    tf.close()
    #save data:
    cs_arr_dataset.append(cs_arr_data)
    cs_err_arr_dataset.append(cs_err_arr_data)

#define:
#cs_arr_dataset       = [cs_arr_data1,cs_arr_data2]
#cs_err_arr_dataset   = [cs_err_arr_data1,cs_err_arr_data2]
#data format: (for both cs and cs err)
#save_arr_cs[0,:]        = ap_AU
#save_arr_cs[1,:]        = cs_arr_AU[:,vi,2,0]                                       #NS-NS exchange
#save_arr_cs[2,:]        = cs_arr_AU[:,vi,0,2]                                       #collisions 12
#save_arr_cs[3,:]        = cs_arr_AU[:,vi,1,2]                                       #collisions 13
#save_arr_cs[4,:]        = cs_arr_AU[:,vi,2,2]                                       #collisions 23
#save_arr_cs[5,:]        = cs_arr_AU[:,vi,0,3]                                       #inspiral collision 12
#save_arr_cs[6,:]        = cs_arr_AU[:,vi,1,3]                                       #inspiral collision 13
#save_arr_cs[7,:]        = cs_arr_AU[:,vi,2,3]                                       #inspiral collision 23
#save_arr_cs[8,:]        = cs_arr_AU[:,vi,0,7]                                       #inspiral 12
#save_arr_cs[9,:]        = cs_arr_AU[:,vi,1,7]                                       #inspiral 13
#save_arr_cs[10,:]       = cs_arr_AU[:,vi,2,7]                                       #inspiral 23
#save_arr_cs[11,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_LT_tlife_limit_arr[:]              #NS-NS exchange with t_life<10**10 years
#save_arr_cs[12,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[0,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
#save_arr_cs[13,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[1,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
#save_arr_cs[14,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[2,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
#save_arr_cs[15,:]       = cs_arr_AU[:,vi,0,8]                                       #TDE 12
#save_arr_cs[16,:]       = cs_arr_AU[:,vi,1,8]                                       #TDE 13
#save_arr_cs[17,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[0,:]                #inspiral 23 with rmin cut
#save_arr_cs[18,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[1,:]                #inspiral 23 with rmin cut
#save_arr_cs[19,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[2,:]                #inspiral 23 with rmin cut
#--------------------------------------------------------------------------------------------------





#--------------------------------------------------------------------------------------------------
#INTERACTION: [[O(1)-O(2)]-O(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 111):

    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)
    #-------------------------------------------------
    
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                        = len(cs_arr_dataset[0][0,:])
    sma_arr_AU                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)

    cs_12_insp                  = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_13_insp                  = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_23_insp                  = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_1213_insp                = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_121323_insp              = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_121323_coll              = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_23_coll                  = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_23bin_insp               = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    cs_1213_collinsp            = np.zeros((nr_dataset,nr_a), dtype=np.float64) 
    
    #-------------------------------------------------

    #-------------------------------------------------
    #get cs and err for each dataset:
    #-------------------------------------------------
    for nd in range(0,nr_dataset):
        #sma:
        sma_arr_AU[nd,:]        = cs_arr_dataset[nd][0,:]
        #cs:
        cs_23bin_insp[nd,:]     = cs_fac*(cs_arr_dataset[nd][1,:])
        
        cs_12_insp[nd,:]        = cs_fac*(cs_arr_dataset[nd][8,:])  + cs_fac*(cs_arr_dataset[nd][2,:])
        cs_13_insp[nd,:]        = cs_fac*(cs_arr_dataset[nd][9,:])  + cs_fac*(cs_arr_dataset[nd][3,:])
        cs_23_insp[nd,:]        = cs_fac*(cs_arr_dataset[nd][10,:]) + cs_fac*(cs_arr_dataset[nd][4,:])
        cs_1213_insp[nd,:]      = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:])# + cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:])
        cs_1213_collinsp[nd,:]  = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:]) + cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:])
                
        cs_121323_insp[nd,:]    = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:] + cs_arr_dataset[nd][10,:]) + cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:]+cs_arr_dataset[nd][4,:])
        cs_121323_coll[nd,:]    = cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:]+cs_arr_dataset[nd][4,:])        
        cs_23_coll[nd,:]        = cs_fac*(cs_arr_dataset[nd][4,:])        
        
    #-------------------------------------------------  

    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------
    fig, ax1 = plt.subplots(figsize=(5, 4))
         
    #plot cross sections:
    ax1.plot(sma_arr_AU[0,:],      cs_12_insp[0,:], marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='blue')
    ax1.plot(sma_arr_AU[0,:],      cs_13_insp[0,:], marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='red')
    ax1.plot(sma_arr_AU[0,:],      cs_23_insp[0,:], marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='green')
    ax1.plot(sma_arr_AU[0,:],      cs_1213_insp[0,:], marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='pink')

    ax1.plot(sma_arr_AU[0,:],      cs_121323_insp[0,:], marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='black')
    ax1.plot(sma_arr_AU[0,:],      cs_121323_coll[0,:], marker='^', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='yellow')
    ax1.plot(sma_arr_AU[0,:],      cs_23_coll[0,:], marker='^', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='orange')
    ax1.plot(sma_arr_AU[0,:],      cs_23bin_insp[0,:], marker='^', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='purple')
    
    ax1.plot(sma_arr_AU[0,:],      cs_1213_collinsp[0,:], marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, color='grey')
    
    cs_23bin_insp
    
    res             = 1000
    a_arr           = 10.**(np.linspace(-5, 5, num=res)) 
    cs_fi           = 0.02*a_arr**(2./7.)
    ax1.plot(a_arr, cs_fi, linestyle=':', alpha=1.0, markersize=7.5, linewidth=1.5)

    res             = 1000
    a_arr           = 10.**(np.linspace(-5, 5, num=res)) 
    cs_fi           = 0.25*a_arr**(1./6.)
    ax1.plot(a_arr, cs_fi, linestyle=':', alpha=1.0, markersize=7.5, linewidth=1.5)
    
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_xlim(1e-4, 1e2)
    ax1.set_ylim(1e-5, 1e2)
    
    plt.show()
    
    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------------------
#INTERACTION: [[BH(1)-BH(2)]-BH(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 50):

    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    Mb1 = 10.
    Mb2 = 20.
    Ms3 = 20.
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)
    #-------------------------------------------------
    
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                        = len(cs_arr_dataset[0][0,:])
    sma_arr_AU                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)

    cs_23_bin                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_23_insp                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)    
    cs_23_tlife                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)    
    
    cs_23_bin_err               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_23_insp_err              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_23_tlife_err             = np.zeros((nr_dataset,nr_a), dtype=np.float64)    
    #-------------------------------------------------

    #-------------------------------------------------
    #get cs and err for each dataset:
    #-------------------------------------------------
    for nd in range(0,nr_dataset):
        #sma:
        sma_arr_AU[nd,:]        = cs_arr_dataset[nd][0,:]
        #cs:
        cs_23_bin[nd,:]         = cs_fac*cs_arr_dataset[nd][1,:]
        cs_23_insp[nd,:]        = cs_fac*cs_arr_dataset[nd][10,:]     
        cs_23_tlife[nd,:]       = cs_fac*cs_arr_dataset[nd][11,:]        
        
        #cs err:
        cs_23_bin_err[nd,:]     = cs_fac*cs_err_arr_dataset[nd][1,:]
        cs_23_insp_err[nd,:]    = cs_fac*cs_err_arr_dataset[nd][10,:]   
        cs_23_tlife_err[nd,:]   = cs_fac*cs_err_arr_dataset[nd][11,:]        
    #-------------------------------------------------  

    #define:
    m1_SI   = M_sun_SI*Mb1
    m2_SI   = M_sun_SI*Mb2
    m3_SI   = M_sun_SI*Ms3
    mi_SI   = M_sun_SI*Mb2
    mj_SI   = M_sun_SI*Ms3  
    mk_SI   = M_sun_SI*Mb1      
    mu_ij   = mi_SI*mj_SI/(mi_SI+mj_SI)
    mu_12   = m1_SI*m2_SI/(m1_SI+m2_SI)
    m_ij    = mi_SI+mj_SI
    m_bs_SI = m1_SI + m2_SI + m3_SI
    Mfac    = mu_ij
    eps         = 85.*np.pi/96.
    Rsch        = 2.*G_new_SI*m_ij/(c_SI**2.)    
    betaGW      = 7./2.

    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------
    f   = plt.figure(figsize=(5,5))
    gs  = gridspec.GridSpec(2, 1,height_ratios=[3,1], hspace=0.0)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
     
    #TOP PLOT:
    
    #plot cross sections:
    ax1.errorbar(sma_arr_AU[0,:],  cs_23_bin[0,:],     yerr=cs_23_bin_err[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],      cs_23_bin[0,:],     marker='s', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, label=r'[BH, BH] binary', mew=1, color = 'white')

    ax1.errorbar(sma_arr_AU[0,:],  cs_23_tlife[0,:],   yerr=cs_23_tlife_err[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],      cs_23_tlife[0,:],   marker='^', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, label=r'[BH, BH] binary, $t_{\rm M}<10^{10}$ years', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:],  cs_23_insp[0,:],    yerr=cs_23_insp_err[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],      cs_23_insp[0,:],    marker='o', linestyle='', alpha=1.0, markersize=7.5, linewidth=1.5, label=r'[BH, BH] GW inspiral', mew=1, color = 'black') 


    #oplot analytical GWinsp scaling:
    res             = 1000
    a_arr           = 10.**(np.linspace(-3, 5, num=res)) 
    cs_23_insp_fit  = 2.1*a_arr**(2./7.)
    Mpfac   = (((mu_ij/m_bs_SI)**2.)*((m_bs_SI/mu_ij)**(3./2.))*((mk_SI*mk_SI)/(m1_SI*m2_SI))*((m_ij/mk_SI)**(1./2.)))**(1./betaGW)
    f_tid   = 0.5
    apu     = ((f_tid/2.)**(1./3.))*((mk_SI/mu_ij)**(2./3.)) + 1.
    insp_Ip = ((1.05/(1.-1.7/betaGW))*np.log(apu)*(apu-1.)**(-(1./(betaGW+1.))))/np.log(apu)
    Nfac    = 20.0  #for EM N \approx 5.
    cs_23_R_AU2     = (Nfac*((((m3_SI/mu_12)**(1./3.))*(2.*np.pi*G_new_SI*m_bs_SI*Rsch)/((1000*vel_cs)**2.))*((m1_SI*m2_SI)/(mi_SI*mj_SI)))/(AU_SI**2.))*np.log(apu)
    cs_23_insp_cal  = cs_23_R_AU2*(eps**(1./betaGW))*insp_Ip*Mpfac*((a_arr*AU_SI/Rsch)**(1./betaGW))
    ax1.plot(a_arr,                cs_23_insp_fit, linestyle='--', linewidth=1.5, label=r'analytical sol., $\sigma_{\rm I} \propto a_{0}^{2/7}$', color = 'black')
    
    #post-interaction mergers:
    a_PIM_arr   = 10.**(np.linspace(-1, 5, num=1000)) 
    PBSij       = 0.9
    tau_sec     = (10.**(10.))*(31536000.)
    sigma_fac_AU    = ((((m3_SI/mu_12)**(1./3.))*(2.*np.pi*G_new_SI*(m1_SI+m2_SI+m3_SI)*Rsch)/((1000*vel_cs)**2.))*((m1_SI*m2_SI)/(mi_SI*mj_SI)))/(AU_SI**2.)
    sigma_PMinsp_ij = sigma_fac_AU*((eps**(2./7.))*PBSij*((c_SI*tau_sec/(np.pi*Rsch))**(2./7.))*((4.*mu_ij/m_ij)**(2./7.))*((m1_SI*m2_SI/(mi_SI*mj_SI))**(1./7.))*((a_PIM_arr*AU_SI/Rsch)**(-1./7.)))
    ax1.plot(a_PIM_arr,                sigma_PMinsp_ij, linestyle='-.', linewidth=1.5, label=r'analytical sol., $\sigma_{\rm M} \propto a_{0}^{-1/7}$', color = 'black')    
    
    #plot settings/labels etc:
    ax1.set_title(r'GW Mergers and Inspirals')
    ax1.set_ylabel(r'cross section $\sigma$ [AU$^{2}$] (10 kms$^{-1}$)')
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.set_xlim(5e-4, 1e2)
    ax1.set_ylim(5e-2, 5e5)
    
    ax1.legend(loc='upper left', numpoints = 1, fontsize = 9, frameon = False, ncol=1)
    
    #final axis settings:
    ax1.tick_params(
        axis='x',           # changes apply to the x,y-axis
        which='both',       # both major and minor ticks are affected
        bottom='on',        # ticks along the bottom edge are off
        top='off',          # ticks along the top edge are off
        labelbottom='off',  # labels along the bottom edge are off
        right='off',
        left='off',
        labelleft='off')
    
    
    #BOTTOM PLOT:
    ax2.plot(sma_arr_AU[0,:],      cs_23_insp[0,:]/cs_23_tlife[0,:],    marker='+', linestyle='', alpha=1.0, markersize=8.0, linewidth=1.5, label=r'test', mew=1, color = 'black')
    ax2.plot(sma_arr_AU[0,:],      cs_23_insp[0,:]/cs_23_tlife[0,:],    marker='', linestyle=':', linewidth=1.0, color = 'grey')

    ax2.set_xlabel(r'$a_{0}$ [AU]')
    ax2.set_ylabel(r'$\sigma_{\rm I}/\sigma_{\rm M}$')
    ax2.set_xscale('log')
    ax2.set_yscale('log')
    ax2.set_xlim(5e-4, 1e2)    
      
    #save and show fig:
    #plt.savefig('BHBHBH_cs_ill_1.eps', bbox_inches='tight')
    plt.show()
    
    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------






#--------------------------------------------------------------------------------------------------
#INTERACTION: [[WD(1)-NS(2)]-NS(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 10):
    
    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)    
    WDradius_arr    = 0.013*((1.43/WDmass_arr)**(1./3.))*((1.-WDmass_arr/1.43)**(0.447))
    #print  0.013*((1.43/0.4)**(1./3.))*((1.-0.4/1.43)**(0.447))
    #print (2.9*10**6)*(0.5**(-1./3.))/R_sun_SI
    #-------------------------------------------------
    
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                        = len(cs_arr_dataset[0][0,:])
    sma_arr_AU                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)

    cs_exchange                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_coll_1213                = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspcoll_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_1213                = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_coll_23                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspcoll_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_23                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_12                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_13                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_1213                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_23_cut_rmin1        = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_23_cut_rmin2        = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_23_cut_rmin3        = np.zeros((nr_dataset,nr_a), dtype=np.float64)

    cs_err_exchange             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_coll_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspcoll_1213        = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_coll_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspcoll_23          = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_12               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_13               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_1213             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_23_cut_rmin1    = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_23_cut_rmin2    = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_23_cut_rmin3    = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    #-------------------------------------------------


    #-------------------------------------------------
    #calc defined cs for each dataset:
    #-------------------------------------------------
    for nd in range(0,nr_dataset):
    
        sma_arr_AU[nd,:]            = cs_arr_dataset[nd][0,:]
   
        #DATA: cs
        #exchange:
        cs_exchange[nd,:]           = cs_fac*cs_arr_dataset[nd][1,:]
        #WD-NS:
        cs_coll_1213[nd,:]          = cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:])
        cs_inspcoll_1213[nd,:]      = cs_fac*(cs_arr_dataset[nd][5,:] + cs_arr_dataset[nd][6,:])
        cs_insp_1213[nd,:]          = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:])
        #NS-NS:
        cs_coll_23[nd,:]            = cs_fac*cs_arr_dataset[nd][4,:]
        cs_inspcoll_23[nd,:]        = cs_fac*cs_arr_dataset[nd][7,:]
        cs_insp_23[nd,:]            = cs_fac*cs_arr_dataset[nd][10,:]
        #TDE
        cs_TDE_12[nd,:]             = cs_fac*cs_arr_dataset[nd][15,:]        
        cs_TDE_13[nd,:]             = cs_fac*cs_arr_dataset[nd][16,:]        
        cs_TDE_1213[nd,:]           = cs_TDE_12[nd,:] + cs_TDE_13[nd,:]
        #GW insp cut: rmin
        cs_insp_23_cut_rmin1[nd,:]  = cs_fac*cs_arr_dataset[nd][17,:]  
        cs_insp_23_cut_rmin2[nd,:]  = cs_fac*cs_arr_dataset[nd][18,:]  
        cs_insp_23_cut_rmin3[nd,:]  = cs_fac*cs_arr_dataset[nd][19,:]  
        
        #DATA: cs-err
        #exchange:
        cs_err_exchange[nd,:]       = cs_fac*cs_err_arr_dataset[nd][1,:]
        #WD-NS:
        cs_err_coll_1213[nd,:]      = cs_fac*np.sqrt(cs_err_arr_dataset[nd][2,:]**2. + cs_err_arr_dataset[nd][3,:]**2.)
        cs_err_inspcoll_1213[nd,:]  = cs_fac*np.sqrt(cs_err_arr_dataset[nd][5,:]**2. + cs_err_arr_dataset[nd][6,:]**2.)
        cs_err_insp_1213[nd,:]      = cs_fac*np.sqrt(cs_err_arr_dataset[nd][8,:]**2. + cs_err_arr_dataset[nd][9,:]**2.)
        #NS-NS:
        cs_err_coll_23[nd,:]        = cs_fac*cs_err_arr_dataset[nd][4,:]
        cs_err_inspcoll_23[nd,:]    = cs_fac*cs_err_arr_dataset[nd][7,:]
        cs_err_insp_23[nd,:]        = cs_fac*cs_err_arr_dataset[nd][10,:]
        #TDE
        cs_err_TDE_12[nd,:]         = cs_fac*cs_err_arr_dataset[nd][15,:]        
        cs_err_TDE_13[nd,:]         = cs_fac*cs_err_arr_dataset[nd][16,:]  
        cs_err_TDE_1213[nd,:]       = np.sqrt(cs_err_TDE_12[nd,:]**2. + cs_err_TDE_13[nd,:]**2.)  
        #GW insp cut: rmin
        cs_err_insp_23_cut_rmin1[nd,:]  = cs_fac*cs_err_arr_dataset[nd][17,:]  
        cs_err_insp_23_cut_rmin2[nd,:]  = cs_fac*cs_err_arr_dataset[nd][18,:]  
        cs_err_insp_23_cut_rmin3[nd,:]  = cs_fac*cs_err_arr_dataset[nd][19,:]   
    #-------------------------------------------------
    
    print cs_insp_23
    print cs_err_insp_23


    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------
    #SETTINGS:
    #GW insp fitting formula params:
    aplot       = np.logspace(-3,10,num=100) 
    fp_beta     = 2./7.         #GW asymptotic slope
    fp_g1       = -2.           #a_0 (transition slop)
    fp_g2       = 0.5-fp_beta   #a_0 (correciton slope at low a)
    fp_r1       = 10.0          #nr res high e (transition slope)
    fp_r2       = -0.5          #nr res high e (slope when m_{WD}> m_{BR})
    fp_ab       = 1e-2          #a break
    fp_mb       = 1.0           #m_{WD} break
    fp_csnorm   = 1.65*(1e-2)   #cs norm    
    
    #def func_GWinsp_cs_fit(a, mWD, mNS):    #AU, Msun
    #def func_aHB_AU(m1, m2, m3, vinf):      #Msun, km/sec
    #def func_aRL_AU(mWD, mNS, RWD):      #Msun, Rsun
    #def func_RWD(mWD):      #Msun
    #def func_GWinsp_rate_fit(a, mWD, mNS):    #AU, Msun
        
    #PLOT cs mass correction factor:
    if (plot_type == 4):
        fig = plt.figure(figsize=(5,4))
        nrmWD           = 100
        rate_arr_norm   = np.zeros(nrmWD, dtype=np.float64) 
        adum            = 1.0   #dummy val for a
        mWD_arr         = np.linspace(0.2,1.4,num=nrmWD) 
        mWDnorm         = 0.6
        #cross section:
        cs_norm         = func_GWinsp_cs_fit(adum,mWDnorm,1.4) 
        cs_arr_norm     = func_GWinsp_cs_fit(adum,mWD_arr,1.4)/cs_norm
        #rates:
        amin = func_aRL_AU(mWDnorm, 1.4, func_RWD(mWDnorm)) #Roche lobe
        amax = func_aHB_AU(mWDnorm, 1.4, 1.4, vel_cs)   #HB limit
        rate_norm   = integrate.quad(func_GWinsp_rate_fit,amin,amax,args=(mWDnorm,1.4))[0]        
        for nw in range(0,nrmWD):
            mWD     = mWD_arr[nw]
            amin    = func_aRL_AU(mWD, 1.4, func_RWD(mWD)) #Roche lobe
            amax    = func_aHB_AU(mWD, 1.4, 1.4, vel_cs)   #HB limit
            rate_arr_norm[nw] = integrate.quad(func_GWinsp_rate_fit,amin,amax,args=(mWD,1.4))[0]/rate_norm
        #plot:
        fig.add_subplot(111).plot(mWD_arr, cs_arr_norm, marker='', linestyle='-', linewidth=2.0, color='black', label=r'cross section $\sigma/\sigma(m_{WD}=0.6M_{\odot})$')
        fig.add_subplot(111).plot(mWD_arr, rate_arr_norm, marker='', linestyle='--', linewidth=2.0, color='black', label=r'rate $\Gamma/\Gamma(m_{WD}=0.6M_{\odot})$')
        #test:
        GWinsp_rate_simpleana_fit   = func_GWinsp_rate_simpleana_fit(adum,mWD_arr,1.4)/func_GWinsp_rate_simpleana_fit(adum,mWDnorm,1.4)
        fig.add_subplot(111).plot(mWD_arr, GWinsp_rate_simpleana_fit, marker='', linestyle=':', linewidth=2.0, color='black', label=r'rate $\Gamma/\Gamma(m_{WD}=0.6M_{\odot})$ (approximation)')

        fig.add_subplot(111).plot([mWDnorm,mWDnorm], [-10,10], marker='', linestyle=':', linewidth=0.5, color='black')
        fig.add_subplot(111).plot([-10,10], [1,1], marker='', linestyle=':', linewidth=0.5, color='black')
        #add labels etc:
        fig.add_subplot(111).set_ylim(0.0,2.25)
        fig.add_subplot(111).set_xlim(0.2,1.4)
        plt.xlabel(r'$m_{WD}$')
        fig.add_subplot(111).legend(loc='lower right', numpoints = 1, fontsize = 7.5, frameon = False, ncol=1)
        #save and show fig:
        plt.savefig('rel_GWinsp_cs.eps', bbox_inches='tight')
        plt.show()
        exit()
    
    #PLOT insp rate wighted with WD mass dist:
    if (plot_type == 5):
        fig = plt.figure(figsize=(5,4))
        #rate per bin:
        mWD_arr = np.linspace(0.2,1.4,num=1000) 
        Dfac    = (1.+(mWD_arr/fp_mb)**fp_r1)**(fp_r2/fp_r1)
        rel_cs  = Dfac*((mWD_arr*(mWD_arr+1.4+1.4))/(mWD_arr+1.4))**(1.+2./7.)  #(mWD_arr/1.4)*(1.4/(mWD_arr+1.4))*(mWD_arr+1.4+1.4)*(1.+(mWD_arr/fp_mb)**fp_r1)**(fp_r2/fp_r1)         
        #WD mass dist (http://iopscience.iop.org/article/10.1088/0067-0049/204/1/5/pdf)
        #gauss 1:
        g_norm     = 0.56 
        g_mean     = 0.589
        g_sigma    = 0.06
        distWDmass_gauss_1  = g_norm*(1./(g_sigma*np.sqrt(2.*np.pi))*np.exp(-((mWD_arr-g_mean)**2.)/(2.*g_sigma**2.)))  
        #gauss 2:
        g_norm     = 0.25 
        g_mean     = 0.587
        g_sigma    = 0.02
        distWDmass_gauss_2  = g_norm*(1./(g_sigma*np.sqrt(2.*np.pi))*np.exp(-((mWD_arr-g_mean)**2.)/(2.*g_sigma**2.)))  
        #gauss 3:
        g_norm     = 0.13 
        g_mean     = 0.822
        g_sigma    = 0.1
        distWDmass_gauss_3  = g_norm*(1./(g_sigma*np.sqrt(2.*np.pi))*np.exp(-((mWD_arr-g_mean)**2.)/(2.*g_sigma**2.)))  
        #gauss 4:
        g_norm     = 0.06 
        g_mean     = 0.389
        g_sigma    = 0.03
        distWDmass_gauss_4  = g_norm*(1./(g_sigma*np.sqrt(2.*np.pi))*np.exp(-((mWD_arr-g_mean)**2.)/(2.*g_sigma**2.)))  
        #final distribution:
        distWDmass_gauss_final      = distWDmass_gauss_1 + distWDmass_gauss_2 + distWDmass_gauss_3 + distWDmass_gauss_4
        norm_WDmass_dist            = max(distWDmass_gauss_final)
        #rate weighted distribution:
        distWDmass_rate_weighted    = rel_cs*distWDmass_gauss_final
        norm_WDmass_rate_weighted   = max(distWDmass_rate_weighted)
        #PLOT:
        #each gauss component:
        fig.add_subplot(111).plot(mWD_arr, distWDmass_gauss_1/norm_WDmass_dist, marker='', linestyle=':', linewidth=1.0)
        fig.add_subplot(111).plot(mWD_arr, distWDmass_gauss_2/norm_WDmass_dist, marker='', linestyle=':', linewidth=1.0)
        fig.add_subplot(111).plot(mWD_arr, distWDmass_gauss_3/norm_WDmass_dist, marker='', linestyle=':', linewidth=1.0)
        fig.add_subplot(111).plot(mWD_arr, distWDmass_gauss_4/norm_WDmass_dist, marker='', linestyle=':', linewidth=1.0)
        #rel GW rate:
        #fig.add_subplot(111).plot(mWD_arr, rel_cs/max(rel_cs), marker='', linestyle='--', linewidth=1.0, color='black', label='relative rate')
        #final gauss:
        fig.add_subplot(111).plot(mWD_arr, distWDmass_gauss_final/norm_WDmass_dist, marker='', linestyle='-', linewidth=2.0, color='black', label='WD mass distribution')
        #rate weighted:
        fig.add_subplot(111).plot(mWD_arr, distWDmass_rate_weighted/norm_WDmass_rate_weighted, marker='', linestyle='-', linewidth=2.0, color='red', label='GW inspiral rate distribution')
        #add labels etc:
        plt.xlabel(r'$m_{WD}$')
        plt.ylabel(r'normalized counts')
        fig.add_subplot(111).legend(loc='upper right', numpoints = 1, fontsize = 8., frameon = False, ncol=1)
        #save and show fig:
        plt.savefig('mWD_dist_and_GWrates.eps', bbox_inches='tight')
        plt.show()
        exit()
        
        
    
    
    
    
    #CROSS SECTION PLOTs:
    if (plot_type == 1 or plot_type == 2 or plot_type == 3):
        fig, ax1 = plt.subplots(figsize=(5, 2.5))
        colornr = [0, 50, 100, 150, 200, 220]
        
    if (plot_type == 31):
        fig, ax1 = plt.subplots(figsize=(5, 4))
        colornr = [0,0, 50,50, 100,100, 150,150, 200,200, 220,220]
    
            
    for nd in range(0,nr_dataset):
        
        #-------------------------------
        #analytical calculations:
        #-------------------------------
        vinf    = 10.   #km/sec (all simulations are (must be) done with 1km/sec!!!!!)
        a0      = sma_arr_AU[nd,:]*AU_U  #Rsun
        a0AU    = sma_arr_AU[nd,:]       #AU
        m1      = WDmass_arr[nd]                #Msun
        m2      = 1.4                           #Msun
        m3      = 1.4                           #Msun
        mi      = m2
        mj      = m3
        mk      = m1
        m_bs    = m1+m2+m3
        m_ij    = mi+mj
        mu_ij   = mi*mj/m_ij
        mu_12   = m1*m2/(m1+m2)
        #for GW inspirals: 
        Es      = np.pi*(85./96.)
        beta    = 7./2.
        Ms      = mu_ij
        Mfac    = ((m1*m2)/(mi*mj))*(((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
        RsMsun  = (2.*G_new_SI*(1.0*M_sun_SI)/(c_SI**2.))           #in m for 1Msun
        C_GW    = np.pi*((85.*np.pi/96.)**(2./7.))*((RsMsun/AU_SI)**(12./7.))*((c_SI/1000.)**2.)    #works for units of Msun, AU, km/s.
        f_tid   = 0.5
        apu     = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
        insp_I  = (1.05/(1.-1.7/beta))*np.log(apu)*(apu-1.)**(-(1./(beta+1.)))
        #Insp cs - ANALYTICAL:    
        #calc the 'N term':
        #Nfac_arr            = np.array([4.37991448, 3.51842139, 2.5029299, 1.95])   #for [0.4, 0.6, 1.0, 1.4]
        #Nfac_23             = Nfac_arr[nd]
        ##calc GW insp cs:
        #cs_insp_23_CALC = Nfac_23*C_GW*insp_I*Mfac*(m_bs*(m_ij**(5./7.))*(a0AU**(2./7.))/(vinf**2.))*((m3/mu_12)**(1./3.))  #*(mu_12/m1)
        #print 'CALC, SIM 23 GW insp:    '
        #print cs_insp_23_CALC[:]
        #print cs_insp_23[nd,:]
        #print cs_insp_23[nd,:]/cs_insp_23_CALC[:]
        #-------------------------------
        
        
            
        #-------------------------------
        #PLOT 1: Cross sections
        #-------------------------------
        if (plot_type == 1):

            #GW inspirals:
            ax1.errorbar(sma_arr_AU[nd,:],  cs_insp_23[nd,:],     yerr=cs_err_insp_23[nd,:],   marker='', linestyle='-', linewidth=2.0, color=plt.cm.terrain(colornr[nd]))
            ax1.plot(sma_arr_AU[nd,:],      cs_insp_23[nd,:],     marker='', linestyle='-', linewidth=1.5, color=plt.cm.terrain(colornr[nd]), label=r' $m_{WD}$ ='+str('%.1f' % (WDmass_arr[nd])))
            ax1.plot(sma_arr_AU[nd,:],      cs_insp_23[nd,:],     marker='x', linestyle='', markersize=4.0, color=plt.cm.terrain(colornr[nd]))
        
            #GW inspirals with rmin cut:
            #ax1.errorbar(sma_arr_AU[nd,:],  cs_insp_23_cut_rmin1[nd,:],     yerr=cs_err_insp_23_cut_rmin1[nd,:],   marker='', linestyle=':', linewidth=2.0, color=plt.cm.terrain(colornr[nd]))
            #ax1.legend(loc='upper right', numpoints = 1, fontsize = 7.5, frameon = False, ncol=1)
                
            #background exchange:
            #ax1.errorbar(sma_arr_AU[nd,:], cs_exchange[nd,:],     yerr=cs_err_exchange[nd,:],  marker='', linestyle='-', linewidth=0.5, color = 'lightgrey')
            #background collisions:
            #ax1.errorbar(sma_arr_AU[nd,:], cs_coll_1213[nd,:],    yerr=cs_err_coll_1213[nd,:], marker='', linestyle='-', linewidth=0.5, color = 'lightgrey')
                    
            ax1.text(0.95, 0.9, 'NS-NS GW inspirals',
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax1.transAxes,
                    color='black', fontsize=10)

            #plot settings:
            ax1.set_xlim([1e-5,1e2])    
            ax1.set_ylim([1e-4,1e0])
            ax1.set_title(r'WD: Point Mass')
        #-------------------------------
            
        #-------------------------------
        #PLOT 2: Cross sections
        #-------------------------------
        if (plot_type == 2):

            #GW inspirals:
            ax1.errorbar(sma_arr_AU[nd,:],  cs_insp_23[nd,:],     yerr=cs_err_insp_23[nd,:],   marker='', linestyle='-', linewidth=2.0, color=plt.cm.terrain(colornr[nd]))
            ax1.plot(sma_arr_AU[nd,:],      cs_insp_23[nd,:],     marker='', linestyle='-', linewidth=1.5, color=plt.cm.terrain(colornr[nd]), label=r' $m_{WD}$ ='+str('%.1f' % (WDmass_arr[nd])))
            ax1.plot(sma_arr_AU[nd,:],      cs_insp_23[nd,:],     marker='x', linestyle='', markersize=4.0, color=plt.cm.terrain(colornr[nd]))
        
            #GW inspirals with rmin cut:
            #ax1.errorbar(sma_arr_AU[nd,:],  cs_insp_23_cut_rmin1[nd,:],     yerr=cs_err_insp_23_cut_rmin1[nd,:],   marker='', linestyle=':', linewidth=2.0, color=plt.cm.terrain(colornr[nd]))
            #ax1.legend(loc='upper right', numpoints = 1, fontsize = 7.5, frameon = False, ncol=1)
        
            #background exchange:
            #ax1.errorbar(sma_arr_AU[nd,:], cs_exchange[nd,:],     yerr=cs_err_exchange[nd,:],  marker='', linestyle='-', linewidth=0.5, color = 'lightgrey')
            #background collisions:
            #ax1.errorbar(sma_arr_AU[nd,:], cs_coll_1213[nd,:],    yerr=cs_err_coll_1213[nd,:], marker='', linestyle='-', linewidth=0.5, color = 'lightgrey')
                    
            ax1.text(0.95, 0.9, 'NS-NS GW inspirals',
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax1.transAxes,
                    color='black', fontsize=10)
       
            #plot settings:
            ax1.set_xlim([1e-5,1e2])    
            ax1.set_ylim([1e-4,1e0])
            ax1.set_title(r'WD: Solid Sphere')
        #-------------------------------        

        #-------------------------------
        #PLOT 3: Cross sections
        #-------------------------------
        if (plot_type == 3):

            #GW inspirals:
            ax1.errorbar(sma_arr_AU[nd,:],  cs_insp_23[nd,:],     yerr=cs_err_insp_23[nd,:],   marker='', linestyle='-', linewidth=2.0, color=plt.cm.terrain(colornr[nd]))
            ax1.plot(sma_arr_AU[nd,:],      cs_insp_23[nd,:],     marker='', linestyle='-', linewidth=1.5, color=plt.cm.terrain(colornr[nd]), label=r' $m_{WD}$ ='+str('%.1f' % (WDmass_arr[nd])))
            ax1.plot(sma_arr_AU[nd,:],      cs_insp_23[nd,:],     marker='x', linestyle='', markersize=4.0, color=plt.cm.terrain(colornr[nd]))
        
            #GW inspirals with rmin cut:
            #ax1.errorbar(sma_arr_AU[nd,:],  cs_insp_23_cut_rmin1[nd,:],     yerr=cs_err_insp_23_cut_rmin1[nd,:],   marker='', linestyle=':', linewidth=2.0, color=plt.cm.terrain(colornr[nd]))
            #ax1.legend(loc='upper right', numpoints = 1, fontsize = 7.5, frameon = False, ncol=1)
        
            #background exchange:
            #ax1.errorbar(sma_arr_AU[nd,:], cs_exchange[nd,:],     yerr=cs_err_exchange[nd,:],  marker='', linestyle='-', linewidth=0.5, color = 'lightgrey')
            #background collisions:
            #ax1.errorbar(sma_arr_AU[nd,:], cs_coll_1213[nd,:],    yerr=cs_err_coll_1213[nd,:], marker='', linestyle='-', linewidth=0.5, color = 'lightgrey')
                    
            ax1.text(0.95, 0.9, 'NS-NS GW Inspirals',
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax1.transAxes,
                    color='black', fontsize=10)
       
            #plot settings:
            ax1.set_xlim([1e-5,1e2])    
            ax1.set_ylim([1e-4,1e0])
            ax1.set_title(r'WD: $\gamma = 5/3$, $n=3$ Polytrope')
        #-------------------------------        
  
  
        
        #-------------------------------
        #PLOT 3: Cross sections
        #-------------------------------
        if (plot_type == 31):
                                 
            #with tides:
            if (tidesyesno_dataset[nd] == 1):                
                a       = sma_arr_AU[nd,:]
                cs      = cs_insp_1213[nd,:] + cs_inspcoll_1213[nd,:] + cs_coll_1213[nd,:] + cs_TDE_1213[nd,:] 
                cs_err  = np.sqrt(cs_err_insp_1213[nd,:]**2. + cs_err_inspcoll_1213[nd,:]**2. + cs_err_coll_1213[nd,:]**2. + cs_err_TDE_1213[nd,:]**2.)
                ax1.errorbar(a, cs,     yerr=cs_err,   marker='', linestyle='-', linewidth=2.0, color='black')
                ax1.plot(a, cs,     marker='', linestyle='-', linewidth=1.5, color='black', label=r' WD-NS mergers $(m_{\rm WD}$ ='+str('%.1f' % (WDmass_arr[nd])) + ', $+$ tides)')
                ax1.plot(a, cs,     marker='x', linestyle='', markersize=4.0, color='black')
            
            #no tides:    
            if (tidesyesno_dataset[nd] == 0):
                a       = sma_arr_AU[nd,:]
                cs      = cs_coll_1213[nd,:]
                cs_err  = np.sqrt(cs_err_insp_1213[nd,:]**2. + cs_err_inspcoll_1213[nd,:]**2. + cs_err_coll_1213[nd,:]**2. + cs_err_TDE_1213[nd,:]**2.)
                ax1.errorbar(a, cs,     yerr=cs_err,   marker='', linestyle=':', linewidth=2.0, color='black')
                ax1.plot(a, cs,     marker='', linestyle=':', linewidth=1.5, color='black', label=r' WD-NS mergers $(m_{\rm WD}$ ='+str('%.1f' % (WDmass_arr[nd])) + ', $-$ tides)')
                ax1.plot(a, cs,     marker='x', linestyle='', markersize=4.0, color='black')                
            
            #RL radius:
            qfac        = m1/m2
            WDrad_AU    = (WDradius_arr[nd])/AU_U
            RL_WD_AU    = WDrad_AU*((0.6*qfac**(2./3.) + np.log(1. + qfac**(1./3.)))/(0.49*qfac**(2./3.)))
            ax1.plot([RL_WD_AU, RL_WD_AU], [1e-10,1e10],     marker='', linestyle=':', linewidth=0.25, color='black')
            
            #oplot cross section scaling:
            sma_AU_arr      = 10.**(np.linspace(np.log10(RL_WD_AU),10,1000))
            cs_ex_arr       = 0.54*sma_AU_arr**(1./6.)
            ax1.plot(sma_AU_arr, cs_ex_arr,     marker='', linestyle='--', linewidth=0.5, color='black')
            ax1.text(0.9, 0.52, r'$\sigma_{\rm insp} \propto a_{0}^{1/6}$ $(\beta = 6)$',
                    verticalalignment='bottom', horizontalalignment='right',
                    transform=ax1.transAxes,
                    color='black', fontsize=14, rotation=20.0)            
            
            ax1.fill([1e-10, RL_WD_AU, RL_WD_AU, 1e-10], [-1e10, 1e-10, 1e10, 1e10], fill=False, hatch='\\')    
            
            #ax1.text(0.95, 0.9, 'NS-WD Mergers',
            #        verticalalignment='bottom', horizontalalignment='right',
            #        transform=ax1.transAxes,
            #        color='black', fontsize=10)
       
            #plot settings:
            ax1.set_xlim([1e-4,1e0])    
            ax1.set_ylim([4e-2,1e0])
            ax1.set_title(r'Neutron Star - White Dwarf Mergers')
        #-------------------------------        

  
  
            
        if (plot_type == 1 or plot_type == 2 or plot_type == 3):
           #oplot a^2/7 scaling (fit by eye):    
           cs_scale_arr    = 0.010*np.array([1.0, 1.7, 2.9, 3.5])
           sma_AU_arr      = 10.**(np.linspace(-10,10,1000))
           cs_nd           = cs_scale_arr[nd]*(sma_AU_arr**(2./7.))
           #cs_nd           = cs_insp_23_CALC
           ax1.plot(sma_AU_arr, cs_nd, marker='', linestyle='--', linewidth=0.5, color=plt.cm.terrain(colornr[nd]))
           
       
        
    plt.yticks([0.5,1.0,1.5])
    
    #Titles/labels:
    ax1.set_ylabel(r'$\sigma$ [AU$^2$] (10 km/sec)')
    ax1.set_xlabel(r'$a_{0}$ [AU]')
    
    #x-scales:
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    #Legends:
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False, ncol=1)
        
    #Save figure:
    if (plot_type == 1):    plt.savefig('insp_cs_WDpp.eps', bbox_inches='tight')
    if (plot_type == 2):    plt.savefig('insp_cs_WDss.eps', bbox_inches='tight')
    if (plot_type == 3):    plt.savefig('insp_cs_WDpt.eps', bbox_inches='tight')
    if (plot_type == 31):   plt.savefig('WD_ss_vs_dt.eps', bbox_inches='tight')
    #if (plot_type == 2): plt.savefig('TEST_WDNS_coll.eps', bbox_inches='tight')
    
    plt.show()


    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------















#--------------------------------------------------------------------------------------------------
#INTERACTION: [[WD(1)-NS(2)]-NS(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 1):

    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    #mass of one of the (equal mass) objects
    MO              = 1.2   #in units of Msun
    #radius of tidal object
    RTO             = 0.013*((1.43/MO)**(1./3.))*((1.-MO/1.43)**(0.447))
    print 'RTO: ', RTO
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)
    #-------------------------------------------------
    
    
    
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                        = len(cs_arr_dataset[0][0,:])
    sma_arr_AU                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)

    cs_exchange                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_coll_1213                = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspcoll_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_1213                = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_coll_23                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspcoll_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_23                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_12                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_13                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_1213                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspTDE_1213             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    
    cs_err_exchange             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_coll_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspcoll_1213        = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_coll_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspcoll_23          = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_12               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_13               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_1213             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspTDE_1213         = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    
    
    #-------------------------------------------------


    #-------------------------------------------------
    #calc defined cs for each dataset:
    #-------------------------------------------------
    for nd in range(0,nr_dataset):
    
        sma_arr_AU[nd,:]            = cs_arr_dataset[nd][0,:]
   
        #DATA: cs
        #exchange:
        cs_exchange[nd,:]           = cs_fac*cs_arr_dataset[nd][1,:]
        #WD-NS:
        cs_coll_1213[nd,:]          = cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:])
        cs_inspcoll_1213[nd,:]      = cs_fac*(cs_arr_dataset[nd][5,:] + cs_arr_dataset[nd][6,:])
        cs_insp_1213[nd,:]          = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:])
        #NS-NS:
        cs_coll_23[nd,:]            = cs_fac*cs_arr_dataset[nd][4,:]
        cs_inspcoll_23[nd,:]        = cs_fac*cs_arr_dataset[nd][7,:]
        cs_insp_23[nd,:]            = cs_fac*cs_arr_dataset[nd][10,:]
        #TDE
        cs_TDE_12[nd,:]             = cs_fac*cs_arr_dataset[nd][15,:]        
        cs_TDE_13[nd,:]             = cs_fac*cs_arr_dataset[nd][16,:]        
        cs_TDE_1213[nd,:]           = cs_TDE_12[nd,:] + cs_TDE_13[nd,:]
        
        cs_inspTDE_1213[nd,:]       = cs_insp_1213[nd,:] + cs_TDE_1213[nd,:]    #WE NEED TO INCL the TDE sample due to gamma=4/3 which results in large osc. This produces more TDEs. Might just be a numerical issue. its ok.
        
        #DATA: cs-err
        #exchange:
        cs_err_exchange[nd,:]       = cs_fac*cs_err_arr_dataset[nd][1,:]
        #WD-NS:
        cs_err_coll_1213[nd,:]      = cs_fac*np.sqrt(cs_err_arr_dataset[nd][2,:]**2. + cs_err_arr_dataset[nd][3,:]**2.)
        cs_err_inspcoll_1213[nd,:]  = cs_fac*np.sqrt(cs_err_arr_dataset[nd][5,:]**2. + cs_err_arr_dataset[nd][6,:]**2.)
        cs_err_insp_1213[nd,:]      = cs_fac*np.sqrt(cs_err_arr_dataset[nd][8,:]**2. + cs_err_arr_dataset[nd][9,:]**2.)
        #NS-NS:
        cs_err_coll_23[nd,:]        = cs_fac*cs_err_arr_dataset[nd][4,:]
        cs_err_inspcoll_23[nd,:]    = cs_fac*cs_err_arr_dataset[nd][7,:]
        cs_err_insp_23[nd,:]        = cs_fac*cs_err_arr_dataset[nd][10,:]
        #TDE
        cs_err_TDE_12[nd,:]         = cs_fac*cs_err_arr_dataset[nd][15,:]        
        cs_err_TDE_13[nd,:]         = cs_fac*cs_err_arr_dataset[nd][16,:]  
        cs_err_TDE_1213[nd,:]       = np.sqrt(cs_err_TDE_12[nd,:]**2. + cs_err_TDE_13[nd,:]**2.)  
    
        cs_err_inspTDE_1213[nd,:]   = np.sqrt(cs_err_insp_1213[nd,:]**2. + cs_err_TDE_1213[nd,:]**2.)
    #-------------------------------------------------
    


    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------

    f = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1], hspace=0.0)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    
    #fig 1: Cross sections
    ax1.errorbar(sma_arr_AU[0,:],   cs_exchange[0,:],     yerr=cs_err_exchange[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],       cs_exchange[0,:],     marker='^', linestyle='', markersize=7.5, linewidth=1.5, label=r'Exchange', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:],   cs_coll_1213[0,:],    yerr=cs_err_coll_1213[0,:], marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],       cs_coll_1213[0,:],    marker='o', linestyle='', markersize=7.5, linewidth=1.5, label=r'WD-CO Collision', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:],   cs_inspTDE_1213[0,:],    yerr=cs_err_inspTDE_1213[0,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    ax1.plot(sma_arr_AU[0,:],       cs_inspTDE_1213[0,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'WD-CO Tidal Inspiral', mew=2, color = 'black')

    ax1.errorbar(sma_arr_AU[0,:],   cs_insp_23[0,:],      yerr=cs_err_insp_23[0,:],   marker='', linestyle=':', linewidth=1.0, color = 'black')
    ax1.plot(sma_arr_AU[0,:],       cs_insp_23[0,:],      marker='x', linestyle='', markersize=7.5, linewidth=5.0, label=r'CO-CO GW Inspiral', mew=2, color = 'black')

    #ax1.plot(sma_arr_AU[0,:],       cs_TDE_1213[0,:],      marker='x', linestyle='', markersize=7.5, linewidth=5.0, label=r'CO-CO GW Inspiral', mew=2, color = 'black')

    #oplot analytical guidelines:
    #hard binary limit:
    v_HB        = 10.0 #km/sec
    m_HB        = MO   #mass
    a_HB_limit  = ((36.5**2.)*(m_HB/v_HB**2.))  #HB limit SMA in units of AU
    ax1.plot([a_HB_limit,a_HB_limit], [1e-10,1e10],    linestyle=':', linewidth=0.5, color = 'black')
    ax1.text(1.5*a_HB_limit, 1e-2, r'HB-limit (10 km/sec)', size = 10,
            horizontalalignment='center',
            verticalalignment='center',
            rotation=-90)    
    ap = np.arange(min(sma_arr_AU[0,:]),a_HB_limit,0.001)
    tf = open('CS_analytical_test.txt',  "r")
    CS_analytical_test_arr = np.loadtxt(tf, dtype=float)
    tf.close()
    #savearr[0,:] = a0_over_RA_arr
    #savearr[1,:] = csp_insp_COTO_CORR
    #savearr[2,:] = csp_insp_GW
    #savearr[3,:] = csp_coll
    fac_tid     = ((MO)*(RTO)/(vel_cs**2))
    x_a0_tid    = (RTO/AU_U)*CS_analytical_test_arr[0,:]
    #collisions:
    #ana_cs_coll_1213    = 2.*fac_tid*CS_analytical_test_arr[3,:]
    #ax1.plot(x_a0_tid, ana_cs_coll_1213, linestyle=':', linewidth=1.5, color = 'black')
    #WD-NS inspirals:
    ana_cs_insp_1213    = 2.*fac_tid*CS_analytical_test_arr[1,:]
    ax1.plot(x_a0_tid, ana_cs_insp_1213,    linestyle='--', linewidth=1.0, color = 'black')    
    #NS-NS inspirals:
    Rs          = MO*(2950.0/R_sun_SI) #Schwarzschild radius Rs of MO in units of Rsun
    fac_gr      = ((MO)*(Rs)/(vel_cs**2))
    x_a0_gr    = (Rs/AU_U)*CS_analytical_test_arr[0,:]
    ana_cs_insp_23 = 1.*fac_gr*CS_analytical_test_arr[2,:]
    ax1.plot(x_a0_gr, ana_cs_insp_23,    linestyle='-.', linewidth=1.0, color = 'black')  
    
    #Fig 2: Ratios (+- tides)
    vA = cs_exchange[0,:]
    vB = cs_exchange[1,:]
    sA = cs_err_exchange[0,:]
    sB = cs_err_exchange[1,:]
    frac_err_exch   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.) 
    vA = cs_coll_1213[0,:]
    vB = cs_coll_1213[1,:]
    sA = cs_err_coll_1213[0,:]
    sB = cs_err_coll_1213[1,:]
    frac_err_coll   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.)
    #exchange:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          yerr=frac_err_exch[:],  marker='',  linestyle='--', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='',  linestyle='--', linewidth=1.0, label=r"Exchange", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='s', markersize=5.0, linestyle='', color='black')
    #coll:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_coll_1213[0,:]/cs_coll_1213[1,:],        yerr=frac_err_coll[:],  marker='',  linestyle='-', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll_1213[0,:]/cs_coll_1213[1,:],        marker='',  linestyle='-', linewidth=1.0, label=r"WD-CO Collision", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll_1213[0,:]/cs_coll_1213[1,:],        marker='s', markersize=5.0, linestyle='', color='black')

    plt.yticks([0.5,1.0,1.5])
    
    #Titles/labels:
    #ax1.set_title('WD,NS,NS')
    ax1.set_ylabel(r'cross section $\sigma$ [AU$^2$] (10 km/sec)')
    ax2.set_xlabel(r'semi-major axis $a_{0}$ [AU]')
    ax2.set_ylabel(r'$\sigma_{+}$/$\sigma_{-}$')

    ax2.plot([a_HB_limit,a_HB_limit], [-10,10],    linestyle=':', linewidth=0.5, color = 'black')
    
    #x-scales:
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    xlim = np.array([5e-4,5e1])
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)    
    #make upper x-axis showing log a/R:
    ax3 = ax1.twiny()
    xscalefac       = AU_U/RTO
    a0_over_RTO_arr = xscalefac*sma_arr_AU[0,:]
    xp = np.log10(a0_over_RTO_arr)
    yp = xp
    ax3.plot(xp, yp,    linewidth=0.0)  #dummy plot
    ax3.set_xlim(np.log10(xscalefac*xlim))
    ax3.set_xlabel(r'Compactness log$(a_{0}/R_{WD})$')

    #y-scales:
    ax1.set_yscale('log')
    ax1.set_ylim([1e-4,1e2])
    ax2.set_ylim([0.5,2.])

    #Legends:
    ax1.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)
    ax2.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False, ncol=2)

    #final axis settings:
    ax1.tick_params(
        axis='x',           # changes apply to the x,y-axis
        which='both',       # both major and minor ticks are affected
        bottom='on',        # ticks along the bottom edge are off
        top='off',           # ticks along the top edge are off
        labelbottom='off',  # labels along the bottom edge are off
        right='off',
        left='off',
        labelleft='off')
    
    ax1.annotate(r'CO$\rightarrow$[WD,CO]', xy=(.7, .925), fontsize = 15.0, xycoords='axes fraction', horizontalalignment='center', verticalalignment='center')
    
    #Save figure:
    plt.savefig('compare_cross_WDNSNS.eps', bbox_inches='tight')
    
    plt.show()


    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------------------
#INTERACTION: [[MS(1)-MS(2)]-MS(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 2):


    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    #mass of one of the (equal mass) objects
    MO              = 1.0   #in units of Msun
    #radius of idal object
    RTO             = 1.0   #in units of Rsun
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)
    #-------------------------------------------------

    #define:
    #cs_arr_dataset       = [cs_arr_data1,cs_arr_data2]
    #cs_err_arr_dataset   = [cs_err_arr_data1,cs_err_arr_data2]
    #data format: (for both cs and cs err)
    #save_arr_cs[0,:]    = ap_AU
    #save_arr_cs[1,:]    = cs_arr_AU[:,vi,1,0] + cs_arr_AU[:,vi,2,0]     #exchange
    #save_arr_cs[2,:]    = cs_arr_AU[:,vi,0,2]                           #collisions 12
    #save_arr_cs[3,:]    = cs_arr_AU[:,vi,1,2]                           #collisions 13
    #save_arr_cs[4,:]    = cs_arr_AU[:,vi,2,2]                           #collisions 23
    #save_arr_cs[5,:]    = cs_arr_AU[:,vi,0,3]                           #inspiral collision 12
    #save_arr_cs[6,:]    = cs_arr_AU[:,vi,1,3]                           #inspiral collision 13
    #save_arr_cs[7,:]    = cs_arr_AU[:,vi,2,3]                           #inspiral collision 23
    #save_arr_cs[8,:]    = cs_arr_AU[:,vi,0,7]                           #inspiral 12
    #save_arr_cs[9,:]    = cs_arr_AU[:,vi,1,7]                           #inspiral 13
    #save_arr_cs[10,:]   = cs_arr_AU[:,vi,2,7]                           #inspiral 23
    
     
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                    = len(cs_arr_dataset[0][0,:])
    sma_arr_AU              = np.zeros((2,nr_a), dtype=np.float64)

    cs_exchange             = np.zeros((2,nr_a), dtype=np.float64)
    cs_coll                 = np.zeros((2,nr_a), dtype=np.float64)
    cs_inspcoll             = np.zeros((2,nr_a), dtype=np.float64)
    cs_insp                 = np.zeros((2,nr_a), dtype=np.float64)
    
    cs_err_exchange         = np.zeros((2,nr_a), dtype=np.float64)
    cs_err_coll             = np.zeros((2,nr_a), dtype=np.float64)
    cs_err_inspcoll         = np.zeros((2,nr_a), dtype=np.float64)
    cs_err_insp             = np.zeros((2,nr_a), dtype=np.float64)
    #-------------------------------------------------


    #-------------------------------------------------
    #calc defined cs for each dataset:
    #-------------------------------------------------
    for nd in range(0,2):
    
        sma_arr_AU[nd,:]            = cs_arr_dataset[nd][0,:]
   
        #DATA: cs
        cs_exchange[nd,:]           = cs_fac*cs_arr_dataset[nd][1,:]
        cs_coll[nd,:]               = cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:] + cs_arr_dataset[nd][4,:])
        cs_inspcoll[nd,:]           = cs_fac*(cs_arr_dataset[nd][5,:] + cs_arr_dataset[nd][6,:] + cs_arr_dataset[nd][7,:])
        cs_insp[nd,:]               = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:] + cs_arr_dataset[nd][10,:])
    
        #DATA: cs-err
        cs_err_exchange[nd,:]       = cs_fac*cs_err_arr_dataset[nd][1,:]
        cs_err_coll[nd,:]           = cs_fac*np.sqrt(cs_err_arr_dataset[nd][2,:]**2. + cs_err_arr_dataset[nd][3,:]**2. + cs_err_arr_dataset[nd][4,:]**2.)
        cs_err_inspcoll[nd,:]       = cs_fac*np.sqrt(cs_err_arr_dataset[nd][5,:]**2. + cs_err_arr_dataset[nd][6,:]**2. + cs_err_arr_dataset[nd][7,:]**2.)
        cs_err_insp[nd,:]           = cs_fac*np.sqrt(cs_err_arr_dataset[nd][8,:]**2. + cs_err_arr_dataset[nd][9,:]**2. + cs_err_arr_dataset[nd][10,:]**2.)
    #-------------------------------------------------


    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------

    f = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1], hspace=0.0)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])


    #fig 1: Cross sections
    ax1.errorbar(sma_arr_AU[0,:], cs_exchange[0,:],     yerr=cs_err_exchange[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:], cs_exchange[0,:],     marker='^', linestyle='', markersize=7.5, linewidth=1.5, label=r'Exchange', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:], cs_coll[0,:],    yerr=cs_err_coll[0,:], marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:], cs_coll[0,:],    marker='o', linestyle='', markersize=7.5, linewidth=1.5, label=r'Collision', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:], cs_insp[0,:],    yerr=cs_err_insp[0,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    ax1.plot(sma_arr_AU[0,:], cs_insp[0,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'Tidal Inspiral', mew=2, color = 'black')

    #ax1.errorbar(sma_arr_AU[0,:], cs_inspcoll[0,:],    yerr=cs_err_inspcoll[0,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    #ax1.plot(sma_arr_AU[0,:], cs_inspcoll[0,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'Tidal Inspiral', mew=2, color = 'black')
    #ax1.errorbar(sma_arr_AU[1,:], cs_inspcoll[1,:],    yerr=cs_err_inspcoll[1,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    #ax1.plot(sma_arr_AU[1,:], cs_inspcoll[1,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'Tidal Inspiral', mew=2, color = 'black')

    #oplot analytical guidelines:
    #hard binary limit:
    v_HB        = 10.0 #km/sec
    m_HB        = MO   #mass
    a_HB_limit  = ((36.5**2.)*(m_HB/v_HB**2.))  #HB limit SMA in units of AU
    ax1.plot([a_HB_limit,a_HB_limit], [1e-10,1e10],    linestyle=':', linewidth=0.5, color = 'black')
    ax1.text(1.15*a_HB_limit, 1e0, r'HB-limit (10 km/sec)', size = 10,
            horizontalalignment='center',
            verticalalignment='center',
            rotation=-90)    
    #ap = np.arange(min(sma_arr_AU[0,:]),1e1,0.001)
    #collisions:
    #ana_cs_coll_1213 = 8.5 + 0.0*ap
    #ax1.plot(ap, ana_cs_coll_1213, linestyle=':', linewidth=1.5, color = 'black')
    #WD-NS inspirals:
    #ana_cs_insp_1213 = 21.0*ap**(1./3.)
    #ax1.plot(ap, ana_cs_insp_1213,  linestyle='--', linewidth=0.5, color = 'black')
    #NS-NS inspirals:
    #ana_cs_insp_23 = 1.75*ap**(2./7.)
    #ax1.plot(ap, ana_cs_insp_23,    linestyle='--', linewidth=0.5, color = 'black')

    #Fig 2: Ratios (+- tides)
    vA = cs_exchange[0,:]
    vB = cs_exchange[1,:]
    sA = cs_err_exchange[0,:]
    sB = cs_err_exchange[1,:]
    frac_err_exch   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.)   
    vA = cs_coll[0,:]
    vB = cs_coll[1,:]
    sA = cs_err_coll[0,:]
    sB = cs_err_coll[1,:]
    frac_err_coll   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.)   
    #exchange:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          yerr=frac_err_exch[:],  marker='',  linestyle='--', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='',  linestyle='--', linewidth=1.0, label=r"Exchange", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='s', markersize=5.0, linestyle='', color='black')
    #coll:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_coll[0,:]/cs_coll[1,:],        yerr=frac_err_coll[:],  marker='',  linestyle='-', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll[0,:]/cs_coll[1,:],        marker='',  linestyle='-', linewidth=1.0, label=r"MS-MS Collision", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll[0,:]/cs_coll[1,:],        marker='s', markersize=5.0, linestyle='', color='black')

    plt.yticks([0.5,1.0,1.5])
    
    #Titles/labels:
    #ax1.set_title('MS,MS,MS')
    ax1.set_ylabel(r'cross section $\sigma$ [AU$^2$] (10 km/sec)')
    ax2.set_xlabel(r'semi-major axis $a_{0}$ [AU]')
    ax2.set_ylabel(r'$\sigma_{+}$/$\sigma_{-}$')

    ax2.plot([a_HB_limit,a_HB_limit], [-10,10],    linestyle=':', linewidth=0.5, color = 'black')
    
    #x-scales:
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    xlim = np.array([0.4,2e1])
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)    
    #make upper x-axis showing log a/R:
    ax3 = ax1.twiny()
    xscalefac       = AU_U/RTO
    a0_over_RTO_arr = xscalefac*sma_arr_AU[0,:]
    xp = np.log10(a0_over_RTO_arr)
    yp = xp
    ax3.plot(xp, yp,    linewidth=0.0)  #dummy plot
    ax3.set_xlim(np.log10(xscalefac*xlim))
    ax3.set_xlabel(r'Compactness log$(a_{0}/R_{MS})$')

    #y-scales:
    ax1.set_yscale('log')
    ax1.set_ylim([1e-2,1e4])
    ax2.set_ylim([0.5,2.])

    #Legends:
    ax1.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)
    ax2.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False, ncol=2)

    #final axis settings:
    ax1.tick_params(
        axis='x',           # changes apply to the x,y-axis
        which='both',       # both major and minor ticks are affected
        bottom='on',        # ticks along the bottom edge are off
        top='off',           # ticks along the top edge are off
        labelbottom='off',  # labels along the bottom edge are off
        right='off',
        left='off',
        labelleft='off')
    
    ax1.annotate(r'MS$\rightarrow$[MS,MS]', xy=(.7, .925), fontsize = 15.0, xycoords='axes fraction', horizontalalignment='center', verticalalignment='center')
    
    #Save figure:
    plt.savefig('compare_cross_MSMSMS.eps', bbox_inches='tight')

    plt.show()


    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------




#--------------------------------------------------------------------------------------------------
#INTERACTION: [[MS(1)-NS(2)]-NS(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 3):


    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    #mass of one of the (equal mass) objects
    MO              = 1.0   #in units of Msun
    #radius of tidal object
    RTO             = 1.0
    print 'RTO: ', RTO
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)
    #-------------------------------------------------
    
    
    
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                        = len(cs_arr_dataset[0][0,:])
    sma_arr_AU                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)

    cs_exchange                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_coll_1213                = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspcoll_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_1213                = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_coll_23                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspcoll_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_insp_23                  = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_12                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_13                   = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_TDE_1213                 = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_inspTDE_1213             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    
    cs_err_exchange             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_coll_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspcoll_1213        = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_1213            = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_coll_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspcoll_23          = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_insp_23              = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_12               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_13               = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_TDE_1213             = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    cs_err_inspTDE_1213         = np.zeros((nr_dataset,nr_a), dtype=np.float64)
    #-------------------------------------------------


    #-------------------------------------------------
    #calc defined cs for each dataset:
    #-------------------------------------------------
    for nd in range(0,nr_dataset):
    
        sma_arr_AU[nd,:]            = cs_arr_dataset[nd][0,:]
   
        #DATA: cs
        #exchange:
        cs_exchange[nd,:]           = cs_fac*cs_arr_dataset[nd][1,:]
        #MS-NS:
        cs_coll_1213[nd,:]          = cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:])
        cs_inspcoll_1213[nd,:]      = cs_fac*(cs_arr_dataset[nd][5,:] + cs_arr_dataset[nd][6,:])
        cs_insp_1213[nd,:]          = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:])
        #NS-NS:
        cs_coll_23[nd,:]            = cs_fac*cs_arr_dataset[nd][4,:]
        cs_inspcoll_23[nd,:]        = cs_fac*cs_arr_dataset[nd][7,:]
        cs_insp_23[nd,:]            = cs_fac*cs_arr_dataset[nd][10,:]
        #TDE
        cs_TDE_12[nd,:]             = cs_fac*cs_arr_dataset[nd][15,:]        
        cs_TDE_13[nd,:]             = cs_fac*cs_arr_dataset[nd][16,:]        
        cs_TDE_1213[nd,:]           = cs_TDE_12[nd,:] + cs_TDE_13[nd,:]
        
        cs_inspTDE_1213[nd,:]       = cs_insp_1213[nd,:] + cs_TDE_1213[nd,:]
        
        #DATA: cs-err
        #exchange:
        cs_err_exchange[nd,:]       = cs_fac*cs_err_arr_dataset[nd][1,:]
        #MS-NS:
        cs_err_coll_1213[nd,:]      = cs_fac*np.sqrt(cs_err_arr_dataset[nd][2,:]**2. + cs_err_arr_dataset[nd][3,:]**2.)
        cs_err_inspcoll_1213[nd,:]  = cs_fac*np.sqrt(cs_err_arr_dataset[nd][5,:]**2. + cs_err_arr_dataset[nd][6,:]**2.)
        cs_err_insp_1213[nd,:]      = cs_fac*np.sqrt(cs_err_arr_dataset[nd][8,:]**2. + cs_err_arr_dataset[nd][9,:]**2.)
        #NS-NS:
        cs_err_coll_23[nd,:]        = cs_fac*cs_err_arr_dataset[nd][4,:]
        cs_err_inspcoll_23[nd,:]    = cs_fac*cs_err_arr_dataset[nd][7,:]
        cs_err_insp_23[nd,:]        = cs_fac*cs_err_arr_dataset[nd][10,:]
        #TDE
        cs_err_TDE_12[nd,:]         = cs_fac*cs_err_arr_dataset[nd][15,:]        
        cs_err_TDE_13[nd,:]         = cs_fac*cs_err_arr_dataset[nd][16,:]  
        cs_err_TDE_1213[nd,:]       = np.sqrt(cs_err_TDE_12[nd,:]**2. + cs_err_TDE_13[nd,:]**2.)  
    
        cs_err_inspTDE_1213[nd,:]   = np.sqrt(cs_err_insp_1213[nd,:]**2. + cs_err_TDE_1213[nd,:]**2.)
    #-------------------------------------------------


    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------
    
    f = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1], hspace=0.0)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    #fig 1: Cross sections
    ax1.errorbar(sma_arr_AU[0,:],   cs_exchange[0,:],     yerr=cs_err_exchange[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],       cs_exchange[0,:],     marker='^', linestyle='', markersize=7.5, linewidth=1.5, label=r'Exchange', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:],   cs_coll_1213[0,:],    yerr=cs_err_coll_1213[0,:], marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:],       cs_coll_1213[0,:],    marker='o', linestyle='', markersize=7.5, linewidth=1.5, label=r'MS-CO Collision', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:],   cs_inspTDE_1213[0,:],    yerr=cs_err_inspTDE_1213[0,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    ax1.plot(sma_arr_AU[0,:],       cs_inspTDE_1213[0,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'MS-CO Tidal Inspiral', mew=2, color = 'black')

    ax1.errorbar(sma_arr_AU[0,:],   cs_insp_23[0,:],      yerr=cs_err_insp_23[0,:],   marker='', linestyle=':', linewidth=1.0, color = 'black')
    ax1.plot(sma_arr_AU[0,:],       cs_insp_23[0,:],      marker='x', linestyle='', markersize=7.5, linewidth=5.0, label=r'CO-CO GW Inspiral', mew=2, color = 'black')

    #oplot analytical guidelines:
    #hard binary limit:
    v_HB        = 10.0 #km/sec
    m_HB        = MO   #mass
    a_HB_limit  = ((36.5**2.)*(m_HB/v_HB**2.))  #HB limit SMA in units of AU
    ax1.plot([a_HB_limit,a_HB_limit], [1e-10,1e10],    linestyle=':', linewidth=0.5, color = 'black')
    ax1.text(1.15*a_HB_limit, 1e0, r'HB-limit (10 km/sec)', size = 10,
            horizontalalignment='center',
            verticalalignment='center',
            rotation=-90)    
    ap = np.arange(min(sma_arr_AU[0,:]),a_HB_limit,0.001)
    tf = open('CS_analytical_test.txt',  "r")
    CS_analytical_test_arr = np.loadtxt(tf, dtype=float)
    tf.close()
    #savearr[0,:] = a0_over_RA_arr
    #savearr[1,:] = csp_insp_COTO_CORR
    #savearr[2,:] = csp_insp_GW
    #savearr[3,:] = csp_coll
    fac_tid     = ((MO)*(RTO)/(vel_cs**2))
    x_a0_tid    = (RTO/AU_U)*CS_analytical_test_arr[0,:]
    #collisions:
    #ana_cs_coll_1213    = 2.*fac_tid*CS_analytical_test_arr[3,:]
    #ax1.plot(x_a0_tid, ana_cs_coll_1213, linestyle=':', linewidth=1.5, color = 'black')
    #WD-NS inspirals:
    ana_cs_insp_1213    = 2.*fac_tid*CS_analytical_test_arr[1,:]
    ax1.plot(x_a0_tid, ana_cs_insp_1213,    linestyle='--', linewidth=1.0, color = 'black')    
    #NS-NS inspirals:
    Rs          = MO*(2950.0/R_sun_SI) #Schwarzschild radius Rs of MO in units of Rsun
    fac_gr      = ((MO)*(Rs)/(vel_cs**2))
    x_a0_gr    = (Rs/AU_U)*CS_analytical_test_arr[0,:]
    ana_cs_insp_23 = 1.*fac_gr*CS_analytical_test_arr[2,:]
    ax1.plot(x_a0_gr, ana_cs_insp_23,    linestyle='-.', linewidth=1.0, color = 'black')  
    

    #Fig 2: Ratios (+- tides)
    vA = cs_exchange[0,:]
    vB = cs_exchange[1,:]
    sA = cs_err_exchange[0,:]
    sB = cs_err_exchange[1,:]
    frac_err_exch   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.) 
    vA = cs_coll_1213[0,:]
    vB = cs_coll_1213[1,:]
    sA = cs_err_coll_1213[0,:]
    sB = cs_err_coll_1213[1,:]
    frac_err_coll   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.)
    #exchange:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          yerr=frac_err_exch[:],  marker='',  linestyle='--', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='',  linestyle='--', linewidth=1.0, label=r"Exchange", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='s', markersize=5.0, linestyle='', color='black')
    #coll:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_coll_1213[0,:]/cs_coll_1213[1,:],        yerr=frac_err_coll[:],  marker='',  linestyle='-', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll_1213[0,:]/cs_coll_1213[1,:],        marker='',  linestyle='-', linewidth=1.0, label=r"MS-CO Collision", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll_1213[0,:]/cs_coll_1213[1,:],        marker='s', markersize=5.0, linestyle='', color='black')

    plt.yticks([0.5,1.0,1.5])

    #Titles/labels:
    #ax1.set_title('WD,NS,NS')
    ax1.set_ylabel(r'cross section $\sigma$ [AU$^2$] (10 km/sec)')
    ax2.set_xlabel(r'semi-major axis $a_{0}$ [AU]')
    ax2.set_ylabel(r'$\sigma_{+}$/$\sigma_{-}$')

    ax2.plot([a_HB_limit,a_HB_limit], [-10,10],    linestyle=':', linewidth=0.5, color = 'black')
    
    #x-scales:
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    xlim = np.array([0.4,2e1])
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)    
    #make upper x-axis showing log a/R:
    ax3 = ax1.twiny()
    xscalefac       = AU_U/RTO
    a0_over_RTO_arr = xscalefac*sma_arr_AU[0,:]
    xp = np.log10(a0_over_RTO_arr)
    yp = xp
    ax3.plot(xp, yp,    linewidth=0.0)  #dummy plot
    ax3.set_xlim(np.log10(xscalefac*xlim))
    ax3.set_xlabel(r'Compactness log$(a_{0}/R_{MS})$')

    #y-scales:
    ax1.set_yscale('log')
    ax1.set_ylim([1e-2,1e4])
    ax2.set_ylim([0.5,2.])

    #Legends:
    ax1.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)
    ax2.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False, ncol=2)

    #final axis settings:
    ax1.tick_params(
        axis='x',           # changes apply to the x,y-axis
        which='both',       # both major and minor ticks are affected
        bottom='on',        # ticks along the bottom edge are off
        top='off',           # ticks along the top edge are off
        labelbottom='off',  # labels along the bottom edge are off
        right='off',
        left='off',
        labelleft='off')

    ax1.annotate(r'CO$\rightarrow$[MS,CO]', xy=(.7, .925), fontsize = 15.0, xycoords='axes fraction', horizontalalignment='center', verticalalignment='center')
                    
    #Save figure:
    plt.savefig('compare_cross_MSNSNS.eps', bbox_inches='tight')

    plt.show()


    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------







#--------------------------------------------------------------------------------------------------
#INTERACTION: [[WD(1)-WD(2)]-WD(3)]:
#--------------------------------------------------------------------------------------------------
if (DATASET == 4):


    #-------------------------------------------------
    #SET:
    #-------------------------------------------------
    #mass of one of the (equal mass) objects
    MO              = 1.2   #in units of Msun
    #radius of tidal object
    RTO             = 0.013*((1.43/MO)**(1./3.))*((1.-MO/1.43)**(0.447))
    #scale cross sections to this vel:
    vel_cs          = 10.0  #in km/sec
    cs_fac          = (1./vel_cs**2.)
    print 0.0041176640*(MO/1.0)**(-1./3.), 0.013*((1.43/MO)**(1./3.))*((1.-MO/1.43)**(0.447))
    #-------------------------------------------------
    
    #define:
    #cs_arr_dataset       = [cs_arr_data1,cs_arr_data2]
    #cs_err_arr_dataset   = [cs_err_arr_data1,cs_err_arr_data2]
    #data format: (for both cs and cs err)
    #save_arr_cs[0,:]        = ap_AU
    #save_arr_cs[1,:]        = cs_arr_AU[:,vi,2,0]                                       #NS-NS exchange
    #save_arr_cs[2,:]        = cs_arr_AU[:,vi,0,2]                                       #collisions 12
    #save_arr_cs[3,:]        = cs_arr_AU[:,vi,1,2]                                       #collisions 13
    #save_arr_cs[4,:]        = cs_arr_AU[:,vi,2,2]                                       #collisions 23
    #save_arr_cs[5,:]        = cs_arr_AU[:,vi,0,3]                                       #inspiral collision 12
    #save_arr_cs[6,:]        = cs_arr_AU[:,vi,1,3]                                       #inspiral collision 13
    #save_arr_cs[7,:]        = cs_arr_AU[:,vi,2,3]                                       #inspiral collision 23
    #save_arr_cs[8,:]        = cs_arr_AU[:,vi,0,7]                                       #inspiral 12
    #save_arr_cs[9,:]        = cs_arr_AU[:,vi,1,7]                                       #inspiral 13
    #save_arr_cs[10,:]       = cs_arr_AU[:,vi,2,7]                                       #inspiral 23
    #save_arr_cs[11,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_LT_tlife_limit_arr[:]              #NS-NS exchange with t_life<10**10 years
    #save_arr_cs[12,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[0,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
    #save_arr_cs[13,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[1,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
    #save_arr_cs[14,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[2,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
    #save_arr_cs[15,:]       = cs_arr_AU[:,vi,0,8]                                       #TDE 12
    #save_arr_cs[16,:]       = cs_arr_AU[:,vi,1,8]                                       #TDE 13
    #save_arr_cs[17,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[0,:]                #inspiral 23 with rmin cut
    #save_arr_cs[18,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[1,:]                #inspiral 23 with rmin cut
    #save_arr_cs[19,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[2,:]                #inspiral 23 with rmin cut  
     
    #-------------------------------------------------
    #define:
    #-------------------------------------------------
    nr_a                    = len(cs_arr_dataset[0][0,:])
    sma_arr_AU              = np.zeros((2,nr_a), dtype=np.float64)

    cs_exchange             = np.zeros((2,nr_a), dtype=np.float64)
    cs_coll                 = np.zeros((2,nr_a), dtype=np.float64)
    cs_inspcoll             = np.zeros((2,nr_a), dtype=np.float64)
    cs_insp                 = np.zeros((2,nr_a), dtype=np.float64)
    cs_TDE                  = np.zeros((2,nr_a), dtype=np.float64)
    
    cs_err_exchange         = np.zeros((2,nr_a), dtype=np.float64)
    cs_err_coll             = np.zeros((2,nr_a), dtype=np.float64)
    cs_err_inspcoll         = np.zeros((2,nr_a), dtype=np.float64)
    cs_err_insp             = np.zeros((2,nr_a), dtype=np.float64)
    #-------------------------------------------------


    #-------------------------------------------------
    #calc defined cs for each dataset:
    #-------------------------------------------------
    for nd in range(0,2):
    
        sma_arr_AU[nd,:]            = cs_arr_dataset[nd][0,:]
   
        #DATA: cs
        cs_exchange[nd,:]           = cs_fac*cs_arr_dataset[nd][1,:]
        cs_coll[nd,:]               = cs_fac*(cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:] + cs_arr_dataset[nd][4,:])
        cs_inspcoll[nd,:]           = cs_fac*(cs_arr_dataset[nd][5,:] + cs_arr_dataset[nd][6,:] + cs_arr_dataset[nd][7,:])
        cs_insp[nd,:]               = cs_fac*(cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:] + cs_arr_dataset[nd][10,:])
        cs_TDE[nd,:]                = cs_fac*(cs_arr_dataset[nd][15,:]+cs_arr_dataset[nd][16,:])        
        #seems like there are no TDEs in this case. So we dont need to show them. its ok.
        
        #DATA: cs-err
        cs_err_exchange[nd,:]       = cs_fac*cs_err_arr_dataset[nd][1,:]
        cs_err_coll[nd,:]           = cs_fac*np.sqrt(cs_err_arr_dataset[nd][2,:]**2. + cs_err_arr_dataset[nd][3,:]**2. + cs_err_arr_dataset[nd][4,:]**2.)
        cs_err_inspcoll[nd,:]       = cs_fac*np.sqrt(cs_err_arr_dataset[nd][5,:]**2. + cs_err_arr_dataset[nd][6,:]**2. + cs_err_arr_dataset[nd][7,:]**2.)
        cs_err_insp[nd,:]           = cs_fac*np.sqrt(cs_err_arr_dataset[nd][8,:]**2. + cs_err_arr_dataset[nd][9,:]**2. + cs_err_arr_dataset[nd][10,:]**2.)
    #-------------------------------------------------


    #-------------------------------------------------    
    #PLOT:
    #-------------------------------------------------

    f = plt.figure(figsize=(5,5))
    gs = gridspec.GridSpec(2, 1,height_ratios=[3,1], hspace=0.0)
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    print sma_arr_AU[0,:]
    #fig 1: Cross sections
    ax1.errorbar(sma_arr_AU[0,:], cs_exchange[0,:],     yerr=cs_err_exchange[0,:],  marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:], cs_exchange[0,:],     marker='^', linestyle='', markersize=7.5, linewidth=1.5, label=r'Exchange', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:], cs_coll[0,:],    yerr=cs_err_coll[0,:], marker='', linestyle=':', linewidth=1.0, color = 'grey')
    ax1.plot(sma_arr_AU[0,:], cs_coll[0,:],    marker='o', linestyle='', markersize=7.5, linewidth=1.5, label=r'Collision', mew=1, color = 'lightgrey')

    ax1.errorbar(sma_arr_AU[0,:], cs_insp[0,:],    yerr=cs_err_insp[0,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    ax1.plot(sma_arr_AU[0,:], cs_insp[0,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'Tidal Inspiral', mew=2, color = 'black')

    #ax1.errorbar(sma_arr_AU[0,:], cs_inspcoll[0,:],    yerr=cs_err_inspcoll[0,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    #ax1.plot(sma_arr_AU[0,:], cs_inspcoll[0,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'Tidal Inspiral', mew=2, color = 'black')
    #ax1.errorbar(sma_arr_AU[1,:], cs_inspcoll[1,:],    yerr=cs_err_inspcoll[1,:], marker='', linestyle=':', linewidth=1.0, color = 'black')
    #ax1.plot(sma_arr_AU[1,:], cs_inspcoll[1,:],    marker='+', linestyle='', markersize=7.5, linewidth=5.0, label=r'Tidal Inspiral', mew=2, color = 'black')

    #oplot analytical guidelines:
    #hard binary limit:
    v_HB        = 10.0 #km/sec
    m_HB        = MO   #mass
    a_HB_limit  = ((36.5**2.)*(m_HB/v_HB**2.))  #HB limit SMA in units of AU
    ax1.plot([a_HB_limit,a_HB_limit], [1e-10,1e10],    linestyle=':', linewidth=0.5, color = 'black')
    ax1.text(1.5*a_HB_limit, 1e-2, r'HB-limit (10 km/sec)', size = 10,
            horizontalalignment='center',
            verticalalignment='center',
            rotation=-90)
    #ap = np.arange(min(sma_arr_AU[0,:]),1e1,0.001)
    #collisions:
    #ana_cs_coll_1213 = 8.5 + 0.0*ap
    #ax1.plot(ap, ana_cs_coll_1213, linestyle=':', linewidth=1.5, color = 'black')
    #WD-NS inspirals:
    #ana_cs_insp_1213 = 21.0*ap**(1./3.)
    #ax1.plot(ap, ana_cs_insp_1213,  linestyle='--', linewidth=0.5, color = 'black')
    #NS-NS inspirals:
    #ana_cs_insp_23 = 1.75*ap**(2./7.)
    #ax1.plot(ap, ana_cs_insp_23,    linestyle='--', linewidth=0.5, color = 'black')

    #Fig 2: Ratios (+- tides)
    vA = cs_exchange[0,:]
    vB = cs_exchange[1,:]
    sA = cs_err_exchange[0,:]
    sB = cs_err_exchange[1,:]
    frac_err_exch   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.)   
    vA = cs_coll[0,:]
    vB = cs_coll[1,:]
    sA = cs_err_coll[0,:]
    sB = cs_err_coll[1,:]
    frac_err_coll   = np.sqrt((sA/vB)**2. + (vA*(sB/(vB**2.)))**2.)   
    #exchange:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          yerr=frac_err_exch[:],  marker='',  linestyle='--', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='',  linestyle='--', linewidth=1.0, label=r"Exchange", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_exchange[0,:]/cs_exchange[1,:],          marker='s', markersize=5.0, linestyle='', color='black')
    #coll:
    ax2.errorbar(   sma_arr_AU[0,:],    cs_coll[0,:]/cs_coll[1,:],        yerr=frac_err_coll[:],  marker='',  linestyle='-', linewidth=1.0, color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll[0,:]/cs_coll[1,:],        marker='',  linestyle='-', linewidth=1.0, label=r"WD-WD Collision", color='black')
    ax2.plot(       sma_arr_AU[0,:],    cs_coll[0,:]/cs_coll[1,:],        marker='s', markersize=5.0, linestyle='', color='black')

    plt.yticks([0.5,1.0,1.5])
    
    #Titles/labels:
    #ax1.set_title('MS,MS,MS')
    ax1.set_ylabel(r'cross section $\sigma$ [AU$^2$] (10 km/sec)')
    ax2.set_xlabel(r'semi-major axis $a_{0}$ [AU]')
    ax2.set_ylabel(r'$\sigma_{+}$/$\sigma_{-}$')

    ax2.plot([a_HB_limit,a_HB_limit], [-10,10],    linestyle=':', linewidth=0.5, color = 'black')

    #x-scales:
    ax1.set_xscale('log')
    ax2.set_xscale('log')
    xlim = np.array([5e-4,5e1])
    ax1.set_xlim(xlim)
    ax2.set_xlim(xlim)    
    #make upper x-axis showing log a/R:
    ax3 = ax1.twiny()
    xscalefac       = AU_U/RTO
    a0_over_RTO_arr = xscalefac*sma_arr_AU[0,:]
    xp = np.log10(a0_over_RTO_arr)
    yp = xp
    ax3.plot(xp, yp,    linewidth=0.0)  #dummy plot
    ax3.set_xlim(np.log10(xscalefac*xlim))
    ax3.set_xlabel(r'Compactness log$(a_{0}/R_{WD})$')

    #y-scales:
    ax1.set_yscale('log')
    ax1.set_ylim([1e-4,1e2])
    ax2.set_ylim([0.5,2.])

    #Legends:
    ax1.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)
    ax2.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False, ncol=2)

    #final axis settings:
    ax1.tick_params(
        axis='x',           # changes apply to the x,y-axis
        which='both',       # both major and minor ticks are affected
        bottom='on',        # ticks along the bottom edge are off
        top='off',           # ticks along the top edge are off
        labelbottom='off',  # labels along the bottom edge are off
        right='off',
        left='off',
        labelleft='off')
    
    ax1.annotate(r'WD$\rightarrow$[WD,WD]', xy=(.7, .925), fontsize = 15.0, xycoords='axes fraction', horizontalalignment='center', verticalalignment='center')
    
    #Save figure:
    plt.savefig('compare_cross_WDWDWD.eps', bbox_inches='tight')

    plt.show()


    exit()
    #-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------















#--------------------------------------------------------------------------------------------------
#GRINDLAY ANALYSIS: [[TO(1)-NS(2)]-NS(3)]
#--------------------------------------------------------------------------------------------------

#-------------------------------------------------
#SET:
#-------------------------------------------------
if (DATASET == 5):
    mTO     = 0.2
    RTO     = 0.018
    mNS     = 1.4
    mNSbin  = mNS+mNS
    mINbin  = mTO+mNS

if (DATASET == 8):
    mTO     = 0.8
    RTO     = mTO**(0.8)
    mNS     = 1.4
    mNSbin  = mNS+mNS
    mINbin  = mTO+mNS
#-------------------------------------------------

#-------------------------------------------------
#define:
#-------------------------------------------------
nr_a                    = len(cs_arr_dataset[0][0,:])
sma_arr_AU              = np.zeros((2,nr_a), dtype=np.float64)

cs_exchange             = np.zeros((2,nr_a), dtype=np.float64)
cs_coll_1213            = np.zeros((2,nr_a), dtype=np.float64)
cs_inspcoll_1213        = np.zeros((2,nr_a), dtype=np.float64)
cs_insp_1213            = np.zeros((2,nr_a), dtype=np.float64)
cs_coll_23              = np.zeros((2,nr_a), dtype=np.float64)
cs_inspcoll_23          = np.zeros((2,nr_a), dtype=np.float64)
cs_insp_23              = np.zeros((2,nr_a), dtype=np.float64)
cs_exchange_tlife       = np.zeros((2,nr_a), dtype=np.float64)
cs_exchange_tlife_rmin_1        = np.zeros((2,nr_a), dtype=np.float64)
cs_exchange_tlife_rmin_2        = np.zeros((2,nr_a), dtype=np.float64)
cs_exchange_tlife_rmin_3        = np.zeros((2,nr_a), dtype=np.float64)
cs_TDE_12               = np.zeros((2,nr_a), dtype=np.float64)
cs_TDE_13               = np.zeros((2,nr_a), dtype=np.float64)
cs_TDE_1213             = np.zeros((2,nr_a), dtype=np.float64)

cs_err_exchange         = np.zeros((2,nr_a), dtype=np.float64)
cs_err_coll_1213        = np.zeros((2,nr_a), dtype=np.float64)
cs_err_inspcoll_1213    = np.zeros((2,nr_a), dtype=np.float64)
cs_err_insp_1213        = np.zeros((2,nr_a), dtype=np.float64)
cs_err_coll_23          = np.zeros((2,nr_a), dtype=np.float64)
cs_err_inspcoll_23      = np.zeros((2,nr_a), dtype=np.float64)
cs_err_insp_23          = np.zeros((2,nr_a), dtype=np.float64)
cs_err_exchange_tlife   = np.zeros((2,nr_a), dtype=np.float64)
cs_err_exchange_tlife_rmin_1    = np.zeros((2,nr_a), dtype=np.float64)
cs_err_exchange_tlife_rmin_2    = np.zeros((2,nr_a), dtype=np.float64)
cs_err_exchange_tlife_rmin_3    = np.zeros((2,nr_a), dtype=np.float64)
cs_err_TDE_12           = np.zeros((2,nr_a), dtype=np.float64)
cs_err_TDE_13           = np.zeros((2,nr_a), dtype=np.float64)
cs_err_TDE_1213         = np.zeros((2,nr_a), dtype=np.float64)
#-------------------------------------------------


#-------------------------------------------------
#calc defined cs for each dataset:
#-------------------------------------------------
for nd in range(0,2):

    sma_arr_AU[nd,:]            = cs_arr_dataset[nd][0,:]

    #DATA: cs
    #exchange NS-NS:
    cs_exchange[nd,:]           = cs_arr_dataset[nd][1,:]
    #WD-NS:
    cs_coll_1213[nd,:]          = cs_arr_dataset[nd][2,:] + cs_arr_dataset[nd][3,:]
    cs_inspcoll_1213[nd,:]      = cs_arr_dataset[nd][5,:] + cs_arr_dataset[nd][6,:]
    cs_insp_1213[nd,:]          = cs_arr_dataset[nd][8,:] + cs_arr_dataset[nd][9,:]
    #NS-NS:
    cs_coll_23[nd,:]            = cs_arr_dataset[nd][4,:]
    cs_inspcoll_23[nd,:]        = cs_arr_dataset[nd][7,:]
    cs_insp_23[nd,:]            = cs_arr_dataset[nd][10,:]
    cs_exchange_tlife[nd,:]     = cs_arr_dataset[nd][11,:]        
    cs_exchange_tlife_rmin_1[nd,:]      = cs_arr_dataset[nd][12,:]
    cs_exchange_tlife_rmin_2[nd,:]      = cs_arr_dataset[nd][13,:]
    cs_exchange_tlife_rmin_3[nd,:]      = cs_arr_dataset[nd][14,:]
    cs_TDE_12[nd,:]             = cs_arr_dataset[nd][15,:]        
    cs_TDE_13[nd,:]             = cs_arr_dataset[nd][16,:]        
    cs_TDE_1213[nd,:]           = cs_TDE_12[nd,:] + cs_TDE_13[nd,:]
    
    #DATA: cs-err
    #exchange NS-NS:
    cs_err_exchange[nd,:]       = cs_err_arr_dataset[nd][1,:]
    #WD-NS:
    cs_err_coll_1213[nd,:]      = np.sqrt(cs_err_arr_dataset[nd][2,:]**2. + cs_err_arr_dataset[nd][3,:]**2.)
    cs_err_inspcoll_1213[nd,:]  = np.sqrt(cs_err_arr_dataset[nd][5,:]**2. + cs_err_arr_dataset[nd][6,:]**2.)
    cs_err_insp_1213[nd,:]      = np.sqrt(cs_err_arr_dataset[nd][8,:]**2. + cs_err_arr_dataset[nd][9,:]**2.)
    #NS-NS:
    cs_err_coll_23[nd,:]        = cs_err_arr_dataset[nd][4,:]
    cs_err_inspcoll_23[nd,:]    = cs_err_arr_dataset[nd][7,:]
    cs_err_insp_23[nd,:]        = cs_err_arr_dataset[nd][10,:]
    cs_err_exchange_tlife[nd,:] = cs_err_arr_dataset[nd][11,:]
    cs_err_exchange_tlife_rmin_1[nd,:]  = cs_err_arr_dataset[nd][12,:]
    cs_err_exchange_tlife_rmin_2[nd,:]  = cs_err_arr_dataset[nd][13,:]
    cs_err_exchange_tlife_rmin_3[nd,:]  = cs_err_arr_dataset[nd][14,:]
    cs_err_TDE_12[nd,:]         = cs_err_arr_dataset[nd][15,:]        
    cs_err_TDE_13[nd,:]         = cs_err_arr_dataset[nd][16,:]  
    cs_err_TDE_1213[nd,:]       = np.sqrt(cs_err_TDE_12[nd,:]**2. + cs_err_TDE_13[nd,:]**2.)
#-------------------------------------------------


#-------------------------------------------------    
#PLOT:
#-------------------------------------------------
#find range in a:
posa = np.where(sma_arr_AU[0,:]*AU_U/RTO > 2.0)[0]
#convert SMA a to orbital period P:
Pdays_arr   = Pdays_in_aAU_mbin(sma_arr_AU[0,:], mINbin)    #def Pdays_in_aAU_mbin(a_arr, mbin)
#print info:
print sma_arr_AU[0,:]*AU_U/RTO

#print Pdays_arr
#exit()
plot_paper  = 1
plot_test   = 0

    
if (plot_test == 1):
    
    fig = plt.figure(figsize=(5, 4))

    #rmin cut array from PV_analyze...:
    #rdist_limit_arr = np.array([2.0, 3.5, 5.0])   #in units R_TO - ONLY PUT 3 VALUES!!!

    #DATASET 1: (NO TIDES)
    dsi = 0 #specify dataset index 0 or 1 
    #NS-NS exchange:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange[dsi, posa],                  yerr=cs_err_exchange[dsi, posa],  marker='', linestyle=':', linewidth=0.5, alpha= 0.75, color = 'red')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange[dsi, posa],                      marker='^', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-NS', mew=1, color = 'red')
    #exchange, cut: t_life
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],            yerr=cs_err_exchange_tlife[dsi, posa],  marker='', linestyle=':', linewidth=1.0, alpha= 0.75, color = 'black')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],                marker='^', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-NS ($t<t_{*}$)', mew=1, color = 'black')
    #exchange, cut: t_life, rmin
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange_tlife_rmin_1[dsi, posa],     yerr=cs_err_exchange_tlife_rmin_1[dsi, posa],  marker='', linestyle=':', linewidth=1.0, alpha= 0.75, color = 'blue')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange_tlife_rmin_1[dsi, posa],         marker='^', linestyle='', markersize=2.5, linewidth=1.5, alpha= 0.75, label=r'NS-NS ($t<t_{*}$, $r<r_{*}$)', mew=1, color = 'blue')
    #NS-NS high ecc insp:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_insp_23[dsi, posa],                   yerr=cs_err_insp_23[dsi, posa],  marker='', linestyle=':', linewidth=0.5, alpha= 0.75, color = 'green')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_insp_23[dsi, posa],                       marker='^', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-NS GW inspiral', mew=1, color = 'green')
    #tidal insp:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_insp_1213[dsi, posa],                 yerr=cs_err_insp_1213[dsi, posa],  marker='', linestyle=':', linewidth=0.5, alpha= 0.75, color = 'purple')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_insp_1213[dsi, posa],                     marker='^', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-TO tidal inspiral', mew=1, color = 'purple')
    #collision:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_coll_1213[dsi, posa],                 yerr=cs_err_coll_1213[dsi, posa],  marker='', linestyle=':', linewidth=0.5, alpha= 0.75, color = 'pink')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_coll_1213[dsi, posa],                     marker='^', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-TO collision', mew=1, color = 'pink')

    #DATASET 2: (WITH TIDES)
    dsi = 1 #specify dataset index 0 or 1 
    #NS-NS exchange:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange[dsi, posa],                  yerr=cs_err_exchange[dsi, posa],  marker='', linestyle='-', linewidth=0.5, color = 'red')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange[dsi, posa],                      marker='s', linestyle='', markersize=7.5, linewidth=1.5, label=r'NS-NS', mew=1, color = 'red')
    #exchange, cut: t_life
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],            yerr=cs_err_exchange_tlife[dsi, posa],  marker='', linestyle='-', linewidth=1.0, alpha= 0.75, color = 'black')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],                marker='s', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-NS ($t<t_{*}$)', mew=1, color = 'black')
    #NS-NS high ecc insp:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_insp_23[dsi, posa],                   yerr=cs_err_insp_23[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha= 0.75, color = 'green')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_insp_23[dsi, posa],                       marker='s', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-NS GW inspiral', mew=1, color = 'green')
    #tidal insp:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_insp_1213[dsi, posa],                 yerr=cs_err_insp_1213[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha= 0.75, color = 'purple')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_insp_1213[dsi, posa],                     marker='s', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-TO tidal inspiral', mew=1, color = 'purple')
    #collision:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_coll_1213[dsi, posa],                 yerr=cs_err_coll_1213[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha= 0.75, color = 'pink')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_coll_1213[dsi, posa],                     marker='s', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.5, label=r'NS-TO collision', mew=1, color = 'pink')
    #TDE:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_TDE_1213[dsi, posa],                  yerr=cs_err_TDE_1213[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha= 0.75, color = 'grey')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_TDE_1213[dsi, posa],                      marker='s', linestyle='', markersize=7.5, linewidth=1.5, alpha= 0.75, label=r'NS-TO TDE', mew=1, color = 'grey')

    #axis settings:
    fig.add_subplot(111).set_xlim(5e-3,2e2)
    fig.add_subplot(111).set_ylim(1e-3,1e2)
    plt.xscale('log')
    plt.yscale('log')
    plt.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)
    plt.xlabel(r'P(days)')
    plt.ylabel(r'$\sigma$(AU$^2$)')




if (plot_paper == 1):
    
    fig = plt.figure(figsize=(5, 4))

    #rmin cut array from PV_analyze...:
    #rdist_limit_arr = np.array([2.0, 3.5, 5.0])   #in units R_TO - ONLY PUT 3 VALUES!!!

    #DATASET 1: (NO TIDES)
    dsi = 0 #specify dataset index 0 or 1 
    #NS-NS exchange:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange[dsi, posa],                  yerr=cs_err_exchange[dsi, posa],  marker='', linestyle=':', linewidth=0.5, alpha=0.5, color = 'red')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange[dsi, posa],                      marker='o', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-NS', mew=1, color = 'red')
    #exchange, cut: t_life
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],            yerr=cs_err_exchange_tlife[dsi, posa],  marker='', linestyle=':', linewidth=1.0, alpha=0.5, color = 'black')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],                marker='o', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-NS ($t<t_{*}$)', mew=1, color = 'black')
    #exchange, cut: t_life, rmin
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange_tlife_rmin_1[dsi, posa],     yerr=cs_err_exchange_tlife_rmin_1[dsi, posa],  marker='', linestyle=':', linewidth=1.0, alpha=0.5, color = 'blue')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange_tlife_rmin_1[dsi, posa],         marker='o', linestyle='', markersize=5.0, linewidth=1.0, alpha=0.5, label=r'NS-NS ($t<t_{*}$, $r<r_{*}$)', mew=1, color = 'blue')
    #collision:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_coll_1213[dsi, posa],                 yerr=cs_err_coll_1213[dsi, posa],  marker='', linestyle=':', linewidth=0.5, alpha=0.5, color = 'pink')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_coll_1213[dsi, posa],                     marker='o', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-TO collision', mew=1, color = 'pink')

    #DATASET 2: (WITH TIDES)
    dsi = 1 #specify dataset index 0 or 1 
    #NS-NS exchange:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange[dsi, posa],                  yerr=cs_err_exchange[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha=0.5, color = 'red')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange[dsi, posa],                      marker='*', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-NS', mew=1, color = 'red')
    #exchange, cut: t_life
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],            yerr=cs_err_exchange_tlife[dsi, posa],  marker='', linestyle='-', linewidth=1.0, alpha=0.5, color = 'black')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_exchange_tlife[dsi, posa],                marker='*', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-NS ($t<t_{*}$)', mew=1, color = 'black')
    #collision:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_coll_1213[dsi, posa],                 yerr=cs_err_coll_1213[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha=0.5, color = 'pink')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_coll_1213[dsi, posa],                     marker='*', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-TO collision', mew=1, color = 'pink')
    #TDE:
    fig.add_subplot(111).errorbar(Pdays_arr[posa], cs_TDE_1213[dsi, posa],                  yerr=cs_err_TDE_1213[dsi, posa],  marker='', linestyle='-', linewidth=0.5, alpha=0.5, color = 'grey')
    fig.add_subplot(111).plot(Pdays_arr[posa], cs_TDE_1213[dsi, posa],                      marker='*', linestyle='', markersize=10.0, linewidth=1.0, alpha=0.5, label=r'NS-TO disruption', mew=1, color = 'grey')
    
    #print Pdays_arr[posa]
    #exit()
    
    #oplot limits:
    limit_a_AU  = 2*RTO/AU_U 
    #convert SMA a to orbital period P:
    limit_P_days    = Pdays_in_aAU_mbin(limit_a_AU, mINbin)    #def Pdays_in_aAU_mbin(a_arr, mbin)
    fig.add_subplot(111).plot([limit_P_days,limit_P_days],[1e-10,1e10], marker='', linestyle=':', linewidth=0.25, color='black', label="2R$_{TO}$")
    
    #axis settings:
    if (DATASET == 5):
        plt.title(r'NS $\rightarrow$ [WD(0.2M$_{\odot}$), NS]')
        fig.add_subplot(111).set_xlim(5e-3,2e2)
        fig.add_subplot(111).set_ylim(1e-2,1e2)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)
        plt.xlabel(r'P(days)')
        plt.ylabel(r'$\sigma$(AU$^2$)')
    if (DATASET == 8):
        plt.title(r'NS $\rightarrow$ [MS(0.8M$_{\odot}$), NS]')
        fig.add_subplot(111).set_xlim(5e-3,2e2)
        fig.add_subplot(111).set_ylim(1e-1,1e2)
        plt.xscale('log')
        plt.yscale('log')
        plt.legend(loc='upper left', numpoints = 1, fontsize = 8.0, frameon = False)
        plt.xlabel(r'P(days)')
        plt.ylabel(r'$\sigma$(AU$^2$)')











#save and show:
if (DATASET == 5): save_name = 'compare_cross_JG_WD02NSNS.pdf'
if (DATASET == 8): save_name = 'compare_cross_JG_MS08NSNS.pdf'

#save figure:
plt.savefig(save_name, bbox_inches='tight')

plt.show()


exit()
#-------------------------------------------------
#--------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------


