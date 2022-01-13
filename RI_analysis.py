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
import matplotlib as mpl
from matplotlib import cm


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
#Data file names and folders:
#-----------------------------------------------------------------
data_folder                 = '/Users/jsamsing/Desktop/TIDES_PROJ/TEST_MC_OUT/'

#dataname:
#filenames_arr   = ['testRI081e5WGR']
#filenames_arr   = ['testRI141e5WGR']
#filenames_arr   = ['testRI141e5WGR_rptid10']
#filenames_arr   = ['RIana_WD02_T1', 'RIana_WD04_T1', 'RIana_WD06_T1', 'RIana_WD08_T1', 'RIana_WD10_T1','RIana_WD12_T1', 'RIana_WD14_T1']
#filenames_arr   = ['RIana_pp04_T1_rptid10', 'RIana_pp06_T1_rptid10', 'RIana_pp08_T1_rptid10', 'RIana_pp10_T1_rptid10', 'RIana_pp12_T1_rptid10','RIana_pp14_T1_rptid10']
#filenames_arr   = ['testRI141e1NGR_rptid10_2000']
#filenames_arr   = ['testRI141e1NGR_rptid10_2000_T1500_IMS10']
#filenames_arr   = ['testRI081e5WGR_rptid10']

filenames_arr   = ['NLRI_1']#['RI_BH102020_10000_WGR_1e4_1000Torb_TT01']#['Test1_BH102020_1e4']

#set:
#WD_mass_arr     = np.array([0.8])
#WD_mass_arr     = np.array([0.8])
#WD_mass_arr     = np.array([0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
#WD_mass_arr     = np.array([0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
#define:
nr_RI_files     = len(filenames_arr)
tf                          = open(data_folder+filenames_arr[0]+'_MCout_arr_endstate_info_INT.txt', "r")
output_Nbody_endstate_INT   = np.loadtxt(tf, dtype=float)
tf.close()
nrsim                       = 1.0*len(output_Nbody_endstate_INT)    #nr sim per bin-sin configuration

#RI_filenames                = ['TWD02_MCout_arr_resonance_analysis.txt', 'TWD04_MCout_arr_resonance_analysis.txt', 'TWD06_MCout_arr_resonance_analysis.txt', 'TWD08_MCout_arr_resonance_analysis.txt', 'TWD10_MCout_arr_resonance_analysis.txt', 'TWD12_MCout_arr_resonance_analysis.txt', 'TWD14_MCout_arr_resonance_analysis.txt', 'TO16_MCout_arr_resonance_analysis.txt', 'TO18_MCout_arr_resonance_analysis.txt', 'TO20_MCout_arr_resonance_analysis.txt', 'TO22_MCout_arr_resonance_analysis.txt', 'TO50_MCout_arr_resonance_analysis.txt']
#scatteringINFO_filenames    = ['TWD02_MCout_arr_scattering_params.txt', 'TWD04_MCout_arr_scattering_params.txt', 'TWD06_MCout_arr_scattering_params.txt', 'TWD08_MCout_arr_scattering_params.txt', 'TWD10_MCout_arr_scattering_params.txt', 'TWD12_MCout_arr_scattering_params.txt', 'TWD14_MCout_arr_scattering_params.txt', 'TO16_MCout_arr_scattering_params.txt', 'TO18_MCout_arr_scattering_params.txt', 'TO20_MCout_arr_scattering_params.txt', 'TO22_MCout_arr_scattering_params.txt', 'TO50_MCout_arr_scattering_params.txt']
#WD_mass_arr         = np.array([0.2,0.4,0.6,0.8,1.0,1.2,1.4,1.6,1.8,2.0,2.2,5.0])
#-----------------------------------------------------------------



#fig = plt.figure(figsize=(5, 4.5))
#a   = np.array([1e-5, 1e-4, 1e-3, 1e-3, 1e-2, 1e-1])
#cs = np.array([8.07933608e-05,   3.70665956e-04,   9.14686436e-04, 0.00076977,  0.00202638,  0.00369657])
#cscalc = cs[5]*(a[:]/a[5])**(2./7.)
#fig.add_subplot(111).plot(a, cs, marker="o", linestyle='none', markersize=8.0, alpha=1.0, color='red', label=r'Inspirals')                  
#fig.add_subplot(111).plot(a, cscalc, marker="o", linestyle='none', markersize=8.0, alpha=1.0, color='blue', label=r'Inspirals')                  
#plt.xscale('log')
#plt.yscale('log')
#plt.show()
#exit()




#-----------------------------------------------------------------
#EXAMPLE: [[WD(1),NS(2)], NS(3)]:
#-----------------------------------------------------------------

cut_val_e   = 0.0
cut_Etot    = 1.0

save_RI_calc_arr    = np.zeros((nr_RI_files,10), dtype=np.float64)

for fc in range(0,nr_RI_files):
    #--------------------------------
    #open data:      
    #--------------------------------
    filename = filenames_arr[fc]
    
    #open RI info:
    tf                          = open(data_folder+filename+'_MCout_arr_resonance_analysis.txt', "r")
    RI_info                     = np.loadtxt(tf, dtype=float)
    #[(0)pc, (1)sc, (2)rc, (3)0/1, (4)out_bin_i, (5)out_bin_j, (6)out_sin_k, (7)E_tot(ij), (8)a_bin(ij), (9)e_bin(ij), (10)E_tot(ijk), (11)a_bin(ijk), (12)e_bin(ijk)]
    tf.close()
	
    #open scattering info params:
    tf                          = open(data_folder+filename+'_MCout_arr_scattering_params.txt', "r")
    scatteringINFO_info         = np.loadtxt(tf, dtype=float)
    #[SMA_bin,vinf_sin,b_max,0,0,0,0,0,0,0]
    tf.close()
	
    #open endstate info:
    tf                          = open(data_folder+filename+'_MCout_arr_endstate_info_INT.txt', "r")
    output_Nbody_endstate_INT   = np.loadtxt(tf, dtype=float)
    #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
    tf.close()
    tf                          = open(data_folder+filename+'_MCout_arr_endstate_info_REAL.txt', "r")
    output_Nbody_endstate_REAL   = np.loadtxt(tf, dtype=float)
    #[E_kin(ij), E_pot(ij), E_tot(ij), a_bin(ij), e_bin(ij), E_kin(ijk), E_pot(ijk), E_tot(ijk), a_bin(ijk), e_bin(ijk)]
    tf.close()
    #--------------------------------
        
    
    #--------------------------------
    #sort data:
    #--------------------------------
    #RI info data:
    #[(0)pc, (1)sc, (2)rc, (3)0/1, (4)out_bin_i, (5)out_bin_j, (6)out_sin_k, (7)E_tot(ij), (8)a_bin(ij), (9)e_bin(ij), (10)E_tot(ijk), (11)a_bin(ijk), (12)e_bin(ijk)]
    #pos all formed (=1) IMSbin:
    pos_IMSbin              = np.where(RI_info[:,3] == 1)[0]        # = 1 indicates an IMS has formed (if detroyed then = 0)
    #pos all IMSbin with ecc > ecc_cut:
    pos_IMSbin_ecc_cut      = np.where(RI_info[:,9] > cut_val_e)[0] #ecc cut on all saved IMS binaries (0 and 1)
    #pos 12,13,23 formed (=1) IMSbin:        
    pos_IMSbin_12           = list(set(pos_IMSbin).intersection(np.where(RI_info[:,6] == 3)[0]))    #12 formed IMSbin
    pos_IMSbin_13           = list(set(pos_IMSbin).intersection(np.where(RI_info[:,6] == 2)[0]))    #13 formed IMSbin
    pos_IMSbin_23           = list(set(pos_IMSbin).intersection(np.where(RI_info[:,6] == 1)[0]))    #23 formed IMSbin
    pos_IMSbin_123          = pos_IMSbin                                                            #12,13,23 formed IMSbin
    #pos 12,13,23 formed (=1) IMSbin with ecc > ecc_cut:
    pos_IMSbin_12_cut_e     = list(set(pos_IMSbin_ecc_cut).intersection(pos_IMSbin_12))
    pos_IMSbin_13_cut_e     = list(set(pos_IMSbin_ecc_cut).intersection(pos_IMSbin_13))
    pos_IMSbin_23_cut_e     = list(set(pos_IMSbin_ecc_cut).intersection(pos_IMSbin_23))
    pos_IMSbin_123_cut_e    = list(set(pos_IMSbin_ecc_cut).intersection(pos_IMSbin_123))   
    
    #endstate output data:
    #[(0)out_end_state_flag, (1)out_bin_i, (2)out_bin_j, (3)out_sin_k, (4)out_IMS_bin_yesno, (5)IMS_rp_counter, (6)IMS_binsin_counter, ...]
    #pos id 10(no id), 5(inspiral), 2(coll):
    pos_id10            = np.where(output_Nbody_endstate_INT[:,0] == 10)[0] #did not finish
    pos_id5             = np.where(output_Nbody_endstate_INT[:,0] == 5)[0]  #inspirals
    pos_id2             = np.where(output_Nbody_endstate_INT[:,0] == 2)[0]  #collisions
    #pos endstate 12,13,23 (we cant use sin_k to sort with here):
    pos_endbin_i1       = np.where(output_Nbody_endstate_INT[:,1] == 1)[0]
    pos_endbin_i2       = np.where(output_Nbody_endstate_INT[:,1] == 2)[0]
    pos_endbin_j2       = np.where(output_Nbody_endstate_INT[:,2] == 2)[0]
    pos_endbin_j3       = np.where(output_Nbody_endstate_INT[:,2] == 3)[0]
    pos_endbin_12       = list(set(pos_endbin_i1).intersection(pos_endbin_j2))
    pos_endbin_13       = list(set(pos_endbin_i1).intersection(pos_endbin_j3))
    pos_endbin_23       = list(set(pos_endbin_i2).intersection(pos_endbin_j3))
    #pos 12,13,23 GW inspirals:
    pos_endbin_12_id5   = list(set(pos_id5).intersection(pos_endbin_12))
    pos_endbin_13_id5   = list(set(pos_id5).intersection(pos_endbin_13))
    pos_endbin_23_id5   = list(set(pos_id5).intersection(pos_endbin_23))
    print len(pos_endbin_12_id5), len(pos_endbin_13_id5), len(pos_endbin_23_id5), len(pos_id2)    
    #--------------------------------
    
        
    #--------------------------------
    #ASSUMING THIS HIERARCHY for further calculations:
    #--------------------------------
    vinf    = 10.   #km/sec (all simulations are (must be) done with 1km/sec!!!!!)
    a0      = scatteringINFO_info[0,0]      #Rsun
    a0AU    = scatteringINFO_info[0,0]/AU_U #AU
    m1      = 20.0 #10.#WD_mass_arr[fc]               #Msun
    m2      = 20.0 #20.#1.4                           #Msun
    m3      = 20.0 #20.#1.4                           #Msun
    mi      = m2
    mj      = m3
    mk      = m1
    m_bs    = m1+m2+m3
    m_ij    = mi+mj
    mu_ij   = mi*mj/m_ij
    mu_12   = m1*m2/(m1+m2)
    E_tot_0 = m1*m2/(2.*a0)
    acAU    = a0AU*((mi*mj)/(m1*m2)) 
    acRsun  = a0AU*((mi*mj)/(m1*m2))*AU_U
    #for GW inspirals: 
    Es      = np.pi*(85./96.)
    beta    = 7./2.
    Ms      = mu_ij
    Mfac    = ((m1*m2)/(mi*mj))*(((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
    RsMsun  = (2.*G_new_SI*(1.0*M_sun_SI)/(c_SI**2.))           #in m for 1Msun
    Rsm_ij  = (2.*G_new_SI*(m_ij*M_sun_SI)/(c_SI**2.))/AU_SI    #in AU for m_ij
    C_GW    = np.pi*((85.*np.pi/96.)**(2./7.))*((RsMsun/AU_SI)**(12./7.))*((c_SI/1000.)**2.)
    f_tid   = 0.5
    apu     = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
    insp_I  = (1.05/(1.-1.7/beta))*np.log(apu)*(apu-1.)**(-(1./(beta+1.)))
    Decc    = (1.-cut_val_e)    
    print C_GW, Mfac, apu, insp_I
    #--------------------------------
    
    
    #--------------------------------
    #GW INSP ANALYSIS: BIN-23
    #--------------------------------
    #Insp cs - ANALYTICAL:    
    #calc the 'N term':
    nr_IMSbin_23        = 1.0*len(pos_IMSbin_23)
    nr_IMSbin_23_cut_e  = 1.0*len(pos_IMSbin_23_cut_e)
    mean_N23            = nr_IMSbin_23/nrsim
    PR_uni23            = nr_IMSbin_23_cut_e/nr_IMSbin_23
    AR_uni23            = Decc*(apu-1.) 
    Nfac_23             = mean_N23*(PR_uni23/AR_uni23)
    #calc GW insp cs:
    cs_insp_23_CALC = C_GW*Nfac_23*insp_I*Mfac*(m_bs*(m_ij**(5./7.))*(a0AU**(2./7.))/(vinf**2.))*((m3/mu_12)**(1./3.))
    print 'cs_insp_23_CALC: (rptid) ', cs_insp_23_CALC
    cs_insp_23_CALC = C_GW*Nfac_23*insp_I*Mfac*(m_bs*(m_ij**(5./7.))*(a0AU**(2./7.))/(vinf**2.))*(mu_12/m1)
    print 'cs_insp_23_CALC: (rplm)  ', cs_insp_23_CALC
    
    #Insp cs - SIMULATION:
    nr_endbin_23_id5    = 1.*len(pos_endbin_23_id5) 
    bmax                = scatteringINFO_info[0,2]/AU_U
    cs_fac              = np.pi*(bmax**2.)/(vinf**2.)   #have to scale by /v^2 since bmax has been calc for 1km/sec.
    cs_insp_23_SIM      = (nr_endbin_23_id5/nrsim)*cs_fac
    print 'cs_insp_23_SIM:  ', cs_insp_23_SIM    
    
    print len(pos_IMSbin_12), len(pos_IMSbin_13), len(pos_IMSbin_23), len(pos_IMSbin)/1000.
    print len(pos_endbin_23_id5)
    print len(pos_id10)
    print 'mean_N23:    ', mean_N23
    print 'Nfac_23:     ', Nfac_23
    #exit()
    #SAVE info:
    save_RI_calc_arr[fc,0:3] = [mean_N23, (PR_uni23/AR_uni23), Nfac_23]#, Nfac_23*insp_I*Mfac, apu, cs_insp_23_CALC]
    #--------------------------------
    
    
    
    
    #--------------------------------
    #PLOTS:
    #--------------------------------
    plot_yesno_1    = 1
    if (plot_yesno_1 == 1):
        
        #ALL formed IMS:
        #energy cut: (used for plots)
        E_tot_IMS       = abs(RI_info[:,10])+abs(RI_info[:,7])
        E_tot_diff      = abs((E_tot_0-E_tot_IMS)/E_tot_0)
        pos_cut_Etot    = np.where(E_tot_diff[:] < cut_Etot)[0]
        pos_IMSbin_12_cut_Etot  = list(set(pos_cut_Etot).intersection(pos_IMSbin_12))    #12 formed IMSbin with dE/E0<...
        pos_IMSbin_13_cut_Etot  = list(set(pos_cut_Etot).intersection(pos_IMSbin_13))    #13 formed IMSbin with dE/E0<...
        pos_IMSbin_23_cut_Etot  = list(set(pos_cut_Etot).intersection(pos_IMSbin_23))    #23 formed IMSbin with dE/E0<...
    
        #INSPIRAL final IMS:
        insp_IMSbin_a  = []
        insp_IMSbin_e  = []
        #choose i,j set:
        pos_insp_list   = pos_endbin_23_id5 #pos wrt: output_Nbody_endstate_INT
        nrinsp = len(pos_insp_list)
        for ic in range(0,nrinsp):
            simid           = pos_insp_list[ic]                         #wrt: output_Nbody_endstate_INT
            pos_insp        = max(np.where(RI_info[:,1] == simid)[0])   #wrt: RI_info  
            #energy cut: (used for plots)
            insp_IMS_E_tot      = abs(RI_info[pos_insp,10])+abs(RI_info[pos_insp,7])
            insp_IMS_E_tot_diff = abs((E_tot_0-insp_IMS_E_tot)/E_tot_0)
            if (insp_IMS_E_tot_diff < cut_Etot):
                insp_a          = RI_info[pos_insp,8]
                insp_e          = RI_info[pos_insp,9]
                insp_IMSbin_a.append(insp_a)
                insp_IMSbin_e.append(insp_e)
    
        #TEST:
        pos         = list(pos_IMSbin_13) + list(pos_IMSbin_23)
        IMS_a       = RI_info[pos,8]/acRsun
        IMS_e       = RI_info[pos,9]
        pos_acut_1  = np.where(IMS_a < 2.0)[0] 
        IMS_a_acut  = IMS_a[pos_acut_1]
        IMS_e_acut  = IMS_e[pos_acut_1]
        
        #plt.hist(IMS_e_acut, bins=20, normed = True)  # arguments are passed to np.histogram
        #binsin_e   = list(output_Nbody_endstate_REAL[pos_endbin_13,4]) + list(output_Nbody_endstate_REAL[pos_endbin_23,4])
        #plt.hist(binsin_e, bins=20, normed = True)  # arguments are passed to np.histogram        
        #plt.show()
        
        
        #xarr    = (list(output_Nbody_endstate_REAL[pos_endbin_13,3]) + list(output_Nbody_endstate_REAL[pos_endbin_23,3]))/acRsun
        #pos     = np.where(xarr[:] > 1.0)[0]
        #RI info data:
        #[(0)pc, (1)sc, (2)rc, (3)0/1, (4)out_bin_i, (5)out_bin_j, (6)out_sin_k, (7)E_tot(ij), (8)a_bin(ij), (9)e_bin(ij), (10)E_tot(ijk), (11)a_bin(ijk), (12)e_bin(ijk)]
        #[E_kin(ij), E_pot(ij), E_tot(ij), a_bin(ij), e_bin(ij), E_kin(ijk), E_pot(ijk), E_tot(ijk), a_bin(ijk), e_bin(ijk)]
        xarr    = RI_info[pos_IMSbin,8]/acRsun
        xarr    = np.log10(1./xarr)
        plt.hist(xarr, range = [-2, 1], bins=30)  # arguments are passed to np.histogram        

        xarr    = output_Nbody_endstate_REAL[:,3]/acRsun
        xarr    = np.log10(1./xarr)
        plt.hist(xarr, range = [-2, 1], histtype = 'step', bins=30)  # arguments are passed to np.histogram        
        
        print acRsun
        
        plt.yscale('log')
        plt.show()
        
        
        exit()
        
        
        #PLOT:
        fig = plt.figure(figsize=(5, 4.0))


        #formed IMS:
        pos     = pos_IMSbin_23_cut_Etot
        IMS_a   = RI_info[pos,8]/acRsun
        IMS_e   = RI_info[pos,9]
        fig.add_subplot(111).plot(IMS_a, IMS_e, marker="+", linestyle='none', markersize=1.5, alpha=1.0, color='grey', label=r'all IMS binaries')                  

        #oplot insp line:
        ap_arr  = np.linspace(1.0, 3.0, num=10000)
        ecc     = 1. - (Es**(1./beta))*Mfac*((a0AU/Rsm_ij)**(1./beta - 1.))*((ap_arr**(1./beta-1.))*(ap_arr-1.)**(-3./(2.*beta)))
        fig.add_subplot(111).plot(ap_arr, ecc, linestyle='-', linewidth=1.0, color='black', label=r'inspiral boundary')                  
        arr_1   = 1.0+0.0*ap_arr
        fig.add_subplot(111).fill_between(ap_arr, ecc, arr_1, where=(ecc < arr_1),                   facecolor='lightgrey',     interpolate=True)

        #insp IMS:
        IMS_a   = insp_IMSbin_a[:]/acRsun
        IMS_e   = insp_IMSbin_e[:]
        fig.add_subplot(111).plot(IMS_a, IMS_e, marker="o", linestyle='none', markersize=3.5, alpha=1.0, color='black', label=r'GW inspirals')                  

        #insp IMS:        
        IMS_a   = output_Nbody_endstate_REAL[pos_id5,3]/acRsun
        IMS_e   = output_Nbody_endstate_REAL[pos_id5,4]
        fig.add_subplot(111).plot(IMS_a, IMS_e, marker="o", linestyle='none', markersize=3.5, alpha=1.0, color='red', label=r'test')                  

        
                
        #examples insp line::
        #ex_a0AU = 1e-6
        #ecc     = 1. - (Es**(1./beta))*Mfac*((ex_a0AU/Rsm_ij)**(1./beta - 1.))*((ap_arr**(1./beta-1.))*(ap_arr-1.)**(-3./(2.*beta)))
        #fig.add_subplot(111).plot(ap_arr, ecc, linestyle=':', linewidth=1.5, color='blue', label=r'example: Inspiral boundary ($10^{-4}$ AU)')                  
        #ex_a0AU = 1e-5
        #ecc     = 1. - (Es**(1./beta))*Mfac*((ex_a0AU/Rsm_ij)**(1./beta - 1.))*((ap_arr**(1./beta-1.))*(ap_arr-1.)**(-3./(2.*beta)))
        #fig.add_subplot(111).plot(ap_arr, ecc, linestyle=':', linewidth=1.5, color='grey', label=r'example: Inspiral boundary ($10^{-5}$ AU)')
        #collision:
        #...
        #labels, legends, etc:
        fig.add_subplot(111).set_xlabel(r'$a/a_{c}$')
        fig.add_subplot(111).set_ylabel(r'$e$')
        fig.add_subplot(111).set_title(r'GW inspirals in orbital phase space')
        fig.add_subplot(111).set_xlim(1.0,2.0)
        fig.add_subplot(111).set_ylim(0.0,1.0)
        fig.add_subplot(111).legend(loc='lower right', numpoints = 1, fontsize = 10.0, frameon = True)
        #save and show:
        #plt.savefig('ae_sim_allIMS_insp_ill.eps', bbox_inches='tight')        
        plt.show()
        exit()
    
    #test:
    #fig = plt.figure(figsize=(5, 4.5))
    #pos = pos_WD_sin#pos_IMSbin
    #fig.add_subplot(111).plot(RI_info[pos,8], RI_info[pos,9], marker="o", linestyle='none', markersize=5.0, alpha=0.5, color='blue')                  
    #plt.show()
    #exit()
    #ascale = 0.2*(1.4*1.4)/(WD_mass_arr[fc]*1.4) 
    #if (fc == 6):
    #    print 2.*((1.4*1.4/(WD_mass_arr[fc]+1.4+1.4))**(1./3.))*WD_mass_arr[fc]/1.4 + 1.
    #    fig.add_subplot(111).plot(RI_info[pos_WD_sin,8]/ascale, RI_info[pos_WD_sin,9], marker="o", linestyle='none', markersize=5.0, alpha=0.5, color='red')                  
    #    plt.show()
    
    plot_yesno_2    = 0
    if (plot_yesno_2 == 1):
    
        if (fc == 0):
            #[0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0, 2.2, 5.0]
            #define and calc:
            m1  = WD_mass_arr[fc]
            m2  = 1.4
            m3  = 1.4
            mi  = m2
            mj  = m3
            mk  = m1
            a0_Rsun = (5e-6)*AU_U   #for insp area illustration
            a_c     = a0_Rsun*(mi*mj)/(m1*m2) 
            m_bs    = m1+m2+m3
            m_ij    = mi+mj
            mu_ij   = mi*mj/m_ij
            #for GW inspirals: 
            beta    = 7./2.
            Ms      = mu_ij
            Rs      = (2.*G_new_SI*(m_ij*M_sun_SI)/(c_SI**2.))/R_sun_SI
            Es      = np.pi*(85./96.)
            Mfac    = ((m1*m2)/(mi*mj))*(((Ms/m_bs)**2.)*((m_bs/mu_ij)**(3./2.))*((mk*mk)/(m1*m2))*((m_ij/mk)**(1./2.)))**(1./beta)
            f_tid   = 0.5
            apu     = ((f_tid/2.)**(1./3.))*((mk/mu_ij)**(2./3.)) + 1.
            ap_arr  = np.linspace(1.0, apu, num=10000)
            ecc     = 1. - Mfac*(Es**(1./beta))*((a0_Rsun/Rs)**(1./beta - 1.))*((ap_arr**(1./beta-1.))*(ap_arr-1.)**(-3./(2.*beta)))
        
            #PLOT:
            fig = plt.figure(figsize=(5, 4.5))
            #a,e dist:
            ascale_plot = scatteringINFO_info[0,0]*(mi*mj)/(m1*m2) 
            ap_dist     = RI_info[pos_WD_sin,8]/ascale_plot
            ecc_dist    = RI_info[pos_WD_sin,9]
            pos_lt_apu  = np.where(ap_dist[:] < apu)[0]
            pos_gt_apu  = np.where(ap_dist[:] > apu)[0]
            fig.add_subplot(111).plot(ap_dist[pos_lt_apu], ecc_dist[pos_lt_apu], marker="o", linestyle='none', markersize=2.0, color='black')                  
            fig.add_subplot(111).plot(ap_dist[pos_gt_apu], ecc_dist[pos_gt_apu], marker="o", linestyle='none', markersize=1.0, color='grey')                  
            #insp a,e:        
            fig.add_subplot(111).plot(ap_arr, ecc, linestyle='--', linewidth=2.5, color='black')                  
            #ecc limit line:
            ecc_limit   = 0.7
            fig.add_subplot(111).plot([1,apu], [ecc_limit,ecc_limit], linestyle='-', linewidth=2.0, color='black')                  
            fig.add_subplot(111).plot([apu,apu], [0,1], linestyle='-', linewidth=2.0, color='black')                          
            #shade areas:
            fig.add_subplot(111).fill([1, apu, apu, 1], [ecc_limit, ecc_limit, 1, 1], fill=False, linestyle='-', linewidth=0.0, hatch='\\')
            fig.add_subplot(111).fill_between(ap_arr, ecc, 1, color='grey')
            fig.add_subplot(111).fill_between([apu,10], 0, 1, color='lightgrey')
            #labels etc:
            fig.add_subplot(111).set_xlabel(r'$a/a_{c}$')
            fig.add_subplot(111).set_ylabel(r'$e$')
            fig.add_subplot(111).set_ylim(0.0,1.0)
            fig.add_subplot(111).set_xlim(1.0,apu)
            #save and show:
            #plt.savefig('a_e_sim_insparea_ill.eps', bbox_inches='tight')        
            plt.show()
            exit()
    #--------------------------------    
    
    
    
exit()  
    
    
    
#--------------------------------
#PLOTS:
#--------------------------------
#save_RI_calc_arr[fc,0:6] = [mean_N23, (PR_uni23/AR_uni23), Nfac_23]#...

fig = plt.figure(figsize=(5, 3.0))

fig.add_subplot(111).plot(WD_mass_arr[:], save_RI_calc_arr[:,0], marker="o", linestyle='-', markersize=3.0, linewidth=2.0, color='black', label=r'$\langle N \rangle$')                  
fig.add_subplot(111).plot(WD_mass_arr[:], save_RI_calc_arr[:,1], marker="o", linestyle='-', markersize=3.0, linewidth=2.0, color='royalblue', label=r'$P(R_{uni})/A_{uni}$')                  
fig.add_subplot(111).plot(WD_mass_arr[:], save_RI_calc_arr[:,2], marker="o", linestyle='-', markersize=3.0, linewidth=2.0, color='pink', label=r'$\mathcal{N} = \langle N \rangle P(R_{uni})/A_{uni}$')                  
fig.add_subplot(111).plot([-10,10], [1,1], marker="", linestyle=':', markersize=2.0, alpha=1.0, color='black')                  

#labels, legends, etc:
fig.add_subplot(111).set_title(r'Weight factor $\mathcal{N}$(NS,NS) - ([WD,NS],NS)')
fig.add_subplot(111).set_xlabel(r'white dwarf mass $m_{WD}/M_{\odot}$')
fig.add_subplot(111).set_ylabel(r'value')
fig.add_subplot(111).set_xlim(0.0, 1.6)
fig.add_subplot(111).set_ylim(0.5, 25.)
plt.yscale('log')
#plt.xscale('log')
fig.add_subplot(111).legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = True)
#save and show:
#plt.savefig('WDNS_weightfac_N.eps', bbox_inches='tight')        
plt.show()

print save_RI_calc_arr[:,2]

exit()
#--------------------------------    
    
    
    
#-----------------------------------------------------------------
#-----------------------------------------------------------------
    
   
   
   
   
   
   
   
   
   
   



   











    
    
    
    
