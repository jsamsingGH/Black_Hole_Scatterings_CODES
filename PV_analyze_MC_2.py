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


#----------------------------------------------------------
#functions:
#----------------------------------------------------------
def tlife_years_in_aAU_e(a_arr, e_arr, m1, m2):
    m = m1
    tlife_years = (1.6*(10.**(17.)))*(a_arr**4.)*(m**(-3.))*(768./425.)*(1.-e_arr**2)**(7./2.) 
    return tlife_years
    
def Pdays_in_aAU_mbin(a_arr, mbin):
    Pdays = 2.*np.pi*np.sqrt((a_arr*AU_SI)**3./(G_new_SI*mbin*M_sun_SI))/(24.*3600.)
    return Pdays


def f_GW_Hz_high_ecc(rp_Rsun_arr, mtot_Msun):   #ONLY FOR ecc ~ 1 !!!!
    f_GW    = ((2.0**(1.1954))/(2.0**(1.5)))*(1./np.pi)*np.sqrt(G_new_SI*mtot_Msun*M_sun_SI)*(1./((rp_Rsun_arr*R_sun_SI)**(3./2.)))
    return f_GW


#GW inspiral life time:
tf = open('a_e_tlife_yrs_interpoltable.txt', "r")
a_e_tlife_yrs_arr   = np.loadtxt(tf, dtype=float)
tf.close()
#read data and make interpolation function:
tdata_e_arr                 = a_e_tlife_yrs_arr[:,1]
tdata_t_arr                 = a_e_tlife_yrs_arr[:,2]
ecc_to_tdata_interpol_func  = interpolate.interp1d(tdata_e_arr, tdata_t_arr)  
#----------------------------------------------------------   

#----------------------------------------------------------
#Units and conversions:
#----------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI       = 299792458.0        #m/s
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#----------------------------------------------------------
#Read in data from files:
#----------------------------------------------------------
#data name:

data_name = 'NIC_WD06NS14NS14_A5'#'Test_NIC_NS14NS14NS14_PN1'#'Test_NIC_NS14NS14NS14_C3'#'test_NI_3NS'#'SIM3_WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a_at6'
#data_name = 'GR1_WD05_NSNS_WT_NGR_001to100d_10kmsec_10000'#'SIM1_BH102020_NT_WGR_0001_10_1kmsec_T2500' #'GR1_WD05_NSNS_NT_NGR_001to100d_10kmsec_10000'

#GR:
#'GR1_WD02_NSNS_WT_WGR_001to100d_10kmsec_10000'
#'GR1_WD02_NSNS_NT_WGR_001to100d_10kmsec_10000'

#'GR1_MS08_NSNS_WT_WGR_001to100d_10kmsec_10000'
#'GR1_MS08_NSNS_NT_WGR_001to100d_10kmsec_10000'

#'GR1_MS02_NSNS_WT_WGR_001to100d_10kmsec_10000'
#'GR1_MS02_NSNS_NT_WGR_001to100d_10kmsec_10000'

#'GR1_MS04_NSNS_WT_WGR_001to100d_10kmsec_10000'
#'GR1_MS04_NSNS_NT_WGR_001to100d_10kmsec_10000'

#'GR1_WD005_NSNS_WT_WGR_001to100d_10kmsec_10000'
#'GR1_WD005_NSNS_NT_WGR_001to100d_10kmsec_10000'

#For new (J) GR sim:
#    use_25PN                = 1
#    evolvetides_threshold   = 0.1
#    tidaldisrup_threshold   = 3.0
#    insp_threshold          = 4.0

#for paper:
#'SIM3_WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a_at6'
#'SIM3_WD12pp12pp12_NT_0001to01AU_1kmsec_50000_rp1a_at6'

#'SIM3_WDWDWD12_WT_NGR_0001_01_1kmsec_rp1a_at6'
#'SIM3_WDWDWD12_NT_NGR_0001_01_1kmsec_rp1a_at6'

#'SIM1_MS10pp10pp10_WT_05to5AU_1kmsec_50000_rp1a_at6'
#'SIM1_MS10pp10pp10_NT_05to5AU_1kmsec_50000_rp1a_at6'


#'SIM2_WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a'
#'SIM2_WDWDWD12_WT_NGR_0001_01_1kmsec_rp1a'

#STATUS:

#WAITING:
#SIM1_WD04_NSNS_WT_WGR_001_1_1kmsec_T1000 #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#SIM1_WD06_NSNS_WT_WGR_001_1_1kmsec_T1000 #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#SIM1_WD10_NSNS_WT_WGR_001_1_1kmsec_T1000 #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#WAITING:
#GR1_WD04_WD06_14NS_NT_NGR_001to100d_10kmsec_10000
#GR1_WD05_NSNS_WT_NGR_001to100d_10kmsec_10000




#ANALYZE:!!!!!
#GR1_WD05_NSNS_NT_NGR_001to100d_10kmsec_10000






#'SIM1_MS10MS10MS10_WT_05to5AU_1kmsec_50000_rp1a_at6'
#'SIM2_MS10pp10pp10_WT_05to5AU_1kmsec_50000_rp1a_at6'
#WAITING: 'SIM1_MS10MS10MS10_NT_05to5AU_1kmsec_50000_rp1a_at6'

#TESTS:
#'SIM4_WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a_at6_64C' - ran out. try with more CPUs
#WAITING: 'SIM4_WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a_at6_256C'

#Paper 2:
#'SIM2_WD10_NSNS_WT_WGR_0001_01_1kmsec_T1000' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#'SIM2_WD06_NSNS_WT_WGR_0001_01_1kmsec_T1000' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#'SIM2_WD04_NSNS_WT_WGR_0001_01_1kmsec_T1000' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.

#'SIM3_WD10_NSNS_WT_WGR_0001_01_1kmsec_T1000' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#'SIM3_WD06_NSNS_WT_WGR_0001_01_1kmsec_T1000' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.
#'SIM3_WD04_NSNS_WT_WGR_0001_01_1kmsec_T1000' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 1000.









#'WD12pp12pp12_WT_0001to01AU_1kmsec_50000_rp1a' #evolvetides_threshold   = 0.1, Tsim_max_unit_Torb  = 350, WGR
#'SIM1_WDWDWD12_WT_NGR_0001_01_1kmsec_rp1a' #evolvetides_threshold   = 0.05, Tsim_max_unit_Torb  = 350 






#'SIM1_WD10_NSNS_WT_WGR_0001_01_1kmsec' #350 orb
#'SIM1_WD06_NSNS_WT_WGR_0001_01_1kmsec' #350 orb
#'SIM1_WD04_NSNS_WT_WGR_0001_01_1kmsec' #350 orb

#for no tides: Torb is 2500
#'SIM2_WD04_NSNS_NT_WGR_0001_01_1kmsec'#'SIM1_pp04_NSNS_NT_WGR_0001_01_1kmsec'
#'SIM2_WD14_NSNS_NT_WGR_0001_01_1kmsec'#'SIM1_WD08_NSNS_NT_WGR_0001_01_1kmsec'




#'R1_WD14_NSNS_NT_WGR_0001_01_1kmsec'#'R4_WD04_NSNS_NT_WGR_0001_01_1kmsec'#'WD04_NSNS_NT_WGR_0001_01_1kmsec'#'TEST_WD02_NSNS_NT_WGR_00003_01_1kmsec'#'TEST2_WD06_NSNS_NT_WGR_00003_01_1kmsec'
#'WD02_NSNS_NT_WGR_00003_01_1kmsec'#'WD02_NSNS_NT_WGR_00003_01_10kmsec'#'WD02_NSNS_WT_WGR_00003_0001_01_1kmsec'
#'WD06_NSNS_WT_WGR_00003_0001_01_1kmsec'#'WD14_NSNS_NT_WGR_m35_m1_10kmsec_ND1'#'WD02_NSNS_WT_WGR_00003_0001_01_1kmsec'
#'WD14_NSNS_WT_WGR_00003_0001_01_10kmsec'#'WD14_NSNS_WT_WGR_00003_0001_01_1kmsec'#'WD14_NSNS_NT_WGR_00003_0001_01_1kmsec'#'WD14_NSNS_WT_WGR_00003_0001_01_1kmsec'

#range of WD masses:
#'WD14_NSNS_WT_WGR_00003_0001_01_1kmsec'#'WD10_NSNS_WT_WGR_00003_0001_01_1kmsec'#'WD06_NSNS_WT_WGR_00003_0001_01_1kmsec'#'WD02_NSNS_WT_WGR_00003_0001_01_1kmsec'

#Grindlay analysis:
#data_name = 'JG_08MS_WGR_NT_tde3_insp2'#'WD08_NSNS_WT_WGR_0001_01_T1'#'JG_02WD_WGR_NT_tde3_insp2'
#'JG_08MS_WGR_WT_tde3_insp2'#'JG_04MS_WGR_WT_tde3_insp2'#'JG_005WD_WGR_NT_tde3_insp2'#'JG_02WD_WGR_WT_tde3_insp2'

#'JG_005WD_WGR_WT'#'JG_02WD_WGR_NT'#'JG_005WD_WT'#'JG_02WD_WT'#'JG_02WD_NT_T2'#'JG_02WD_NT'#'JG_02WD_WT'#'JG_02WD_NT'#'JG_02WD_NT_T1'#'BH_T1'#'T1'#'JG_02WD_WT_T1'#'JG_02WD_NT'#'T1'#'JG_04_T3'
#'JG_T1'#'BHI_T3'#'T1'#'CEL_T1'#'NT1'#'WD_NT_T1'

#'WD10WD10WD10_NTWGR_0001to01AU_1kmsec_50000'#'MS10NS10NS10_NTWGR_05to5AU_1kmsec_50000'
#data_name = 'MS10NS10NS10_WTWGR_05to5AU_1kmsec_50000'
#data_name = '1MS1MS1MS_WTWGR_05to5AU_1kmsec_50000'
#'WD14NSNS_NT_0001to01AU_1kmsec_50000'#'1MS1MS1MS_NTNGR_05to5AU_1kmsec_50000'#'WD14NSNS_WT_0001to01AU_1kmsec'#'1MS1MS1MS_WTWGR_05to5AU_1kmsec_100000'#'1MS1MS1MS_WTWGR_05to5AU_1kmsec_50000'
#'1MS1MS1MS_WTWGR_05to5AU_1kmsec_50000'#'1MS1MS1MS_WTWGR_05to5AU_1kmsec_50000'#'testSIM3'#'testSIM1'#'WD14NSNS_WTWGR_0001to1AU_1kmsec'
#'test'#'WD14NSNS_WT_0001to01AU_1kmsec'#'WD14NSNS_WT_0001to01AU_1kmsec'
#NS14MS10NS14_WT' 'MSMS10NS_005_5000_NT'#'MSMS10NS_005_5000_WT'#'3MS_005_5000_NT'#'TT6_pointmass'
#'TT6_pointmass'#'TT6_NT'#'GR_NSWD_1_5000'#'pointmass_NSWD_1_5000'#'GR_NSWD_1_5000'#'GRT_NSWD05'
#'GRT_NSWD05'#'GRT_NSWD_1_5000'#'GRT_NSWD05'#'GR_NSWD_1_5000'#'PH1'#'T1'#'TT5'#'TT6'#'test'#'TT2'#'notides1_'#'WT_'

#WD14NSNS_NT_0001to01AU_1kmsec_50000
#WD14NSNS_WT_0001to01AU_1kmsec

#input data folder:
data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'
#read data:
tf = open(data_folder+data_name+'MC_settings_list_INT.txt', "r")
MC_settings_list_INT        = np.loadtxt(tf, dtype=int)         #[nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 0, 0]
tf.close()
tf = open(data_folder+data_name+'MC_settings_list_REAL.txt', "r")
MC_settings_list_REAL       = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'outlist_MC_info_INT.txt', "r")
outlist_MC_info_INT         = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'outlist_MC_info_REAL.txt', "r")
outlist_MC_info_REAL        = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_endstate_INT.txt', "r")
output_Nbody_endstate_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_endstate_REAL.txt', "r")
output_Nbody_endstate_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_usrstate_INT.txt', "r")
output_Nbody_usrstate_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_usrstate_REAL.txt', "r")
output_Nbody_usrstate_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_info_INT.txt', "r")
output_Nbody_xtra_info_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_info_REAL.txt', "r")
output_Nbody_xtra_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
#----------------------------------------------------------
#----------------------------------------------------------
#Define:
#----------------------------------------------------------
nr_SMA                          = MC_settings_list_INT[0]
nr_vinf                         = MC_settings_list_INT[1]
nr_scatterings_per_paramcomb    = MC_settings_list_INT[2]
#----------------------------------------------------------
#reshape input arrays for easy analysis:
#----------------------------------------------------------
RS_avs_outlist_MC_info_INT          = outlist_MC_info_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)                #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
RS_avs_outlist_MC_info_REAL         = outlist_MC_info_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)               #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]
RS_avs_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)          #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
RS_avs_output_Nbody_endstate_REAL   = output_Nbody_endstate_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)         #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
RS_avs_output_Nbody_usrstate_INT    = output_Nbody_usrstate_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10, 10)      #[id, bin_i, bin_j, IMS_rp_counter, IMS_binsin_counter]
RS_avs_output_Nbody_usrstate_REAL   = output_Nbody_usrstate_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10, 10)     #DONT USE ANYMORE
RS_avs_output_Nbody_xtra_info_INT   = output_Nbody_xtra_info_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)
RS_avs_output_Nbody_xtra_info_REAL  = output_Nbody_xtra_info_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)        #[rmin_12, rmin_13, rmin_23]
#From data define:
MCvalues_SMA_arr                    = RS_avs_outlist_MC_info_REAL[:,0,0,0]
MCvalues_vinf_arr                   = RS_avs_outlist_MC_info_REAL[0,:,0,1]
#----------------------------------------------------------

#----------------------------------------------------------
#Allocate arrays:
#----------------------------------------------------------
#below: 'sim' is the direct count from the sim. This can be slightly biased since some interactions may not have been completed wihtin the time limit. 
nr_eachstate_sim_DI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
nr_eachstate_sim_RI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
#below: 'fin' is our (best try) corrected nr of endstates, where we try to correct for the set of non-completed interactions.
nr_eachstate_fin_DI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
nr_eachstate_fin_RI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
#probability and cross sections:
#DI:
prob_DI_arr                         = np.zeros((nr_SMA, nr_vinf), dtype=np.float64) 
prob_eachstate_DI                   = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
err_prob_eachstate_DI               = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
cross_section_allstates_DI          = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
err_cross_section_allstates_DI      = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64) 
#RI:
prob_RI_arr                         = np.zeros((nr_SMA, nr_vinf), dtype=np.float64) 
prob_eachstate_RI                   = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
err_prob_eachstate_RI               = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
cross_section_allstates_RI          = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
err_cross_section_allstates_RI      = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
#DI,RI:
prob_eachstate_DIRI                 = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
err_prob_eachstate_DIRI             = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
cross_section_allstates_DIRI        = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64)
err_cross_section_allstates_DIRI    = np.zeros((nr_SMA, nr_vinf, 3,10), dtype=np.float64) 
#make arr of object combinations: 
obj_combs_arr = np.array([[1,2],[1,3],[2,3]])
#recoil kick energy array:
Etot_ij_k_DI_RI_arr = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 2, 3), dtype=np.float64)     #[a, v, nr_scatterings_per_paramcomb, DI/RI,12/13/23] - save E in that 'pos'.
Etot_ij_k_DI_RI_arr[:,:,:,:,:] = -1 #initialize (use as a threshold for organizing later...)
Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  
Efin_Rmin_arr[:,:,:,:]  = -1    #initialize
#----------------------------------------------------------






#----------------------------------------------------------
#Calc cross sections:
#----------------------------------------------------------
#We calc the CS for each combination of SMA(ac),vinf(vc):
for ac in range(0,nr_SMA):                               # ac-loop over SMA
    for vc in range(0,nr_vinf):                          # vc-loop over vinf
        
        #---------------------------------------
        #set threshold between DI and RI:
        #---------------------------------------
        thresholdnr_RIDI    = 3 #This cut is used e.g. to divide the no endstates into different states assuming no endstate is a RI.
        #---------------------------------------
        #Endstate counters (reset):
        #---------------------------------------
        #Direct interactions (DI):
        nrctot_DI                       = 0 #counter
        DI_nrc_123_binsin               = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_ionization           = [0,0,0]   #we only use the first entry
        DI_nrc_123_headon_collisions    = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_insp_collisions      = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_userstate_20         = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_userstate_21         = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_userstate_22         = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_inspiral             = [0,0,0]   #format: [12,13,23]
        DI_nrc_123_tidalevent           = [0,0,0]   #format: [12,13,23]
        #resonance interactions (RI):
        nrctot_RI                       = 0 #counter
        RI_nrc_123_binsin               = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_ionization           = [0,0,0]   #we only use the first entry
        RI_nrc_123_headon_collisions    = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_insp_collisions      = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_userstate_20         = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_userstate_21         = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_userstate_22         = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_inspiral             = [0,0,0]   #format: [12,13,23]
        RI_nrc_123_tidalevent           = [0,0,0]   #format: [12,13,23]
        #no endstate:
        nrc_no_endstate                 = 0 #counting nr of interactions that did not finish in time     
        #---------------------------------------
        
        #---------------------------------------
        #loop over nr scatterings and count endstates:
        #---------------------------------------
        for sc in range(0,nr_scatterings_per_paramcomb):  # sc-loop over nr scatterings per param combination (a(ac),v(vc),...)          
            
            #ENDSTATE CODES: 1 = TDE, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 5 = Inspiral, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...
            #IF ENDSTATE code is < 10 then the code reached a physical OK endstate - if >=10 then the code did not converge towards a physical solution.
			#Return_3body_Info_REAL(1:10)		= [E_kin(ij), E_pot(ij), E_tot(ij), a_bin(ij), e_bin(ij), E_kin(ijk), E_pot(ijk), E_tot(ijk), a_bin(ijk), e_bin(ijk)]
            
            #-------------------------
            #Define:
            #-------------------------
            #RS_avs_output_Nbody_endstate_INT format: [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
            endstate_numseq             = RS_avs_output_Nbody_endstate_INT[ac,vc,sc,0:4] #numseq: [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k]
            endstate_id                 = endstate_numseq[0]
            #nr of pericenter passages from IMS to collision: (if collisions is an endstate)
            nr_peri_passages            = RS_avs_output_Nbody_endstate_INT[ac,vc,sc,5]
            #nr of IMS binary changes from ini to final endstate:
            nr_IMS_binsin               = RS_avs_output_Nbody_endstate_INT[ac,vc,sc,6]
            #userstate 10x10 array block:
            usrstate_INT                = RS_avs_output_Nbody_usrstate_INT[ac,vc,sc,:,:]
            usrstate_REAL               = RS_avs_output_Nbody_usrstate_REAL[ac,vc,sc,:,:]
            #orbital params:
            fin_Etot_ijk                = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,7]
            fin_a_ij                    = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,3]   
            fin_e_ij                    = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,4]
            ini_Etot                    = RS_avs_outlist_MC_info_REAL[ac,vc,sc,4]   #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]
            ini_Ltot                    = RS_avs_outlist_MC_info_REAL[ac,vc,sc,5]   #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]
            #Rmin:
            Rmin_12_13_23               = RS_avs_output_Nbody_xtra_info_REAL[ac,vc,sc,0:3]
            #-------------------------
            
            
            #-------------------------
            #make list of userstates:
            #-------------------------
            pos_userstates  = np.where(usrstate_INT[:,0] != 0)[0]    #need this [0] in the end.
            nr_userstates   = len(pos_userstates)
            if (nr_userstates > 0):
                arr_userstate_ids           = usrstate_INT[pos_userstates,0]
                arr_nrp_passages_IMStoState = usrstate_INT[pos_userstates,3]
                arr_userstate_numseq        = usrstate_INT[pos_userstates,0:3]  #[userstateID(20,21,22,...), bin_i, bin_j]
            #-------------------------
             
                                     
            #-------------------------
            #Direct Interaction (DI) - with ok endstate:
            #-------------------------
            if (nr_IMS_binsin <= thresholdnr_RIDI   and endstate_id < 10):
                #count nr DI:
                nrctot_DI = nrctot_DI+1
                #Endstate: Ionization
                if (endstate_numseq[0] == 4):   DI_nrc_123_ionization[0]    = DI_nrc_123_ionization[0]+1                    
                #Loop over obj comb:
                for noc in range(0,3):
                    #2body object combination:  #[[1,2],[1,3],[2,3]]
                    obj2list = list(obj_combs_arr[noc,:])
                    #-----------------
                    #Endstates: 
                    #-----------------
                    #endstate_numseq: [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k]
                    #Endstate: Binary_single endstate
                    if ((endstate_numseq[0:3] == ([3]+obj2list)).all()): DI_nrc_123_binsin[noc]     = DI_nrc_123_binsin[noc]+1
                    #Endstate: collisions
                    if ((endstate_numseq[0:3] == ([2]+obj2list)).all() and nr_peri_passages == 0): DI_nrc_123_headon_collisions[noc]    = DI_nrc_123_headon_collisions[noc]+1
                    if ((endstate_numseq[0:3] == ([2]+obj2list)).all() and nr_peri_passages >= 1): DI_nrc_123_insp_collisions[noc]      = DI_nrc_123_insp_collisions[noc]+1
                    #Endstate: Inspiral
                    if ((endstate_numseq[0:3] == ([5]+obj2list)).all()): DI_nrc_123_inspiral[noc]   = DI_nrc_123_inspiral[noc]+1
                    #Endstate: Tidal Event:
                    if ((endstate_numseq[0:3] == ([1]+obj2list)).all()): DI_nrc_123_tidalevent[noc] = DI_nrc_123_tidalevent[noc]+1
                    #-----------------
                    #-----------------
                    #Userstates: 
                    #-----------------
                    #userstate_numseq: [X, bin_i, bin_j, IMS_rp_counter, IMS_binsin_counter]
                    #2body object combination:
                    if (nr_userstates > 0):
                        #loop over nr userstates detected:
                        for uc in range(0,nr_userstates):
                            userstate_numseq    = arr_userstate_numseq[uc,:]   #[userstateID(20,21,22,...), bin_i, bin_j]
                            nrp_passages        = arr_nrp_passages_IMStoState[uc]
                            if ((userstate_numseq[:] == ([20]+obj2list)).all() and nrp_passages >= 1): DI_nrc_123_userstate_20[noc]  = DI_nrc_123_userstate_20[noc] + 1
                            if ((userstate_numseq[:] == ([21]+obj2list)).all() and nrp_passages >= 1): DI_nrc_123_userstate_21[noc]  = DI_nrc_123_userstate_21[noc] + 1
                            if ((userstate_numseq[:] == ([22]+obj2list)).all() and nrp_passages >= 1): DI_nrc_123_userstate_22[noc]  = DI_nrc_123_userstate_22[noc] + 1
                    #-----------------
                    #if endstate is a bin-sin:
                    #-----------------
                    if ((endstate_numseq[0:3] == ([3]+obj2list)).all()): Etot_ij_k_DI_RI_arr[ac,vc,sc,0,noc]    = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,7]
                    #-----------------
                    #any endstate:
                    #-----------------
                    if ((endstate_numseq[1:3] == (obj2list)).all()): Efin_Rmin_arr[ac,vc,sc,0:11] = [endstate_id, noc, 0, Rmin_12_13_23[0], Rmin_12_13_23[1], Rmin_12_13_23[2], fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_Etot, ini_Ltot]
                    #-----------------
            #-------------------------
            
            
            #-------------------------
            #Resonance Interaction (RI) - with ok endstate:
            #-------------------------
            if (nr_IMS_binsin > thresholdnr_RIDI    and endstate_id < 10):        
                #count nr RI:
                nrctot_RI = nrctot_RI+1
                #Endstate: Ionization
                if (endstate_numseq[0] == 4):   RI_nrc_123_ionization[0]    = RI_nrc_123_ionization[0]+1                    
                #Loop over obj comb:
                for noc in range(0,3):
                    #2body object combination:   [[1,2],[1,3],[2,3]]
                    obj2list = list(obj_combs_arr[noc,:])
                    #-----------------
                    #Endstates: 
                    #-----------------
                    #endstate_numseq: [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k]
                    #Endstate: Binary_single endstate
                    if ((endstate_numseq[0:3] == ([3]+obj2list)).all()): RI_nrc_123_binsin[noc]    = RI_nrc_123_binsin[noc]+1
                    #Endstate: collisions
                    if ((endstate_numseq[0:3] == ([2]+obj2list)).all() and nr_peri_passages == 0): RI_nrc_123_headon_collisions[noc]    = RI_nrc_123_headon_collisions[noc]+1
                    if ((endstate_numseq[0:3] == ([2]+obj2list)).all() and nr_peri_passages >= 1): RI_nrc_123_insp_collisions[noc]      = RI_nrc_123_insp_collisions[noc]+1
                    #Endstate: Inspiral
                    if ((endstate_numseq[0:3] == ([5]+obj2list)).all()): RI_nrc_123_inspiral[noc]   = RI_nrc_123_inspiral[noc]+1
                    if ((endstate_numseq[0:3] == ([5]+obj2list)).all()): print nr_peri_passages
                    #Endstate: Tidal Event:
                    if ((endstate_numseq[0:3] == ([1]+obj2list)).all()): RI_nrc_123_tidalevent[noc] = RI_nrc_123_tidalevent[noc]+1
                    #-----------------
                    #-----------------
                    #Userstates: 
                    #-----------------
                    #userstate_numseq: [X, bin_i, bin_j, IMS_rp_counter, IMS_binsin_counter]
                    #2body object combination:
                    if (nr_userstates > 0):
                        #loop over nr userstates detected:
                        for uc in range(0,nr_userstates):
                            userstate_numseq    = arr_userstate_numseq[uc,:]   #[userstateID(20,21,22,...), bin_i, bin_j]
                            nrp_passages        = arr_nrp_passages_IMStoState[uc]
                            if ((userstate_numseq[:] == ([20]+obj2list)).all() and nrp_passages >= 1): RI_nrc_123_userstate_20[noc]  = RI_nrc_123_userstate_20[noc] + 1
                            if ((userstate_numseq[:] == ([21]+obj2list)).all() and nrp_passages >= 1): RI_nrc_123_userstate_21[noc]  = RI_nrc_123_userstate_21[noc] + 1
                            if ((userstate_numseq[:] == ([22]+obj2list)).all() and nrp_passages >= 1): RI_nrc_123_userstate_22[noc]  = RI_nrc_123_userstate_22[noc] + 1
                    #-----------------
                    #if endstate is a bin-sin:
                    #-----------------
                    if ((endstate_numseq[0:3] == ([3]+obj2list)).all()): Etot_ij_k_DI_RI_arr[ac,vc,sc,1,noc]    = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,7]
                    #-----------------     
                    #-----------------
                    #any endstate:
                    #-----------------
                    if ((endstate_numseq[1:3] == (obj2list)).all()): Efin_Rmin_arr[ac,vc,sc,0:11] = [endstate_id, noc, 1, Rmin_12_13_23[0], Rmin_12_13_23[1], Rmin_12_13_23[2], fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_Etot, ini_Ltot]
                    #-----------------
            #-------------------------
            
            
            #-------------------------
            #No identified endstate:
            #-------------------------
            if (endstate_id >= 10): nrc_no_endstate = nrc_no_endstate+1     #counting nr of interactions that did not finish in time 
            #-------------------------
            
                      
        #---------------------------------------
        #---------------------------------------
        
                
        
        #---------------------------------------
        #Number and Probability of different endstates:
        #---------------------------------------
        #Change arr format of our sim state counts:        
        for noc in range(0,3):
            #2body object combination: [[1,2],[1,3],[2,3]]
            nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]
            nr_eachstate_sim_RI[noc,:]  = [RI_nrc_123_binsin[noc], RI_nrc_123_ionization[noc], RI_nrc_123_headon_collisions[noc], RI_nrc_123_insp_collisions[noc], RI_nrc_123_userstate_20[noc], RI_nrc_123_userstate_21[noc], RI_nrc_123_userstate_22[noc], RI_nrc_123_inspiral[noc], RI_nrc_123_tidalevent[noc], 0.0]        
        #ADJUST nr endstates ASSUMING the no-endstate-set is DOMINATED BY RESONANCE INTERACTIONS (RIs):
        if (nrctot_RI > 0):     #correct if we have RI sample
            nrctot_fin_RI = nrctot_RI + nrc_no_endstate
            nr_eachstate_fin_DI[:,:]  = nr_eachstate_sim_DI[:,:]
            nr_eachstate_fin_RI[:,:]  = nrctot_fin_RI*(1.*nr_eachstate_sim_RI[:,:]/(1.*nrctot_RI))  
        if (nrctot_RI == 0):    #no correction if nr RI is 0
            nr_eachstate_fin_DI[:,:]  = nr_eachstate_sim_DI[:,:]
            nr_eachstate_fin_RI[:,:]  = nr_eachstate_sim_RI[:,:]
            if (nrc_no_endstate > 0): print 'WARNING: no RI sample to recon no-id sample with.'  
        #corrected final total number of each endstate from both RI,DIs: 
        nr_eachstate_fin_RIDI   = nr_eachstate_fin_DI[:,:] + nr_eachstate_fin_RI[:,:]           
                
        #probabilities:
        #DI: (ONLY USING DI SAMPLE)
        if (nrctot_DI > 0): 
            prob_eachstate_DI[ac,vc,:,:]        = (1./(1.*nr_scatterings_per_paramcomb))*nr_eachstate_fin_DI[:,:]
            err_prob_eachstate_DI[ac,vc,:,:]    = (1./(1.*nr_scatterings_per_paramcomb))*np.sqrt(nr_eachstate_fin_DI[:,:])
            prob_DI_arr[ac,vc]                  = (1./(1.*nr_scatterings_per_paramcomb))*(1.*nrctot_DI)
        #RI: (ONLY USING RI SAMPLE)
        if (nrctot_RI > 0):
            prob_eachstate_RI[ac,vc,:,:]        = (1./(1.*nr_scatterings_per_paramcomb))*nr_eachstate_fin_RI[:,:]
            err_prob_eachstate_RI[ac,vc,:,:]    = (1./(1.*nr_scatterings_per_paramcomb))*np.sqrt(nr_eachstate_fin_RI[:,:])
            prob_RI_arr[ac,vc]                  = (1./(1.*nr_scatterings_per_paramcomb))*(1.*nrctot_fin_RI)
        #DI,RI:    
        prob_eachstate_DIRI[ac,vc,:,:]          = (1./(1.*nr_scatterings_per_paramcomb))*nr_eachstate_fin_RIDI[:,:]          #since we also incluse user-states in the array is the total prob not 1. Everything ok.
        err_prob_eachstate_DIRI[ac,vc,:,:]      = (1./(1.*nr_scatterings_per_paramcomb))*np.sqrt(nr_eachstate_fin_RIDI[:,:])        
        #---------------------------------------
        #---------------------------------------

        #---------------------------------------
        #cross sections:
        #---------------------------------------    
        bmax_ij = RS_avs_outlist_MC_info_REAL[ac,vc,0,2]
        CSfac   = np.pi*(bmax_ij**2.) 
        #DI:        
        cross_section_allstates_DI[ac,vc,:,:]         = CSfac*prob_eachstate_DI[ac,vc,:,:]
        err_cross_section_allstates_DI[ac,vc,:,:]     = CSfac*err_prob_eachstate_DI[ac,vc,:,:]        
        #RI:
        cross_section_allstates_RI[ac,vc,:,:]         = CSfac*prob_eachstate_RI[ac,vc,:,:]
        err_cross_section_allstates_RI[ac,vc,:,:]     = CSfac*err_prob_eachstate_RI[ac,vc,:,:]
        #DI,RI:       
        cross_section_allstates_DIRI[ac,vc,:,:]       = CSfac*prob_eachstate_DIRI[ac,vc,:,:]
        err_cross_section_allstates_DIRI[ac,vc,:,:]   = CSfac*err_prob_eachstate_DIRI[ac,vc,:,:]
        #---------------------------------------
        #---------------------------------------
        
        #---------------------------------------
        #print test info:
        #---------------------------------------
        print 'TEST INFO:'
        print 'CSfac', bmax_ij, CSfac 
        print ac,vc
        print nr_eachstate_sim_DI
        print nr_eachstate_sim_RI
        #print prob_eachstate_RI[ac,vc,:,:]
        #print nr_eachstate_sim_DI
        print 'DI FIN:'
        print nr_eachstate_fin_DI
        print 'RI FIN:'
        print nr_eachstate_fin_RI
        print 'RIDI SIM:'
        print nr_eachstate_sim_DI+nr_eachstate_sim_RI
        print 'RIDI FIN:'
        print nr_eachstate_fin_RIDI
        #print prob_eachstate_DIRI[ac,vc,:,:]
        print nrctot_DI, nrctot_RI, nrc_no_endstate#, (1.*nrc_no_endstate)/(1.*nrctot_RI), prob_DI_arr[ac,vc], prob_RI_arr[ac,vc] 
        print cross_section_allstates_DIRI[ac,vc,:,:] 
        print 'fraction:', (1.*nrc_no_endstate)/(1.*nr_scatterings_per_paramcomb)
        #---------------------------------------
        
#----------------------------------------------------------        
#----------------------------------------------------------







#RS_avs_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)       #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
#print (sum(RS_avs_output_Nbody_endstate_INT[0,0,:,6])/50000.)/2.
#print (sum(RS_avs_output_Nbody_endstate_INT[1,0,:,6])/50000.)/2.
#print (sum(RS_avs_output_Nbody_endstate_INT[2,0,:,6])/50000.)/2.
#exit()





#--------------------------------------------------------------------------------------------------------------------
#Analysis:
#--------------------------------------------------------------------------------------------------------------------













#----------------------------------------------------------
#Analysis:
#----------------------------------------------------------
fig = plt.figure(figsize=(12, 8))

MTO     = 1.4       #m1
RTO     = 1e-5
mNS     = 1.4       #m2,m3
mNSbin  = mNS+mNS
mINbin  = MTO+mNS


#MTO     = 10.       #m1
#RTO     = 1e-5
#mNS     = 20.       #m2,m3
#mNSbin  = mNS+mNS
#mINbin  = MTO+mNS

#MTO     = 0.5
#RTO     = 0.013*((1.43/MTO)**(1./3.))*((1.-MTO/1.43)**(0.447))  #1e-5
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = MTO+mNS
#print 'RTO: ', RTO

#MTO     = 1.0
#RTO     = 1.0**(0.8)
#mNS     = 1.0
#mNSbin  = mNS+mNS
#mINbin  = MTO+mNS
#print 'RTO: ', RTO

#MTO     = 0.4
#RTO     = MTO**(0.8)
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = MTO+mNS
#print 'RTO: ', RTO

#RTO     = 0.4**(0.8)
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.4

#RTO     = 0.018
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.2

#RTO     = 0.03
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.05

#RTO     = 0.8**(0.8)
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.8

#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  

#in this plot we have 'v' as a constant and vary 'a':
vi = 0
#calc frac with t<tlife:
tlife_limit     = 10.**10.          #in years
rdist_limit_arr = np.array([2.0, 3.0, 5.0])   #in units R_TO - ONLY PUT 3 VALUES!!!

frac_tlife_LT_tlife_limit_arr   = np.zeros(nr_SMA, dtype=np.float64)
frac_tlife_rmin_cut_arr         = np.zeros((3,nr_SMA), dtype=np.float64)
frac_rmin_cut_arr               = np.zeros((3,nr_SMA), dtype=np.float64)
info_arr_NSNS_bin               = np.zeros((nr_SMA,nr_scatterings_per_paramcomb,10), dtype=np.float64)
info_arr_NSNS_bin[:,:,:]        = -1.
#----------------------------------
for ac in range(0,nr_SMA):
#----------------------------------    
    #------------------------
    #EXCHANGE NS-NS ENDSTATE:
    #------------------------
    #find pos for exch to NS-NS (2,3)
    pos_id_3        = np.where(Efin_Rmin_arr[ac,vi,:,0] == 3)[0]
    pos_bin_23      = np.where(Efin_Rmin_arr[ac,vi,:,1] == 2)[0]
    pos_NSNS        = list(set(pos_id_3).intersection(pos_bin_23))
    nr_NSNS         = len(pos_NSNS)
    if (nr_NSNS > 0):
        #get orbital a,e for this set:
        a_Rsun_arr      = Efin_Rmin_arr[ac,vi,pos_NSNS,7]
        a_AU_arr        = a_Rsun_arr/AU_U
        e_arr           = Efin_Rmin_arr[ac,vi,pos_NSNS,8]
        #calc corresponding t_life times (see tlife_calc.py):
        tfac        = (a_Rsun_arr**(4.))/(mNS*mNS*(mNS+mNS))
        tbin_yr_arr = np.zeros(nr_NSNS, dtype=np.float64)
        for ec in range(0,nr_NSNS):
            #read ecc:
            e_bin = e_arr[ec]
            #calc tlife from interpolation:
            if (e_bin <  max(tdata_e_arr)):
                tbin_yr = tfac[ec]*ecc_to_tdata_interpol_func(e_bin)
            if (e_bin >= max(tdata_e_arr)):
                tbin_yr_circular = tfac[ec]*ecc_to_tdata_interpol_func(0.0)
                tbin_yr = (768./425.)*tbin_yr_circular*((1.-e_bin**2.)**(7./2.))
            #save tlife in array:    
            tbin_yr_arr[ec] = tbin_yr    
        #apply tlife threshold:
        pos_tlife       = np.where(tbin_yr_arr[:] < tlife_limit)[0]
        frac_tlife_LT_tlife_limit = (1.*len(pos_tlife)/(1.*nr_NSNS))
        frac_tlife_LT_tlife_limit_arr[ac] = frac_tlife_LT_tlife_limit
        #calc and save info about (all) the NSNS bins:
        rmin_12     = Efin_Rmin_arr[ac,vi,pos_NSNS,3] 
        rmin_13     = Efin_Rmin_arr[ac,vi,pos_NSNS,4]
        rmin_arr    = np.array([min([rmin_12[i],rmin_13[i]]) for i in range(0,nr_NSNS)])    #min(rmin) between 1,2 or 1,3 (TO-CO)
        Pdays_arr   = Pdays_in_aAU_mbin(a_AU_arr, mNSbin)
        EfinNSNSbs  = Efin_Rmin_arr[ac,vi,pos_NSNS,6]
        info_arr_NSNS_bin[ac,0:nr_NSNS,0] = a_AU_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,1] = e_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,2] = tbin_yr_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,3] = rmin_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,4] = Pdays_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,5] = EfinNSNSbs
        #apply further cuts (rmin, ...):
        for rc in range(0,3):
            pos_rmin        = np.where(rmin_arr[:] > RTO*rdist_limit_arr[rc])[0]
            pos_tr_cut      = list(set(pos_tlife).intersection(pos_rmin))
            frac_tr_cut     = (1.*len(pos_tr_cut)/(1.*nr_NSNS))
            frac_tlife_rmin_cut_arr[rc,ac] = frac_tr_cut
            print ac, frac_tlife_LT_tlife_limit, frac_tr_cut, (1.*len(pos_rmin)/(1.*nr_NSNS)), RTO*rdist_limit_arr[rc]
    #------------------------
    #GW INSPIRAL ENDSTATE:
    #------------------------
    #Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  
    #find pos for GW inspiral NS-NS (2,3)
    pos_id_insp     = np.where(Efin_Rmin_arr[ac,vi,:,0] == 5)[0]
    pos_bin_23      = np.where(Efin_Rmin_arr[ac,vi,:,1] == 2)[0]
    pos_NSNS        = list(set(pos_id_insp).intersection(pos_bin_23))
    nr_NSNS         = len(pos_NSNS) # = nr NS-NS inspirals
    if (nr_NSNS > 0):
        #find rmin (to TO) for each NS-NS insp:
        rmin_12     = Efin_Rmin_arr[ac,vi,pos_NSNS,3] 
        rmin_13     = Efin_Rmin_arr[ac,vi,pos_NSNS,4]
        rmin_arr    = np.array([min([rmin_12[i],rmin_13[i]]) for i in range(0,nr_NSNS)])
        #apply rmin cut on all NS-NS inspirals:       
        for rc in range(0,3):   #loop over 3 different rmin cuts
            pos_rmin                    = np.where(rmin_arr[:] > RTO*rdist_limit_arr[rc])[0]
            frac_rmin_cut               = (1.*len(pos_rmin)/(1.*nr_NSNS))
            frac_rmin_cut_arr[rc,ac]    = frac_rmin_cut
            print 'INSPIRALS INFO:  ', 'frac_rmin_cut = ', frac_rmin_cut
    #------------------------    
#----------------------------------    
    
    
#scalings:
ap_AU           = MCvalues_SMA_arr[:]/AU_U
cs_AU           = (R_sun_SI/AU_SI)**2.
cs_arr_AU       = cs_AU*cross_section_allstates_DIRI[:,:,:,:]   
cs_err_arr_AU   = cs_AU*err_cross_section_allstates_DIRI[:,:,:,:]
Pdays_arr       = Pdays_in_aAU_mbin(ap_AU, mINbin)

#Save data:
nr_a        = len(ap_AU)
nr_cs       = 20    #dont change. We just always keep space for ...
save_arr_cs     = np.zeros((nr_cs,nr_a), dtype=np.float64)
save_arr_cs_err = np.zeros((nr_cs,nr_a), dtype=np.float64)
#decide which cs to save:
#2body object combination: [[1,2],[1,3],[2,3]]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]
#format: cs_arr_AU[a,v,obj comb, cs]
#cross sections:
save_arr_cs[0,:]        = ap_AU
save_arr_cs[1,:]        = cs_arr_AU[:,vi,2,0]                                       #NS-NS exchange
save_arr_cs[2,:]        = cs_arr_AU[:,vi,0,2]                                       #collisions 12
save_arr_cs[3,:]        = cs_arr_AU[:,vi,1,2]                                       #collisions 13
save_arr_cs[4,:]        = cs_arr_AU[:,vi,2,2]                                       #collisions 23
save_arr_cs[5,:]        = cs_arr_AU[:,vi,0,3]                                       #inspiral collision 12
save_arr_cs[6,:]        = cs_arr_AU[:,vi,1,3]                                       #inspiral collision 13
save_arr_cs[7,:]        = cs_arr_AU[:,vi,2,3]                                       #inspiral collision 23
save_arr_cs[8,:]        = cs_arr_AU[:,vi,0,7]                                       #inspiral 12
save_arr_cs[9,:]        = cs_arr_AU[:,vi,1,7]                                       #inspiral 13
save_arr_cs[10,:]       = cs_arr_AU[:,vi,2,7]                                       #inspiral 23
save_arr_cs[11,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_LT_tlife_limit_arr[:]              #NS-NS exchange with t_life<10**10 years
save_arr_cs[12,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[0,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs[13,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[1,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs[14,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[2,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs[15,:]       = cs_arr_AU[:,vi,0,8]                                       #TDE 12
save_arr_cs[16,:]       = cs_arr_AU[:,vi,1,8]                                       #TDE 13
save_arr_cs[17,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[0,:]                #inspiral 23 with rmin cut
save_arr_cs[18,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[1,:]                #inspiral 23 with rmin cut
save_arr_cs[19,:]       = cs_arr_AU[:,vi,2,7]*frac_rmin_cut_arr[2,:]                #inspiral 23 with rmin cut

#cross sections err:
save_arr_cs_err[0,:]    = ap_AU
save_arr_cs_err[1,:]    = cs_err_arr_AU[:,vi,2,0]                                   #NS-NS exchange
save_arr_cs_err[2,:]    = cs_err_arr_AU[:,vi,0,2]                                   #collisions 12
save_arr_cs_err[3,:]    = cs_err_arr_AU[:,vi,1,2]                                   #collisions 13
save_arr_cs_err[4,:]    = cs_err_arr_AU[:,vi,2,2]                                   #collisions 23
save_arr_cs_err[5,:]    = cs_err_arr_AU[:,vi,0,3]                                   #inspiral collision 12
save_arr_cs_err[6,:]    = cs_err_arr_AU[:,vi,1,3]                                   #inspiral collision 13
save_arr_cs_err[7,:]    = cs_err_arr_AU[:,vi,2,3]                                   #inspiral collision 23
save_arr_cs_err[8,:]    = cs_err_arr_AU[:,vi,0,7]                                   #inspiral 12
save_arr_cs_err[9,:]    = cs_err_arr_AU[:,vi,1,7]                                   #inspiral 13
save_arr_cs_err[10,:]   = cs_err_arr_AU[:,vi,2,7]                                   #inspiral 23
save_arr_cs_err[11,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_LT_tlife_limit_arr[:]) #NS-NS exchange with t_life<10**10 years
save_arr_cs_err[12,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_rmin_cut_arr[0,:])     #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs_err[13,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_rmin_cut_arr[1,:])     #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs_err[14,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_rmin_cut_arr[2,:])     #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs_err[15,:]   = cs_err_arr_AU[:,vi,0,8]                                   #TDE 12
save_arr_cs_err[16,:]   = cs_err_arr_AU[:,vi,1,8]                                   #TDE 13
save_arr_cs_err[17,:]   = cs_err_arr_AU[:,vi,2,7]*np.sqrt(frac_rmin_cut_arr[0,:])   #inspiral 23 with rmin cut
save_arr_cs_err[18,:]   = cs_err_arr_AU[:,vi,2,7]*np.sqrt(frac_rmin_cut_arr[1,:])   #inspiral 23 with rmin cut
save_arr_cs_err[19,:]   = cs_err_arr_AU[:,vi,2,7]*np.sqrt(frac_rmin_cut_arr[2,:])   #inspiral 23 with rmin cut

#save cs:
tf = open('cs_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs,   fmt='%5f')
tf.close()
#save cs err:
tf = open('cs_err_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs_err,   fmt='%5f')
tf.close()







print 'ORBIT INFO:'
print ap_AU, Pdays_arr, MCvalues_SMA_arr[:]/RTO
##2body object combination: [[1,2],[1,3],[2,3]]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]

#plot cs:
csp = cs_arr_AU[:,vi,2,0]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)
csp = cs_arr_AU[:,vi,2,0]*frac_tlife_LT_tlife_limit_arr[:]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=10.0, linewidth=1.5, alpha=0.75)

csp = cs_arr_AU[:,vi,0,2] + cs_arr_AU[:,vi,1,2]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

csp = cs_arr_AU[:,vi,0,8] + cs_arr_AU[:,vi,1,8]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

fig.add_subplot(221).set_ylim(0.0,5.0)
plt.xscale('log')


#plot distributions:
#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  


#INPUT SAVE INDEX:
ai = 0

x = info_arr_NSNS_bin[ai,:,3]/RTO
y = np.log10(info_arr_NSNS_bin[ai,:,2])
fig.add_subplot(222).plot(x,y, marker='o', linestyle='none', markersize=2.0, linewidth=1.5, alpha=0.75)

x = info_arr_NSNS_bin[ai,:,5]
y = info_arr_NSNS_bin[ai,:,2]
fig.add_subplot(223).plot(x,y, marker='o', linestyle='none', markersize=2.0, linewidth=1.5, alpha=0.75)
plt.yscale('log')

x = info_arr_NSNS_bin[ai,:,3]/RTO
y = info_arr_NSNS_bin[ai,:,5]
fig.add_subplot(224).plot(x,y, marker='o', linestyle='none', markersize=2.0, linewidth=1.5, alpha=0.75)

plt.show()
#exit()

print 'cs info:'
print ap_AU
print cs_arr_AU[:,0,2,7]/(10.**2.)

#save NS-NS info:
#info_arr_NSNS_bin               = np.zeros((nr_SMA,nr_scatterings_per_paramcomb,10), dtype=np.float64)
#info_arr_NSNS_bin[:,:,:]        = -1.
#info_arr_NSNS_bin[ac,0:nr_NSNS,0] = a_AU_arr
#info_arr_NSNS_bin[ac,0:nr_NSNS,1] = e_arr
#info_arr_NSNS_bin[ac,0:nr_NSNS,2] = tbin_yr_arr
#info_arr_NSNS_bin[ac,0:nr_NSNS,3] = rmin_arr
#info_arr_NSNS_bin[ac,0:nr_NSNS,4] = Pdays_arr
#info_arr_NSNS_bin[ac,0:nr_NSNS,5] = EfinNSNSbs
tf = open('NSNS_info_' + data_name, "w")
np.savetxt(tf, info_arr_NSNS_bin[ai,:,:],   fmt='%5f')
tf.close()


exit()
#----------------------------------------------------------
















#----------------------------------------------------------
#Grindlay analysis:
#----------------------------------------------------------
fig = plt.figure(figsize=(12, 8))

MTO     = 0.2
RTO     = 0.0041176640*(MTO/1.0)**(-1./3.)
mNS     = 1.4
mNSbin  = mNS+mNS
mINbin  = MTO+mNS

#RTO     = 0.4**(0.8)
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.4

#RTO     = 0.018
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.2

#RTO     = 0.03
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.05

#RTO     = 0.8**(0.8)
#mNS     = 1.4
#mNS     = 1.4
#mNSbin  = mNS+mNS
#mINbin  = 1.4+0.8

#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  

#in this plot we have 'v' as a constant and vary 'a':
vi = 0
#calc frac with t<tlife:
tlife_limit     = 10.**10.          #in years
rdist_limit_arr = np.array([2.0, 3.5, 5.0])   #in units R_TO - ONLY PUT 3 VALUES!!!
frac_tlife_LT_tlife_limit_arr   = np.zeros(nr_SMA, dtype=np.float64)
frac_tlife_rmin_cut_arr         = np.zeros((3,nr_SMA), dtype=np.float64)
info_arr_NSNS_bin               = np.zeros((nr_SMA,nr_scatterings_per_paramcomb,10), dtype=np.float64)
info_arr_NSNS_bin[:,:,:]        = -1.
for ac in range(0,nr_SMA):
    #find pos for exch to NS-NS (2,3)
    pos_id_3        = np.where(Efin_Rmin_arr[ac,vi,:,0] == 3)[0]
    pos_bin_23      = np.where(Efin_Rmin_arr[ac,vi,:,1] == 2)[0]
    pos_NSNS        = list(set(pos_id_3).intersection(pos_bin_23))
    nr_NSNS         = len(pos_NSNS)
    if (nr_NSNS > 0):
        #get orbital a,e for this set:
        a_Rsun_arr      = Efin_Rmin_arr[ac,vi,pos_NSNS,7]
        a_AU_arr        = a_Rsun_arr/AU_U
        e_arr           = Efin_Rmin_arr[ac,vi,pos_NSNS,8]
        #calc corresponding t_life times (see tlife_calc.py):
        tfac        = (a_Rsun_arr**(4.))/(mNS*mNS*(mNS+mNS))
        tbin_yr_arr = np.zeros(nr_NSNS, dtype=np.float64)
        for ec in range(0,nr_NSNS):
            #read ecc:
            e_bin = e_arr[ec]
            #calc tlife from interpolation:
            if (e_bin <  max(tdata_e_arr)):
                tbin_yr = tfac[ec]*ecc_to_tdata_interpol_func(e_bin)
            if (e_bin >= max(tdata_e_arr)):
                tbin_yr_circular = tfac[ec]*ecc_to_tdata_interpol_func(0.0)
                tbin_yr = (768./425.)*tbin_yr_circular*((1.-e_bin**2.)**(7./2.))
            #save tlife in array:    
            tbin_yr_arr[ec] = tbin_yr    
        #apply tlife threshold:
        pos_tlife       = np.where(tbin_yr_arr[:] < tlife_limit)[0]
        frac_tlife_LT_tlife_limit = (1.*len(pos_tlife)/(1.*nr_NSNS))
        frac_tlife_LT_tlife_limit_arr[ac] = frac_tlife_LT_tlife_limit
        #calc and save info about (all) the NSNS bins:
        rmin_12     = Efin_Rmin_arr[ac,vi,pos_NSNS,3] 
        rmin_13     = Efin_Rmin_arr[ac,vi,pos_NSNS,4]
        rmin_arr    = np.array([min([rmin_12[i],rmin_13[i]]) for i in range(0,nr_NSNS)])
        Pdays_arr   = Pdays_in_aAU_mbin(a_AU_arr, mNSbin)
        EfinNSNSbs  = Efin_Rmin_arr[ac,vi,pos_NSNS,6]
        info_arr_NSNS_bin[ac,0:nr_NSNS,0] = a_AU_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,1] = e_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,2] = tbin_yr_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,3] = rmin_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,4] = Pdays_arr
        info_arr_NSNS_bin[ac,0:nr_NSNS,5] = EfinNSNSbs
        #apply further cuts (rmin, ...):
        for rc in range(0,3):
            pos_rmin        = np.where(rmin_arr[:] > RTO*rdist_limit_arr[rc])[0]
            pos_tr_cut      = list(set(pos_tlife).intersection(pos_rmin))
            frac_tr_cut     = (1.*len(pos_tr_cut)/(1.*nr_NSNS))
            frac_tlife_rmin_cut_arr[rc,ac] = frac_tr_cut
            print ac, frac_tlife_LT_tlife_limit, frac_tr_cut, (1.*len(pos_rmin)/(1.*nr_NSNS)), RTO*rdist_limit_arr[rc]
        
#scalings:
ap_AU           = MCvalues_SMA_arr[:]/AU_U
cs_AU           = (R_sun_SI/AU_SI)**2.
cs_arr_AU       = cs_AU*cross_section_allstates_DIRI[:,:,:,:]   
cs_err_arr_AU   = cs_AU*err_cross_section_allstates_DIRI[:,:,:,:]
Pdays_arr       = Pdays_in_aAU_mbin(ap_AU, mINbin)

print 'ORBIT INFO:'
print ap_AU, Pdays_arr, MCvalues_SMA_arr[:]/RTO
##2body object combination: [[1,2],[1,3],[2,3]]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]

#plot cs:
csp = cs_arr_AU[:,vi,2,0]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)
csp = cs_arr_AU[:,vi,2,0]*frac_tlife_LT_tlife_limit_arr[:]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=10.0, linewidth=1.5, alpha=0.75)

csp = cs_arr_AU[:,vi,0,2] + cs_arr_AU[:,vi,1,2]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

csp = cs_arr_AU[:,vi,0,8] + cs_arr_AU[:,vi,1,8]
fig.add_subplot(221).plot(Pdays_arr,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

fig.add_subplot(221).set_ylim(0.0,5.0)
plt.xscale('log')


#plot distributions:
#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  

ai = 3

x = info_arr_NSNS_bin[ai,:,3]/RTO
y = np.log10(info_arr_NSNS_bin[ai,:,2])
fig.add_subplot(222).plot(x,y, marker='o', linestyle='none', markersize=2.0, linewidth=1.5, alpha=0.75)

x = info_arr_NSNS_bin[ai,:,5]
y = info_arr_NSNS_bin[ai,:,2]
fig.add_subplot(223).plot(x,y, marker='o', linestyle='none', markersize=2.0, linewidth=1.5, alpha=0.75)
plt.yscale('log')

x = info_arr_NSNS_bin[ai,:,3]/RTO
y = info_arr_NSNS_bin[ai,:,5]
fig.add_subplot(224).plot(x,y, marker='o', linestyle='none', markersize=2.0, linewidth=1.5, alpha=0.75)

plt.show()
#exit()


#Save data:
nr_a        = len(ap_AU)
nr_cs       = 20    #dont change. We just always keep space for ...
save_arr_cs     = np.zeros((nr_cs,nr_a), dtype=np.float64)
save_arr_cs_err = np.zeros((nr_cs,nr_a), dtype=np.float64)
#decide which cs to save:
#2body object combination: [[1,2],[1,3],[2,3]]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]
#format: cs_arr_AU[a,v,obj comb, cs]
#cross sections:
save_arr_cs[0,:]        = ap_AU
save_arr_cs[1,:]        = cs_arr_AU[:,vi,2,0]                                       #NS-NS exchange
save_arr_cs[2,:]        = cs_arr_AU[:,vi,0,2]                                       #collisions 12
save_arr_cs[3,:]        = cs_arr_AU[:,vi,1,2]                                       #collisions 13
save_arr_cs[4,:]        = cs_arr_AU[:,vi,2,2]                                       #collisions 23
save_arr_cs[5,:]        = cs_arr_AU[:,vi,0,3]                                       #inspiral collision 12
save_arr_cs[6,:]        = cs_arr_AU[:,vi,1,3]                                       #inspiral collision 13
save_arr_cs[7,:]        = cs_arr_AU[:,vi,2,3]                                       #inspiral collision 23
save_arr_cs[8,:]        = cs_arr_AU[:,vi,0,7]                                       #inspiral 12
save_arr_cs[9,:]        = cs_arr_AU[:,vi,1,7]                                       #inspiral 13
save_arr_cs[10,:]       = cs_arr_AU[:,vi,2,7]                                       #inspiral 23
save_arr_cs[11,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_LT_tlife_limit_arr[:]              #NS-NS exchange with t_life<10**10 years
save_arr_cs[12,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[0,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs[13,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[1,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs[14,:]       = cs_arr_AU[:,vi,2,0]*frac_tlife_rmin_cut_arr[2,:]                  #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs[15,:]       = cs_arr_AU[:,vi,0,8]                                       #TDE 12
save_arr_cs[16,:]       = cs_arr_AU[:,vi,1,8]                                       #TDE 13
#cross sections err:
save_arr_cs_err[0,:]    = ap_AU
save_arr_cs_err[1,:]    = cs_err_arr_AU[:,vi,2,0]                                   #NS-NS exchange
save_arr_cs_err[2,:]    = cs_err_arr_AU[:,vi,0,2]                                   #collisions 12
save_arr_cs_err[3,:]    = cs_err_arr_AU[:,vi,1,2]                                   #collisions 13
save_arr_cs_err[4,:]    = cs_err_arr_AU[:,vi,2,2]                                   #collisions 23
save_arr_cs_err[5,:]    = cs_err_arr_AU[:,vi,0,3]                                   #inspiral collision 12
save_arr_cs_err[6,:]    = cs_err_arr_AU[:,vi,1,3]                                   #inspiral collision 13
save_arr_cs_err[7,:]    = cs_err_arr_AU[:,vi,2,3]                                   #inspiral collision 23
save_arr_cs_err[8,:]    = cs_err_arr_AU[:,vi,0,7]                                   #inspiral 12
save_arr_cs_err[9,:]    = cs_err_arr_AU[:,vi,1,7]                                   #inspiral 13
save_arr_cs_err[10,:]   = cs_err_arr_AU[:,vi,2,7]                                   #inspiral 23
save_arr_cs_err[11,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_LT_tlife_limit_arr[:]) #NS-NS exchange with t_life<10**10 years
save_arr_cs_err[12,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_rmin_cut_arr[0,:])     #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs_err[13,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_rmin_cut_arr[1,:])     #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs_err[14,:]   = cs_err_arr_AU[:,vi,2,0]*np.sqrt(frac_tlife_rmin_cut_arr[2,:])     #NS-NS exchange with t_life<10**10 years and r_min > r_cut
save_arr_cs_err[15,:]   = cs_err_arr_AU[:,vi,0,8]                                   #TDE 12
save_arr_cs_err[16,:]   = cs_err_arr_AU[:,vi,1,8]                                   #TDE 13

#save cs:
tf = open('cs_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs,   fmt='%5f')
tf.close()
#save cs err:
tf = open('cs_err_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs_err,   fmt='%5f')
tf.close()


exit()
#----------------------------------------------------------






#----------------------------------------------------------
#BH mergers:
#----------------------------------------------------------
fig = plt.figure(figsize=(10, 10))

#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [(0) endstate_id, (1) 12(0)/13(1)/23(2), (2) DI(0)/RI(1), (3) Rmin_12, (4) Rmin_13, (5) Rmin_23,
#(6) fin_Etot_ijk, (7) fin_a_ij, (8) fin_e_ij, (9) ini_E, (10) ini_L]]  

#initial is: a=1R_sun, m=1.

#IT IS IMPORTANT THAT WE USE CORRECT NORMALIZATIONS FOR IC!!!!
#INCLUDER EVT MCvalues_SMA_arr[0]!!!! FOR NORMALIZATION!!!!!!

#select [1,2](BH,BH) bin-sin outcomes:
ai = 0
vi = 0
pos_id_3        = np.where(Efin_Rmin_arr[ai,vi,:,0] == 3)[0]
pos_bin_12      = np.where(Efin_Rmin_arr[ai,vi,:,1] == 0)[0]
pos_fin         = list(set(pos_id_3).intersection(pos_bin_12))
nr_BHBH         = len(pos_fin) 
#define: (use only pos_fin)
Efin_arr    = Efin_Rmin_arr[ai,vi,pos_fin,6]
a_Rsun_arr  = Efin_Rmin_arr[ai,vi,pos_fin,7]
e_arr       = Efin_Rmin_arr[ai,vi,pos_fin,8]

mBH         = 1.

#calc corresponding t_life times (see tlife_calc.py):
tfac        = (a_Rsun_arr**(4.))/(mBH*mBH*(mBH+mBH))
tbin_yr_arr = np.zeros(nr_BHBH, dtype=np.float64)
for ec in range(0,nr_BHBH):
    #read ecc:
    e_bin = e_arr[ec]
    #calc tlife from interpolation:
    if (e_bin <  max(tdata_e_arr)):
        tbin_yr = tfac[ec]*ecc_to_tdata_interpol_func(e_bin)
    if (e_bin >= max(tdata_e_arr)):
        tbin_yr_circular = tfac[ec]*ecc_to_tdata_interpol_func(0.0)
        tbin_yr = (768./425.)*tbin_yr_circular*((1.-e_bin**2.)**(7./2.))
    #save tlife in array:    
    tbin_yr_arr[ec] = tbin_yr    


fig.add_subplot(221).plot(1.-e_arr[:], tbin_yr_arr[:], marker='o', linestyle='none', markersize=5.0, linewidth=1.5, alpha=0.75)
plt.xscale('log')
plt.yscale('log')

fig.add_subplot(222).hist(np.log10(tbin_yr_arr), bins=100, normed=True, histtype='step', cumulative=1)
plt.yscale('log')
fig.add_subplot(222).set_ylim(1e-10,1e10)

fig.add_subplot(223).hist(np.log10(tbin_yr_arr), bins=100, normed=True, histtype='step')
plt.yscale('log')
fig.add_subplot(223).set_ylim(1e-10,1e10)

#CHANGE REST OF SCRIPT! TLIFE ETC!!!!!
plt.show()
exit()


#initial is: a=1, m=1.
cs_11   = cross_section_allstates_DI[ai,vi,0,0]#         = CSfac*prob_eachstate_DI[ac,vc,:,:]
#vary m:
nr_m    = 10
m_arr   = np.linspace(1.0, 50, num=nr_m)
R_arr   = np.zeros(nr_m, dtype=np.float64)
for mc in range(0,nr_m):
    m       = m_arr[mc]
    a_HB    = 1.0*m*(36.5/50.)**2.
    cs_m    = a_HB*(m**(-1./3.))#cs_11
    d_a_au  = a_HB*Efin_Rmin_arr[ai,vi,pos_fin,7]/AU_U  
    d_e     = Efin_Rmin_arr[ai,vi,pos_fin,8]

    t_life  = (1.6*(10.**(17.)))*(d_a_au**4.)*(m**(-3.))*(768./425.)*(1.-d_e**2)**(7./2.)   #in years
    
    fig.add_subplot(122).hist(np.log10(t_life), bins=100, normed=True, histtype='step', cumulative=1)
    
    pos_t   = np.where(t_life < 1e4)[0]
    print len(pos_t)
    R_arr[mc]   = cs_m*len(pos_t)
    

fig.add_subplot(122).step(m_arr,R_arr)


#fig.add_subplot(122).hist(np.log10(t_life), bins=100, normed=True, histtype='step')
plt.yscale('log')
fig.add_subplot(122).set_ylim(1e-10,1e10)

plt.show()
exit()
#----------------------------------------------------------





#----------------------------------------------------------
#HILLS:
#----------------------------------------------------------
fig = plt.figure(figsize=(10, 10))

#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  

RTO = 1.0

ai = 0
vi = 0
pos_id_3        = np.where(Efin_Rmin_arr[ai,vi,:,0] == 3)[0]
pos_bin_13      = np.where(Efin_Rmin_arr[ai,vi,:,1] == 1)[0]
pos_bin_23      = np.where(Efin_Rmin_arr[ai,vi,:,1] == 2)[0]
pos_bin_1323    = np.hstack((pos_bin_13,pos_bin_23))
pos_fin         = list(set(pos_id_3).intersection(pos_bin_1323))

print len(pos_bin_13)
print len(pos_bin_23)
print Efin_Rmin_arr[ai,vi,pos_bin_13,0]

x = Efin_Rmin_arr[ai,vi,pos_fin,3]/RTO
y = Efin_Rmin_arr[ai,vi,pos_fin,6]
fig.add_subplot(111).plot(x,y, marker="o", linestyle='none', markersize=2.0, alpha=0.5)

plt.show()
exit()
#----------------------------------------------------------




#stellar bin-sin interactions:

fig = plt.figure(figsize=(10, 10))

#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  

RTO = 0.0041176640*(1.4/1.0)**(-1./3.)

ai = 0
vi = 0
pos_id_3    = np.where(Efin_Rmin_arr[ai,vi,:,0] == 3)[0]
pos_bin_23  = np.where(Efin_Rmin_arr[ai,vi,:,1] == 2)[0]
pos_fin     = list(set(pos_id_3).intersection(pos_bin_23))

x = Efin_Rmin_arr[ai,vi,pos_fin,3]/RTO
y = Efin_Rmin_arr[ai,vi,pos_fin,6]
fig.add_subplot(221).plot(x,y, marker="o", linestyle='none', markersize=2.0, alpha=0.5)
x = Efin_Rmin_arr[ai,vi,pos_fin,4]/RTO
y = Efin_Rmin_arr[ai,vi,pos_fin,6]
fig.add_subplot(221).plot(x,y, marker="o", linestyle='none', markersize=2.0, alpha=0.5)
fig.add_subplot(221).set_ylim(0,90)

fig.add_subplot(222).hist(y, bins=100, normed=True, histtype='step')
fig.add_subplot(222).set_xlim(0,50)

x = Efin_Rmin_arr[ai,vi,pos_fin,7]**4.
y = Efin_Rmin_arr[ai,vi,pos_fin,8]
fig.add_subplot(223).plot(x,y, marker="o", linestyle='none', markersize=2.0, alpha=0.5)

plt.show()
exit()
#----------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------
#--------------------------------------------------------------------------------------------------------------------


















#----------------------------------------------------------
#PLOT Histograms:
#----------------------------------------------------------

fig = plt.figure(figsize=(12, 10))

#plot for some a,v (indices):
ai = 0
vi = 0
#loop over bin-sin combinations:
for noc in range(0,3):
    #2body object combination:   [[1,2],[1,3],[2,3]]
    obj2list = list(obj_combs_arr[noc,:])
    #-------------------------
    #Direct:
    #-------------------------
    DIRIi = 0
    pos_E_arr   = np.where(Etot_ij_k_DI_RI_arr[ai,vi,:,DIRIi,noc] > 0.0)[0]  #need this [0] in the end.
    nr_E_arr    = len(pos_E_arr)
    #if values exist, plot:
    if (nr_E_arr > 0):
        E_arr       = Etot_ij_k_DI_RI_arr[ai,vi,pos_E_arr,DIRIi,noc]
        print obj2list, np.mean(E_arr)
        E_arr_norm  = E_arr/(5e-3)
        #plot histogram:
        fig.add_subplot(221).hist(E_arr_norm, bins=50, normed=True, histtype='step', label=str(obj2list))          
        fig.add_subplot(221).set_xlim(1e-1,50)
        fig.add_subplot(221).set_ylim(1e-5,2)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        #plot cumulative histogram:
        fig.add_subplot(222).hist(E_arr_norm, bins=50, normed=True, histtype='step', label=str(obj2list), cumulative=True)          
        fig.add_subplot(222).set_xlim(1e-1,15)
        fig.add_subplot(222).set_ylim(1e-5,1.1)
        plt.legend()
            
    #-------------------------      
    #Resonance:
    #-------------------------
    DIRIi = 1    
    pos_E_arr   = np.where(Etot_ij_k_DI_RI_arr[ai,vi,:,DIRIi,noc] > 0.0)[0]  #need this [0] in the end.
    nr_E_arr    = len(pos_E_arr)
    #if values exist, plot:
    if (nr_E_arr > 0):
        E_arr       = Etot_ij_k_DI_RI_arr[ai,vi,pos_E_arr,DIRIi,noc]
        print obj2list, np.mean(E_arr)
        E_arr_norm  = E_arr/(5e-3)
        #plot histogram:
        fig.add_subplot(223).hist(E_arr_norm, bins=50, normed=True, histtype='step', label=str(obj2list))          
        fig.add_subplot(223).set_xlim(1e-1,50)
        fig.add_subplot(223).set_ylim(1e-5,2)
        plt.yscale('log')
        plt.xscale('log')
        plt.legend()
        #plot cumulative histogram:
        fig.add_subplot(224).hist(E_arr_norm, bins=50, normed=True, histtype='step', label=str(obj2list), cumulative=True)                  
        fig.add_subplot(224).set_xlim(1e-1,15)
        fig.add_subplot(224).set_ylim(1e-5,1.1)
        plt.legend()
        
    #-------------------------    
    #Direct and Resonance:
    #-------------------------
    


    #-------------------------

plt.show()
#----------------------------------------------------------
#----------------------------------------------------------










#E,L and collisions:

#MAKE SURE WE COMPARE SAME DATASET!! VERY SENSITIVE!!!

#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  

fig = plt.figure(figsize=(10, 10))

#Constant E - Vary L:
#set E(a)
ai = 5
vi = 0
#all L from experiment:
pos_okdata  = np.where(Efin_Rmin_arr[ai,vi,:,0] > -1)[0]
pos_RI      = np.where(Efin_Rmin_arr[ai,vi,:,2] == 1)[0]
pos_fin     = list(set(pos_okdata).intersection(pos_RI))
data        = Efin_Rmin_arr[ai,vi,pos_fin,10]
min_L       = min(data)
max_L       = max(data)
N_bins      = 10.
HIST_L_ALL      = np.histogram(data, range=(min_L,max_L), bins=N_bins, density=True)[0]
#find coll from RI:
pos_id_2    = np.where(Efin_Rmin_arr[ai,vi,:,0] == 2)[0]
pos_fin     = list(set(pos_id_2).intersection(pos_RI))
data        = Efin_Rmin_arr[ai,vi,pos_fin,10]
HIST_L_coll_RI  = np.histogram(data, range=(min_L,max_L), bins=N_bins, density=True)[0]
#weighted histogram:
P_HIST          = HIST_L_coll_RI/HIST_L_ALL
#make x-plot arr:
L_arr       = np.linspace(min_L,max_L, num=N_bins)
#plot:
fig.add_subplot(131).plot(L_arr, HIST_L_ALL)
fig.add_subplot(131).plot(L_arr, HIST_L_coll_RI)
fig.add_subplot(131).plot(L_arr, P_HIST)

#Constant L - Vary E:
#choose L bin:
Lbin_min = 10.
Lbin_max = 12.
HT = np.zeros(nr_SMA, dtype=np.float64) 

for ac in range(0,nr_SMA):                               # ac-loop over SMA

    pos_RI      = np.where(Efin_Rmin_arr[ac,vi,:,2] == 1)[0]

    pos_GT_Lmin     = np.where(Efin_Rmin_arr[ac,vi,:,10] > Lbin_min)[0]
    pos_LT_Lmax     = np.where(Efin_Rmin_arr[ac,vi,:,10] < Lbin_max)[0]
    pos_GL_LMAX     = list(set(pos_GT_Lmin).intersection(pos_LT_Lmax))
    pos_GL_LMAX_RI  = list(set(pos_GL_LMAX).intersection(pos_RI))
    nr_all      = len(pos_GL_LMAX_RI)
    
    pos_id_2    = np.where(Efin_Rmin_arr[ac,vi,:,0] == 2)[0]
    pos_coll_RI_GL_LMAX = list(set(pos_id_2).intersection(pos_GL_LMAX_RI))
    nr_coll     = len(pos_coll_RI_GL_LMAX)    
    #calc coll frac:
    frac_coll_all = (1.*nr_coll)/(1.*nr_all)
    
    if (nr_all > 0): HT[ac] = frac_coll_all
    

x = 1./MCvalues_SMA_arr
y = HT
fig.add_subplot(132).plot(x, y)


for ac in range(0,nr_SMA):                               # ac-loop over SMA
    ai = ac
    vi = 0
    pos_id_2    = np.where(Efin_Rmin_arr[ai,vi,:,0] == 2)[0]
    pos_RI      = np.where(Efin_Rmin_arr[ai,vi,:,2] == 1)[0]
    pos_fin     = list(set(pos_id_2).intersection(pos_RI))
    x = Efin_Rmin_arr[ai,vi,pos_fin,9]
    y = Efin_Rmin_arr[ai,vi,pos_fin,10]

    fig.add_subplot(133).plot(x,y, marker="o", linestyle='none', markersize=5.0, alpha=0.5)

plt.show()
exit()












#----------------------------------------------------------
#Distribution of endstate params:
#----------------------------------------------------------
vi = 0
#ENDSTATE CODES: 1 = TDE, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 5 = Inspiral, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...
#Efin_Rmin_arr           = np.zeros((nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 20), dtype=np.float64)       #[a, v, nr_scatterings_per_paramcomb, [endstate_id, 12(0)/13(1)/23(2), DI(0)/RI(1), Rmin_12, Rmin_13, Rmin_23, fin_Etot_ijk, fin_a_ij, fin_e_ij, ini_E, ini_L]]  
#All (bin-sin (TO-CO), tidal insp, GW insp):

#Distribution if r_p do not depend on the masses,
#but we might use M and R for other things.

MTO = 0.6                                                   #SET THIS!
RTO = 0.013*((1.43/MTO)**(1./3.))*((1.-MTO/1.43)**(0.447))  #SET THIS!

MCO = 1.4
RCO = 1.7246e-5

MO  = MTO   #used only in equal mass case for f_GW axis.

vel_cs  = 10.0  #in km/sec

fig, ax1 = plt.subplots(figsize=(5, 4))

ax1_rp_min  = 1e-5
ax1_rp_max  = 1e2

colornr = [50,150,200]
#colornr  = np.linspace(0, 225, num=nr_SMA, dtype=int)


for ac in range(0,nr_SMA):
    
    #calc cs fac:
    bmax_AU     = RS_avs_outlist_MC_info_REAL[ac,0,0,2]/AU_U     #outlist_MC_info_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)            #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]
    CSfac_AU    = (1./nr_scatterings_per_paramcomb)*(np.pi*(bmax_AU**2.))*(1./vel_cs**2.) 
    
    pos_bs      = np.where(Efin_Rmin_arr[ac,vi,:,0] == 3)[0]
    pos_insp    = np.where(Efin_Rmin_arr[ac,vi,:,0] == 5)[0]
    pos_TDE     = np.where(Efin_Rmin_arr[ac,vi,:,0] == 1)[0]
    pos_12      = np.where(Efin_Rmin_arr[ac,vi,:,1] == 0)[0]
    pos_13      = np.where(Efin_Rmin_arr[ac,vi,:,1] == 1)[0]
    pos_23      = np.where(Efin_Rmin_arr[ac,vi,:,1] == 2)[0]
    
    pos_bs_12   = list(set(pos_bs).intersection(pos_12))
    pos_bs_13   = list(set(pos_bs).intersection(pos_13))
    pos_bs_23   = list(set(pos_bs).intersection(pos_23))
    
    pos_insp_12 = list(set(pos_insp).intersection(pos_12))
    pos_insp_13 = list(set(pos_insp).intersection(pos_13))
    pos_insp_23 = list(set(pos_insp).intersection(pos_23))
    
    pos_TDE_12 = list(set(pos_TDE).intersection(pos_12))
    pos_TDE_13 = list(set(pos_TDE).intersection(pos_13))
    pos_TDE_23 = list(set(pos_TDE).intersection(pos_23))
    
    pos_bs_1213         = pos_bs_12 + pos_bs_13
    pos_inspTDE_121323  = (pos_insp_12 + pos_insp_13 + pos_insp_23) + (pos_TDE_12 + pos_TDE_13 + pos_TDE_23)
    
    pos_inspTDE_1213    = pos_insp_12 + pos_insp_13 + pos_TDE_12 + pos_TDE_13
    pos_insp_23         = pos_insp_23
    
    #plot: bin-sin
    if (len(pos_bs_1213) > 1):
        pos_arr     = pos_bs_1213
        a_arr       = Efin_Rmin_arr[ac,vi,pos_arr,7]
        e_arr       = Efin_Rmin_arr[ac,vi,pos_arr,8] 
        rp_arr      = (1.-e_arr)*a_arr          #calc of r_p here assumes that no circulation or dissipation has taken place. Which is ok for binary-unbound-single endstates.
        pos_rp_GT_RTO   = np.where(rp_arr > RTO)[0] 
        if (len(pos_rp_GT_RTO) > 1):
            rp_arr          = rp_arr[pos_rp_GT_RTO]
            warr            = CSfac_AU*(1.+ np.zeros(len(rp_arr), dtype=np.float64))
            ax1.hist(np.log10(rp_arr), bins=100, range = [np.log10(ax1_rp_min),np.log10(ax1_rp_max)], weights = warr, normed=False, linewidth=0.25, histtype='step', color=plt.cm.terrain(colornr[ac]))
            
    #plot: tidal insp + TDE
    if (len(pos_inspTDE_1213) > 1):
        pos_arr     = pos_inspTDE_1213
        a_arr       = Efin_Rmin_arr[ac,vi,pos_arr,7]
        e_arr       = Efin_Rmin_arr[ac,vi,pos_arr,8] 
        rp_arr      = (1.-e_arr**2.)*a_arr/2.   #calc of r_p here assumes that the initial e was very close to 1 and L is concerved.
        pos_rp_GT_RTO   = np.where(rp_arr > RTO)[0] 
        if (len(pos_rp_GT_RTO) > 1):
            rp_arr          = rp_arr[pos_rp_GT_RTO]    
            warr            = CSfac_AU*(1.+ np.zeros(len(rp_arr), dtype=np.float64))
            ax1.hist(np.log10(rp_arr), bins=100, range = [np.log10(ax1_rp_min),np.log10(ax1_rp_max)], weights = warr, normed=False, linewidth=1.0, histtype='step', color=plt.cm.terrain(colornr[ac]))

    #plot: GW insp
    if (len(pos_insp_23) > 1):
        pos_arr     = pos_insp_23
        a_arr       = Efin_Rmin_arr[ac,vi,pos_arr,7]
        e_arr       = Efin_Rmin_arr[ac,vi,pos_arr,8] 
        rp_arr      = (1.-e_arr**2.)*a_arr/2.   #calc of r_p here assumes that the initial e was very close to 1 and L is concerved.
        pos_rp_GT_RTO   = np.where(rp_arr > 2.*RCO)[0] 
        if (len(pos_rp_GT_RTO) > 1):
            rp_arr          = rp_arr[pos_rp_GT_RTO]
            warr            = CSfac_AU*(1.+ np.zeros(len(rp_arr), dtype=np.float64))  
            ax1.hist(np.log10(rp_arr), bins=100, range = [np.log10(ax1_rp_min),np.log10(ax1_rp_max)], weights = warr, normed=False, linewidth=1.0, histtype='step', color=plt.cm.terrain(colornr[ac]))
                
    #for legend:
    ax1.plot(-1, -1, linewidth=1.0, linestyle='-', color=plt.cm.terrain(colornr[ac]), label=r'$a_{0}$ [AU]='+str('%.0e' % (MCvalues_SMA_arr[ac]/AU_U)))


#def f_GW_Hz_high_ecc(rp_Rsun_arr, mtot_Msun):   #ONLY FOR ecc ~ 1 !!!!
#x-axis limits:
ax1_xlim    = [np.log10(ax1_rp_min),                            np.log10(ax1_rp_max)]
ax3_xlim    = [np.log10(f_GW_Hz_high_ecc(ax1_rp_min,(MO+MO))),  np.log10(f_GW_Hz_high_ecc(ax1_rp_max,(MO+MO)))]
#make upper x-axis:
ax3 = ax1.twiny()
ax3.plot(-1, -1,    linewidth=0.0)  #dummy plot
#set x-axis limits:
ax1.set_xlim(ax1_xlim)
ax3.set_xlim(ax3_xlim)
#set y-axis limit:
ax1.set_ylim(1e-5,1e1)
#axis settings:
plt.yscale('log')
#text:
ax1.text(-4, 2e-2, r'NS-NS GW inspirals', size = 6,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=0)
ax1.text(-2.5, 9e-2, r'WD-NS tidal inspirals', size = 6,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=0)
ax1.text(0, 3e0, r'WD-NS bin-sin endstates', size = 6,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=0)
#labels:  
ax1.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)
ax1.set_xlabel(r'log $r_{p}$ $[R_{\odot}]$')
ax3.set_xlabel(r'log $f_{GW}$ $[Hz]$ ($m=1.2M_{\odot}$, $e \sim 1)$')
ax1.set_ylabel(r'cross section $\sigma$ [AU$^{2}$] (10 km s$^{-1}$) ')

#Save figure:
#plt.savefig('rp_hist_WDNSNS.eps', bbox_inches='tight')
plt.show()

exit()
#----------------------------------------------------------













#----------------------------------------------------------
#PLOT (Piet Hut):
#----------------------------------------------------------
#cross_section_allstates_DIRI[ac,vc,:,:]       = CSfac*prob_eachstate_DIRI[ac,vc,:,:]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]
#nr_eachstate_sim_RI[noc,:]  = [RI_nrc_123_binsin[noc], RI_nrc_123_ionization[noc], RI_nrc_123_headon_collisions[noc], RI_nrc_123_insp_collisions[noc], RI_nrc_123_userstate_20[noc], RI_nrc_123_userstate_21[noc], RI_nrc_123_userstate_22[noc], RI_nrc_123_inspiral[noc], RI_nrc_123_tidalevent[noc], 0.0]       

fig = plt.figure(figsize=(12, 8))

v_vc        = (0.264134504066*np.sqrt(1./10))
v_vscaled   = MCvalues_vinf_arr[:]/v_vc
cs_scalefac = np.pi*(MCvalues_SMA_arr[0]**2)
#---------------------
#plot: exchange:
#---------------------
csp_exchange            = cross_section_allstates_DIRI[0,:,0,0]
cs      = csp_exchange/cs_scalefac
fig.add_subplot(111).plot(v_vscaled,cs, marker='o', markersize=5.0, linewidth=2.5, alpha=0.75, color = 'green')

csp_exchange            = cross_section_allstates_DIRI[0,:,1,0]
cs      = csp_exchange/cs_scalefac
fig.add_subplot(111).plot(v_vscaled,cs, marker='o', markersize=5.0, linewidth=2.5, alpha=0.75, color = 'blue')

csp_exchange            = cross_section_allstates_DIRI[0,:,2,0]
cs      = csp_exchange/cs_scalefac
fig.add_subplot(111).plot(v_vscaled,cs, marker='o', markersize=5.0, linewidth=2.5, alpha=0.75, color = 'red')

#---------------------
#plot: ionization:
#---------------------
csp_ionization          = cross_section_allstates_DIRI[0,:,0,1]
cs      = csp_ionization/cs_scalefac
fig.add_subplot(111).plot(v_vscaled,cs, marker='o', markersize=5.0, linewidth=2.5, alpha=0.75, color = 'orange')
#---------------------
#plot: collisions:
#---------------------
#csp_headon_collisions   = cross_section_allstates_DIRI[0,:,0,2]+cross_section_allstates_DIRI[0,:,1,2]+cross_section_allstates_DIRI[0,:,2,2]
#cs      = csp_headon_collisions/cs_scalefac
#fig.add_subplot(111).plot(v_vscaled,cs, marker='o', markersize=5.0, linewidth=2.5, alpha=0.75)
#---------------------

#plt.xlim(0.4,   10)
#plt.ylim(0.01,  10)

plt.yscale('log')
plt.xscale('log')
plt.show()
exit()
#----------------------------------------------------------








#----------------------------------------------------------
#PLOT Cross sections:
#----------------------------------------------------------

fig = plt.figure(figsize=(6, 5))

#DEFINE:
#in this plot we have 'v' as a constant and vary 'a':
vi = 0
#scale from R_sun to AU:
ap_AU           = MCvalues_SMA_arr[:]/AU_U
cs_AU           = (R_sun_SI/AU_SI)**2.
cs_arr_AU       = cs_AU*cross_section_allstates_DIRI[:,:,:,:]   
cs_err_arr_AU   = cs_AU*err_cross_section_allstates_DIRI[:,:,:,:]

##2body object combination: [[1,2],[1,3],[2,3]]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]

#plot cs:
csp = (cs_arr_AU[:,vi,1,0] + cs_arr_AU[:,vi,2,0])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle=':', markersize=10.0, linewidth=1.5, alpha=0.75, label=r'Exchange')

#plot cs:
csp = (cs_arr_AU[:,vi,0,8])+(cs_arr_AU[:,vi,1,8])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, color='green', label=r'TDE')

#plot cs:
csp = (cs_arr_AU[:,vi,0,2]+cs_arr_AU[:,vi,1,2])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle=':', markersize=10.0, linewidth=1.5, alpha=0.75, label=r'Head-on Coll [12,13]')

#plot cs:
csp = (cs_arr_AU[:,vi,0,3]+cs_arr_AU[:,vi,1,3]+cs_arr_AU[:,vi,2,3])
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, label=r'Insp coll')

#plot cs:
csp = (cs_arr_AU[:,vi,0,7]+cs_arr_AU[:,vi,1,7])
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle=':', markersize=10.0, linewidth=1.5, alpha=0.75, label=r'Insp [12,13]')

#plot cs:
csp = (cs_arr_AU[:,vi,2,7])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle=':', markersize=10.0, linewidth=1.5, alpha=0.75, label=r'Insp [23]')

#plot cs:
csp = (cs_arr_AU[:,vi,0,2]+cs_arr_AU[:,vi,1,2]+cs_arr_AU[:,vi,2,2] + cs_arr_AU[:,vi,0,3]+cs_arr_AU[:,vi,1,3]+cs_arr_AU[:,vi,2,3] + cs_arr_AU[:,vi,0,7]+cs_arr_AU[:,vi,1,7]+cs_arr_AU[:,vi,2,7] + cs_arr_AU[:,vi,0,8]+cs_arr_AU[:,vi,1,8]+cs_arr_AU[:,vi,2,8])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75, color='black', label=r'All coll [12,13,23]')

#plot cs:
csp = (cs_arr_AU[:,vi,0,8]+cs_arr_AU[:,vi,1,8]+cs_arr_AU[:,vi,2,8])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, label=r'TEST 3')

csp = (cs_arr_AU[:,vi,0,3]+cs_arr_AU[:,vi,1,3]+cs_arr_AU[:,vi,2,3]) + (cs_arr_AU[:,vi,0,7]+cs_arr_AU[:,vi,1,7])+(cs_arr_AU[:,vi,2,7])
print csp
fig.add_subplot(111).plot(ap_AU,csp, marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, label=r'ALL')


#plot settings:
csp_exchange = (cs_arr_AU[:,vi,1,0] + cs_arr_AU[:,vi,2,0])
fig.add_subplot(111).set_xlim(0.75*min(ap_AU), 5.*max(ap_AU))
#fig.add_subplot(111).set_ylim(0.01*min([1e-6]+csp_exchange), 10*max(csp_exchange))
fig.add_subplot(111).set_ylim(1e-5, 1e8)
plt.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)
plt.xlabel(r'initial SMA $a_{0}$ [AU]')
plt.ylabel(r'cross section $\sigma$ [AU$^2$] (x km/sec)')
plt.title(data_name)
plt.yscale('log')
plt.xscale('log')

#Save figure:
plt.savefig('crossections_' + data_name + '.pdf', bbox_inches='tight')



#ap_AU           = MCvalues_SMA_arr[:]/AU_U
#cs_AU           = (R_sun_SI/AU_SI)**2.
#cs_arr_AU       = cs_AU*cross_section_allstates_DIRI[:,:,:,:]   
#cs_err_arr_AU   = cs_AU*err_cross_section_allstates_DIRI[:,:,:,:]



#Save data:
nr_a        = len(ap_AU)
nr_cs       = 20    #dont change. We just always keep space for ...
save_arr_cs     = np.zeros((nr_cs,nr_a), dtype=np.float64)
save_arr_cs_err = np.zeros((nr_cs,nr_a), dtype=np.float64)
#decide which cs to save:
#2body object combination: [[1,2],[1,3],[2,3]]
#nr_eachstate_sim_DI[noc,:]  = [DI_nrc_123_binsin[noc], DI_nrc_123_ionization[noc], DI_nrc_123_headon_collisions[noc], DI_nrc_123_insp_collisions[noc], DI_nrc_123_userstate_20[noc], DI_nrc_123_userstate_21[noc], DI_nrc_123_userstate_22[noc], DI_nrc_123_inspiral[noc], DI_nrc_123_tidalevent[noc], 0.0]
#format: cs_arr_AU[a,v,obj comb, cs]
#cross sections:
save_arr_cs[0,:]        = ap_AU
save_arr_cs[1,:]        = cs_arr_AU[:,vi,1,0] + cs_arr_AU[:,vi,2,0]         #exchange
save_arr_cs[2,:]        = cs_arr_AU[:,vi,0,2]                               #collisions 12
save_arr_cs[3,:]        = cs_arr_AU[:,vi,1,2]                               #collisions 13
save_arr_cs[4,:]        = cs_arr_AU[:,vi,2,2]                               #collisions 23
save_arr_cs[5,:]        = cs_arr_AU[:,vi,0,3]                               #inspiral collision 12
save_arr_cs[6,:]        = cs_arr_AU[:,vi,1,3]                               #inspiral collision 13
save_arr_cs[7,:]        = cs_arr_AU[:,vi,2,3]                               #inspiral collision 23
save_arr_cs[8,:]        = cs_arr_AU[:,vi,0,7]                               #inspiral 12
save_arr_cs[9,:]        = cs_arr_AU[:,vi,1,7]                               #inspiral 13
save_arr_cs[10,:]       = cs_arr_AU[:,vi,2,7]                               #inspiral 23
save_arr_cs[15,:]       = cs_arr_AU[:,vi,0,8]                                       #TDE 12
save_arr_cs[16,:]       = cs_arr_AU[:,vi,1,8]                                       #TDE 13
#cross sections err:
save_arr_cs_err[0,:]    = ap_AU
save_arr_cs_err[1,:]    = np.sqrt(cs_err_arr_AU[:,vi,1,0]**2 + cs_err_arr_AU[:,vi,2,0]**2) #exchange
save_arr_cs_err[2,:]    = cs_err_arr_AU[:,vi,0,2]                           #collisions 12
save_arr_cs_err[3,:]    = cs_err_arr_AU[:,vi,1,2]                           #collisions 13
save_arr_cs_err[4,:]    = cs_err_arr_AU[:,vi,2,2]                           #collisions 23
save_arr_cs_err[5,:]    = cs_err_arr_AU[:,vi,0,3]                           #inspiral collision 12
save_arr_cs_err[6,:]    = cs_err_arr_AU[:,vi,1,3]                           #inspiral collision 13
save_arr_cs_err[7,:]    = cs_err_arr_AU[:,vi,2,3]                           #inspiral collision 23
save_arr_cs_err[8,:]    = cs_err_arr_AU[:,vi,0,7]                           #inspiral 12
save_arr_cs_err[9,:]    = cs_err_arr_AU[:,vi,1,7]                           #inspiral 13
save_arr_cs_err[10,:]   = cs_err_arr_AU[:,vi,2,7]                           #inspiral 23
save_arr_cs_err[15,:]   = cs_err_arr_AU[:,vi,0,8]                                   #TDE 12
save_arr_cs_err[16,:]   = cs_err_arr_AU[:,vi,1,8]                                   #TDE 13
#save cs:
tf = open('cs_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs,   fmt='%5f')
tf.close()
#save cs err:
tf = open('cs_err_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs_err,   fmt='%5f')
tf.close()

#show figures:
plt.show()

exit()

#----------------------------------------------------------
#----------------------------------------------------------




