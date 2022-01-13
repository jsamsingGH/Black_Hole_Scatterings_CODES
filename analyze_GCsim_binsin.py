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
import csv



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


#----------------------------------------------------------
#GW inspiral life time:
#----------------------------------------------------------
tf = open('a_e_tlife_yrs_interpoltable.txt', "r")
a_e_tlife_yrs_arr   = np.loadtxt(tf, dtype=float)
tf.close()
#read data and make interpolation function:
tdata_e_arr                 = a_e_tlife_yrs_arr[:,1]
tdata_t_arr                 = a_e_tlife_yrs_arr[:,2]
ecc_to_tdata_interpol_func  = interpolate.interp1d(tdata_e_arr, tdata_t_arr)  
#----------------------------------------------------------


#----------------------------------------------------------
#Save analyzed data:
#----------------------------------------------------------
save_data_name  = 'MOCCA_test_1_'
#----------------------------------------------------------
#----------------------------------------------------------
#Read in SIM ICs data:
#----------------------------------------------------------
#data name:
data_name       = 'MOCCA_test_1_'    
#input data folder:
data_folder     = '/Users/jsamsing/Desktop/TIDES_PROJ/'
#read data:
tf = open(data_folder+data_name+'nbody_params_INT.txt', "r")
nbody_params_INT_all        = np.loadtxt(tf, dtype=int)         #[0, use_25PN, 1, 1000000000, 0, 0, nr_few_body, 0,0,0]
tf.close()
tf = open(data_folder+data_name+'nbody_params_REAL.txt', "r")
nbody_params_REAL_all       = np.loadtxt(tf, dtype=float)       #[0.01, simtime_U, 0.0, 0.01, simtim_secs, 0.5, 0.0, insp_threshold, 0,0]
tf.close()
tf = open(data_folder+data_name+'obj_info.txt', "r")
obj_info_all                = np.loadtxt(tf, dtype=float)       #[b1_const_arr, b2_const_arr, b3_const_arr]
tf.close()
tf = open(data_folder+data_name+'pos_vel.txt', "r")
pos_vel_all                 = np.loadtxt(tf, dtype=float)       #[posvel(1), posvel(2), posvel(3)]
tf.close()
#----------------------------------------------------------
#----------------------------------------------------------
#Read in SIM results data:
#----------------------------------------------------------
#data name:
data_name = 'MOCCA_test_1_'
#input data folder:
data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'
#read data:
tf = open(data_folder+data_name+'MC_settings_list_INT.txt', "r")
MC_settings_list_INT        = np.loadtxt(tf, dtype=int)                     #[0, 0, 0, nr_tot_scatterings, 0]
tf.close()
tf = open(data_folder+data_name+'MC_settings_list_REAL.txt', "r")
MC_settings_list_REAL       = np.loadtxt(tf, dtype=float)                   #[0.0, 0.0, 0.0, 0.0, 0.0]    
tf.close()
tf = open(data_folder+data_name+'outlist_MC_info_INT.txt', "r")
outlist_MC_info_INT         = np.loadtxt(tf, dtype=int)                     #[icidc, 0, 0, 0, 0, 0, 0, 0, 0, 0]
tf.close()
tf = open(data_folder+data_name+'outlist_MC_info_REAL.txt', "r")
outlist_MC_info_REAL        = np.loadtxt(tf, dtype=float)                   #[0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
tf.close()
tf = open(data_folder+data_name+'output_Nbody_endstate_INT.txt', "r")
output_Nbody_endstate_INT   = np.loadtxt(tf, dtype=int)                     #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
tf.close()
tf = open(data_folder+data_name+'output_Nbody_endstate_REAL.txt', "r")
output_Nbody_endstate_REAL  = np.loadtxt(tf, dtype=float)                   #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_info_INT.txt', "r")
output_Nbody_xtra_info_INT   = np.loadtxt(tf, dtype=int)                    #[0, 0, 0, 0, 0, 0, 0, 0, 0, 0]
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_info_REAL.txt', "r")
output_Nbody_xtra_info_REAL  = np.loadtxt(tf, dtype=float)                  #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e, 0.0, 0.0, 0.0, 0.0, 0.0]
tf.close()
#1 = TDE, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 5 = Inspiral, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...
#----------------------------------------------------------


#----------------------------------------------------------
#Define:
#----------------------------------------------------------
nr_tot_scatterings  = len(output_Nbody_endstate_INT)
#----------------------------------------------------------
#Allocate arrays:
#----------------------------------------------------------
simresults_arr_INT                  = np.zeros((nr_tot_scatterings,10), dtype='i')  # ... 
simresults_arr_REAL                 = np.zeros((nr_tot_scatterings,10), dtype='d')  #
#----------------------------------------------------------


#----------------------------------------------------------
#loop over all scatterings:
#----------------------------------------------------------

for nc in range(0,nr_tot_scatterings):
           
    #-----------------------------------
    #read info:
    #-----------------------------------
    #endstate:
    endstate_id     = output_Nbody_endstate_INT[nc,0]
    endstate_bini   = output_Nbody_endstate_INT[nc,1]
    endstate_binj   = output_Nbody_endstate_INT[nc,2]
    endstate_bink   = 6 - (endstate_bini + endstate_binj)
    #mass:
    m1_Msun     = obj_info_all[nc,0]
    m2_Msun     = obj_info_all[nc,10]
    m3_Msun     = obj_info_all[nc,20]
    m123_Msun   = np.array([m1_Msun, m2_Msun, m3_Msun])
    mi_Msun     = m123_Msun[endstate_bini-1]
    mj_Msun     = m123_Msun[endstate_binj-1]
    mk_Msun     = m123_Msun[endstate_bink-1]
    #-----------------------------------
    
    #-----------------------------------
    #input/define:
    #-----------------------------------
    #GW threshold:
    f_GW    = 10.0 #in Hz
    #mass in SI:
    m1_SI   = m1_Msun*M_sun_SI
    m2_SI   = m2_Msun*M_sun_SI
    m3_SI   = m3_Msun*M_sun_SI
    mi_SI   = mi_Msun*M_sun_SI
    mj_SI   = mj_Msun*M_sun_SI
    mk_SI   = mk_Msun*M_sun_SI
    #-----------------------------------
    
    #-----------------------------------
    if (endstate_id == 3):
    #-----------------------------------     
        #initial a,e after BS endstate:
        a_BS_SI     = output_Nbody_endstate_REAL[nc,3]*R_sun_SI        
        e_BS        = output_Nbody_endstate_REAL[nc,4]
    
        #calc tlife: (fast, assumes e>>0 limit)
        tbin_yr_circular    = ((a_BS_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*mi_SI*mj_SI*(mi_SI+mj_SI)/(c_SI**5.)))/sec_year
        tbin_yr             = ((768./425.)*tbin_yr_circular*((1.-e_BS**2.)**(7./2.)))            
        log10tbin_yr        = np.log10(tbin_yr)
        #calc v_kick:
        Etot_binsin         = output_Nbody_endstate_REAL[nc,7]
        ESI_fac             = G_new_SI*(M_sun_SI**2.)/R_sun_SI
        Etot_SI             = Etot_binsin*ESI_fac
        mbin                = mi_SI+mj_SI
        msin                = mk_SI
        v_inf_bin           = np.sqrt(2.*Etot_SI*((mbin*(1.+mbin/msin))**(-1.)))   #SI units
        v_kick_kms          = v_inf_bin/1000.    
    
        #propagating ini a,e to a,e(f_GW)
        c0      = a_BS_SI/((e_BS**(12./19.)/(1.-e_BS**2.))*((1.+(121./304.)*e_BS**2.)**(870./2299.)))
        #a(ecc)  = (c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))
        func    = lambda ecc : f_GW - (1./np.pi)*np.sqrt(G_new_SI*(mi_SI+mj_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
        ecc_initial_guess   = 1e-8
        fGW_BS_e            = fsolve(func, ecc_initial_guess)   #value of bin ecc when the bin is at f_GW   
        log10fGW_BS_e       = np.log10(fGW_BS_e)
        
        #save results:
        simresults_arr_INT[nc,:]    = 0
        simresults_arr_REAL[nc,:]   = [log10tbin_yr, v_kick_kms, log10fGW_BS_e, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        
        print np.log10(tbin_yr), a_BS_SI/AU_SI, e_BS, v_kick_kms, m123_Msun, mi_Msun, mj_Msun, mk_Msun         
    #-----------------------------------           
    
    #-----------------------------------
    if (endstate_id == 5):
    #-----------------------------------    
        #initial a,e at formation of IMS
        a_IMS_SI    = output_Nbody_xtra_info_REAL[nc,3]*R_sun_SI
        e_IMS       = output_Nbody_xtra_info_REAL[nc,4]
        #propagating ini a,e to a,e(f_GW)
        c0      = a_IMS_SI/((e_IMS**(12./19.)/(1.-e_IMS**2.))*((1.+(121./304.)*e_IMS**2.)**(870./2299.)))
        #a(ecc)  = (c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))
        func    = lambda ecc : f_GW - (1./np.pi)*np.sqrt(G_new_SI*(mi_SI+mj_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
        ecc_initial_guess   = 0.001
        fGW_IMS_e           = fsolve(func, ecc_initial_guess)   #value of bin ecc when the bin is at f_GW
        log10fGW_IMS_e      = np.log10(fGW_IMS_e)
        
        #save results:
        simresults_arr_INT[nc,:]    = 0
        simresults_arr_REAL[nc,:]   = [0.0, 0.0, log10fGW_IMS_e, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #maybe also save:
        # - nr pericenter passages in the LIGO band.
        # - ...
    #-----------------------------------

#----------------------------------------------------------   
#----------------------------------------------------------


#----------------------------------------------------------
#save output results:
#----------------------------------------------------------
tf = open(save_data_name+'simresults_arr_REAL.txt', "w")
np.savetxt(tf, simresults_arr_REAL,   fmt='%5f')
tf.close()
#----------------------------------------------------------




exit()













