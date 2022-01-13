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
from matplotlib import colors as ccolor
from matplotlib import cm
import re    
import matplotlib.gridspec as gridspec

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
yr_sec     = 31536000.0    
#----------------------------------------------------------


#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


#----------------------------------------------------------
#input data folders and names:
#----------------------------------------------------------
#data folders:
input_data_folder   = '/Users/jsamsing/Desktop/TIDES_PROJ/MOCCA_ICs/MC_1/'  #initial conditions and input files
output_data_folder  = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'       #output data from N-body simulations

#data names:
#oridata_name        = 'MT1000_inter-binsin-bhs-13gyr-m100.dat'#'MOCCATEST_inter-binsin-bhs-13gyr-m100.dat' #            #name of original data file (all lines before any filters are applied)
#icdata_name         = 'MT1000_case1_'#'MOCCATEST_case1A_'#'MOCCATEST_case0_' #                                          #name of the ICs generated for this case (case1, case2, ..., etc.).
#simdata_name        = 'OUT_MT1000_case1_'#'OUTSIM_MOCCATEST_case1A_'#'OUTSIM_MOCCATEST_case0_' #                        #name of the N-body sim output files resulting from the ICs with name "icdata_name".

oridata_name        = 'MOCCATEST_inter-binsin-bhs-13gyr-m100.dat' #            #name of original data file (all lines before any filters are applied)
icdata_name         = 'MOCCATEST_case1A_'#'MOCCATEST_case0_' #                                          #name of the ICs generated for this case (case1, case2, ..., etc.).
simdata_name        = 'OUTSIM_MOCCATEST_case1A_'#'OUTSIM_MOCCATEST_case0_' #                        #name of the N-body sim output files resulting from the ICs with name "icdata_name".
#----------------------------------------------------------

#----------------------------------------------------------
#calc info and save data:
#----------------------------------------------------------
calcinfo_and_save_data   = 0
#----------------------------------------------------------

#----------------------------------------------------------
#get info:
#----------------------------------------------------------
#nr ori ICs:
tf              = open(input_data_folder+oridata_name,  "r")
nr_ori_ICs      = len(tf.readlines())
#nr tot sim:
tf              = open(output_data_folder+simdata_name+'output_Nbody_endstate_INT.txt',   "r")
nr_tot_sims     = len(tf.readlines())
#nr sim per IC:
nr_sim_per_IC   = 5 #5 #10
#nr sim ICs:
nr_sim_ICs      = nr_tot_sims/nr_sim_per_IC

print nr_ori_ICs, nr_tot_sims, nr_sim_per_IC, nr_sim_ICs
#----------------------------------------------------------


#----------------------------------------------------------
#....
#----------------------------------------------------------
#outlist_MC_info_INT                #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
#outlist_MC_info_REAL               #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, R_val, b_val, f_val, g_val]
#output_Nbody_endstate_INT          #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
#output_Nbody_endstate_REAL         #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#output_Nbody_xtra_info_REAL        #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
#output_Nbody_xtra_2_info_REAL      #[out_pos_CMij_wrt_sink(x,y,z), out_vel_CMij_wrt_sink(x,y,z), Lorb_ij(x,y,z)]
#----------------------------------------------------------


fGW_limit   = 10.0     #in Hz
efGW_limit  = 0.1
N_IMS       = 20.0


#--------------------------------------------------------------
if (calcinfo_and_save_data == 1):
#--------------------------------------------------------------

    #----------------------------------------------------------
    #Open files:
    #----------------------------------------------------------
    tf_oridata                          = open(input_data_folder+oridata_name,  "r")
    tf_infofile                         = open(input_data_folder+icdata_name+'sim_info_1.txt',  "r")

    tf_output_Nbody_endstate_INT        = open(output_data_folder+simdata_name+'output_Nbody_endstate_INT.txt',  "r")
    tf_output_Nbody_endstate_REAL       = open(output_data_folder+simdata_name+'output_Nbody_endstate_REAL.txt',  "r")
    tf_output_Nbody_xtra_info_REAL      = open(output_data_folder+simdata_name+'output_Nbody_xtra_info_REAL.txt',  "r")
    tf_output_Nbody_xtra_2_info_REAL    = open(output_data_folder+simdata_name+'output_Nbody_xtra_2_info_REAL.txt',  "r")
    #----------------------------------------------------------

    #----------------------------------------------------------
    #define and initialize:
    #----------------------------------------------------------
    #arrays:
    save_calc_arr_INT   = np.zeros((nr_tot_sims, 30), dtype=int)
    save_calc_arr_REAL  = np.zeros((nr_tot_sims, 30), dtype=float)
    
    rndnum_arr          = np.random.random((nr_tot_sims,10))
    
    ier = 0
    
    #counters:
    tc = 0
    #----------------------------------------------------------

    #----------------------------------------------------------
    #loop over all ICs:
    #----------------------------------------------------------
    for ni in range(0, nr_ori_ICs):
    
        #---------------------------------------
        #Read line at this 'ni':
        #---------------------------------------
        #tf_oridata:
        splitline_oridata   = re.split(r'\s{1,}', tf_oridata.readline())
        #tf_infofile:
        splitline_infofile  = re.split(r'\s{1,}', tf_infofile.readline())
    
        #sim yes/no::
        sim_yesno           = int(float(splitline_infofile[1]))    
    
        #---------------------------
        if (sim_yesno == 1):
        #---------------------------
            #-----------------------
            #read and define:
            #-----------------------
            line0 = 1
            #mass:
            mb1         = float(splitline_oridata[30-line0])*1.0                    #Msun
            mb2         = float(splitline_oridata[31-line0])*1.0                    #Msun
            mb3         = float(splitline_oridata[32-line0])*1.0                    #Msun
            mb123_arr   = np.array([mb1,mb2,mb3])
            #initial orbital parameters:
            a_0         = float(splitline_oridata[39-line0])*1.0                    #Rsun 
            e_0         = float(splitline_oridata[40-line0])*1.0
            #interaction time:
            tint_yr     = float(splitline_oridata[29-line0])*(10.**(6.0))           #in years
            #Initial central velocity of the cluster model in units of km/s:        
            vc0_kms     = float(splitline_oridata[16-line0])*1.0                    #km/sec
            #primordial binary yes(1)/no(0):
            binobjs_id_diff = abs(int(float(splitline_oridata[24-line0]))-int(float(splitline_oridata[25-line0])))
            if (binobjs_id_diff == 1): primorbin_yesno  = 1
            if (binobjs_id_diff != 1): primorbin_yesno  = 0
            #escape velocity of the cluster at the time of interaction:
            vesc_kms    = np.sqrt(2.0*abs(float(splitline_oridata[112-line0])*1.0)) #km/sec
            #initial cluster mass:
            m0_cluster      = float(splitline_oridata[10-line0])*1.0
            #initial cluster concentration:
            rho0_cluster    = float(splitline_oridata[17-line0])*1.0
            #initial cluster binary fraction:
            binfrac_cluster = float(splitline_oridata[4-line0])*1.0
            #-----------------------
            #loop over sim per IC:
            #-----------------------
            for ns in range(0, nr_sim_per_IC):
                #-------------------
                #read data:
                #-------------------
                #tf_output_Nbody_endstate_INT:
                splitline = re.split(r'\s{1,}', tf_output_Nbody_endstate_INT.readline())
                endstate_id         = int(float(splitline[1]))
                out_bin_i           = int(float(splitline[2]))
                out_bin_j           = int(float(splitline[3]))
                out_bin_k           = 6 - (out_bin_i+out_bin_j)
                nr_IMS              = int(float(splitline[7]))
                #tf_output_Nbody_endstate_REAL:
                splitline = re.split(r'\s{1,}', tf_output_Nbody_endstate_REAL.readline())
                a_ij                = float(splitline[3])
                e_ij                = float(splitline[4])
                Etot_ijk            = float(splitline[7])
                #output_Nbody_xtra_info_REAL:
                splitline = re.split(r'\s{1,}', tf_output_Nbody_xtra_info_REAL.readline())
                MRini_IMS_a         = float(splitline[3])
                MRini_IMS_e         = float(splitline[4])
                #output_Nbody_xtra_2_info_REAL:
                splitline = re.split(r'\s{1,}', tf_output_Nbody_xtra_2_info_REAL.readline())
                pos_CMijk_x         = float(splitline[0])
                pos_CMijk_y         = float(splitline[1])
                pos_CMijk_z         = float(splitline[2])
                pos_CMijk_xyz       = np.array([pos_CMijk_x, pos_CMijk_y, pos_CMijk_z])
                vel_CMijk_x         = float(splitline[3])
                vel_CMijk_y         = float(splitline[4])
                vel_CMijk_z         = float(splitline[5])
                vel_CMijk_xyz       = np.array([vel_CMijk_x, vel_CMijk_y, vel_CMijk_z])
                Lorb_ij_x           = float(splitline[6])
                Lorb_ij_y           = float(splitline[7])
                Lorb_ij_z           = float(splitline[8])
                #-------------------
                #calc/define:
                #-------------------
                #masses:
                mbi = mb123_arr[out_bin_i-1]
                mbj = mb123_arr[out_bin_j-1]
                mbk = mb123_arr[out_bin_k-1]
                mbi_SI  = mbi*M_sun_SI
                mbj_SI  = mbj*M_sun_SI
                mbk_SI  = mbk*M_sun_SI
                #-------------------
                #calc: ecc(X Hz)
                #-------------------
                #initialize:
                ecc_fGW               = 0.0
                ier = 1 #for solver info
                #coll:
                if (endstate_id == 2):
                    ecc_fGW           = 1.0                               #well-defined collision endstates.
                #GW insp:
                if (endstate_id == 5 and nr_IMS > 0):
                    #initial a,e:
                    ini_a_SI    = MRini_IMS_a*R_sun_SI
                    ini_e_SI    = abs(MRini_IMS_e - 1e-20)   #make sure 0 < ecc < 1 !!!
                    #check if its born in band:
                    ini_rp_SI   = ini_a_SI*(1.-ini_e_SI)
                    ini_fGW     = (1./np.pi)*(G_new_SI*(mbi_SI+mbj_SI)/(ini_rp_SI**3.))**(1./2.)
                    if (ini_fGW > fGW_limit):
                        ecc_fGW = ini_e_SI
                    if (ini_fGW < fGW_limit):                    
                        #propagating ini a,e to a,e(fGW_limit)
                        c0      = ini_a_SI/((ini_e_SI**(12./19.)/(1.-ini_e_SI**2.))*((1.+(121./304.)*ini_e_SI**2.)**(870./2299.)))
                        #a(ecc)  = (c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))
                        func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mbi_SI+mbj_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
                        ecc_initial_guess   = 1e-8
                        ecc_fGW, info, ier, msg = fsolve(func, ecc_initial_guess, full_output=1)   #value of bin ecc when the bin is at fGW_limit
                        if (ier == 1):  ecc_fGW = ecc_fGW[0]
                        if (ier != 1):  ecc_fGW = -1
                #bin-sin:
                if (endstate_id == 3):
                    #initial a,e:
                    ini_a_SI    = a_ij*R_sun_SI
                    ini_e_SI    = abs(e_ij - 1e-20)          #make sure 0 < ecc < 1 !!!
                    #check if its born in band:
                    ini_rp_SI   = ini_a_SI*(1.-ini_e_SI)
                    ini_fGW     = (1./np.pi)*(G_new_SI*(mbi_SI+mbj_SI)/(ini_rp_SI**3.))**(1./2.)
                    if (ini_fGW > fGW_limit):
                        ecc_fGW = ini_e_SI
                    if (ini_fGW < fGW_limit):                    
                        #propagating ini a,e to a,e(fGW_limit)
                        c0      = ini_a_SI/((ini_e_SI**(12./19.)/(1.-ini_e_SI**2.))*((1.+(121./304.)*ini_e_SI**2.)**(870./2299.)))
                        #a(ecc)  = (c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))
                        func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mbi_SI+mbj_SI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
                        ecc_initial_guess   = 1e-8
                        ecc_fGW, info, ier, msg = fsolve(func, ecc_initial_guess, full_output=1)   #value of bin ecc when the bin is at fGW_limit
                        if (ier == 1):  ecc_fGW = ecc_fGW[0]
                        if (ier != 1):  ecc_fGW = -1
                #------------------- 
                #save a,e at formation:
                #-------------------
                a_form = -1.0
                e_form = -1.0
                if (endstate_id == 5):
                    a_form    = MRini_IMS_a
                    e_form    = MRini_IMS_e
                if (endstate_id == 3):
                    a_form    = a_ij
                    e_form    = e_ij
                #------------------- 
                #calc: v_kick
                #-------------------
                #initialize:
                v_kick_kms          = -1.0
                #bin-sin:
                if (endstate_id == 3):
                    ESI_fac         = G_new_SI*(M_sun_SI**2.)/R_sun_SI
                    Etot_SI         = Etot_ijk*ESI_fac
                    mbin_SI         = mbi_SI+mbj_SI
                    msin_SI         = mbk_SI
                    v_inf_bin       = np.sqrt(2.*Etot_SI*((mbin_SI*(1.+mbin_SI/msin_SI))**(-1.)))   #SI units
                    v_kick_kms      = v_inf_bin/1000.   #vel kick of binary at infinity.
                #escape of binary yes/no:
                if (v_kick_kms < vesc_kms): bin_esc_yesno = 0
                if (v_kick_kms > vesc_kms): bin_esc_yesno = 1
                #-------------------
                #calc: t_life_ij
                #-------------------
                #initialize:
                tinsp_yr            = -1.0
                #bin-sin:
                if (endstate_id == 3):
                    #initial a,e:
                    ini_a_SI    = a_ij*R_sun_SI
                    ini_e_SI    = e_ij
                    #calc tlife: fast procedure (assumes e>>0 limit):
                    tinsp_yr_circular   = ((ini_a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*mbi_SI*mbj_SI*(mbi_SI+mbj_SI)/(c_SI**5.)))/yr_sec
                    tinsp_yr            = ((768./425.)*tinsp_yr_circular*((1.-ini_e_SI**2.)**(7./2.)))            
                #-------------------
                #calc: Lorb angle
                #-------------------
                #initialize:
                theta_Lorb_zvec     = -1.0
                #bin-sin:
                if (endstate_id == 3):
                    unit_z_vec  = np.array([0.0, 0.0, 1.0])
                    Lorb_ij     = np.array([Lorb_ij_x, Lorb_ij_y, Lorb_ij_z])
                    len_unit_z_vec  = 1.0
                    len_Lorb_ij     = np.sqrt(sum(Lorb_ij[:]**2.))
                    theta_Lorb_zvec = np.arccos(np.dot(Lorb_ij, unit_z_vec)/(len_Lorb_ij*len_unit_z_vec))
                #-------------------
                #calc eccentric GW probability:
                #-------------------                
                m_average   = ((mb1+mb2+mb3)/3.)    # average mass in Msun
                a_0_SI      = a_0*R_sun_SI
                r_fGW_SI    = (2.*G_new_SI*(m_average*M_sun_SI)/((fGW_limit**2.)*(np.pi**2.)))**(1./3.)
                F_efGW      = ((efGW_limit**(12./19.))/(1.+efGW_limit))*((1. + (121./304.)*(efGW_limit**2.))**(870./2299.))
                r_ec_SI     = (r_fGW_SI/F_efGW)*(1./2.)*((1. + (121./304.))**(870./2299.))
                P_rp_LT_rec         = (2.*r_ec_SI/a_0_SI)*N_IMS #P for forming during the resonating state (N_IMS tries)
                P_rp_LT_rec_finbin  = (2.*r_ec_SI/a_0_SI)*1.0   #P for the final ejected bin (only 1 try per interaction)
                #-------------------         
                #SPIN: calc chi_eff:
                #-------------------
                #initialize:
                chi_eff     = -1.0
                                
                if (endstate_id == 3):
                    #L orbital unit vec:
                    Lorb_ij         = np.array([Lorb_ij_x, Lorb_ij_y, Lorb_ij_z])
                    len_Lorb_ij     = np.sqrt(sum(Lorb_ij[:]**2.))
                    Lorb_ij_unitvec = Lorb_ij/len_Lorb_ij 
                    
                    #BH spin vectors:
                    #BH1:
                    theta   = np.arccos(2.*rndnum_arr[tc,0]-1.)         #P flat in cos(theta)	- pos on unitsphere
                    phi     = ((2.*np.pi)*rndnum_arr[tc,1])			    #P flat in phi			- pos on unitsphere
                    S1_rnd_unit = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
                    #BH2:
                    theta   = np.arccos(2.*rndnum_arr[tc,2]-1.)         #P flat in cos(theta)	- pos on unitsphere
                    phi     = ((2.*np.pi)*rndnum_arr[tc,3])			    #P flat in phi			- pos on unitsphere
                    S2_rnd_unit = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
                    #BH3:
                    theta   = np.arccos(2.*rndnum_arr[tc,4]-1.)         #P flat in cos(theta)	- pos on unitsphere
                    phi     = ((2.*np.pi)*rndnum_arr[tc,5])			    #P flat in phi			- pos on unitsphere
                    S3_rnd_unit = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
                    #group together:
                    S123_rnd_unit   = [S1_rnd_unit, S2_rnd_unit, S3_rnd_unit]
                    
                    #Primordial bin NO:
                    if (primorbin_yesno == 0):
                        #spin vecs for BH 123:
                        S1_unit = S123_rnd_unit[0]
                        S2_unit = S123_rnd_unit[1]
                        S3_unit = S123_rnd_unit[2]
                        S1_a    = 0.5
                        S2_a    = 0.5
                        S3_a    = 0.5
                        S123_unit_vec   = np.array([S1_unit,    S2_unit,    S3_unit])
                        S123_a          = np.array([S1_a,       S2_a,       S3_a])  
                        #spin vecs for BH ij:
                        #BH i:
                        Si_unit = S123_unit_vec[out_bin_i-1]    
                        Si_a    = S123_a[out_bin_i-1]    
                        #BH j:
                        Sj_unit = S123_unit_vec[out_bin_j-1]    
                        Sj_a    = S123_a[out_bin_j-1]    
                    
                    #Primordial bin YES:
                    if (primorbin_yesno == 1):
                        #spin vecs for BH 123:
                        S1_unit = np.array([0.0, 0.0, 1.0])
                        S2_unit = np.array([0.0, 0.0, 1.0])
                        S3_unit = S123_rnd_unit[2]
                        S1_a    = 0.5
                        S2_a    = 0.5
                        S3_a    = 0.5
                        S123_unit_vec   = np.array([S1_unit,    S2_unit,    S3_unit])
                        S123_a          = np.array([S1_a,       S2_a,       S3_a])  
                        #spin vecs for BH ij:
                        #BH i:
                        Si_unit = S123_unit_vec[out_bin_i-1]    
                        Si_a    = S123_a[out_bin_i-1]    
                        #BH j:
                        Sj_unit = S123_unit_vec[out_bin_j-1]    
                        Sj_a    = S123_a[out_bin_j-1]    
                    
                    #calc of chi effective:                    
                    chi_eff = np.dot(((1./(mbi+mbj))*(Si_a*mbi*Si_unit + Sj_a*mbj*Sj_unit)), Lorb_ij_unitvec)                
                #-------------------
                #chirp mass:
                #-------------------
                M_chirp = ((mbi*mbj)**(3./5.))/((mbi + mbj)**(1./5.))
                #-------------------
                #t_insp of second merger:
                #-------------------
                #define:
                m_m     = mbi_SI+mbj_SI
                m_s     = mbk_SI
                m_sm    = m_s+m_m
                mu_sm   = m_s*m_m/m_sm
                Rs      = (2.*G_new_SI*m_sm/(c_SI**2.))
                Ms      = mu_sm
                eps     = 85.*np.pi/96.
                beta    = 7./2.
                #calc:
                t_SM_yrs                    = 10e20 #initialize
                if (endstate_id == 2 or endstate_id == 5):  #endstate must be coll(2) or GWinsp(5)
                    if (Etot_ijk < 0.0):                    #bin-sin must be bound    
                        #before kick (no kick NK):
                        NK_rel_pos_vec_SI   = pos_CMijk_xyz*R_sun_SI
                        NK_rel_vel_vec_SI   = vel_CMijk_xyz*(1000./kmsec_U)
                        #kick velocity vec...INCLUDE LATER:       
                        vel_vec_kick_SI     = np.array([0.0, 0.0, 0.0])
                        #after kick:
                        rel_pos_vec_SI = NK_rel_pos_vec_SI
                        rel_vel_vec_SI = NK_rel_vel_vec_SI + vel_vec_kick_SI
                        #calc orb params:
                        rel_L_vec_SI    = mu_sm*np.cross(rel_pos_vec_SI, rel_vel_vec_SI)
                        rel_r_SI        = np.sqrt(sum(rel_pos_vec_SI**2.))
                        rel_v_SI        = np.sqrt(sum(rel_vel_vec_SI**2.))
                        rel_L_SI        = np.sqrt(sum(rel_L_vec_SI**2.))
                        asm_SI  = 1./((2./rel_r_SI) - (rel_v_SI**2.)/(G_new_SI*m_sm))
                        esm_SI  = np.sqrt(1. - ((rel_L_SI**2.)/(G_new_SI*m_sm*asm_SI*(mu_sm**2.))))
                        rpsm_SI = asm_SI*(1.-esm_SI)
                        #check if bin is still bound:
                        if (asm_SI > 0.0):
                            t_SM_yrs        = (2.*np.pi*(rpsm_SI**(beta))*np.sqrt(asm_SI)*(m_sm*mu_sm*(Rs**(1.-beta))/(eps*(Ms**(2.))*np.sqrt(G_new_SI*m_sm))))/yr_sec
                #-------------------
            
    
                #-------------------
                #save calc info data:
                #-------------------
                save_calc_arr_INT[tc,0:7]     = [endstate_id, out_bin_i, out_bin_j, out_bin_k, primorbin_yesno, bin_esc_yesno, ier]       
                save_calc_arr_REAL[tc,0:24]   = [mbi, mbj, mbk, tint_yr, tinsp_yr, vc0_kms, ecc_fGW, v_kick_kms, theta_Lorb_zvec, P_rp_LT_rec, a_0, e_0, a_ij, e_ij, chi_eff, M_chirp, P_rp_LT_rec_finbin, vesc_kms, a_form, e_form, m0_cluster, rho0_cluster, t_SM_yrs, binfrac_cluster]
                #-------------------
                                
            
                #-------------------
                #print info:
                #-------------------
                print tc+1, endstate_id, ier, P_rp_LT_rec, v_kick_kms, chi_eff, vc0_kms, ecc_fGW, vesc_kms, t_SM_yrs, binfrac_cluster
                #-------------------
            
                #-------------------
                #update total counter:
                #-------------------
                tc = tc+1
                #-------------------
        #---------------------------        
         
    #----------------------------------------------------------
    #save calc info data:
    #----------------------------------------------------------
    tf = open(output_data_folder+simdata_name+'save_calc_arr_INT.txt', "w")
    np.savetxt(tf, save_calc_arr_INT,   fmt='%5i')
    tf.close()

    tf = open(output_data_folder+simdata_name+'save_calc_arr_REAL.txt', "w")
    np.savetxt(tf, save_calc_arr_REAL,   fmt='%5f')
    tf.close()
    #----------------------------------------------------------
    #----------------------------------------------------------

#--------------------------------------------------------------
#--------------------------------------------------------------










#--------------------------------------------------------------
#read calc info data:
#--------------------------------------------------------------
tf = open(output_data_folder+simdata_name+'save_calc_arr_INT.txt', "r")
save_calc_arr_INT   = np.loadtxt(tf, dtype=int)         #[(0)endstate_id, (1)out_bin_i, (2)out_bin_j, (3)out_bin_k, (4)primorbin_yesno, (5)bin_esc_yesno, (6)ier]  
tf.close()

tf = open(output_data_folder+simdata_name+'save_calc_arr_REAL.txt', "r")
save_calc_arr_REAL  = np.loadtxt(tf, dtype=float)       #[(0)mbi, (1)mbj, (2)mbk, (3)tint_yr, (4)tinsp_yr, (5)vc0_kms, (6)ecc_fGW, (7)v_kick_kms, (8)theta_Lorb_zvec, (9)P_rp_LT_rec, (10)a_0, (11)e_0, (12)a_ij, (13)e_ij, (14)chi_eff, (15)M_chirp, (16)P_rp_LT_rec_finbin, (17)vesc_kms, (18)a_form, (19)e_form, (20)m0_cluster, (21)rho0_cluster, (22)t_SM_yrs, (23)binfrac_cluster]
tf.close()
#--------------------------------------------------------------
#--------------------------------------------------------------


#--------------------------------------------------------------
#analyze data:
#--------------------------------------------------------------
#define:
mij_all         = save_calc_arr_REAL[:,0] + save_calc_arr_REAL[:,1]
tint_yrs_all    = save_calc_arr_REAL[:,3]
tinsp_yrs_all   = save_calc_arr_REAL[:,4]
P_ecc_cap_all   = save_calc_arr_REAL[:,9]
P_rp_LT_rec_finbin_all = save_calc_arr_REAL[:,16]
chi_eff_all     = save_calc_arr_REAL[:,14]
Mchirp_all      = save_calc_arr_REAL[:,15]
ecc_fGW_all     = save_calc_arr_REAL[:,6]
vc0_kms_all     = save_calc_arr_REAL[:,5]
t_SM_yrs_all    = save_calc_arr_REAL[:,22]
binfrac_all     = save_calc_arr_REAL[:,23]
theta_Lorb_zvec_all = save_calc_arr_REAL[:,8]


m0_cluster_all      = save_calc_arr_REAL[:,20]
rho0_cluster_all    = save_calc_arr_REAL[:,21]

a_atform_all    = save_calc_arr_REAL[:,18]
e_atform_all    = save_calc_arr_REAL[:,19]
rp_atform_all   = a_atform_all*(1.0 - e_atform_all)
fGW_atform_all  = (1./np.pi)*np.sqrt(G_new_SI*(mij_all*M_sun_SI)/((rp_atform_all*R_sun_SI)**3.))

tmerge_yrs_all  = tint_yrs_all + tinsp_yrs_all

#pos:
pos_id3             = np.where(save_calc_arr_INT[:,0] == 3)[0]
pos_id5             = np.where(save_calc_arr_INT[:,0] == 5)[0]
pos_id10            = np.where(save_calc_arr_INT[:,0] == 10)[0]
pos_bin12           = np.where(save_calc_arr_INT[:,3] == 3)[0]
pos_bin13           = np.where(save_calc_arr_INT[:,3] == 2)[0]
pos_bin23           = np.where(save_calc_arr_INT[:,3] == 1)[0]
pos_pmbin           = np.where(save_calc_arr_INT[:,4] == 1)[0]
pos_escGC           = np.where(save_calc_arr_INT[:,5] == 1)[0]

pos_ecc             = np.where(ecc_fGW_all[:] > 0.1)[0]

pos_pmbinYES        = np.where(save_calc_arr_INT[:,4] == 1)[0]
pos_pmbinNO         = np.where(save_calc_arr_INT[:,4] == 0)[0]

print len(pos_pmbinYES)
print len(pos_pmbinNO)

pos_MchirpW1l       = np.where(Mchirp_all[:] > 10.0)[0]
pos_MchirpW1u       = np.where(Mchirp_all[:] < 20.0)[0]
pos_MchirpW1        = list(set(pos_MchirpW1l).intersection(pos_MchirpW1u))

pos_mij             = np.where(mij_all[:] < 60.0)[0]

pos_tmergeLTtH      = np.where(tmerge_yrs_all[:]    < 10.**(10.))[0]
pos_tintLTtH        = np.where(tint_yrs_all[:]      < 10.**(10.))[0]

pos_id3_escGC       = list(set(pos_id3).intersection(pos_escGC))

pos_binfrac_GTv1    = np.where(binfrac_all[:] >= 0.8)[0]

t1 = 1.*(10.**(9.))
t2 = 1.*(10.**(10.))
pos_tmergeGTt1      = np.where(tmerge_yrs_all[:] > t1)[0]
pos_tmergeLTt2      = np.where(tmerge_yrs_all[:] < t2)[0]
pos_tmerge12        = list(set(pos_tmergeGTt1).intersection(pos_tmergeLTt2))
pos_tintGTt1        = np.where(tint_yrs_all[:] > t1)[0]
pos_tintLTt2        = np.where(tint_yrs_all[:] < t2)[0]
pos_tint12          = list(set(pos_tintGTt1).intersection(pos_tintLTt2))




pos_id3_escGC_tmergeLTtH            = list(set(pos_id3_escGC).intersection(pos_tmergeLTtH))
pos_id5_tintLTtH                    = list(set(pos_id5).intersection(pos_tintLTtH))

pos_id3_escGC_tmerge12              = list(set(pos_id3_escGC).intersection(pos_tmerge12))
pos_id5_tint12                      = list(set(pos_id5).intersection(pos_tint12))





pos_id3_escGC_tmergeLTtH_pmbinYES       = list(set(pos_id3_escGC_tmergeLTtH).intersection(pos_pmbinYES))
pos_id3_escGC_tmergeLTtH_pmbinNO        = list(set(pos_id3_escGC_tmergeLTtH).intersection(pos_pmbinNO))
pos_id3_escGC_tmergeLTtH_pmbinYES_bin12 = list(set(pos_id3_escGC_tmergeLTtH_pmbinYES).intersection(pos_bin12))

pos_id3_escGC_tmergeLTtH_MchirpW1           = list(set(pos_id3_escGC_tmergeLTtH).intersection(pos_MchirpW1))
pos_id3_escGC_tmergeLTtH_pmbinYES_MchirpW1  = list(set(pos_id3_escGC_tmergeLTtH_pmbinYES).intersection(pos_MchirpW1))
pos_id3_escGC_tmergeLTtH_pmbinNO_MchirpW1   = list(set(pos_id3_escGC_tmergeLTtH_pmbinNO).intersection(pos_MchirpW1))


pos_id3_escGC_bin12             = list(set(pos_id3_escGC).intersection(pos_bin12))

pos_id3_escGC_pmbin             = list(set(pos_id3_escGC).intersection(pos_pmbin))

pos_id3_escGC_pmbin_bin12       = list(set(pos_id3_escGC_pmbin).intersection(pos_bin12))

pos_id3_escGC_mij               = list(set(pos_id3_escGC).intersection(pos_mij))
#--------------------------------------------------------------

#completeness fraction:
print 'pos10/all', (1.*len(pos_id10))/(1.*nr_tot_sims), 1.*nr_tot_sims, 1.*nr_tot_sims/5.












fig, ax1 = plt.subplots(figsize=(8, 6))

pos = list(set(pos_id3_escGC_tmergeLTtH).intersection(pos_binfrac_GTv1))
plotarr = theta_Lorb_zvec_all[pos]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [0.0, 4.0], normed=False, linewidth=3.0, histtype='step', color = 'red')

pos = list(set(pos_id3_escGC_tmergeLTtH_pmbinYES).intersection(pos_binfrac_GTv1))
plotarr = theta_Lorb_zvec_all[pos]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [0.0, 4.0], normed=False, linewidth=3.0, histtype='step', color = 'blue')

pos = list(set(pos_id3_escGC_tmergeLTtH_pmbinNO).intersection(pos_binfrac_GTv1))
plotarr = theta_Lorb_zvec_all[pos]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [0.0, 4.0], normed=False, linewidth=3.0, histtype='step', color = 'orange')

pos = list(set(pos_id3_escGC_tmergeLTtH_pmbinYES_bin12).intersection(pos_binfrac_GTv1))
plotarr = theta_Lorb_zvec_all[pos]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [0.0, 4.0], normed=False, linewidth=3.0, histtype='step', color = 'black')

plt.show()



exit()












#--------------------------------------------------------------
#LISA+LIGO figure:
#--------------------------------------------------------------

fig, ax = plt.subplots(figsize=(6, 5))

#LISA+LIGO regions:
#LISA:
xl  = 1e-3
xu  = 1e-1
rect    = mpl.patches.Rectangle((xl,1e0), (xu-xl), 1e10, color='lightblue')
ax.add_patch(rect)
ax.text(1e-2, 6, 'LISA', color='black', fontsize = 12, horizontalalignment='center')
#LIGO:
xl  = 1e1
xu  = 1e3
rect    = mpl.patches.Rectangle((xl,1e0), (xu-xl), 1e10, color='lightgrey')
ax.add_patch(rect)
ax.text(1e2, 6, 'LIGO', color='black', fontsize = 12, horizontalalignment='center')

#explanations:
ax.text(1e-1, 40, r'BBH evolution $\rightarrow$ $\rightarrow$ $\rightarrow$', color='black', fontsize = 10, horizontalalignment='left')


#plot:
#GW captures:
#pp  = pos_id5_tintLTtH
pp  = pos_id5_tint12
xp  = fGW_atform_all[pp]
yp  = Mchirp_all[pp]
plt.plot(xp[:], yp[:], marker='^', linestyle = '', fillstyle='full', markeredgewidth=0.0, markersize=5, alpha=0.75, c='red',    label='binary-single GW mergers')
fGW_p5 = xp
#Ejections:
#pp  = pos_id3_escGC_tmergeLTtH
pp  = pos_id3_escGC_tmerge12
xp  = fGW_atform_all[pp]
yp  = Mchirp_all[pp]
plt.plot(xp[:], yp[:], marker='o', linestyle = '', fillstyle='full', markeredgewidth=0.0, markersize=5, alpha=0.75, c='black',  label='ejected GW mergers')
fGW_p3 = xp

#axis settings:
ax.set_xlim([1e-5,  1e4])
ax.set_ylim([5,75])

ax.set_xlabel(r'GW frequency at formation $f_{fm}$ [Hz]')
ax.set_ylabel(r'Chirp mass $\mathscr{M}_{c}$ [$M_{\odot}$]')
ax.set_title(r'BBH mergers in LISA & LIGO')

ax.set_xscale('log')
ax.set_yscale('log')

ax.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = True)

plt.savefig(icdata_name+'LISALIGO.pdf', bbox_inches='tight')     
plt.show()

#print info:
print len(pos_id3_escGC_tmergeLTtH), len(pos_id5_tintLTtH)

fGW_p53 = np.array(list(fGW_p5)+list(fGW_p3))
print 1.0*len(np.where(fGW_p3[:] > 0.1)[0])/(1.0*len(fGW_p3[:]))
print 1.0*len(np.where(fGW_p53[:] > 0.1)[0])/(1.0*len(fGW_p53[:]))

print 1.0*len(np.where(fGW_p53[:] > 0.1)[0])
print 1.0*len(fGW_p53[:])


exit()
#--------------------------------------------------------------













print t_SM_yrs_all[pos_id5]
#TEST:
fig, ax1 = plt.subplots(figsize=(6, 3))
input_arr       = np.log10(t_SM_yrs_all[pos_id5])
print input_arr
ax1.hist(input_arr, bins=20, range = [0, 20], linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
plt.show()

exit()







#--------------------------------------------------------------
#FIGURE 1:
#--------------------------------------------------------------


#------------------------------
#settings/calc:
#------------------------------
tint_MEGAyrs_all    = tint_yrs_all/(10.**(6.))
tinsp_MEGAyrs_all   = tinsp_yrs_all/(10.**(6.))
tmerge_MEGAyrs_all  = tint_MEGAyrs_all + tinsp_MEGAyrs_all


#lin bins:
#bins_arr    = np.arange(0.0, 10.**(4.), 100.) 
#log bins:
nr_bins     = 100.
logmin      = -2.0
logmax      = 4.0
log_deltarange  = logmax - logmin
log_deltax      = log_deltarange/nr_bins
bins_arr            = 10.**(np.arange(logmin, logmax+log_deltax, log_deltax))
bins_centers_arr    = 10.**(np.log10(bins_arr[0:nr_bins]) + 0.5*log_deltax)


#PLOT SETTINGS:
f   = plt.figure(figsize=(5,8))
gs  = gridspec.GridSpec(2, 1,height_ratios=[2,1], hspace=0.0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])


xmin = 7.*10.**(1.0)
xmax = 1.*10.**(4.0)
#------------------------------


#------------------------------
#FIGURE ax1:
#------------------------------
#plot:

plotarr     = tmerge_MEGAyrs_all[pos_id3_escGC]
input_arr   = plotarr
hist_vals_setA, bins_info_arr, patches =  ax1.hist(input_arr, bins=bins_arr, normed=False, linewidth=3.0, histtype='step', color= 'black',      alpha = 1.0, label = r'escaper mergers (EM)')

plotarr     = tint_MEGAyrs_all[pos_id5]
input_arr   = plotarr
hist_vals_setB, bins_info_arr, patches =  ax1.hist(input_arr, bins=bins_arr, normed=False, linewidth=3.0, histtype='step', color= 'crimson',    alpha = 1.0, label = r'$3$-body capture mergers (CM)')

plotarr     = tint_MEGAyrs_all[pos_ecc]
input_arr   = plotarr
hist_vals_setC, bins_info_arr, patches =  ax1.hist(input_arr, bins=bins_arr, normed=False, linewidth=3.0, histtype='step', color= 'steelblue',  alpha = 1.0, label = r'$e>0.1$, $10$ Hz (eCM)')

plotarr     = tint_MEGAyrs_all[:]
weight_arr  = P_ecc_cap_all[:]
input_arr   = plotarr
hist_vals_setD =  np.histogram(input_arr, bins=bins_arr, weights = weight_arr, normed=False)[0]
ax1.plot(bins_centers_arr, hist_vals_setD, linestyle='--', linewidth=1.5, color = 'steelblue', alpha = 1.0, label = r'$e>0.1$, $10$ Hz (eCM, analytical)')

plotarr     = tint_MEGAyrs_all[pos_id3_escGC]
weight_arr  = P_rp_LT_rec_finbin_all[pos_id3_escGC]
input_arr   = plotarr
hist_vals_setE =  np.histogram(input_arr, bins=bins_arr, weights = weight_arr, normed=False)[0]
ax1.plot(bins_centers_arr, hist_vals_setE, linestyle='--', linewidth=1.5, color = 'orange', alpha = 1.0, label = r'$e>0.1$, $10$ Hz (eEM, analytical)')

#plotarr     = tint_MEGAyrs_all[pos_id3_escGC]
#input_arr   = plotarr
#hist_vals_setF, bins_info_arr, patches =  ax1.hist(input_arr, bins=bins_arr, normed=False, linewidth=3.0, histtype='step', color= 'orange',  alpha = 1.0, label = r'test..')

#plotarr     = tint_MEGAyrs_all[:]
#input_arr   = plotarr
#hist_vals_setG, bins_info_arr, patches =  ax1.hist(input_arr, bins=bins_arr, normed=False, linewidth=3.0, histtype='step', color= 'purple',  alpha = 1.0, label = r'test..')


#axis settings:
ax1.set_xlim([xmin,  xmax])
ax1.set_ylim([1e-2,1e4])

ax1.set_ylabel(r'number events $N$')

ax1.set_xscale('log')
ax1.set_yscale('log')

ax1.set_title('Black Hole Mergers')

ax1.legend(loc='best', numpoints = 1, fontsize = 10.0, frameon = False)

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
#------------------------------


#------------------------------
#FIGURE ax2:
#------------------------------
#plot:

ax2.plot([1e-10,1e10], [1.0,  1.0],     linestyle=':', linewidth=1.5, color = 'grey')
ax2.plot([1e-10,1e10], [0.1,  0.1],     linestyle=':', linewidth=1.5, color = 'grey')
ax2.plot([1e-10,1e10], [0.01, 0.01],    linestyle=':', linewidth=1.5, color = 'grey')

input_arr   = (1.*hist_vals_setB)/(1.*hist_vals_setA)
ax2.plot(bins_centers_arr, input_arr,     marker='o', linestyle=':', markersize=5.0, linewidth=1.5, markeredgewidth=0, color = 'crimson',  label = r'$N_{\rm CM}/N_{\rm EM}$')

input_arr   = (1.*hist_vals_setC)/(1.*hist_vals_setA)
ax2.plot(bins_centers_arr, input_arr,     marker='^', linestyle=':', markersize=5.0, linewidth=1.5, markeredgewidth=0, color = 'steelblue', label = r'$N_{\rm eCM}/N_{\rm EM}$')

input_arr   = (1.*hist_vals_setE)/(1.*hist_vals_setA)
ax2.plot(bins_centers_arr, input_arr, linestyle='--', linewidth=1.5, color = 'orange', label = r'$N_{\rm eEM}/N_{\rm EM}$')

#axis settings:
ax2.set_xlim([xmin,  xmax])
ax2.set_ylim([5e-5, 9e0])

ax2.set_xlabel(r'Time $t$ [Myr]')
ax2.set_ylabel(r'$N/N_{\rm EM}$')

ax2.set_xscale('log')
ax2.set_yscale('log')

ax2.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False, ncol=3)
#------------------------------


#------------------------------
#save and show:
#------------------------------
plt.savefig(icdata_name+'HIST_events.pdf', bbox_inches='tight')     
plt.show()
#------------------------------

exit()
#--------------------------------------------------------------











#TEST:

fig, ax1 = plt.subplots(figsize=(6, 3))

input_arr       = Mchirp_all[pos_id3_escGC_tmerge12]
ax1.hist(input_arr, bins=20, range = [0.0, 40.], normed=True, linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
input_arr       = Mchirp_all[pos_id5_tint12]
ax1.hist(input_arr, bins=20, range = [0.0, 40.], normed=True, linewidth=2.0, histtype='step', linestyle = '-', color = 'red')

input_arr       = Mchirp_all[pos_id3_escGC_tmergeLTtH]
ax1.hist(input_arr, bins=20, range = [0.0, 40.], normed=True, linewidth=2.0, histtype='step', linestyle = ':', color = 'black')
input_arr       = Mchirp_all[pos_id5_tintLTtH]
ax1.hist(input_arr, bins=20, range = [0.0, 40.], normed=True, linewidth=2.0, histtype='step', linestyle = ':', color = 'red')

ax1.set_yscale('log')

#------------------------------
#save and show:
#------------------------------
plt.savefig('test_mchirp.pdf', bbox_inches='tight')     
plt.show()

exit()




#TEST:
fig, ax1 = plt.subplots(figsize=(6, 5))
input_arr       = np.log10(m0_cluster_all[:])
ax1.hist(input_arr, bins=100)#, range = [0.0, 200.], normed=True, linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
#ax1.set_xscale('log')
ax1.set_yscale('log')
#------------------------------
#save and show:
#------------------------------
plt.show()
#exit()



#TEST:
fig, ax1 = plt.subplots(figsize=(6, 5))
input_arr       = np.log10(rho0_cluster_all[:])
ax1.hist(input_arr, bins=100)#, range = [0.0, 200.], normed=True, linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
#ax1.set_xscale('log')
ax1.set_yscale('log')
#------------------------------
#save and show:
#------------------------------
plt.show()
#exit()


exit()








#TEST:

fig, ax1 = plt.subplots(figsize=(6, 5))

#input_arr       = Mchirp_all[pos_id3_escGC_tmergeLTtH]
input_arr       = Mchirp_all[pos_id3_escGC_tmerge12]
ax1.hist(input_arr, bins=500, range = [0.0, 200.], normed=True, linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
#input_arr       = Mchirp_all[pos_id5]
#ax1.hist(input_arr, bins=500, range = [0.0, 200.], normed=True, linewidth=2.0, histtype='step', linestyle = '-', color = 'red')
ax1.set_yscale('log')

#------------------------------
#save and show:
#------------------------------
plt.savefig('TEST.pdf', bbox_inches='tight')     
plt.show()

exit()











#--------------------------------------------------------------
#FIGURE 3:
#--------------------------------------------------------------

#PROBLEMS WITH THE ECC CALCULAITON: MANY 0.00000000e+00 and MANY 1e-6: IS NOT WORKING FOR MANY VALUES!!!

print np.log10(ecc_fGW_all[pos_id5[0:1000]])
print (ecc_fGW_all[pos_id5[0:1000]])

fig, ax1 = plt.subplots(figsize=(6, 5))
ecc_fGW_arr     = ecc_fGW_all[pos_id3_escGC_tmergeLTtH]
input_arr       = np.log10(ecc_fGW_arr)
ax1.hist(input_arr, bins=200, range = [-8, 0.1], normed=False, linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
ecc_fGW_arr     = ecc_fGW_all[pos_id5]
input_arr       = np.log10(ecc_fGW_arr)
ax1.hist(input_arr, bins=200, range = [-8, 0.1], normed=False, linewidth=2.0, histtype='step', linestyle = '-', color = 'red')

print save_calc_arr_INT[pos_id5,6]
print len(pos_id5)

#------------------------------
#save and show:
#------------------------------
plt.savefig('TEST.pdf', bbox_inches='tight')     
plt.show()

exit()
#--------------------------------------------------------------











#--------------------------------------------------------------
#FIGURE 2:
#--------------------------------------------------------------
#THERE ARE TO MANY LINES: JUST A MESS CANT USE IT!
fig = plt.figure(figsize=(6, 5))


posplot = pos_id3_escGC
fig.add_subplot(111).plot(np.log10(e_atform_all[posplot]), np.log10(fGW_atform_all[posplot]), markeredgewidth=0, marker='o', linewidth = 0, markersize=5.0, alpha=0.5, color='black')
#draw f(e) evolution lines:
nrp = len(posplot)
for nc in range(0, nrp):
    #get pos:
    pid = posplot[nc]
    #define/calc:
    ain     = a_atform_all[pid]*R_sun_SI    #SI
    ein     = e_atform_all[pid]             #SI
    Mtot    = mij_all[pid]*M_sun_SI         #SI
    c0      = ain*(((ein**(12./19.)/(1.-ein**2.))*((1.+(121./304.)*(ein**2.))**(870./2299.)))**(-1.0))
    eparr   = 10.**(np.arange(-5.0, np.log10(ein + 1e-5), 1e-1))
    Fearr   = (eparr**(12./19.)/(1.+eparr))*((1.+(121./304.)*(eparr**2.))**(870./2299.))
    fparr   = (1./np.pi)*((G_new_SI*Mtot/((c0**3.)*(Fearr**3.)))**(1./2.))
    #plot:
    fig.add_subplot(111).plot(np.log10(eparr), np.log10(fparr), linewidth = 1, linestyle = '-', alpha=0.5, color='black')
    

posplot = pos_id5
fig.add_subplot(111).plot(np.log10(e_atform_all[posplot]), np.log10(fGW_atform_all[posplot]), markeredgewidth=0, marker='o', linewidth = 0, markersize=5.0, alpha=0.5, color='red')
#draw f(e) evolution lines:
nrp = len(posplot)
for nc in range(0, nrp):
    #get pos:
    pid = posplot[nc]
    #define/calc:
    ain     = a_atform_all[pid]*R_sun_SI    #SI
    ein     = e_atform_all[pid]             #SI
    Mtot    = mij_all[pid]*M_sun_SI         #SI
    c0      = ain*(((ein**(12./19.)/(1.-ein**2.))*((1.+(121./304.)*(ein**2.))**(870./2299.)))**(-1.0))
    eparr   = 10.**(np.arange(-5.0, np.log10(ein + 1e-5), 1e-1))
    Fearr   = (eparr**(12./19.)/(1.+eparr))*((1.+(121./304.)*(eparr**2.))**(870./2299.))
    fparr   = (1./np.pi)*((G_new_SI*Mtot/((c0**3.)*(Fearr**3.)))**(1./2.))
    #plot:
    fig.add_subplot(111).plot(np.log10(eparr), np.log10(fparr), linewidth = 1, linestyle = '-', alpha=0.5, color='red')


plt.show()

exit()

#--------------------------------------------------------------







#plot:

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

tint_MEGAyrs_all    = tint_yrs_all/(10.**(6.))
tinsp_MEGAyrs_all   = tinsp_yrs_all/(10.**(6.))
tmerge_MEGAyrs_all  = tint_MEGAyrs_all + tinsp_MEGAyrs_all

fig, ax1 = plt.subplots(figsize=(10, 6))

bins_arr    = np.arange(0.0, 10.**(4.) + 100., 100.) 

plotarr     = tint_MEGAyrs_all[:]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', color= 'black', alpha = 0.7, label = 'all int')

plotarr     = tint_MEGAyrs_all[pos_id3_escGC]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', color= 'blue', alpha = 0.7, label = 'eject bin int')



plotarr     = tmerge_MEGAyrs_all[pos_id3_escGC]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', color= 'green', alpha = 0.7, label = 'eject bin mergers')

plotarr     = tint_MEGAyrs_all[pos_id5]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', color= 'purple', alpha = 0.7, label = 'GW capture mergers')



plotarr     = tint_MEGAyrs_all[pos_ecc]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', color= 'red', alpha = 0.7, label = 'eccentric mergers (EM)')

plotarr     = tint_MEGAyrs_all[:]
weight_arr  = P_ecc_cap_all[:]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, weights = weight_arr,  linewidth=3.0, histtype='step', color= 'red', alpha = 0.7, linestyle = '--', label = 'EM analytical (all)')

plotarr     = tint_MEGAyrs_all[pos_id3_escGC]
weight_arr  = P_rp_LT_rec_finbin_all[pos_id3_escGC]
input_arr   = plotarr
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, weights = weight_arr,  linewidth=3.0, histtype='step', color= 'red', alpha = 0.7, linestyle = ':', label = 'EM analytical (eject bin)')


#plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_pmbinYES_bin12]
#input_arr   = plotarr
#ax1.hist(input_arr, bins=50, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')



ax1.legend(loc='lower left', numpoints = 1, fontsize = 10.0, frameon = False)

plt.yscale('log')
plt.xscale('log')

ax1.set_xlabel(r'Time [Myr]')
ax1.set_ylabel(r'Rate')
ax1.set_xlim([1e2,  1e4])
#ax1.set_ylim([min(y_input),  max(y_input)])

#save and show:
plt.savefig(icdata_name+'TEST1.pdf', bbox_inches='tight')     
plt.show()




exit()












#PROBLEMS WITH THE ECC CALCULAITON: MANY 0.00000000e+00 and MANY 1e-6: IS NOT WORKING FOR MANY VALUES!!!

print np.log10(ecc_fGW_all[pos_id5[0:1000]])
print (ecc_fGW_all[pos_id5[0:1000]])
fig, ax1 = plt.subplots(figsize=(6, 5))
ecc_fGW_arr     = ecc_fGW_all[pos_id3_escGC_tmergeLTtH]
input_arr       = np.log10(ecc_fGW_arr)
ax1.hist(input_arr, bins=100, range = [-8, 0.1], normed=False, linewidth=2.0, histtype='step', linestyle = '-', color = 'black')
ecc_fGW_arr     = ecc_fGW_all[pos_id5]
input_arr       = np.log10(ecc_fGW_arr)
ax1.hist(input_arr, bins=100, range = [-8, 0.1], normed=False, linewidth=2.0, histtype='step', linestyle = '-', color = 'red')
plt.show()



exit()


plotarr = tmerge_MEGAyrs_all[pos_id3_escGC_pmbin]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = '2')

plotarr = tmerge_MEGAyrs_all[pos_id3_escGC_pmbin_bin12]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = '3')

plotarr = tmerge_MEGAyrs_all[pos_id3_escGC_mij]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = 'A')

plotarr = tmerge_MEGAyrs_all[pos_id3_escGC_bin12]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = 'B')


plotarr     = tint_MEGAyrs_all[:]
weight_arr  = P_ecc_cap_all[:]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, weights = weight_arr,  linewidth=3.0, histtype='step', label = 'C')

plotarr     = tint_MEGAyrs_all[pos_id3_escGC]
weight_arr  = P_ecc_cap_all[pos_id3_escGC]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, weights = weight_arr,  linewidth=3.0, histtype='step', label = 'D')



plotarr = tint_MEGAyrs_all[pos_id3_escGC_pmbin]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = '4')

plotarr = tint_MEGAyrs_all[pos_id3_escGC]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = '5')

plotarr = tint_MEGAyrs_all[:]
input_arr   = plotarr
bins_arr    = np.arange(0.0, 10.**(6.) + 100., 100.) 
ax1.hist(input_arr, bins=bins_arr, range = [0.0, 10.**(6.)], normed=False, linewidth=3.0, histtype='step', label = '6')


ax1.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)

plt.yscale('log')
plt.xscale('log')



plt.show()







fig = plt.figure(figsize=(6, 5))

v_kick_kms_arr      = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH,7]
a_0_AU_arr          = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH,10]/AU_U
a_ij_AU_arr         = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH,12]/AU_U
fig.add_subplot(111).plot(v_kick_kms_arr[:], np.log10(a_0_AU_arr[:]),       markeredgewidth=0, marker='o', linewidth = 0, markersize=3.0, alpha=0.5, color='black')
#fig.add_subplot(111).plot(v_kick_kms_arr[:], np.log10(a_ij_AU_arr[:]),      markeredgewidth=0, marker='o', linewidth = 0, markersize=3.0, alpha=0.5, color='red')

v_kick_kms_arr      = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH_pmbinYES,7]
a_0_AU_arr          = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH_pmbinYES,10]/AU_U
a_ij_AU_arr         = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH_pmbinYES,12]/AU_U
fig.add_subplot(111).plot(v_kick_kms_arr[:], np.log10(a_0_AU_arr[:]),       markeredgewidth=0, marker='o', linewidth = 0, markersize=3.0, alpha=0.5, color='yellow')
#fig.add_subplot(111).plot(v_kick_kms_arr[:], np.log10(a_ij_AU_arr[:]),      markeredgewidth=0, marker='o', linewidth = 0, markersize=3.0, alpha=0.5, color='green')

plt.show()







fig = plt.figure(figsize=(6, 5))

chirp_mass_arr      = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH,15]
a_0_AU_arr          = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH,10]/AU_U
fig.add_subplot(111).plot(chirp_mass_arr[:], np.log10(a_0_AU_arr[:]),       markeredgewidth=0, marker='o', linewidth = 0, markersize=3.0, alpha=0.5, color='black')

chirp_mass_arr      = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH_pmbinYES,15]
a_0_AU_arr          = save_calc_arr_REAL[pos_id3_escGC_tmergeLTtH_pmbinYES,10]/AU_U
fig.add_subplot(111).plot(chirp_mass_arr[:], np.log10(a_0_AU_arr[:]),       markeredgewidth=0, marker='o', linewidth = 0, markersize=3.0, alpha=0.5, color='yellow')

plt.show()







fig, ax1 = plt.subplots(figsize=(8, 6))

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_pmbinYES]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_pmbinNO]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_pmbinYES_bin12]
input_arr   = plotarr
ax1.hist(input_arr, bins=50, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plt.show()





fig, ax1 = plt.subplots(figsize=(8, 6))

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_MchirpW1]
input_arr   = plotarr
ax1.hist(input_arr, bins=20, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_pmbinNO_MchirpW1]
input_arr   = plotarr
ax1.hist(input_arr, bins=20, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plotarr = chi_eff_all[pos_id3_escGC_tmergeLTtH_pmbinYES_MchirpW1]
input_arr   = plotarr
ax1.hist(input_arr, bins=20, range = [-1.0, 1.0], normed=False, linewidth=3.0, histtype='step')

plt.show()








exit()












#----------------------------------------------------------
#open data: SIM OUTPUT
#----------------------------------------------------------
tf = open(output_data_folder+simdata_name+'MC_settings_list_INT.txt', "r")
MC_settings_list_INT        = np.loadtxt(tf, dtype=int)         #[nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 0, 0]
tf.close()
tf = open(output_data_folder+simdata_name+'MC_settings_list_REAL.txt', "r")
MC_settings_list_REAL       = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(output_data_folder+simdata_name+'outlist_MC_info_INT.txt', "r")
outlist_MC_info_INT         = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(output_data_folder+simdata_name+'outlist_MC_info_REAL.txt', "r")
outlist_MC_info_REAL        = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(output_data_folder+simdata_name+'output_Nbody_endstate_INT.txt', "r")
output_Nbody_endstate_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(output_data_folder+simdata_name+'output_Nbody_endstate_REAL.txt', "r")
output_Nbody_endstate_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(output_data_folder+simdata_name+'output_Nbody_xtra_info_INT.txt', "r")
output_Nbody_xtra_info_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(output_data_folder+simdata_name+'output_Nbody_xtra_info_REAL.txt', "r")
output_Nbody_xtra_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(output_data_folder+simdata_name+'output_Nbody_xtra_2_info_REAL.txt', "r")
output_Nbody_xtra_2_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
#----------------------------------------------------------





















