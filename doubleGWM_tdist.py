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
from matplotlib import rcParams

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
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
m_parsec    = 3.086*(10**16.)   #m
#----------------------------------------------------------

mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'




analyze_data_YN = 0



#----------------------------------------------------------
#dataset:
#----------------------------------------------------------
#data name:
data_name = 'DMproj_WGR_202020Msun_a1e1_v001vc_fbfull_gamma00pi_fbres500'
#BH mass:
M1          = 20.
#GW kick:
vkick_kmsec = 0.0
#GW freq threshold:
fGW_limit   = 10.0  #in Hz: we calc ecc of both 1st and 2nd merg at this fGW limit.
#----------------------------------------------------------


#list of datanames:
#'DMproj_WGR_202020Msun_a1e2_v001vc_fbfull_gamma00pi_fbres500'#'DMproj_WGR_202020Msun_a1e2_v001vc_fbfull_gamma00pi_fbres500'#'test_NO3'#'NEWTOP_WGR_202020Msun_a1e2_v001vc_fbfull_gamma00pi_fbres500'
#or...
#DMproj_WGR_202020Msun_a1e2_v001vc_fbfull_gamma00pi_fbres500
#DMproj_WGR_202020Msun_a1e15_v001vc_fbfull_gamma00pi_fbres500
#DMproj_WGR_202020Msun_a1e1_v001vc_fbfull_gamma00pi_fbres500
#DMproj_WGR_202020Msun_a1e05_v001vc_fbfull_gamma00pi_fbres500
#DMproj_WGR_202020Msun_a1e0_v001vc_fbfull_gamma00pi_fbres500


#----------------------------------------------------------
#Data:
#----------------------------------------------------------
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
tf = open(data_folder+data_name+'output_Nbody_xtra_info_INT.txt', "r")
output_Nbody_xtra_info_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_info_REAL.txt', "r")
output_Nbody_xtra_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_2_info_REAL.txt', "r")
output_Nbody_xtra_2_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
#----------------------------------------------------------
#Define:
#----------------------------------------------------------
#MC_settings_list_INT #[nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 0, 0]
res_pa  = int(np.sqrt(MC_settings_list_INT[2]))
sma_a0  = outlist_MC_info_REAL[0,0]
nrtot   = len(outlist_MC_info_INT[:,0])
#----------------------------------------------------------


#----------------------------------------------------------
#----------------------------------------------------------
#outlist_MC_info_INT                #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
#outlist_MC_info_REAL               #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, R_val, b_val, f_val, g_val]
#output_Nbody_endstate_INT          #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
#output_Nbody_endstate_REAL         #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#output_Nbody_xtra_info_REAL        #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
#output_Nbody_xtra_2_info_REAL      #[out_pos_CMij_wrt_sink(x,y,z), out_vel_CMij_wrt_sink(x,y,z), Lvec_ij(x,y,z)]

#reshape to order: [b/g, f, 10]
#xy_outlist_MC_info_REAL         = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   
#xy_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(res_pa, res_pa, 10)             #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
#xy_output_Nbody_endstate_REAL   = output_Nbody_endstate_REAL.reshape(res_pa, res_pa, 10)            #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#xy_output_Nbody_xtra_info_REAL  = output_Nbody_xtra_info_REAL.reshape(res_pa, res_pa, 10)           #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
#----------------------------------------------------------



#--------------------------------------------------------------
if (analyze_data_YN == 1):
#--------------------------------------------------------------
    
    
    #----------------------------------------------------------
    #define:
    #----------------------------------------------------------
    ecc_1M_XHz_arr          = np.zeros(nrtot, dtype=np.float64)
    
    tGW_1M_yrs_arr          = np.zeros(nrtot, dtype=np.float64)

    tGW_2M_yrs_arr          = np.zeros(nrtot, dtype=np.float64)
    tGW_2M_yrs_arr[:]       = 1e20  #initialize

    ecc_2M_XHz_arr          = np.zeros(nrtot, dtype=np.float64)
    ecc_2M_XHz_arr[:]       = -1    #initialize

    fGW_2nd_bin_SI_arr      = np.zeros(nrtot, dtype=np.float64)
    fGW_2nd_bin_SI_arr[:]   = -1    #initialize

    rel_vel_NK_kmsec_arr    = np.zeros(nrtot, dtype=np.float64)
    rel_vel_NK_kmsec_arr[:] = -1    #initialize
    #----------------------------------------------------------


    #----------------------------------------------------------
    #calc properties of 1st GW merger:
    #----------------------------------------------------------
    #define:
    mSI_1   = M1*M_sun_SI
    mSI_2   = M1*M_sun_SI
    
    
    #loop over all f,b interactions:
    for c in range(0,nrtot):
        
        #endstate id:
        endstate_id     = output_Nbody_endstate_INT[c,0]
        
        #initialize
        ecc_Xhz         = -1
        tGW_1M_yrs      = 1e20
        
        
        #coll:
        if (endstate_id == 2):
            ecc_Xhz     = 1.0
            tGW_1M_yrs  = 0.0
            fGW_SI      = 1e20 
        #3b-merger
        if (endstate_id == 5):
            #initial a,e:
            ini_a_SI    = output_Nbody_xtra_info_REAL[c,3]*R_sun_SI
            ini_e_SI    = output_Nbody_xtra_info_REAL[c,4]
            rp_SI       = ini_a_SI*(1.-ini_e_SI)
            fGW_SI      = (1./np.pi)*np.sqrt(G_new_SI*(mSI_1 + mSI_2)/(rp_SI**3.))
            tGW_1M_yrs  = ((768./425.)*((ini_a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*mSI_1*mSI_2*(mSI_1+mSI_2)/(c_SI**5.)))*((1.-ini_e_SI**2.)**(7./2.)))/yr_sec                     
        #2b-merger
        if (endstate_id == 3):
            #initial a,e:
            ini_a_SI    = output_Nbody_endstate_REAL[c,3]*R_sun_SI
            ini_e_SI    = output_Nbody_endstate_REAL[c,4]
            rp_SI       = ini_a_SI*(1.-ini_e_SI)
            fGW_SI      = (1./np.pi)*np.sqrt(G_new_SI*(mSI_1 + mSI_2)/(rp_SI**3.))
            tGW_1M_yrs  = ((768./425.)*((ini_a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*mSI_1*mSI_2*(mSI_1+mSI_2)/(c_SI**5.)))*((1.-ini_e_SI**2.)**(7./2.)))/yr_sec                     
        
        #check initial fGW and evolv ecc:
        if (fGW_SI > fGW_limit):
            ecc_Xhz     = 1.0        
        if (fGW_SI < fGW_limit):
            #propagating ini a,e to a,e(fGW_limit)
            c0      = ini_a_SI/((ini_e_SI**(12./19.)/(1.-ini_e_SI**2.))*((1.+(121./304.)*ini_e_SI**2.)**(870./2299.)))
            func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mSI_1+mSI_2)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 1e-8
            ecc_Xhz             = fsolve(func, ecc_initial_guess)   #value of bin ecc when the bin is at fGW_limit
            
        #save:
        ecc_1M_XHz_arr[c]       = ecc_Xhz
        tGW_1M_yrs_arr[c]       = tGW_1M_yrs     
    #----------------------------------------------------------


    #----------------------------------------------------------
    #calc properties of 2nd GW merger:
    #----------------------------------------------------------
    #define:
    vkick_SI    = vkick_kmsec*1000.
    rndnum_arr  = np.random.random((nrtot,2))
    rnd_phi     = ((2.*np.pi)*rndnum_arr[:,0])			    #P flat in phi			- pos on unitsphere
    rnd_theta   = np.arccos(2.*rndnum_arr[:,1]-1.)          #P flat in cos(theta)	- pos on unitsphere
    m_i         = (1.*M1)*M_sun_SI
    m_j         = (2.*M1)*M_sun_SI
    m_ij        = m_i+m_j
    mu_ij       = m_i*m_j/m_ij

    for c in range(0,nrtot):
        Etot                = output_Nbody_endstate_REAL[c,7]
        endstate_id         = output_Nbody_endstate_INT[c,0]
        rel_vel_NK_kmsec    = 0.0   #initialize
        if (Etot < 0.0):                                #bin-sin must be bound
            if (endstate_id == 2 or endstate_id == 5):  #endstate must be coll(2) or GWinsp(5)
                #before kick (no kick NK):
                NK_rel_pos_vec_SI   = output_Nbody_xtra_2_info_REAL[c,0:3]*R_sun_SI
                NK_rel_vel_vec_SI   = output_Nbody_xtra_2_info_REAL[c,3:6]*(1000./kmsec_U)
                #kick velocity vec (isotropic vkick dist):       
                phi                 = rnd_phi[c]
                theta               = rnd_theta[c]
                vel_vec_kick_SI     = vkick_SI*np.array([(np.sin(theta)*np.cos(phi)), (np.sin(theta)*np.sin(phi)), (np.cos(theta))])
                #after kick:
                rel_pos_vec_SI = NK_rel_pos_vec_SI
                rel_vel_vec_SI = NK_rel_vel_vec_SI + vel_vec_kick_SI            
                #calc orb params:
                rel_L_vec_SI    = mu_ij*np.cross(rel_pos_vec_SI, rel_vel_vec_SI)
                rel_r_SI        = np.sqrt(sum(rel_pos_vec_SI**2.))
                rel_v_SI        = np.sqrt(sum(rel_vel_vec_SI**2.))
                rel_L_SI        = np.sqrt(sum(rel_L_vec_SI**2.))
                a_SI    = 1./((2./rel_r_SI) - (rel_v_SI**2.)/(G_new_SI*m_ij))
                e_SI    = np.sqrt(1. - ((rel_L_SI**2.)/(G_new_SI*m_ij*a_SI*(mu_ij**2.))))
                rp_SI   = a_SI*(1.-e_SI)
                #define/save info:
                rel_vel_NK_kmsec        = np.sqrt(sum(NK_rel_vel_vec_SI**2.))/1000.
                rel_vel_NK_kmsec_arr[c] = rel_vel_NK_kmsec       
                #if 2nd binary is still bound:            
                if (a_SI > 0.0):
                    ini_a_SI    = a_SI
                    ini_e_SI    = e_SI
                    #calc tGW merge:
                    tGW_2M_yrs              = ((768./425.)*((ini_a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m_i*m_j*(m_i+m_j)/(c_SI**5.)))*((1.-ini_e_SI**2.)**(7./2.)))/yr_sec            
                    tGW_2M_yrs_arr[c]       = tGW_2M_yrs
                    #fGW for 2nd binary at formation:
                    fGW_SI                  = (1./np.pi)*np.sqrt(G_new_SI*(m_i + m_j)/(rp_SI**3.))
                    fGW_2nd_bin_SI_arr[c]   = fGW_SI 
                    #ecc at Xhz for second GW merger:
                    if (fGW_SI > fGW_limit):
                        ecc_Xhz     = 1.0        
                    if (fGW_SI < fGW_limit):
                        #propagating ini a,e to a,e(fGW_limit)
                        c0      = ini_a_SI/((ini_e_SI**(12./19.)/(1.-ini_e_SI**2.))*((1.+(121./304.)*ini_e_SI**2.)**(870./2299.)))
                        func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(m_i+m_j)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
                        ecc_initial_guess       = 1e-8
                        ecc_Xhz                 = fsolve(func, ecc_initial_guess)   #value of bin ecc when the bin is at fGW_limit
                        ecc_2M_XHz_arr[c]       = ecc_Xhz
    #----------------------------------------------------------
    #save:
    #----------------------------------------------------------
    output_all_arrays   = [ecc_1M_XHz_arr, tGW_1M_yrs_arr, tGW_2M_yrs_arr, ecc_2M_XHz_arr, fGW_2nd_bin_SI_arr, rel_vel_NK_kmsec_arr]
    saveoutput_arr      = output_all_arrays  
    np.savez(data_name + '_output_all_arrays', saveoutput_arr)    
#--------------------------------------------------------------
#--------------------------------------------------------------


#----------------------------------------------------------
#LOAD DATA:
#----------------------------------------------------------
#output_all_arrays   = [ecc_1M_XHz_arr, tGW_1M_yrs_arr, tGW_2M_yrs_arr, ecc_2M_XHz_arr, fGW_2nd_bin_SI_arr, rel_vel_NK_kmsec_arr]
output_arr              = np.load(data_name + '_output_all_arrays' + '.npz')['arr_0']
output_all_arrays       = output_arr
#define:
ecc_1M_XHz_arr          = output_all_arrays[0]
tGW_1M_yrs_arr          = output_all_arrays[1]
tGW_2M_yrs_arr          = output_all_arrays[2]
ecc_2M_XHz_arr          = output_all_arrays[3]
fGW_2nd_bin_SI_arr      = output_all_arrays[4]
rel_vel_NK_kmsec_arr    = output_all_arrays[5]
#----------------------------------------------------------


#----------------------------------------------------------
#outlist_MC_info_INT                #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
#outlist_MC_info_REAL               #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, R_val, b_val, f_val, g_val]
#output_Nbody_endstate_INT          #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
#output_Nbody_endstate_REAL         #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#output_Nbody_xtra_info_REAL        #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
#output_Nbody_xtra_2_info_REAL      #[out_pos_CMij_wrt_sink(x,y,z), out_vel_CMij_wrt_sink(x,y,z)]

#reshape to order: [b/g, f, 10]
#xy_outlist_MC_info_REAL         = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   
#xy_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(res_pa, res_pa, 10)             #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
#xy_output_Nbody_endstate_REAL   = output_Nbody_endstate_REAL.reshape(res_pa, res_pa, 10)            #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#xy_output_Nbody_xtra_info_REAL  = output_Nbody_xtra_info_REAL.reshape(res_pa, res_pa, 10)           #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
#----------------------------------------------------------


#----------------------------------------------------------
#sort and define:
#----------------------------------------------------------
#ecc limit:
ecc_limit   = 0.1
#time-limit for 1st/2nd GW merger:
t_limit     = 10.**5.   #in years

pos_b_LT_0          = np.where(outlist_MC_info_REAL[:,7] < 0.0)[0]
pos_b_GT_0          = np.where(outlist_MC_info_REAL[:,7] > 0.0)[0]

pos_id5             = np.where(output_Nbody_endstate_INT[:,0] == 5)[0]
pos_id2             = np.where(output_Nbody_endstate_INT[:,0] == 2)[0]
pos_ELT0            = np.where(output_Nbody_endstate_REAL[:,7] < 0.0)[0]

pos_2M_LT_tlimit    = np.where(tGW_2M_yrs_arr[:] < t_limit)[0]
pos_1M_LT_tlimit    = np.where(tGW_1M_yrs_arr[:] < t_limit)[0]

pos_1M_GT_ecclimit  = np.where(ecc_1M_XHz_arr[:] > ecc_limit)[0]

pos_id3             = np.where(output_Nbody_endstate_INT[:,0] == 3)[0]
pos_s1              = np.where(output_Nbody_endstate_INT[:,3] == 1)[0]
pos_s2              = np.where(output_Nbody_endstate_INT[:,3] == 2)[0]

pos_id25                                    = list(pos_id2) + list(pos_id5)
pos_id25_ELT0                               = list(set(pos_id25).intersection(pos_ELT0))
pos_id25_ELT0_2M_LT_tlimit                  = list(set(pos_id25_ELT0).intersection(pos_2M_LT_tlimit))
pos_id25_ELT0_2M_LT_tlimit_1M_GT_ecclimit   = list(set(pos_id25_ELT0_2M_LT_tlimit).intersection(pos_1M_GT_ecclimit))

pos_s12             = list(pos_s1) + list(pos_s2)
pos_exch            = list(set(pos_id3).intersection(pos_s12))
#----------------------------------------------------------







#----------------------------------------------------------
#plot: f,b EXAMPLE:
#----------------------------------------------------------
pos_2b              = list(set(pos_id3).intersection(pos_1M_LT_tlimit))
pos_2b_GT_ecclimit  = list(set(pos_2b).intersection(pos_1M_GT_ecclimit))
pos_3b              = list(set(pos_id25_ELT0).intersection(pos_1M_LT_tlimit))
pos_3b_GT_ecclimit  = list(set(pos_3b).intersection(pos_1M_GT_ecclimit))
pos_exch_no2b       = list(set(pos_exch).difference(pos_2b))

#2nd bin: LISA
fGW_Ll      = 1e-3
fGW_Ul      = 1e-1
pos_2bin_fGW_LISA   = np.where((fGW_2nd_bin_SI_arr[:] > fGW_Ll) & (fGW_2nd_bin_SI_arr[:] < fGW_Ul))[0]
#2nd bin: LIGO
fGW_Ll      = 10
fGW_Ul      = 1e3
pos_2bin_fGW_LIGO   = np.where((fGW_2nd_bin_SI_arr[:] > fGW_Ll) & (fGW_2nd_bin_SI_arr[:] < fGW_Ul))[0]

#define:
xy_outlist_MC_info_REAL = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   
vel_vc                  = np.sqrt((3./2.)*M1/sma_a0)                                        #EQUAL MASS LIMIT!
vel_v                   = xy_outlist_MC_info_REAL[0,0,1]/vel_vc                             #EQUAL MASS LIMIT!
imp_b_unita0v_arr       = xy_outlist_MC_info_REAL[:,0,7]/(sma_a0/vel_v)
phase_f_unitrad_arr     = xy_outlist_MC_info_REAL[0,:,8]
X, Y                    = np.meshgrid(phase_f_unitrad_arr, imp_b_unita0v_arr)
minf    = min(phase_f_unitrad_arr)
maxf    = max(phase_f_unitrad_arr)
minb    = min(imp_b_unita0v_arr)
maxb    = max(imp_b_unita0v_arr)

#FIG 1:

pf  = 5
fig = plt.figure(figsize=(pf*5, pf*3))
#define:
plot_arr                        = np.zeros(nrtot, dtype=int)
#values for plotting:
plot_arr[:]                     = 0
plot_arr[pos_2b_GT_ecclimit]    = 1
plot_arr[pos_3b_GT_ecclimit]    = 2
plot_arr[pos_exch_no2b]         = 3
#change dim:
plot_arr_xy = plot_arr.reshape(res_pa, res_pa)
#plot:
ax      = plt.subplot(111)
cMap    = ccolor.ListedColormap(['dimgray', 'deeppink', 'deepskyblue', 'black', ])
Z       = plot_arr_xy[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=3.0, alpha=1.0)
#color bar:
bounds  = np.array([0,1,2,3,4])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .95)
cbar.ax.set_yticklabels([r'Remain.', r'$2$-body', r'$3$-body', r'Exch.'], rotation=270, fontsize = pf*10)
ax.plot([-10,10], [0,0], marker = '', linestyle = ':', color = 'white', linewidth = 2*pf)
plt.xticks(fontsize=10*pf)
plt.yticks(fontsize=10*pf)
#axis min max f,b;
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
#labels:
ax.set_xlabel(r'binary phase $f$ [radians]', fontsize=10*pf)
ax.set_ylabel(r'impact parameter $b^{\prime}$', fontsize=10*pf)
ax.text(0.1, -2.95,  r'$e>0.1$ at $10\ Hz$', horizontalalignment='left', verticalalignment='center', transform=ax.transData, fontsize = 10*pf, color = 'white')
ax.text(0.05, 0.95, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$, $a = 0.1$ AU, $e = 0$', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize = pf*10, color='white')
#save and show fig:
save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/'
plt.savefig(save_data_folder + data_name + '_phasespace_2.pdf', bbox_inches='tight')
plt.show()

#exit()

pf  = 5
fig = plt.figure(figsize=(pf*5, pf*3))
#define:
plot_arr                        = np.zeros(nrtot, dtype=int)
#values for plotting:
plot_arr[:]                     = 0
plot_arr[pos_2b]                = 1
plot_arr[pos_3b]                = 2
plot_arr[pos_exch_no2b]         = 3
#change dim:
plot_arr_xy = plot_arr.reshape(res_pa, res_pa)
#plot:
ax      = plt.subplot(111)
cMap    = ccolor.ListedColormap(['dimgray', 'deeppink', 'deepskyblue', 'black', ])
Z       = plot_arr_xy[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=3.0, alpha=1.0)
#color bar:
bounds  = np.array([0,1,2,3,4])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .95)
cbar.ax.set_yticklabels([r'Remain.', r'$2$-body', r'$3$-body', r'Exch.'], rotation=270, fontsize = pf*10)
ax.plot([-10,10], [0,0], marker = '', linestyle = ':', color = 'white', linewidth = 2*pf)
plt.xticks(fontsize=10*pf)
plt.yticks(fontsize=10*pf)
#axis min max f,b;
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
#labels:
ax.set_xlabel(r'binary phase $f$ [radians]', fontsize=10*pf)
ax.set_ylabel(r'impact parameter $b^{\prime}$', fontsize=10*pf)
#ax.text(0.1, -2.95,  r'All', horizontalalignment='left', verticalalignment='center', transform=ax.transData, fontsize = 10*pf, color = 'white')
ax.text(0.05, 0.95, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$, $a = 0.1$ AU, $e = 0$', horizontalalignment='left', verticalalignment='center', transform=ax.transAxes, fontsize = pf*10, color='white')
#save and show fig:
save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/'
plt.savefig(save_data_folder + data_name + '_phasespace_1.pdf', bbox_inches='tight')
plt.show()



exit()

#FIG 3:
fig = plt.figure(figsize=(5, 3))
#define:
plot_arr                        = np.zeros(nrtot, dtype=int)
#values for plotting:
plot_arr[:]                     = 0
plot_arr[pos_2bin_fGW_LISA]     = 1
plot_arr[pos_2bin_fGW_LIGO]     = 2
#change dim:
plot_arr_xy = plot_arr.reshape(res_pa, res_pa)
#plot:
ax      = plt.subplot(111)
cMap    = ccolor.ListedColormap(['black', 'gold', 'lime', ])
Z       = plot_arr_xy[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=2.0, alpha=1.0)
#color bar:
bounds  = np.array([0,1,2,3])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .95)
cbar.ax.set_yticklabels(['rest', '2nd(LISA)', '2nd(LIGO)'], rotation=270, fontsize = 8)
ax.plot([-10,10], [0,0], marker = '', linestyle = ':', color = 'white')
#analytical bounds:
f_LISA  = 1e-3
b0  = -np.sqrt(3.)/2.
rf  = ((G_new_SI*(3.*M1*M_sun_SI)/((f_LISA**2.)*(np.pi**2.)))**(1./3.))/R_sun_SI
Cc  = (16./3.)*(rf/sma_a0)  
db  = b0*np.sqrt(Cc)
b1  = b0+db
b2  = b0-db
ax.plot([-10,10], [b1,b1], marker = '', linestyle = ':', color = 'red')
ax.plot([-10,10], [b2,b2], marker = '', linestyle = ':', color = 'red')
#axis min max f,b;
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
#labels:
ax.set_xlabel(r'binary phase $f$ [radians]')
ax.set_ylabel(r'impact parameter $b^{\prime}$')
ax.text(0.1, -2.95,  r'2nd merger', horizontalalignment='left', verticalalignment='center', transform=ax.transData, fontsize = 10, color = 'white')
#save and show fig:
save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/'
plt.savefig(save_data_folder + data_name + '_phasespace_3.eps', bbox_inches='tight')
plt.show()
#----------------------------------------------------------




exit()

















#----------------------------------------------------------
#plot: f,b EXAMPLE:
#----------------------------------------------------------
#n_nr_SI     = (10.**(4.))/(m_parsec**3.)
#v_SI        = 10.0*(10.**3.)
#m_SI        = M1*M_sun_SI
#a_SI        = sma_a0*R_sun_SI
#sigma_bs_SI = 2.*np.pi*G_new_SI*(3.*m_SI)*a_SI/(v_SI**2.)
#Gamma_bs_SI = n_nr_SI*sigma_bs_SI*v_SI
#tbs_SI      = 1./Gamma_bs_SI
#tbs_year    = tbs_SI/yr_sec
#t_limit     = tbs_year


#define:
plot_arr    = np.zeros(nrtot, dtype=int)
#values for plotting:
plot_arr[:]                                 = 0 #remaining set
plot_arr[pos_exch]                          = 1
plot_arr[pos_1M_GT_ecclimit]                   = 2
plot_arr[pos_id25_ELT0_2M_LT_tlimit_1M_GT_ecclimit]  = 3
#change dim:
plot_arr_xy             = plot_arr.reshape(res_pa, res_pa)
xy_outlist_MC_info_REAL = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   

#define:
vel_vc              = np.sqrt((3./2.)*M1/sma_a0)                #EQUAL MASS LIMIT!
vel_v               = xy_outlist_MC_info_REAL[0,0,1]/vel_vc     #EQUAL MASS LIMIT!
imp_b_unita0v_arr   = xy_outlist_MC_info_REAL[:,0,7]/(sma_a0/vel_v)
phase_f_unitrad_arr = xy_outlist_MC_info_REAL[0,:,8]

#make x,y grid:
X, Y    = np.meshgrid(phase_f_unitrad_arr, imp_b_unita0v_arr)

fig = plt.figure(figsize=(6, 5))#plt.figure(figsize=(15, 12.5))
ax  = plt.subplot(111)
cMap    = ccolor.ListedColormap(['black', '#377eb8', '#ffff99', 'crimson', ])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
Z       = plot_arr_xy[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=3.0)
#color bar:
bounds = np.array([0,1,2,3,4])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .5)
cbar.ax.set_yticklabels([r'Remaining', 'Exchange', r'$e_{1}>0.1$', r'$e_{1}>0.1$, $t_{12}<10$ yrs'], rotation=270, fontsize = 20)

#axis min max f,b;
minf    = min(phase_f_unitrad_arr)
maxf    = max(phase_f_unitrad_arr)
minb    = min(imp_b_unita0v_arr)
maxb    = max(imp_b_unita0v_arr)
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
#labels:
#ax.set_title(r'Topology')
ax.set_xlabel(r'binary phase $f$ [radians]')
ax.set_ylabel(r'impact parameter $b^{\prime}$')

#save and show fig:
#save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/doubleGWmergers/'
#plt.savefig(save_data_folder + data_name + '_3BH_phasespace_ex.eps', bbox_inches='tight')
plt.show()
#----------------------------------------------------------



exit()






#----------------------------------------------------------
#plot rel vel:
#----------------------------------------------------------
fig = plt.figure(figsize=(15, 12.5))
ax  = plt.subplot(111)
#define:
plot_arr    = np.zeros(nrtot, dtype=np.float64)
#values for plotting:
plot_arr[:] = np.log10(rel_vel_NK_kmsec_arr[:])
#change dim:
plot_arr_xy             = plot_arr.reshape(res_pa, res_pa)
xy_outlist_MC_info_REAL = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   
#define:
vel_vc              = np.sqrt((3./2.)*M1/sma_a0)                #EQUAL MASS LIMIT!
vel_v               = xy_outlist_MC_info_REAL[0,0,1]/vel_vc     #EQUAL MASS LIMIT!
imp_b_unita0v_arr   = xy_outlist_MC_info_REAL[:,0,7]/(sma_a0/vel_v)
phase_f_unitrad_arr = xy_outlist_MC_info_REAL[0,:,8]
#make x,y grid:
X, Y    = np.meshgrid(phase_f_unitrad_arr, imp_b_unita0v_arr)
Z       = plot_arr_xy[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap='magma_r', linewidth=0, rasterized=True, vmin=3.0, vmax=4.0)
plt.colorbar()
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)

#save and show fig:
save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/doubleGWmergers/'
plt.savefig(save_data_folder + data_name + '_3BH_phasespace_relvel.eps', bbox_inches='tight')
plt.show()
#----------------------------------------------------------


exit()







#----------------------------------------------------------
#TEST:
#----------------------------------------------------------
m_i     = (1.*M1)*M_sun_SI
m_j     = (2.*M1)*M_sun_SI
m_ij    = m_i+m_j
mu_ij   = m_i*m_j/m_ij
Rs      = (2.*G_new_SI*m_ij/(c_SI**2.))
Ms      = mu_ij
eps     = 85.*np.pi/96.
beta    = 7./2.

pos     = np.where(output_Nbody_endstate_INT[:,0] == 5)[0]
a       = output_Nbody_endstate_REAL[pos,8]
e       = output_Nbody_endstate_REAL[pos,9]
bp      = outlist_MC_info_REAL[pos,7]/(sma_a0/0.01)

#calc rp:
rp_sim  = (a*(1.-e)/sma_a0)
rp_cal  = (3./16.)*((2./np.sqrt(3.))*bp - 1.)**2.

#calc tlife:
a_SI    = a*R_sun_SI
rp_SI   = a_SI*(1.-e)
tl_sim  = ((2.*np.pi*(rp_SI**(beta))*np.sqrt(a_SI)*(m_ij*mu_ij*(Rs**(1.-beta))/(eps*(Ms**(2.))*np.sqrt(G_new_SI*m_ij)))))/yr_sec
tl_cal  = ((4.*np.sqrt(2.)/85.)*np.sqrt(a_SI)*(c_SI**5.)/((G_new_SI**3.)*((M1*M_sun_SI)**3.))*(rp_SI**(7./2.)))/yr_sec

print ((((sma_a0*R_sun_SI)**4.)/(4.*(6.*64./5.)*(G_new_SI**3.)*((M1*M_sun_SI)**3.)/(c_SI**5.))))/yr_sec

print tl_sim
print tl_cal

#plot:
fig = plt.figure(figsize=(10, 8))
ax  = plt.subplot(111)

ax.plot(bp, np.log10(tl_sim), marker='o', markersize=2.0, linestyle='', color='black')
ax.plot(bp, np.log10(tl_cal), marker='o', markersize=2.0, linestyle='', color='red')

#ax.plot( bp, rp_sim, marker='o', markersize=2.0, linestyle='', color='black')
#ax.plot(-bp, rp_cal, marker='o', markersize=2.0, linestyle='', color='red')

plt.show()
exit()
#----------------------------------------------------------



        

