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
font = {'family' : 'serif',
        'size'   : 20}
mpl.rc('font', **font)


#data name:
data_name = 'Test_Lvec_1'#'TEST_SPN_2'#'Test_Lvec_BH20BH20S1'#'Test_Lvec_2_rminlarge'#'Test_Lvec_2_rminsmall'#'Test_Lvec_2'#'Test_Lvec_1'
M1 = 20.


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
tf = open(data_folder+data_name+'output_Nbody_xtra_info_INT.txt', "r")
output_Nbody_xtra_info_INT   = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_info_REAL.txt', "r")
output_Nbody_xtra_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'output_Nbody_xtra_2_info_REAL.txt', "r")
output_Nbody_xtra_2_info_REAL  = np.loadtxt(tf, dtype=float)
tf.close()


nr_SMA                          = MC_settings_list_INT[0]
nr_vinf                         = MC_settings_list_INT[1]
nr_scatterings_per_paramcomb    = MC_settings_list_INT[2]

#----------------------------------------------------------
#Define:
#----------------------------------------------------------
#MC_settings_list_INT #[nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 0, 0]
#res_pa  = int(np.sqrt(MC_settings_list_INT[2]))
sma_a0  = outlist_MC_info_REAL[0,0]
print sma_a0/AU_U
#----------------------------------------------------------
#reshape input arrays for easy analysis:
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
#Analyze: EQUAL MASS EXAMPLE.
#----------------------------------------------------------
#define:
tot_nr_scat     = nr_scatterings_per_paramcomb
info_arr_all    = np.zeros((tot_nr_scat, 10), dtype=np.float64)

#loop over data:
for c in range(0, tot_nr_scat):
    
    #define:
    mSI = M1*M_sun_SI
            
    #calc v_kick:
    Etot_binsin         = output_Nbody_endstate_REAL[c,7]
    ESI_fac             = G_new_SI*(M_sun_SI**2.)/R_sun_SI
    Etot_SI             = Etot_binsin*ESI_fac
    mbin                = mSI+mSI
    msin                = mSI
    v_inf_bin           = np.sqrt(2.*Etot_SI*((mbin*(1.+mbin/msin))**(-1.)))   #SI units
    v_kick_kms          = v_inf_bin/1000.
    
    #Lorb angle:
    unit_z_vec  = np.array([0.0, 0.0, 1.0])
    Lorb_ij     = output_Nbody_xtra_2_info_REAL[c,6:9]
    len_unit_z_vec  = 1.0
    len_Lorb_ij     = np.sqrt(sum(Lorb_ij[:]**2.))
    theta_Lorb_zvec = np.arccos(np.dot(Lorb_ij, unit_z_vec)/(len_Lorb_ij*len_unit_z_vec))
    
    #save info:
    info_arr_all[c,0] = v_kick_kms
    info_arr_all[c,1] = theta_Lorb_zvec
    

##TEST collisons plot:
#pos_1 = np.where(output_Nbody_endstate_INT[:,0] == 2)[0]
#pos_2 = np.where(output_Nbody_endstate_REAL[:,2] < 100.0)[0]
#pos   = list(set(pos_1).intersection(pos_2))
##plot:
#fig, ax1 = plt.subplots(figsize=(8, 6))
#input_arr = info_arr_all[pos,1]
#ax1.hist(input_arr, bins=100, range = [0.0, np.pi], normed=False, linewidth=3.0, histtype='step', color='black')
#angle_arr = np.arange(0.0, np.pi, 0.01)
#ax1.plot(angle_arr, ((1.*len(pos))*(1./2.)*np.sin(angle_arr))*(np.pi/(1.*100.)), linestyle='--', linewidth=2.0, color='green')
#plt.show()
#exit()


#define:
a_binij_all_unit_a0 = output_Nbody_endstate_REAL[:,3]/sma_a0


#organize data:
pos_cut_vkick   = np.where(info_arr_all[:,0] > 50.0)[0]
pos_cut_SMAred  = np.where(a_binij_all_unit_a0 < 0.5)[0]  

pos_cut_IMSc1   = np.where(output_Nbody_endstate_INT[:,6] <= 1)[0]
pos_cut_IMSc2   = np.where(output_Nbody_endstate_INT[:,6] >= 2)[0]

pos_id3         = np.where(output_Nbody_endstate_INT[:,0] == 3)[0]

pos_s1          = np.where(output_Nbody_endstate_INT[:,3] == 1)[0]
pos_s2          = np.where(output_Nbody_endstate_INT[:,3] == 2)[0]
pos_s3          = np.where(output_Nbody_endstate_INT[:,3] == 3)[0]

pos_id3_s1      = list(set(pos_id3).intersection(pos_s1))
pos_id3_s2      = list(set(pos_id3).intersection(pos_s2))
pos_id3_s3      = list(set(pos_id3).intersection(pos_s3))



#TEST PLOT:
pos_intcut = pos_cut_vkick
#pos_intcut = pos_cut_SMAred
pos_1 = list(set(pos_intcut).intersection(pos_id3_s1))
pos_2 = list(set(pos_intcut).intersection(pos_id3_s2))
pos_3 = list(set(pos_intcut).intersection(pos_id3_s3))
#define:
nrHbins_0pi = 250
nr_pos123   = len(pos_1) + len(pos_2) + len(pos_3) 

#print test:
print len(pos_1), len(pos_2), len(pos_3)
print len(list(set(pos_id3).intersection(pos_intcut))), nr_pos123




#plot:
fig, ax1 = plt.subplots(figsize=(8, 6))

input_arr = info_arr_all[pos_3,1]
ax1.hist(input_arr, bins=nrHbins_0pi, range = [0.0, np.pi], normed=False, linewidth=3.0, histtype='step', color='black', alpha=0.5, label=r'bin$[1,2]$ (all)')

pos = list(set(pos_cut_IMSc1).intersection(pos_3))
input_arr = info_arr_all[pos,1]
ax1.hist(input_arr, bins=nrHbins_0pi, range = [0.0, np.pi], normed=False, linewidth=1.0, histtype='step', color = 'grey', alpha=0.5, label=r'bin$[1,2]$ ($N_{\rm IMS} \leq 1$)', linestyle=':')
pos = list(set(pos_cut_IMSc2).intersection(pos_3))
input_arr = info_arr_all[pos,1]
ax1.hist(input_arr, bins=nrHbins_0pi, range = [0.0, np.pi], normed=False, linewidth=1.0, histtype='step', color = 'red', alpha=0.5, label=r'bin$[1,2]$ ($N_{\rm IMS} \geq 2$)', linestyle=':')

input_arr = info_arr_all[pos_2,1]
ax1.hist(input_arr, bins=nrHbins_0pi, range = [0.0, np.pi], normed=False, linewidth=3.0, histtype='step', color='blue', alpha=0.5,  label=r'bin$[1,3]$')

input_arr = info_arr_all[pos_1,1]
ax1.hist(input_arr, bins=nrHbins_0pi, range = [0.0, np.pi], normed=False, linewidth=3.0, histtype='step', color='orange', alpha=0.5,   label=r'bin$[2,3]$')

angle_arr = np.arange(0.0, np.pi, 0.01)
ax1.plot(angle_arr, ((1.*nr_pos123/3.)*(1./2.)*np.sin(angle_arr))*(np.pi/(1.*nrHbins_0pi)), linestyle='--', linewidth=2.0, color='green', label='isotropic')

#labels etc:
ax1.legend(loc='upper right', numpoints = 1, fontsize = 15.0, frameon = False)
ax1.set_xlabel(r'angle $\theta$')
ax1.set_ylabel(r'Counts')
ax1.set_title(r'[BH$_{1}$, BH$_{2}$] $\leftarrow$ BH$_{3}$')
ax1.set_xlim(0.0,np.pi)

#save figure:
plt.savefig(data_name + '_Langle_dist_test.pdf', bbox_inches='tight')

plt.show()

exit()




#plot:
fig, ax1 = plt.subplots(figsize=(8, 6))
input_arr = a_binij_all_unit_a0[:]
ax1.hist(input_arr, bins=100, range = [0.0,2.0], normed=False, linewidth=3.0, histtype='step', color='black')
plt.show()
exit()







#----------------------------------------------------------
































        

