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
from matplotlib import colors as c
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
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif',
        'size'   : 35}
mpl.rc('font', **font)



#data name:
data_name = 'TOP_WGR_101010Msun_a1e0_v001vc_fbfull_gamma00pi_fbres500'
M1 = 10.
M2 = 10.
M3 = 10.



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


#----------------------------------------------------------
#Define:
#----------------------------------------------------------
#MC_settings_list_INT #[nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 0, 0]
res_pa  = int(np.sqrt(MC_settings_list_INT[2]))
sma_a0  = outlist_MC_info_REAL[0,0]
#----------------------------------------------------------
#reshape input arrays for easy analysis:
#----------------------------------------------------------
#outlist_MC_info_INT                #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
#outlist_MC_info_REAL               #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, R_val, b_val, f_val, g_val]
#output_Nbody_endstate_INT          #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
#output_Nbody_endstate_REAL         #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#output_Nbody_xtra_info_REAL        #[rmin_12, rmin_13, rmin_23]

#reshape to order: [b/g, f, 10]
#xy_outlist_MC_info_REAL         = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   
#xy_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(res_pa, res_pa, 10)             #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
#xy_output_Nbody_endstate_REAL   = output_Nbody_endstate_REAL.reshape(res_pa, res_pa, 10)            #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
#xy_output_Nbody_xtra_info_REAL  = output_Nbody_xtra_info_REAL.reshape(res_pa, res_pa, 10)           #[rmin_12, rmin_13, rmin_23]
#----------------------------------------------------------



print 'test'


pos_b_LT_0  = np.where(outlist_MC_info_REAL[:,7] < 0.0)[0]
pos_b_GT_0  = np.where(outlist_MC_info_REAL[:,7] > 0.0)[0]

pos_id5     = np.where(output_Nbody_endstate_INT[:,0] == 5)[0]
pos_id2     = np.where(output_Nbody_endstate_INT[:,0] == 2)[0]
pos_E_LT_0  = np.where(output_Nbody_endstate_REAL[:,7] < 0.0)[0]

pos_id25    = list(set(pos_id5).intersection(pos_id2))

















