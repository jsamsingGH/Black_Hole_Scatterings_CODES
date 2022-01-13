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
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI 
#----------------------------------------------------------


#print ((1./50.)/((2.*30.*Rsch_1Msun_unitRsun)/AU_U))*4.76711397547
#exit()

#------------------------------------------------------------------------------------------
#read data and define:
#------------------------------------------------------------------------------------------
#data name and folder:
data_folder     = '/Users/jsamsing/Desktop/TIDES_PROJ/MOCCA_ICs/MC_1/'
data_name       = 'MOCCATEST_case1_'

#read data:
tf = open(data_folder+data_name+'nbody_params_INT.txt', "r")
nbody_params_INT_all        = np.loadtxt(tf, dtype=int)
tf.close()
tf = open(data_folder+data_name+'nbody_params_REAL.txt', "r")
nbody_params_REAL_all       = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'obj_info.txt', "r")
obj_info_all                = np.loadtxt(tf, dtype=float)
tf.close()
tf = open(data_folder+data_name+'pos_vel.txt', "r")
pos_vel_all                 = np.loadtxt(tf, dtype=float)
tf.close()

#define:
nr_tot_scatterings  = len(nbody_params_INT_all)
#------------------------------------------------------------------------------------------

#3170
#5293
#33898
#34274
#64974
#66009
#86108
#86791
#99698
#100276
#116781
#117538
#134447
#137938

print nr_tot_scatterings

for tc in range(0,nr_tot_scatterings):                                 
    print tc
    if (tc == 3170):
        print 'TEST'
        #---------------------------------------
        #get IC info from input file:
        #---------------------------------------    
        #Nbody params:
        nbody_params_arr_1_INT      = nbody_params_INT_all[tc,:]
        nbody_params_arr_2_REAL     = nbody_params_REAL_all[tc,:]
        #obj 1:
        b1_const_arr    = obj_info_all[tc, 0:10]
        b1_posxyz_CM    = pos_vel_all[tc,  0:3]
        b1_velxyz_CM    = pos_vel_all[tc,  3:6]
        b1_q            = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
        b1_qdot         = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
        #obj 2:
        b2_const_arr    = obj_info_all[tc, 10:20]
        b2_posxyz_CM    = pos_vel_all[tc,  6:9]
        b2_velxyz_CM    = pos_vel_all[tc,  9:12]
        b2_q            = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
        b2_qdot         = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
        #obj 3:
        b3_const_arr    = obj_info_all[tc, 20:30]
        b3_posxyz_CM    = pos_vel_all[tc,  12:15]
        b3_velxyz_CM    = pos_vel_all[tc,  15:18]
        b3_q            = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
        b3_qdot         = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
        #---------------------------------------
        #open file:
        #---------------------------------------
        fn = 'MCinput_Nbody.txt'
        text_file = open(fn, "w")
        #nr particles in the sim   
        text_file.write('3' + '\n')
        #nbody code settings:
       
        aini    = np.sqrt(sum((b1_posxyz_CM-b2_posxyz_CM)**2.))
        
        nbody_params_arr_1_INT[5]   = 1
        nbody_params_arr_2_REAL[8]  = (1./10.)*aini     
        
        np.savetxt(text_file, nbody_params_arr_1_INT[None],   fmt='%10i')
        np.savetxt(text_file, nbody_params_arr_2_REAL[None],  fmt='%10f')
        #body 1:
        np.savetxt(text_file, b1_const_arr[None],  fmt='%10f')
        np.savetxt(text_file, b1_posxyz_CM[None],  fmt='%10f')
        np.savetxt(text_file, b1_velxyz_CM[None],  fmt='%10f')
        np.savetxt(text_file, b1_q,  fmt='%10f')
        np.savetxt(text_file, b1_qdot,  fmt='%10f')
        #body 2:
        np.savetxt(text_file, b2_const_arr[None],  fmt='%10f')
        np.savetxt(text_file, b2_posxyz_CM[None],  fmt='%10f')
        np.savetxt(text_file, b2_velxyz_CM[None],  fmt='%10f')
        np.savetxt(text_file, b2_q,  fmt='%10f')
        np.savetxt(text_file, b2_qdot,  fmt='%10f')
        #body 3:
        np.savetxt(text_file, b3_const_arr[None],  fmt='%10f')
        np.savetxt(text_file, b3_posxyz_CM[None],  fmt='%10f')
        np.savetxt(text_file, b3_velxyz_CM[None],  fmt='%10f')
        np.savetxt(text_file, b3_q,  fmt='%10f')
        np.savetxt(text_file, b3_qdot,  fmt='%10f')
        #close file:
        text_file.close()   
        #run:
        subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
        exit()
        #---------------------------------------
    
    
#------------------------------------------------------------------------------------------
#check ... 
#------------------------------------------------------------------------------------------
output_data_folder  = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'       #output data from N-body simulations
filename            = 'slurm-3662695.out'
#read data:
tf = open(output_data_folder+filename, "r")
slurm_file          = np.loadtxt(tf, dtype=int)
tf.close()

nrtot       = len(slurm_file)
arr_sorted = np.array(sorted(slurm_file))

c = 0
for tc in range(0,nrtot):
    if (arr_sorted[c] == arr_sorted[c+1]):
        c = c + 2
    if (arr_sorted[c] != arr_sorted[c+1]):
        print arr_sorted[c]
        c = c + 1
        
#exit()
#------------------------------------------------------------------------------------------


