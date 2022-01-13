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







#----------------------------------------------------------
#MOCCA proj:
#----------------------------------------------------------
path_data           = '/Users/jsamsing/Desktop/TIDES_PROJ/MOCCA_ICs/MC_1/'  #input and output data folder.
#-------------------------------------------
#INPUT data:
#-------------------------------------------
name_INput_IC_data  = 'MOCCATEST_inter-binsin-bhs-13gyr-m100.dat'
INput_IC_data  = path_data + name_INput_IC_data
#get total number of lines:
tf_indata   = open(INput_IC_data, "r")
nr_ICs      = len(tf_indata.readlines())
print 'nr ICs:', nr_ICs
#-------------------------------------------
#OUTPUT data:
#-------------------------------------------
name_OUTput_IC_data = 'MOCCATEST_case1A_'
#open files:
fn = path_data + name_OUTput_IC_data + 'nbody_params_INT.txt'
tf_nbody_params_INT     = open(fn, "w")

fn = path_data + name_OUTput_IC_data + 'nbody_params_REAL.txt'
tf_nbody_params_REAL    = open(fn, "w")

fn = path_data + name_OUTput_IC_data + 'obj_info.txt'
tf_obj_info             = open(fn, "w")

fn = path_data + name_OUTput_IC_data + 'pos_vel.txt'
tf_pos_vel              = open(fn, "w")

fn = path_data + name_OUTput_IC_data + 'sim_info_1.txt'
tf_sim_info_1           = open(fn, "w")
#-------------------------------------------
#Setting:
#-------------------------------------------
#nr scatterings per IC:
nr_sim_per_IC       = 5
#max_sim time:
Tsim_max_unit_Torb  = 2500.0       #max sim time in nr initial orbital times
#tidal sim-thresholds:
rp_max_FtidFbin     = 1.0  
r_simsurf_FtidFbin  = 0.01         #this must be << rp_max_FtidFbin
#define:
nr_tot_sim  = nr_ICs*nr_sim_per_IC  #this is an upper limit as some of the ICs will be skipped if not passing several thresholds...
rndnum_arr  = np.random.random((nr_tot_sim,5))
totc        = 0
#-------------------------------------------
#SIM CASE OPTIONS:
#-------------------------------------------
#input sim CASE:
sim_case_opt = 1    #CHANGE FILENAME ACCORDING TO THIS CHOICE!

#sim CASE 0:
if (sim_case_opt == 0):
    #sim settings:
    incl_25PN           = 0
    use_true_BH_rad     = 1
    use_eff_BH_rad      = 0
    sim_input_b         = 1
    sim_isotropic_b     = 0

#sim CASE 1:
if (sim_case_opt == 1):
    #sim settings:
    incl_25PN           = 1
    use_true_BH_rad     = 1
    use_eff_BH_rad      = 0
    sim_input_b         = 1
    sim_isotropic_b     = 0

#sim CASE 2:
if (sim_case_opt == 2):
    #sim settings:
    incl_25PN           = 1
    use_true_BH_rad     = 1
    use_eff_BH_rad      = 0
    sim_input_b         = 0
    sim_isotropic_b     = 1

#sim CASE 3:
if (sim_case_opt == 3):
    #sim settings:
    incl_25PN           = 0
    use_true_BH_rad     = 0
    use_eff_BH_rad      = 1
    sim_input_b         = 0
    sim_isotropic_b     = 1
    #sim param vals:
    eff_rad_fracrperi   = (1./20.)
#-------------------------------------------








#-------------------------------------------
#3body-BH_interactions:
#-------------------------------------------
#Open file:
tf_indata   = open(INput_IC_data, "r")
for nIC in range(0, nr_ICs):
    
    #---------------------------------------
    #Read input data line by line:
    #---------------------------------------
    line_indata         = tf_indata.readline()
    splitline_indata    = re.split(r'\s{1,}', line_indata)
    #---------------------------------------
    #IC info from INPUT file:
    #---------------------------------------    
    line0 = 1
    #units:
    MCunits_NbodyU_pos  = float(splitline_indata[109-line0])*1.0            #Units:mc2rsun
    MCunits_NbodyU_vel  = float(splitline_indata[110-line0])*kmsec_U        #Units:mc2kms * kmsec_U
    #mass:
    mb1         = float(splitline_indata[30-line0])*1.0                     #Msun
    mb2         = float(splitline_indata[31-line0])*1.0                     #Msun
    mb3         = float(splitline_indata[32-line0])*1.0                     #Msun
    #radius:
    Rb1         = float(splitline_indata[33-line0])*1.0                     #Rsun
    Rb2         = float(splitline_indata[34-line0])*1.0                     #Rsun
    Rb3         = float(splitline_indata[35-line0])*1.0                     #Rsun
    #orbit params:
    SMA_abin    = float(splitline_indata[39-line0])*1.0                     #Rsun 
    ecc_ebin    = float(splitline_indata[40-line0])*1.0
    impact_b    = float(splitline_indata[44-line0])*MCunits_NbodyU_pos      #Rsun
    relvel_v    = float(splitline_indata[45-line0])*MCunits_NbodyU_vel      #Rsun,Msun,...        
    r_peri      = SMA_abin*(1. - ecc_ebin)                                  #Rsun
    #regulate obj rad based on sim case:
    if (use_true_BH_rad == 1):
        Rb1 = Rb1
        Rb2 = Rb2
        Rb3 = Rb3
    if (use_eff_BH_rad  == 1):
        Rb1 = eff_rad_fracrperi*r_peri
        Rb2 = eff_rad_fracrperi*r_peri
        Rb3 = eff_rad_fracrperi*r_peri
    #from this point everything is in our Nbody code units.    
    #---------------------------------------
    #define:
    #---------------------------------------
    mass_tot    = mb1 + mb2 + mb3               #Msun
    mass_bin    = mb1 + mb2                     #Msun
    mass_sin    = mb3                           #Msun
    mass_red    = mass_bin*mass_sin/mass_tot    #Msun, reduced mass of bin-sin system
    #---------------------------------------
    #calc:
    #---------------------------------------
    #max impact param (b):
    rp_max          = SMA_abin*((1./rp_max_FtidFbin)*((mb1+mb2)*mb3/(mb1*mb2)))**(1./3.)            
    b_max           = rp_max*np.sqrt(1.+(2.*mass_tot)/(rp_max*(relvel_v**2.)))
    bimp_over_bmax  = impact_b/b_max
    #characteristic vel v_c:
    v_binsin_c      = np.sqrt((mb1*mb2*(mb1+mb2+mb3)/(mb3*(mb1+mb2)))*(1./SMA_abin))
    vinf_over_vc    = relvel_v/v_binsin_c
    #---------------------------------------
    #flags:
    #---------------------------------------
    sim_yesno                                           = 1     #(yes) initialize
    if (bimp_over_bmax  > 1.0):             sim_yesno   = 0     #(no)  
    if (vinf_over_vc    > 1.0):             sim_yesno   = 0     #(no)  
    #---------------------------------------
    #write sim info to file:
    #---------------------------------------
    sim_info_1_arr  = np.array([sim_yesno,0,0,0,0,0,0,0,0,0])
    np.savetxt(tf_sim_info_1, sim_info_1_arr[None], fmt='%10f')
    #---------------------------------------    
    #-------------------------------------------
    #Make 'rnd sample' IC for OUTPUT:
    #-------------------------------------------
    if (sim_yesno == 1):
    #-------------------------------------------
        for ns in range(0, nr_sim_per_IC):
            #-----------------------------------
            #setup binary:
            #-----------------------------------
            #find rnd phase f: (see: http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1983ApJ...268..319H&defaultprint=YES&filetype=.pdf)
            rnd02pi             = (2.*np.pi)*rndnum_arr[totc,0]
            func                = lambda eccE : rnd02pi - (eccE - ecc_ebin*np.sin(eccE))    
            eccE_initial_guess  = 1.0
            sol_eccE            = fsolve(func, eccE_initial_guess)
            bin_f               = 2.*np.arctan((((1.+ecc_ebin)/(1.-ecc_ebin))**(1./2.))*np.tan(sol_eccE/2.))[0]
            #calc/define:
            Mb12_SI     = (mb1+mb2)*M_sun_SI
            pfac_SI     = ((SMA_abin*R_sun_SI)*(1.-ecc_ebin**2))
            hfac_SI     = (np.sqrt(G_new_SI*Mb12_SI*pfac_SI))
            rdist_SI    = (pfac_SI/(1.+ecc_ebin*np.cos(bin_f)))
            v_ang_SI    = (hfac_SI/rdist_SI)
            v_rad_SI    = ((ecc_ebin*np.sin(bin_f)/pfac_SI)*hfac_SI)
            #calc pos and vel of reduced-mass object:
            posx_U  = - ((rdist_SI*np.cos(bin_f))/R_sun_SI)
            posy_U  =   ((rdist_SI*np.sin(bin_f))/R_sun_SI)
            velx_U  = - ((v_ang_SI*np.sin(bin_f) - v_rad_SI*np.cos(bin_f))*(kmsec_U/1000.))
            vely_U  = - ((v_rad_SI*np.sin(bin_f) + v_ang_SI*np.cos(bin_f))*(kmsec_U/1000.))
            #define pos/vel vecs:
            posvec_U = np.array([posx_U, posy_U, 0.0])
            velvec_U = np.array([velx_U, vely_U, 0.0])
            #change to binary COM:
            posxyz_1_U    = - (mb2/(mb1+mb2))*posvec_U
            velxyz_1_U    = - (mb2/(mb1+mb2))*velvec_U
            posxyz_2_U    =   (mb1/(mb1+mb2))*posvec_U
            velxyz_2_U    =   (mb1/(mb1+mb2))*velvec_U
            #-----------------------------------
            #setup incoming single:
            #-----------------------------------                
            #define/set:
            r_sampsurf  = 1e20*SMA_abin                 #Rsun, just make sure to put in a very big number.
            r_simsurf   = SMA_abin*((1./r_simsurf_FtidFbin)*((mb1+mb2)*mb3/(mb1*mb2)))**(1./3.)
            #Total energy of bin-sin system (where bin is a single par with mass Mbin)
            Etot_binsin = (1./2.)*mass_red*(relvel_v**2.)
            #rnd sampling parameters:
            psi_sampsurf        = (2.*np.pi)*rndnum_arr[totc,1]			                #P flat in psi			- angle on b disc.
            theta               = np.arccos(2.*rndnum_arr[totc,2]-1.)                   #P flat in cos(theta)	- pos on unitsphere
            phi                 = ((2.*np.pi)*rndnum_arr[totc,3])			            #P flat in phi			- pos on unitsphere
            if (sim_input_b     == 1): b_sampsurf = impact_b
            if (sim_isotropic_b == 1): b_sampsurf = b_max*np.sqrt(rndnum_arr[totc,4])   #P flat in b^2			- rad pos on 0->bmax disc
            #calc:
            rb_sampsurf         = np.sqrt(r_sampsurf**2 + b_sampsurf**2)    #dist to point 'b' at sampling surface.
            v_sampsurf          = np.sqrt((2./mass_red)*(Etot_binsin + mass_tot*mass_red/rb_sampsurf))
            #(b,psi,v) at simulation surface (r_simsurf):
            #By conservation of E,L we now propagate v,b onto the rsim surface:
            v_simsurf           = np.sqrt((2./mass_red)*(Etot_binsin + mass_tot*mass_red/r_simsurf))
            b_simsurf           = b_sampsurf*(v_sampsurf/v_simsurf)
            #angle dist is conserved onto rsim surface:
            psi_simsurf         = psi_sampsurf
            #Send obj3 in from a pos at the r_sim sphere:
            #initial x,y,z pos of body 3 at simplane:
            r_simplane  = np.sqrt(r_simsurf**2 - b_simsurf**2)
            pos_ini     = np.array([b_simsurf*np.cos(psi_simsurf), b_simsurf*np.sin(psi_simsurf), r_simplane])
            #rotate plane into a rnd pos on rsim sphere:
            #Make rotation matrices:
            Ry	= np.array([[np.cos(theta),0,np.sin(theta)], [0,1,0], [-np.sin(theta),0,np.cos(theta)]])
            Rz	= np.array([[np.cos(phi),-np.sin(phi),0],[np.sin(phi),np.cos(phi),0],[0,0,1]])
            #Final pos: rotate ini pos into new pos set by theta, phi:
            pos_Ry          = np.dot(Ry, pos_ini)   #rotate around y
            pos_RyRz        = np.dot(Rz, pos_Ry)    #rotate around z
            posxyz_3_U      = pos_RyRz              #final pos
            #Final vel: we know the vel is by constr parallel with the final pos simplane unitvec:
            unitvec_RyRz_simplane   = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
            velxyz_3_U              = - v_simsurf*unitvec_RyRz_simplane
            #perform check:
            #if (b_simsurf > r_simsurf): print 'STOP! choose other tidal thresholds!'    #check
            #print 'r_simsurf/SMA_abin', r_simsurf/SMA_abin
            #print 'rp_max/SMA_abin', rp_max/SMA_abin
            #print 'b_max/impact_b', b_max/impact_b
            #-----------------------------------
            #final pos,vel of obj 1,2,3:
            #-----------------------------------
            #posxyz_1_U, velxyz_1_U
            #posxyz_2_U, velxyz_2_U
            #posxyz_3_U, velxyz_3_U
            #find pos,vel of CM:
            pos_CM  = (mb1*posxyz_1_U + mb2*posxyz_2_U + mb3*posxyz_3_U)/mass_tot
            vel_CM  = (mb1*velxyz_1_U + mb2*velxyz_2_U + mb3*velxyz_3_U)/mass_tot
            #correct body 1,2,3 pos,vel to CM:
            #pos:
            posxyz_1_U  = posxyz_1_U - pos_CM
            posxyz_2_U  = posxyz_2_U - pos_CM
            posxyz_3_U  = posxyz_3_U - pos_CM
            #vel:
            velxyz_1_U  = velxyz_1_U - vel_CM
            velxyz_2_U  = velxyz_2_U - vel_CM
            velxyz_3_U  = velxyz_3_U - vel_CM
            #collect for output:
            posvelxyz_1_U  = np.append(posxyz_1_U,velxyz_1_U)
            posvelxyz_2_U  = np.append(posxyz_2_U,velxyz_2_U)
            posvelxyz_3_U  = np.append(posxyz_3_U,velxyz_3_U)
            #-----------------------------------
            #prepare info to output:
            #-----------------------------------
            #nr objects:
            nr_few_body = 3
            #Set simulation time:
            Torb_inibin = 2.*np.pi*np.sqrt((SMA_abin**3.)/(mass_bin))        #orbital time of ini target bin in code units.
            #select and define sim time:
            simtime_U   = Tsim_max_unit_Torb*Torb_inibin
            simtim_secs = 100000000 #just put a very high number. Might be set to remaining sim time at Master process. See mpi modules.        
            #Nbody code settings:
            if (incl_25PN == 1):
                use_25PN        = 1
                insp_threshold  = 5.0   # = a*(Ri+Rj).
                insp_SMA        = (1./25.)*SMA_abin
            if (incl_25PN == 0):
                use_25PN        = 0            
                insp_threshold  = 0.0   # = a*(Ri+Rj).
                insp_SMA        = 0.0
            # [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, n_particles=.., ...]
            nbody_params_arr_1_INT  = np.array([0, use_25PN, 1, 1000000000, 10, 0, nr_few_body, 0,0,0], dtype='i')                           #OUTPUT
            # [scale_dt, max_sim_time, evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in mpi code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, insp_SMA, ...]
            nbody_params_arr_2_REAL = np.array([0.01, simtime_U, 0.0, 0.01, simtim_secs, 0.1, 0.0, insp_threshold, insp_SMA, 0], dtype='d') #OUTPUT
            #object 1,2,3 info:
            b1_const_arr = np.array([mb1, Rb1, 0.0, 0.0, 0.0, 0, 1, 0,0,0], dtype='d')
            b2_const_arr = np.array([mb2, Rb2, 0.0, 0.0, 0.0, 0, 1, 0,0,0], dtype='d')
            b3_const_arr = np.array([mb3, Rb3, 0.0, 0.0, 0.0, 0, 1, 0,0,0], dtype='d')
            b123_obj_info   = np.append(b1_const_arr,np.append(b2_const_arr,b3_const_arr))                                                  #OUTPUT
            #object 1,2,3 pos,vel:
            b123_pos_vel    = np.append(posvelxyz_1_U,np.append(posvelxyz_2_U,posvelxyz_3_U))                                               #OUTPUT
            #-----------------------------------
            #-----------------------------------
            #write IC data to files:
            #-----------------------------------
            np.savetxt(tf_nbody_params_INT,     nbody_params_arr_1_INT[None],   fmt='%10i')     #np.array([0, use_25PN, 1, 1000000000, 0, 0, nr_few_body, 0,0,0], dtype='i')
            np.savetxt(tf_nbody_params_REAL,    nbody_params_arr_2_REAL[None],  fmt='%10f')     #np.array([0.01, simtime_U, 0.0, 0.01, simtim_secs, 0.5, 0.0, insp_threshold, insp_SMA, 0], dtype='d')
            np.savetxt(tf_obj_info,             b123_obj_info[None],            fmt='%10f')     #np.append(b1_const_arr,np.append(b2_const_arr,b3_const_arr))
            np.savetxt(tf_pos_vel,              b123_pos_vel[None],             fmt='%10f')     #np.append(posvelxyz_1_U,np.append(posvelxyz_2_U,posvelxyz_3_U))
            #-----------------------------------                
            #-----------------------------------
            #update total counter:
            #-----------------------------------
            totc = totc + 1
            #-----------------------------------
            
            #print nIC, totc, SMA_abin/AU_U, ecc_ebin, mb1, mb2, mb3            
    
#-------------------------------------------
#-------------------------------------------

#-------------------------------------------
#close files:
#-------------------------------------------
tf_nbody_params_INT.close()   
tf_nbody_params_REAL.close()   
tf_obj_info.close()   
tf_pos_vel.close()   
tf_sim_info_1.close()   
#-------------------------------------------
#create and save info file:
#....
#-------------------------------------------

exit()
#----------------------------------------------------------
#----------------------------------------------------------


























#test IC:
#open file:
fn = 'MCinput_Nbody.txt'
text_file = open(fn, "w")
#nr particles in the sim   
text_file.write('2' + '\n')
#nbody code settings:
nbody_params_arr_1_INT[2] = 0
np.savetxt(text_file, nbody_params_arr_1_INT[None],   fmt='%10i')
np.savetxt(text_file, nbody_params_arr_2_REAL[None],  fmt='%10f')
#body 1:
b1_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b1_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
np.savetxt(text_file, b1_const_arr[None],  fmt='%10f')
np.savetxt(text_file, posxyz_1_U[None],  fmt='%10f')
np.savetxt(text_file, velxyz_1_U[None],  fmt='%10f')
np.savetxt(text_file, b1_q,  fmt='%10f')
np.savetxt(text_file, b1_qdot,  fmt='%10f')
#body 2:
b2_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b2_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
np.savetxt(text_file, b2_const_arr[None],  fmt='%10f')
np.savetxt(text_file, posxyz_2_U[None],  fmt='%10f')
np.savetxt(text_file, velxyz_2_U[None],  fmt='%10f')
np.savetxt(text_file, b2_q,  fmt='%10f')
np.savetxt(text_file, b2_qdot,  fmt='%10f')
#close file:
text_file.close()   
#run:
subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
print mb1, mb2, mb3, SMA_abin/AU_U
print posxyz_1_U
print velxyz_1_U
exit()


#test IC:
#open file:
fn = 'MCinput_Nbody.txt'
text_file = open(fn, "w")
#nr particles in the sim   
text_file.write('3' + '\n')
#nbody code settings:
np.savetxt(text_file, nbody_params_arr_1_INT[None],   fmt='%10i')
np.savetxt(text_file, nbody_params_arr_2_REAL[None],  fmt='%10f')
#body 1:
b1_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b1_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
np.savetxt(text_file, b1_const_arr[None],  fmt='%10f')
np.savetxt(text_file, posxyz_1_U[None],  fmt='%10f')
np.savetxt(text_file, velxyz_1_U[None],  fmt='%10f')
np.savetxt(text_file, b1_q,  fmt='%10f')
np.savetxt(text_file, b1_qdot,  fmt='%10f')
#body 2:
b2_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b2_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
np.savetxt(text_file, b2_const_arr[None],  fmt='%10f')
np.savetxt(text_file, posxyz_2_U[None],  fmt='%10f')
np.savetxt(text_file, velxyz_2_U[None],  fmt='%10f')
np.savetxt(text_file, b2_q,  fmt='%10f')
np.savetxt(text_file, b2_qdot,  fmt='%10f')
#body 3:
b3_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b3_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
np.savetxt(text_file, b3_const_arr[None],  fmt='%10f')
np.savetxt(text_file, posxyz_3_U[None],  fmt='%10f')
np.savetxt(text_file, velxyz_3_U[None],  fmt='%10f')
np.savetxt(text_file, b3_q,  fmt='%10f')
np.savetxt(text_file, b3_qdot,  fmt='%10f')
#close file:
text_file.close()   
#run:
subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
print mb1, mb2, mb3, SMA_abin/AU_U
print posxyz_1_U
print velxyz_1_U
exit()






#posxyz_1_U  = np.array([0.2       , -0.34641016,  0.])
#posxyz_2_U  = np.array([-0.15      ,  0.25980762,  0. ])
#velxyz_1_U  = np.array([1.10656667,  4.47213595,  0.])
#velxyz_2_U  = np.array([-0.829925 , -3.35410197,  0.])
posrel12    = posxyz_1_U - posxyz_2_U
velrel12    = velxyz_1_U - velxyz_2_U
mu12        = mb1*mb2/(mb1 + mb2)
lang        = mu12*np.sqrt(np.sum(np.cross(posrel12, velrel12)**2.))
Eenr        = (1./2.)*mu12*np.sum(velrel12**2.) - mb1*mb2/(np.sqrt(np.sum(posrel12**2.)))
ecal        = np.sqrt(1. + 2.*Eenr*lang**2./(mu12*((mb1 + mb2)*mu12)**2.))
acal        = mb1*mb2/(2.*Eenr)
print ecal
print acal
exit()




















