from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
import sys
from itertools import combinations
import itertools
from scipy import linalg
from subprocess import call
from mpl_toolkits.mplot3d import Axes3D
import subprocess
from scipy.optimize import fsolve
from scipy import interpolate
import matplotlib as mpl
from matplotlib import colors as ccolor
from matplotlib import cm
import re    
import csv


#----------------------------------------------------------
#Units and conversions:
#----------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI        = 299792458.0       #m/s
M_sun_SI    = 1.989*(10.**30.)  #kg
R_sun_SI    = 695800000.        #m
AU_SI       = 149597871000.     #m 
G_new_SI    = 6.67*(10.**(-11.))
AU_U        = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U     = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
sec_year    = 31536000.
m_parsec    = 3.086*(10**16.)   #m
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
#----------------------------------------------------------


#------------------------------------------------------------------------------------------
def func_return_binsin_posvel(inputfunc_binsin_posvel):  
#------------------------------------------------------------------------------------------

    #define (CU is for 'code units'):
    m1_CU       = inputfunc_binsin_posvel[0]    #b1_mass
    m2_CU       = inputfunc_binsin_posvel[1]    #b2_mass
    m3_CU       = inputfunc_binsin_posvel[2]    #b3_mass
    a_CU        = inputfunc_binsin_posvel[3]    #SMA_bin
    e_in        = inputfunc_binsin_posvel[4]    #ecc_ini
    fbin_01     = inputfunc_binsin_posvel[5]    #fbin_01
    i_rot       = inputfunc_binsin_posvel[6]    #i_rot
    w_rot       = inputfunc_binsin_posvel[7]    #w_rot
    Omega_rot   = inputfunc_binsin_posvel[8]    #Omega_rot
    v3_CU       = inputfunc_binsin_posvel[9]    #vinf_sin
    bimp_CU     = inputfunc_binsin_posvel[10]   #b_imp
    rsim_CU     = inputfunc_binsin_posvel[11]   #r_simsurf
    #calc:
    m12_CU      = m1_CU + m2_CU
    m123_CU     = m1_CU + m2_CU + m3_CU
    mr123_CU    = m12_CU*m3_CU/m123_CU
    
    #-----------------------------------------------------------------
    #Set up 2-body (1,2) system:
    #-----------------------------------------------------------------
    #find rnd phase f: (see: http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1983ApJ...268..319H&defaultprint=YES&filetype=.pdf)
    rnd02pi             = (2.*np.pi)*fbin_01
    func                = lambda eccE : rnd02pi - (eccE - e_in*np.sin(eccE))    
    eccE_initial_guess  = 1.0
    sol_eccE            = fsolve(func, eccE_initial_guess)
    bin_f               = 2.*np.arctan((((1.+e_in)/(1.-e_in))**(1./2.))*np.tan(sol_eccE/2.))[0]
    #calc/define:
    Mb12_SI     = m12_CU*M_sun_SI
    pfac_SI     = ((a_CU*R_sun_SI)*(1.-e_in**2.))
    hfac_SI     = (np.sqrt(G_new_SI*Mb12_SI*pfac_SI))
    rdist_SI    = (pfac_SI/(1.+e_in*np.cos(bin_f)))
    v_ang_SI    = (hfac_SI/rdist_SI)
    v_rad_SI    = ((e_in*np.sin(bin_f)/pfac_SI)*hfac_SI)
    #calc pos/vel of reduced-mass object:
    posx_U      = - ((rdist_SI*np.cos(bin_f))/R_sun_SI)
    posy_U      =   ((rdist_SI*np.sin(bin_f))/R_sun_SI)
    velx_U      = - ((v_ang_SI*np.sin(bin_f) - v_rad_SI*np.cos(bin_f))*(kmsec_U/1000.))
    vely_U      = - ((v_rad_SI*np.sin(bin_f) + v_ang_SI*np.cos(bin_f))*(kmsec_U/1000.))
    #define pos/vel vecs:
    posvec_U    = np.array([-posx_U, -posy_U, 0.0])    #the '-'signs here change the bin direction such that the 'pericenter vector' points in the pos-x-axis direction.
    velvec_U    = np.array([-velx_U, -vely_U, 0.0])    #the '-'signs here change the bin direction such that the 'pericenter vector' points in the pos-x-axis direction.
    #change to binary COM:
    b1_posxyz_binCM     = - (m2_CU/m12_CU)*posvec_U
    b1_velxyz_binCM     = - (m2_CU/m12_CU)*velvec_U
    b2_posxyz_binCM     =   (m1_CU/m12_CU)*posvec_U
    b2_velxyz_binCM     =   (m1_CU/m12_CU)*velvec_U
    #-----------------------------------------------------------------
    #-----------------------------------------------------------------
    #Set up 3-body (1,2 - 3) system:
    #-----------------------------------------------------------------  
    #---------------------------------------
    #incoming single obj3 pos,vel:
    #---------------------------------------
    #--------------------
    #calc/define: (some eqs taken from this site: http://www.braeunig.us/space/orbmech.htm)
    #--------------------
    E_orb               = (1./2.)*mr123_CU*(v3_CU**2.)      # = E_tot
    L_orb               = abs(mr123_CU*v3_CU*bimp_CU)       # = L tot
    ecc_orb             = np.sqrt(1. + 2.*E_orb*(L_orb**2.)/(mr123_CU*(m123_CU*mr123_CU)**(2.)))
    a_orb               = -m123_CU*mr123_CU/(2.*E_orb)
    theta_orb_A         = np.arccos(-1./ecc_orb)            #theta (asymptotic)    
    #--------------------
    #set up system where v_inf is parallel to the x-axis:
    #--------------------
    r_orb               = rsim_CU                           #constant. Not varied.
    theta_orb_r         = np.arccos((a_orb*(1.-ecc_orb**2.) - r_orb)/(ecc_orb*r_orb))
    phi_orb_r           = theta_orb_A - theta_orb_r    
    rx_orb              = r_orb*np.cos(phi_orb_r)
    ry_orb              = r_orb*np.sin(phi_orb_r)
    v_alpha             = np.arctan(ecc_orb*np.sin(theta_orb_r)/(1.+ecc_orb*np.cos(theta_orb_r)))
    v_beta              = phi_orb_r + v_alpha - np.pi/2.
    v_orb_r             = np.sqrt(2.*E_orb/mr123_CU + 2.*m123_CU/r_orb)
    vx_orb_r            = abs(v_orb_r*np.cos(v_beta))  #vx abs val
    vy_orb_r            = abs(v_orb_r*np.sin(v_beta))  #vy abs val
    #pos/vel before rotations:
    b3_posxyz_binCM     = np.array([rx_orb, ry_orb, 0.0])       #not final! we rotate below...
    b3_velxyz_binCM     = np.array([-vx_orb_r, -vy_orb_r, 0.0]) #not final! we rotate below...
    #--------------------
    #rotate to correction configuration:
    #--------------------
    #initialize:
    rotang = (np.pi - theta_orb_A)
    Rz	    = np.array([[np.cos(rotang),-np.sin(rotang),0],[np.sin(rotang),np.cos(rotang),0],[0,0,1]])
    b3_posxyz_binCM     = np.dot(Rz, b3_posxyz_binCM)
    b3_velxyz_binCM     = np.dot(Rz, b3_velxyz_binCM)               
    rotang = -np.pi/2.
    Rz	    = np.array([[np.cos(rotang),-np.sin(rotang),0],[np.sin(rotang),np.cos(rotang),0],[0,0,1]])
    b3_posxyz_binCM     = np.dot(Rz, b3_posxyz_binCM)
    b3_velxyz_binCM     = np.dot(Rz, b3_velxyz_binCM)               
    #rot: w_rot
    rotang  = w_rot
    Rz	    = np.array([[np.cos(rotang),-np.sin(rotang),0],[np.sin(rotang),np.cos(rotang),0],[0,0,1]])
    b3_posxyz_binCM     = np.dot(Rz, b3_posxyz_binCM)
    b3_velxyz_binCM     = np.dot(Rz, b3_velxyz_binCM)               
    #rot: i_rot
    rotang  = -i_rot
    Ry	    = np.array([[np.cos(-rotang),0,np.sin(-rotang)], [0,1,0], [-np.sin(-rotang),0,np.cos(-rotang)]])
    b3_posxyz_binCM     = np.dot(Ry, b3_posxyz_binCM)
    b3_velxyz_binCM     = np.dot(Ry, b3_velxyz_binCM)             
    #initialize:
    rotang  = -np.pi/2.
    Rz	    = np.array([[np.cos(rotang),-np.sin(rotang),0],[np.sin(rotang),np.cos(rotang),0],[0,0,1]])
    b3_posxyz_binCM     = np.dot(Rz, b3_posxyz_binCM)
    b3_velxyz_binCM     = np.dot(Rz, b3_velxyz_binCM)               
    #rot: Omega_rot
    rotang  = Omega_rot
    Rz	    = np.array([[np.cos(rotang),-np.sin(rotang),0],[np.sin(rotang),np.cos(rotang),0],[0,0,1]])
    b3_posxyz_binCM     = np.dot(Rz, b3_posxyz_binCM)
    b3_velxyz_binCM     = np.dot(Rz, b3_velxyz_binCM)               
    #finish rotations.
    #--------------------
    #final output:  b3_posxyz_binCM, b3_velxyz_binCM
    #---------------------------------------
    #---------------------------------------
    #Change pos,vel to 3body system CM:
    #---------------------------------------
    #find pos,vel of CM:
    pos_CM  = (m1_CU*b1_posxyz_binCM + m2_CU*b2_posxyz_binCM + m3_CU*b3_posxyz_binCM)/m123_CU
    vel_CM  = (m1_CU*b1_velxyz_binCM + m2_CU*b2_velxyz_binCM + m3_CU*b3_velxyz_binCM)/m123_CU
    #correct body 1,2,3 pos,vel to CM:
    #pos:
    b1_posxyz_CM = b1_posxyz_binCM - pos_CM
    b2_posxyz_CM = b2_posxyz_binCM - pos_CM
    b3_posxyz_CM = b3_posxyz_binCM - pos_CM
    #vel:
    b1_velxyz_CM = b1_velxyz_binCM - vel_CM
    b2_velxyz_CM = b2_velxyz_binCM - vel_CM
    b3_velxyz_CM = b3_velxyz_binCM - vel_CM
    #---------------------------------------
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    #return:
    #-----------------------------------------------------------------
    return [b1_posxyz_CM, b2_posxyz_CM, b3_posxyz_CM, b1_velxyz_CM, b2_velxyz_CM, b3_velxyz_CM]
    #-----------------------------------------------------------------
#------------------------------------------------------------------------------------------    
    

#------------------------------------------------------------------------------------------
def func_runcode():  
#------------------------------------------------------------------------------------------
    #-----------------------------------------------------------------
    #Write param file to Nbody code and run:
    #-----------------------------------------------------------------
    #open file:
    fn = 'MCinput_Nbody.txt'
    text_file = open(fn, "w")
    #nr particles in the sim   
    text_file.write('3' + '\n')
    #nbody code settings: 
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
    #Run Nbody_AffineTides_solver.exe with 'MCinput_Nbody.txt' as input:
    subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
    #-----------------------------------------------------------------
#------------------------------------------------------------------------------------------


#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)



#------------------------------------------------------------------------------------------
#INPUT SECTION:
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Specify Properties for the 3 objects:
#------------------------------------------------------------------------------------------
#Order: Body 1,2 will be in a binary and Body 3 is incoming.
#---------------------------------------
#IC Obj 1:  (binary 1)
#---------------------------------------
b1_mass      = 20.0
b1_rad       = Rsch_1Msun_unitRsun*b1_mass
b1_gas_n     = 1.5
b1_gas_gamma = 5./3.
b1_Mqp       = 0.1022999*(b1_mass*(b1_rad**2))
b1_evoTides_yesno = 0
b1_RigidSph_yesno = 1
b1_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b1_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
b1_const_arr = np.array([b1_mass, b1_rad, b1_gas_n, b1_gas_gamma, b1_Mqp, b1_evoTides_yesno, b1_RigidSph_yesno, 0,0,0], dtype='d')
#---------------------------------------
#IC Obj 2:  (binary 2)
#---------------------------------------
b2_mass      = b1_mass
b2_rad       = Rsch_1Msun_unitRsun*b2_mass
b2_gas_n     = 3.0
b2_gas_gamma = 5./3.
b2_Mqp       = 0.0376788*(b2_mass*(b2_rad**2))
b2_evoTides_yesno = 0
b2_RigidSph_yesno = 1
b2_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b2_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
b2_const_arr = np.array([b2_mass, b2_rad, b2_gas_n, b2_gas_gamma, b2_Mqp, b2_evoTides_yesno, b2_RigidSph_yesno, 0,0,0], dtype='d')
#---------------------------------------
#IC Obj 3:  (single 3)
#---------------------------------------
b3_mass      = b1_mass
b3_rad       = Rsch_1Msun_unitRsun*b3_mass
b3_gas_n     = 3.0
b3_gas_gamma = 5./3.
b3_Mqp       = 0.0376788*(b3_mass*(b3_rad**2))
b3_evoTides_yesno = 0
b3_RigidSph_yesno = 1
b3_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
b3_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
b3_const_arr = np.array([b3_mass, b3_rad, b3_gas_n, b3_gas_gamma, b3_Mqp, b3_evoTides_yesno, b3_RigidSph_yesno, 0,0,0], dtype='d')
#---------------------------------------
#------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------
#Nbody code settings:
#------------------------------------------------------------------------------------------
#important parameters:
use_12PN                = 0
use_25PN                = 1         #this also controls all 2.5PN effects in the code below incl. GW evol of the BBH etc. set = 0 for no 2.5PN effects at all.
insp_threshold          = 5.0       #a_min = insp_threshold*(Ri+Rj).
insp_SMA_abin           = (1./10.)  #units of a_bin
tidaldisrup_threshold   = 5.0
evolvetides_threshold   = min([0.05, 32.*(b2_mass/b1_mass)*(insp_threshold*(b1_rad+b2_rad)/b1_rad)**(-3.), 32.*(b3_mass/b1_mass)*(insp_threshold*(b1_rad+b3_rad)/b1_rad)**(-3.)]) #works for: TO is only '1' - OR all identical objects. - the 0.05 in front is just to make sure its small.
outputinfo_screenfiles  = 2 #options: 0,1,2
de_analysis             = 1 #otherwise set to 0.
# [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, n_particles=3, de_analysis, ...]
nbody_params_arr_1_INT  = np.array([use_12PN, use_25PN, 1, 1000000000, 10, outputinfo_screenfiles, 3, de_analysis,0,0], dtype='i')
# [scale_dt, max_sim_time (set later in code!!!), evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, insp_SMA (set later in code!!!), ...]
nbody_params_arr_2_REAL = np.array([0.01, -1.0, evolvetides_threshold, -1.0, -1.0, 0.5, tidaldisrup_threshold, insp_threshold, -1.0, 0], dtype='d')
#------------------------------------------------------------------------------------------




#----------------------------------------------------------------------------------------------
#SIM: simexp1
#----------------------------------------------------------------------------------------------
#filename = raw_input('filename: ')

Nsamp       = 10

vdis_kms    = 10.0      #DISPERSION vel
vesc_kms    = 50.0      #ESCAPE     vel
m_Msun      = b1_mass
n_nr_SI     = (10.**(5.))/(m_parsec**3.)

tH_SI       = (10.**10.)*sec_year
Nims        = 20
delta_bs    = 7./9.

vdis_SI = vdis_kms*1000.
vesc_SI = vesc_kms*1000.
m_SI    = m_Msun*M_sun_SI

#calc a_HB:
a_HB_SI     = (3./2.)*(G_new_SI*m_SI/(vdis_SI**2.))
#calc a_ej:
a_ej_SI     = (1./6.)*(1./delta_bs - 1.)*(G_new_SI*m_SI/(vesc_SI**2.))

#model settings:
only_analytical_YN  = 0
de_12ord_term_input = 12 # '=1' for only using 1. ord. '12' for using both 1+2 ord. terms.

#use our hybrid model with 3-body sims:
if (only_analytical_YN  == 0):
    #Encounter zones:
    rpSE_Ua         = 2.0       #SE = Strong Encounter
    rpWS_Ua_min     = 5.0       #WS = Weak-Strong
    rpWW_Ua         = 50.0      #WW = Weak-Weak (maximum pericenter dist.)
#use only analytical models:
if (only_analytical_YN  == 1):
    #Encounter zones:
    rpSE_Ua         = 2.0       #SE = Strong Encounter
    rpWS_Ua_min     = -1        #must be set to -1
    rpWW_Ua         = 50.0      #WW = Weak-Weak (maximum pericenter dist.)

#define:
m123_SI         = m_SI + m_SI + m_SI

#Plotting/Info/Test settings:
make_ill_figure_YN  = 0
#----------------------------------------------------------------------------------------------



#----------------------------------------------------
#PLOT/TEST section:
#----------------------------------------------------
if (make_ill_figure_YN == 1):
    fig = plt.figure(figsize=(8, 10))
    ax1 = fig.add_subplot(311)
    ax2 = fig.add_subplot(312)
    ax3 = fig.add_subplot(313)
#----------------------------------------------------



#----------------------------------------------------------------------------------------------
#run:
#----------------------------------------------------------------------------------------------
for nrs in range(0,Nsamp):

    print 'NOW RUNNING: ', nrs+1, ' out of: ', Nsamp
    #initialize etc.:
    runloop         = 0
    cc              = 0
    enc_type        = 0 
    out_breakID     = 0    
    #initial BBH a,e values:
    a_bin_0         = 0.5*a_HB_SI                               #we stop the code if a>a_HB, so don't choose a0 to be exactly a_HB.
    e_bin_0         = np.sqrt(np.random.random_sample())
    out_info_array  = np.zeros((10000,13), dtype=np.float64)    #temp array used for store data for saving.
    #we save: [a_bin_0, e_bin_0, dtime, t_insp, a_bin_dt, e_bin_dt, rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN]
    
    #START NEW SIMULATION:
    while (runloop == 0):
    

    
        #STATUS: binary state: a_bin_0, e_bin_0



        #-----------------------------------------------------------------------------
        #Initial state: a0,e0,... no breaks!
        #-----------------------------------------------------------------------------
        #------------------------------------------------
        #initilize:
        #------------------------------------------------
        Lz_sign     = 1.0
        endstate_ID = -1
        #------------------------------------------------
        #calc dt, tinsp, etc.:
        #------------------------------------------------
        #time before next encounter:
        rpmax   = rpWW_Ua*a_bin_0
        b_max   = rpmax*np.sqrt(((2.*G_new_SI*m123_SI)/(rpmax*(vdis_SI**2.))))  #we here assume HB limit! If we change it, change also other calcs below incl. rpENC.
        s_max   = np.pi*(b_max**2.)
        dtime   = 1./(n_nr_SI*s_max*vdis_SI)
        #GW inspiral time:
        rp_0        = a_bin_0*(1.-e_bin_0)
        dErp_0      = (85.*np.pi/12.)*(G_new_SI**(7./2.))*(c_SI**(-5.))*(m_SI**(9./2.))*(rp_0**(-7./2.))
        if (use_25PN == 1): t_insp  = np.pi*np.sqrt(2.*G_new_SI)*(m_SI**(3./2.))*np.sqrt(a_bin_0)/dErp_0    #I have checked: this is equivalent to Peters high ecc limit with (1+e) = 2. so OK!
        if (use_25PN == 0): t_insp  = float('inf')        
        #print a_bin_0/AU_SI, e_bin_0, t_insp/sec_year, dtime/sec_year, t_insp/dtime, float('inf')
        #------------------------------------------------
        #-----------------------------------------------------------------------------
    
    
        #-----------------------------------------------------------------------------
        #save info 0: (binary initial state: no breaks)
        #-----------------------------------------------------------------------------
        saveinfo_arr    = [a_bin_0, e_bin_0, dtime, t_insp]
        out_info_array[cc,0:4] = saveinfo_arr[:]
        #we save: [ [a_bin_0, e_bin_0, dtime, t_insp], [a_bin_dt, e_bin_dt], [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN] ]
        print saveinfo_arr
        #-----------------------------------------------------------------------------
    
    
        #-----------------------------------------------------------------------------
        #EVOLVE binary:
        #-----------------------------------------------------------------------------
        #------------------------------------------------
        #BREAK CHECK: BBH inspiral during dt...
        #------------------------------------------------
        #-------------------------------
        #merge: YES
        #-------------------------------
        if (t_insp < dtime):
            out_breakID = 1                                                         #BREAK!!! endstate
            out_encID   = 0
            out_a       = a_bin_0
            out_e       = e_bin_0
            out_rpENC   = -1
            break
        #-------------------------------
        #merge: NO
        #-------------------------------
        if (t_insp > dtime):
            #evolve a,e:
            a_bin_dt    = a_bin_0*(1. - dtime/t_insp)**2.
            e_bin_dt    = (1.0 - rp_0/a_bin_dt) #assumes const. rp (rp_0 = rp_dt) during evolution.
            if (a_bin_dt < rp_0):
                out_breakID = 6                                                     #BREAK!!! endstate
                out_encID   = 0
                out_a       = a_bin_0
                out_e       = e_bin_0
                out_rpENC   = -1
                break
                #NOTE: if at<r0 then our app. breaks down as the BBH circularizes. what do we do? dont happen often, so now
                #just put a break. study the importance later. OK?
        #-------------------------------
        #print a_bin_dt/a_bin_0, t_insp/dtime, a_bin_0/rp_0, a_bin_dt/rp_0, e_bin_dt
        #------------------------------------------------
        #-----------------------------------------------------------------------------
    
    
        #-----------------------------------------------------------------------------
        #save info 1: (binary after GW evolution)
        #-----------------------------------------------------------------------------
        saveinfo_arr    = [a_bin_dt, e_bin_dt]
        out_info_array[cc,4:6] = saveinfo_arr[:]
        #we save: [ [a_bin_0, e_bin_0, dtime, t_insp], [a_bin_dt, e_bin_dt], [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN] ]
        print saveinfo_arr
        #-----------------------------------------------------------------------------



        #STATUS: binary state: a_bin_dt, e_bin_dt
    
    
    
        #-----------------------------------------------------------------------------
        #PERFORM ENCOUNTER:
        #-----------------------------------------------------------------------------
        #------------------------------------------------    
        #encounter params:    
        #------------------------------------------------
        bENC        = b_max*np.sqrt(np.random.random_sample())
        rpENC       = ((bENC**2.)*(vdis_SI**2.))/(2.*G_new_SI*m123_SI)
        rpENC_Ua    = rpENC/a_bin_dt
        #------------------------------------------------
           
        #------------------------------------------------
        #(1) STRONG encounter:
        #------------------------------------------------
        if (rpENC_Ua < rpSE_Ua):
        
            enc_type = 1        #STRONG enc.        (=1)
        
            #-------------------------------
            #CHECK if BBH merges during the interaction:
            #-------------------------------
            if (use_25PN == 0):
                binsin_merge_YN = 0
            #check for binsin merger if use_25PN = 1:
            if (use_25PN == 1):
                #initialize:
                binsin_merge_YN = 0
                #calc minimum ecc for merger (e_cap):
                Rsch_SI = (2.*G_new_SI*m_SI)/(c_SI**2.)
                rcap_SI = 1.8*Rsch_SI*(a_bin_dt/Rsch_SI)**(2./7.)
                e_cap   = 1.0 - (rcap_SI/a_bin_dt)
                #sample BBH ecc and check for merger:
                for nisc in range(0,Nims):
                    if (binsin_merge_YN == 0):
                        bs_ecc  = np.sqrt(np.random.random_sample())
                        #print bs_ecc, e_cap, a_bin_dt
                        if (bs_ecc > e_cap):
                            binsin_merge_YN = 1
                            bs_merg_a   = a_bin_dt
                            bs_merg_e   = bs_ecc
                            break
            #-------------------------------
            #MERGE during interaction:
            #-------------------------------
            if (binsin_merge_YN == 1):
                  out_breakID = 3                                                   #BREAK!!! endstate
                  out_encID   = enc_type
                  out_a       = bs_merg_a
                  out_e       = bs_merg_e
                  out_rpENC   = rpENC_Ua
                  break
            #-------------------------------
            #SURVIVE interaction:
            #-------------------------------
            if (binsin_merge_YN == 0):
                #update a:
                a_bin_enc   = a_bin_dt*delta_bs
                #update e from P(e)=2e:
                e_bin_enc   = np.sqrt(np.random.random_sample())
            #-------------------------------
        #------------------------------------------------
        #END STRONG encounter. out: [a_bin_enc, e_bin_enc]
        #------------------------------------------------
    
        #------------------------------------------------
        #WEAK encounter:
        #------------------------------------------------
        if (rpENC_Ua > rpSE_Ua):
        
            #-------------------------------
            #define:
            #-------------------------------
            deFAC       = -(15.*np.pi/16.)*(((2.*(m_SI**2.))/(m123_SI*(m_SI+m_SI)))**(1./2.))*(e_bin_dt*np.sqrt(1.0-e_bin_dt**2.))
            #-------------------------------
            #-------------------------------
            #determine rpWS_Ua:
            #-------------------------------
            #we consider the case where the sin terms in eq. for de is = 1. so 'de' used below is the maximum possible.
            rpDEeq1E_Ua     = (abs(deFAC)/(1.0-e_bin_dt))**(2./3.)         #rp where abs(de) = 1-e (important for high ecc)
            rpDEeqE_Ua      = (abs(deFAC)/(e_bin_dt))**(2./3.)             #rp where abs(de) = e   (important for low  ecc)
            if (only_analytical_YN == 0): rpWS_Ua   = max([rpWS_Ua_min, rpDEeq1E_Ua, rpDEeqE_Ua])   #choose correct WS limit.
            if (only_analytical_YN == 1): rpWS_Ua   = -1
            #print e_bin_dt, rpENC_Ua, rpDEeq1E_Ua, rpDEeqE_Ua, rpWS_Ua
            #-------------------------------
                
            #-------------------------------
            #(2) WS Encounter: Weak-Strong (run code)
            #-------------------------------
            if (rpENC_Ua < rpWS_Ua):
            
                enc_type = 2    #WEAK-STRONG enc.   (=2)           
            
                #---------------------------
                #define params for code:
                #---------------------------
                r_simsurf_Ua    = 10.0*rpENC_Ua
                #define input. all in CU (code units):
                SMA_bin         = a_bin_dt/R_sun_SI
                ecc_ini         = e_bin_dt
                vinf_sin        = vdis_kms*kmsec_U
                b_imp           = bENC/R_sun_SI
                r_simsurf       = r_simsurf_Ua*SMA_bin
                #stopping criteria: simulation time
                pfac    = 2.*(rpENC_Ua*SMA_bin)
                Dfac    = np.tan((np.arccos((pfac/r_simsurf) - 1.0))/2.)
                m123_CU     = b1_mass + b2_mass + b3_mass
                Tpara_CU    = np.sqrt((pfac**3.)/(m123_CU))*(Dfac + (1./3.)*(Dfac**3.))
                nbody_params_arr_2_REAL[1]  = Tpara_CU  #MAX TIME in code-units
                nbody_params_arr_2_REAL[4]  = 10.0      #MAX TIME in seconds (just choose some reasonable number)
                #stopping criteria: inspiral
                nbody_params_arr_2_REAL[8]  = insp_SMA_abin*SMA_bin        
                #---------------------------       
            
                #---------------------------
                #simulate encounter:
                #---------------------------
                #initialize:
                sim_OK      = 0
                nr_simtries = 0
                #---------------------------
                while (sim_OK == 0):
                #---------------------------
                    nr_simtries = nr_simtries + 1
                    #-----------------------
                    #orbital orientations:
                    #-----------------------
                    fbin_01     = np.random.random_sample()
                    i_rot       = np.arccos(2.*np.random.random_sample()-1.)    #P flat in cos(theta)
                    w_rot       = np.random.random_sample()*(2.*np.pi)
                    Omega_rot   = np.random.random_sample()*(2.*np.pi)               
                    #-----------------------
                    #-----------------------
                    #run 3-body code:
                    #-----------------------
                    inputfunc_binsin_posvel = [b1_mass, b2_mass, b3_mass, SMA_bin, ecc_ini, fbin_01, i_rot, w_rot, Omega_rot, vinf_sin, b_imp, r_simsurf]
                    RF_binsin_posvel        = func_return_binsin_posvel(inputfunc_binsin_posvel)   #function returns: [b1_posxyz_CM, b2_posxyz_CM, b3_posxyz_CM, b1_velxyz_CM, b2_velxyz_CM, b3_velxyz_CM]        
                    #from output define:
                    [b1_posxyz_CM, b2_posxyz_CM, b3_posxyz_CM, b1_velxyz_CM, b2_velxyz_CM, b3_velxyz_CM] = RF_binsin_posvel[:]
                    func_runcode()  #function writes IC .txt file and run code.
                    #-----------------------
                    #-----------------------
                    #analyze/save output from code:
                    #-----------------------
                    #open data:            
                    tf = open('Nbody_endsim_info_data.txt',  "r")
                    endstate_ID             = int(tf.readline().split()[0])                     #endstate id
                    Return_3body_Info_INT   = np.array((tf.readline().split()), dtype=int)      #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno]
                    Return_3body_Info_REAL  = np.array((tf.readline().split()), dtype=float)    #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
                    tf.close()
                    #define:
                    sim_Etotbin     = Return_3body_Info_REAL[2]  
                    sim_Etotbinsin  = Return_3body_Info_REAL[7]
                    #-----------------------
                    #determine if sim is OK/not OK:
                    #-----------------------
                    #CHECK sim OK:
                    if (sim_Etotbin < 0.0) and (endstate_ID == 10 or endstate_ID == 5):
                        sim_OK  = 1
                    #CHECK sim NOT OK:    
                    if (nr_simtries == 5):  #just choose some reasonable number
                        sim_OK  = -1                    
                    #-----------------------                    
                #---------------------------
 
                #---------------------------
                #if sim = NOT OK:
                #---------------------------
                if (sim_OK == -1):
                    out_breakID = -1                                                #BREAK!!! endstate
                    out_encID   = -1
                    out_a       = a_bin_dt
                    out_e       = e_bin_dt
                    out_rpENC   = rpENC_Ua
                    break
                #---------------------------
                
                #---------------------------     
                #if sim = OK:
                #---------------------------
                #not necessary to run this for outID = 5.... but ok
                if (sim_OK == 1):
                    #-----------------------
                    #define:
                    #-----------------------
                    menc        = b1_mass
                    m12         = b1_mass + b2_mass
                    mu12        = (b1_mass*b2_mass)/m12
                    m123        = b1_mass + b2_mass + b3_mass
                    mu123       = m12*b3_mass/m123
                    Etotbin     = sim_Etotbin
                    Etotbinsin  = sim_Etotbinsin
                    vdis_CU     = (vdis_SI/1000.)*kmsec_U
                    n_nr_CU     = n_nr_SI*(R_sun_SI**3.)
                    bin1_id     = Return_3body_Info_INT[1]
                    bin2_id     = Return_3body_Info_INT[2]
                    #-----------------------
                    #determine if sin is bound (binsin_bound_YN = 1/0):
                    #-----------------------
                    if (Etotbinsin > 0.0):
                        binsin_bound_YN = 0
                    if (Etotbinsin < 0.0):
                        a_binsin    = - m123*mu123/(2.*Etotbinsin)
                        Torb_binsin = 2.*np.pi*np.sqrt(a_binsin**3./m123)                                                
                        Tenc_binsin = 1./(n_nr_CU*(np.pi*(a_binsin**2.)*(1. + 2.*(m123+menc)/(a_binsin*vdis_CU**2)))*vdis_CU)
                        if (Torb_binsin > Tenc_binsin): binsin_bound_YN = 0
                        if (Torb_binsin < Tenc_binsin): binsin_bound_YN = 1
                    #out: binsin_bound_YN = 1/0
                    #-----------------------
                    #Binary-single bound to each other: YES
                    #-----------------------
                    if (binsin_bound_YN == 1):
                        da_fin_SI   = 0.0
                        de_fin      = 0.0   
                    #out: [da_fin_SI, de_fin, ...]                     
                    #-----------------------
                    #Binary-single bound to each other: NO
                    #-----------------------
                    if (binsin_bound_YN == 0):
                        #open data:
                        tf = open('NbodyTides_dataout_fin_posvel3N.txt',  "r")
                        NbodyTides_dataout_fin_posvel3N = np.loadtxt(tf, dtype=np.float64)
                        tf.close()
                        b1_posxyz   = NbodyTides_dataout_fin_posvel3N[0,0:3]
                        b2_posxyz   = NbodyTides_dataout_fin_posvel3N[0,3:6]
                        b3_posxyz   = NbodyTides_dataout_fin_posvel3N[0,6:9]
                        b1_velxyz   = NbodyTides_dataout_fin_posvel3N[1,0:3]
                        b2_velxyz   = NbodyTides_dataout_fin_posvel3N[1,3:6]
                        b3_velxyz   = NbodyTides_dataout_fin_posvel3N[1,6:9]
                        b123_posxyz = [b1_posxyz, b2_posxyz, b3_posxyz]
                        b123_velxyz = [b1_velxyz, b2_velxyz, b3_velxyz]
                        #define:
                        bin1_posxyz = b123_posxyz[bin1_id-1]
                        bin2_posxyz = b123_posxyz[bin2_id-1]
                        bin1_velxyz = b123_velxyz[bin1_id-1]
                        bin2_velxyz = b123_velxyz[bin2_id-1]
                        #analyze:
                        r_vec   = np.array([bin1_posxyz[0]-bin2_posxyz[0], bin1_posxyz[1]-bin2_posxyz[1], bin1_posxyz[2]-bin2_posxyz[2]])
                        v_vec   = np.array([bin1_velxyz[0]-bin2_velxyz[0], bin1_velxyz[1]-bin2_velxyz[1], bin1_velxyz[2]-bin2_velxyz[2]])
                        L_vec   = mu12*(np.cross(r_vec,v_vec))
                        r_len   = np.sqrt(sum(r_vec**2.))
                        v_len   = np.sqrt(sum(v_vec**2.))
                        L_len   = np.sqrt(sum(L_vec**2.))
                        E_orb   = (1./2.)*mu12*(v_len**2.) - m12*mu12/r_len
                        a_fin   = - m12*mu12/(2.*E_orb)
                        e_fin   = np.sqrt(1. + 2.*E_orb*(L_len**2.)/(mu12*((m12*mu12)**2.)))
                        #define:
                        da_fin_SI   = (a_fin - SMA_bin)*R_sun_SI
                        de_fin      = e_fin - ecc_ini
                        Lz_sign     = L_vec[2]/abs(L_vec[2])
                    #out: [da_fin_SI, de_fin, ...]      
                    #-----------------------
                    #TEST: Analytical Predictions
                    #-----------------------
                    #assumptions: equal mass and parabolic (E=1)
                    #define:
                    m12_SI  = m_SI+m_SI
                    m3_SI   = m_SI
                    eps_SA  = ((m3_SI**(2.)/(m12_SI*(m12_SI+m3_SI)))*(rpENC_Ua**(-3.))*(1.+1.)**(-3.))**(1./2.)
                    de_1ord_term    = ((eps_SA**1.)*(15.*np.pi/4.)*e_bin_dt*np.sqrt(1.-e_bin_dt**2.)*np.sin(2.*(-Omega_rot))*(np.sin(i_rot))**(2.))
                    de_2ord_term    = ((eps_SA**2.)*(3.*np.pi*e_bin_dt/512.))*(100.*(1.-e_bin_dt**2.)*np.sin(2.*(-Omega_rot))*((5.*np.cos(i_rot)+3.*np.cos(3.*i_rot))*np.cos(2.*(-w_rot)) + 6.*np.sin(i_rot)*np.sin(2.*i_rot)) + 4.*np.cos(2.*i_rot)*(3.*np.pi*(81.*(e_bin_dt**2.) - 56.) + 200.*(1.-e_bin_dt**2.)*np.cos(2.*(-Omega_rot))*np.sin(2.*(-w_rot))) + 3.*np.pi*(200.*(e_bin_dt**2.)*(np.sin(i_rot)**4.)*np.cos(4.*(-Omega_rot)) + 8.*(16.*(e_bin_dt**2.) + 9.)*(np.sin(2.*i_rot)**2.)*np.cos(2.*(-Omega_rot)) + (39.*(e_bin_dt**2.) + 36.)*np.cos(4.*i_rot) - 299.*(e_bin_dt**2.) + 124.))  
                    #Heggie/Rasio result:
                    de_HR           = de_1ord_term 
                    e_bin_enc_WSHR  = e_bin_dt + de_HR                    
                    #Hamers/Samsing:
                    de_HS           = de_1ord_term + de_2ord_term
                    e_bin_enc_WSHS  = e_bin_dt + de_HS    
                    #-----------------------   
                #---------------------------   
            #out: [da_fin_SI, de_fin, ...]        
            #-------------------------------    
            
            #-------------------------------
            #(3) WW Encounter: Weak-Weak (use analytical calc.)
            #-------------------------------
            if (rpENC_Ua > rpWS_Ua):
            
                enc_type = 3    #WEAK-WEAK enc.     (=3)
            
                #---------------------------
                #orbital orientations:
                #---------------------------
                fbin_01     = np.random.random_sample()
                i_rot       = np.arccos(2.*np.random.random_sample()-1.)    #P flat in cos(theta)
                w_rot       = np.random.random_sample()*(2.*np.pi)
                Omega_rot   = np.random.random_sample()*(2.*np.pi)
                #---------------------------
                #calc de:
                #---------------------------                
                #define:
                m12_SI  = m_SI+m_SI
                m3_SI   = m_SI
                eps_SA  = ((m3_SI**(2.)/(m12_SI*(m12_SI+m3_SI)))*(rpENC_Ua**(-3.))*(1.+1.)**(-3.))**(1./2.)
                de_1ord_term        = ((eps_SA**1.)*(15.*np.pi/4.)*e_bin_dt*np.sqrt(1.-e_bin_dt**2.)*np.sin(2.*(-Omega_rot))*(np.sin(i_rot))**(2.))
                #Heggie/Rasio result:
                if (de_12ord_term_input == 1):
                    de_analyt   =  de_1ord_term 
                #Hamers/Samsing:
                if (de_12ord_term_input == 12):
                    de_2ord_term    = ((eps_SA**2.)*(3.*np.pi*e_bin_dt/512.))*(100.*(1.-e_bin_dt**2.)*np.sin(2.*(-Omega_rot))*((5.*np.cos(i_rot)+3.*np.cos(3.*i_rot))*np.cos(2.*(-w_rot)) + 6.*np.sin(i_rot)*np.sin(2.*i_rot)) + 4.*np.cos(2.*i_rot)*(3.*np.pi*(81.*(e_bin_dt**2.) - 56.) + 200.*(1.-e_bin_dt**2.)*np.cos(2.*(-Omega_rot))*np.sin(2.*(-w_rot))) + 3.*np.pi*(200.*(e_bin_dt**2.)*(np.sin(i_rot)**4.)*np.cos(4.*(-Omega_rot)) + 8.*(16.*(e_bin_dt**2.) + 9.)*(np.sin(2.*i_rot)**2.)*np.cos(2.*(-Omega_rot)) + (39.*(e_bin_dt**2.) + 36.)*np.cos(4.*i_rot) - 299.*(e_bin_dt**2.) + 124.))
                    de_analyt   =  de_1ord_term + de_2ord_term                
                #correct for possible analytical-model breakdown:
                e_fin_analyt    = e_bin_dt + de_analyt
                if (e_fin_analyt >= 1.0):
                    de_analyt   = 0.0
                if (e_fin_analyt <= 0.0):
                    de_analyt   = 0.0               
                #define final output da,de:
                da_fin_SI   = 0.0
                de_fin      = de_analyt
                #---------------------------
            #out: [da_fin_SI, de_fin, ...]
            #-------------------------------
        
            #-------------------------------   
            #update a,e:
            #-------------------------------
            a_bin_enc   = a_bin_dt + da_fin_SI
            e_bin_enc   = e_bin_dt + de_fin
            #-------------------------------
        #------------------------------------------------
        #END WEAK encounter. out: [a_bin_enc, e_bin_enc]
        #------------------------------------------------
        
        #-----------------------------------------------------------------------------   
        #END ENCOUNTER. out: [a_bin_enc, e_bin_enc]
        #-----------------------------------------------------------------------------
    
    
    
        #STATUS: binary state: a_bin_enc, e_bin_enc        
        
        
    
        #-----------------------------------------------------------------------------
        #BREAK CHECK: encounter merger/ejection/...
        #-----------------------------------------------------------------------------
        #-------------------------------
        #Endstates:
        #-------------------------------
        #Ejected from STRONG encounter:
        if (enc_type == 1 and a_bin_dt < a_ej_SI):
            out_breakID = 2                                                         #BREAK!!! endstate
            out_encID   = enc_type
            out_a       = a_bin_enc
            out_e       = e_bin_enc
            out_rpENC   = rpENC_Ua
            break
        #Merge during WEAK encounter: (NOTE: in theory we don't need this, but in practice we do!)
        if (enc_type == 2 and endstate_ID == 5):
            out_breakID = 5                                                         #BREAK!!! endstate
            out_encID   = enc_type
            out_a       = a_bin_enc
            out_e       = e_bin_enc 
            out_rpENC   = rpENC_Ua
            break
        #-------------------------------    
        #Model/breakdown checks:
        #-------------------------------
        if (e_bin_enc <= 0.0):
            out_breakID = 10                                                        #BREAK!!!
            out_encID   = enc_type
            out_a       = a_bin_enc
            out_e       = e_bin_enc
            out_rpENC   = rpENC_Ua
            break
        if (e_bin_enc >= 1.0):
            out_breakID = 11                                                        #BREAK!!!
            out_encID   = enc_type
            out_a       = a_bin_enc
            out_e       = e_bin_enc
            out_rpENC   = rpENC_Ua
            break
        if (a_bin_enc > a_HB_SI):
            out_breakID = 12                                                        #BREAK!!!
            out_encID   = enc_type
            out_a       = a_bin_enc
            out_e       = e_bin_enc
            out_rpENC   = rpENC_Ua
            break
        #-------------------------------
        #-----------------------------------------------------------------------------
    
    
        #-----------------------------------------------------------------------------
        #Calc additional info etc.:
        #-----------------------------------------------------------------------------
        #t_insp:
        rp_bin_enc  = a_bin_enc*(1.-e_bin_enc)
        dE_bin_enc  = (85.*np.pi/12.)*(G_new_SI**(7./2.))*(c_SI**(-5.))*(m_SI**(9./2.))*(rp_bin_enc**(-7./2.))
        if (use_25PN == 1): t_insp_enc  = np.pi*np.sqrt(2.*G_new_SI)*(m_SI**(3./2.))*np.sqrt(a_bin_enc)/dE_bin_enc
        if (use_25PN == 0): t_insp_enc  = float('inf')                
        #t binary-single SE:
        rpSE        = rpSE_Ua*a_bin_enc
        bSE         = rpSE*np.sqrt(((2.*G_new_SI*m123_SI)/(rpSE*(vdis_SI**2.))))
        sSE         = np.pi*(bSE**2.)
        dtSE        = 1./(n_nr_SI*sSE*vdis_SI)  #time between strong-encounters (set by rpSE_Ua).
        if (t_insp_enc < dtSE): BBHmerge_SE_YN  = 1
        if (t_insp_enc > dtSE): BBHmerge_SE_YN  = 0
        #-----------------------------------------------------------------------------        
         
            
        #-----------------------------------------------------------------------------
        #save info 2: (bin survives encounter: no break)
        #-----------------------------------------------------------------------------
        saveinfo_arr    = [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN]
        out_info_array[cc,6:13] = saveinfo_arr[:]
        #we save: [ [a_bin_0, e_bin_0, dtime, t_insp], [a_bin_dt, e_bin_dt], [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN] ]
        print saveinfo_arr
        #-----------------------------------------------------------------------------
    
    
        #----------------------------------------------------
        #PLOT/TEST section:
        #----------------------------------------------------
        if (make_ill_figure_YN == 1):
            #calc:
            rp_enc  = a_bin_enc*(1.-e_bin_enc)                                    
            fGW_enc = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp_enc**3.))        
            #------------------------------------------------
            #PLOT:
            #------------------------------------------------
            #(1) STRONG encounter:
            if (enc_type == 1):
                markertype = 'x'
                markersize = 8
                ax1.plot(cc, e_bin_enc,         marker=markertype, color = 'red',               markeredgewidth=2, markersize=1.0*markersize, alpha = 0.75)
                ax1.plot([cc,cc], [0.0, 1.0], linestyle = ':', color = 'red')
                ax2.plot(cc, rpENC_Ua,          marker=markertype, color = 'red',               markeredgewidth=2, markersize=1.0*markersize, alpha = 0.75)
                ax2.plot([cc,cc], [0.0, 1e3], linestyle = ':', color = 'red')
                ax3.plot(cc, fGW_enc,           marker=markertype, color = 'red',               markeredgewidth=2, markersize=1.0*markersize, alpha = 0.75)
                ax3.plot([cc,cc], [0.0, 1e3], linestyle = ':', color = 'red')
        
        
            #(2) WS Encounter: Weak-Strong    
            if (enc_type == 2):
                #standard outcome:
                markertype = 'o'
                markersize = 6
                ax1.plot(cc, e_bin_enc,         marker=markertype, color = 'black',             markeredgewidth=0, markersize=1.0*markersize, alpha = 0.75)
                ax1.plot(cc, e_bin_enc_WSHR,    marker=markertype, markeredgecolor = 'forestgreen', markeredgewidth=1, markersize=1.0*markersize, fillstyle = 'none', alpha = 0.75)
                ax1.plot(cc, e_bin_enc_WSHS,    marker=markertype, markeredgecolor = 'dodgerblue',  markeredgewidth=1, markersize=1.0*markersize, fillstyle = 'none', alpha = 0.75)
                ax2.plot(cc, rpENC_Ua,          marker=markertype, color = 'black',             markeredgewidth=0, markersize=1.0*markersize, alpha = 0.75)
                ax3.plot(cc, fGW_enc,           marker=markertype, color = 'black',             markeredgewidth=0, markersize=1.0*markersize, alpha = 0.75)
                #special outcomes:
                if (Lz_sign < 0.0):
                    ax1.plot(cc, e_bin_enc,     marker='v', markeredgecolor = 'purple', markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none', alpha = 0.75)
                    ax2.plot(cc, rpENC_Ua,      marker='v', markeredgecolor = 'purple', markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none', alpha = 0.75)
                    ax3.plot(cc, fGW_enc,       marker='v', markeredgecolor = 'purple', markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none', alpha = 0.75)
                #if (binsin_bound_YN == 1):
                    #ax1.plot(cc, e_bin_enc,     marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')
                    #ax2.plot(cc, rpENC_Ua,      marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')
                    #ax3.plot(cc, fGW_enc,       marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')
                
                
            #(3) WW Encounter: Weak-Weak
            if (enc_type == 3):
                markertype = 'o'
                markersize = 3
                ax1.plot(cc, e_bin_enc,         marker=markertype, color = 'grey',              markeredgewidth=0, markersize=1.0*markersize, alpha = 0.5)
                ax2.plot(cc, rpENC_Ua,          marker=markertype, color = 'grey',              markeredgewidth=0, markersize=1.0*markersize, alpha = 0.5)
                ax3.plot(cc, fGW_enc,           marker=markertype, color = 'grey',              markeredgewidth=0, markersize=1.0*markersize, alpha = 0.5)
        
        
            #other:
            if (BBHmerge_SE_YN == 1):
                markertype = 'D'
                markersize = 6
                ax1.plot(cc, e_bin_enc,         marker=markertype, markeredgecolor = 'black',   markeredgewidth=1, markersize=markersize, fillstyle = 'none', alpha = 0.5)
                ax2.plot(cc, rpENC_Ua,          marker=markertype, markeredgecolor = 'black',   markeredgewidth=1, markersize=markersize, fillstyle = 'none', alpha = 0.5)
                ax3.plot(cc, fGW_enc,           marker=markertype, markeredgecolor = 'black',   markeredgewidth=1, markersize=markersize, fillstyle = 'none', alpha = 0.5)  
        #ax1.plot(cc, a_bin_enc/a_HB_SI, marker='x', color = 'blue', markersize=5)
        #----------------------------------------------------
    
    
        #-----------------------------------------------------------------------------
        #initialize/update
        #-----------------------------------------------------------------------------
        a_bin_0 = a_bin_enc
        e_bin_0 = e_bin_enc
        cc      = cc+1
        #-----------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------    
    #------------------------------------------------------------------------------------------
    #save info 3: (endstate info)
    #------------------------------------------------------------------------------------------
    saveinfo_arr    = [out_breakID, out_encID, out_a, out_e, out_rpENC]
    out_info_array[cc+1, 0:5] = saveinfo_arr[:]
    #print endstate info:
    print saveinfo_arr
    #------------------------------------------------------------------------------------------
    #save all output data to file:
    #------------------------------------------------------------------------------------------
    #Trim output array:
    out_info_array  = out_info_array[0:cc+2,:] #RESIZE ARRAY TO SAVE SPACE!
    #this is what we save:
    #in code:       [ [a_bin_0, e_bin_0, dtime, t_insp], [a_bin_dt, e_bin_dt], [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN] ]
    #after break:   [ [out_breakID, out_encID, out_a, out_e, out_rpENC] ]
    #Save info array to file:
    save_file_dir   = '/Users/jsamsing/Desktop/TIDES_PROJ/BBHpert_files/'
    #save_file_dir   = '/scratch/gpfs/jsamsing/' #TIGER CLUSTER
    file_name       = 'testLV' +  str(int(1000000000*np.random.random_sample()))
    #TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_
    tf = open(save_file_dir+file_name, "w")
    np.savetxt(tf, out_info_array,   fmt='%8f')
    tf.close()
    #------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------










#----------------------------------------------------
#PLOT/TEST section:
#----------------------------------------------------
if (make_ill_figure_YN == 1):

    #ax1:
    #final output point:
    markertype = '*'
    markersize = 10
    ax1.plot(cc, out_e, linewidth = 0.0, marker=markertype, color = 'black', markeredgewidth=2, markersize=1.0*markersize, alpha = 1, label='endstate')
    ax1.plot([0.0,10000], [1.0, 1.0], linestyle = '-', color = 'black')
    #dummy plots for legend:
    markertype = 'x'
    markersize = 8
    ax1.plot(-1,-1,     linewidth = 0.0, marker=markertype, color = 'red',                   markeredgewidth=2, markersize=1.0*markersize, alpha = 0.75, label=r'SE')
    markertype = 'o'
    markersize = 6
    ax1.plot(-1, -1,    linewidth = 0.0, marker=markertype, color = 'black',                 markeredgewidth=0, markersize=1.0*markersize, alpha = 0.75, label=r'WS, num.')
    ax1.plot(-1, -1,    linewidth = 0.0, marker=markertype, markeredgecolor = 'forestgreen', markeredgewidth=1, markersize=1.0*markersize, fillstyle = 'none', alpha = 0.75, label=r'WS, HR96')
    ax1.plot(-1, -1,    linewidth = 0.0, marker=markertype, markeredgecolor = 'dodgerblue',  markeredgewidth=1, markersize=1.0*markersize, fillstyle = 'none', alpha = 0.75, label=r'WS, HS19')
    ax1.plot(-1, -1,    linewidth = 0.0, marker='v', markeredgecolor = 'purple',             markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none', alpha = 0.75, label=r'orbit flip')
    markertype = 'o'
    markersize = 3
    ax1.plot(-1, -1,    linewidth = 0.0, marker=markertype, color = 'grey',                  markeredgewidth=0, markersize=1.0*markersize, alpha = 0.5, label=r'WW, analytical.')
    markertype = 'D'
    markersize = 6
    ax1.plot(-1, -1,    linewidth = 0.0, marker=markertype, markeredgecolor = 'black',       markeredgewidth=1, markersize=markersize, fillstyle = 'none', alpha = 0.5, label=r'$\tau < t_{\rm SE}$')
    #axis settings:
    ax1.set_xlabel(r'scattering counter')
    ax1.set_ylabel(r'eccentricity $e$')
    ax1.set_xlim(0, 1.05*cc)
    ax1.set_ylim(0, 1.3)
    ax1.legend(loc='upper left', numpoints = 1, fontsize = 10.0, ncol = 4, frameon = False)

    #ax2:
    markertype = '*'
    markersize = 10
    ax2.plot(cc, out_rpENC, marker=markertype, color = 'black', markeredgewidth=2, markersize=1.0*markersize, alpha = 1)
    ax2.plot([0,10000], [rpSE_Ua,rpSE_Ua],             linestyle = '--', color = 'red')
    ax2.plot([0,10000], [rpWS_Ua_min,rpWS_Ua_min],     linestyle = '--', color = 'black')
    rect1 = mpl.patches.Rectangle((0.0,1e-3), 10000, rpSE_Ua-1e-3, color='red', alpha = 0.25)
    ax2.add_patch(rect1)
    rect1 = mpl.patches.Rectangle((0.0,rpSE_Ua), 10000, rpWW_Ua-rpSE_Ua, color='grey', alpha = 0.25)
    ax2.add_patch(rect1)
    ax2.text(10, 0.4, r'STRONG', verticalalignment='center', fontsize=12)    
    ax2.text(10, 10, r'WEAK', verticalalignment='center', fontsize=12)    
    ax2.set_xlabel(r'scattering counter')
    ax2.set_ylabel(r'$r_{\rm p}/a$')
    ax2.set_xlim(0, 1.05*cc)
    ax2.set_ylim(1e-1, 2*rpWW_Ua)
    ax2.set_yscale("log")

    #ax3:
    markertype = '*'
    markersize = 10
    rp_out  = out_a*(1.-out_e)                                    
    fGW_out = (1./np.pi)*np.sqrt(2.*G_new_SI*m_SI/(rp_out**3.))        
    ax3.plot(cc, fGW_out,   marker=markertype, color = 'black', markeredgewidth=2, markersize=1.0*markersize, alpha = 1)
    rect1 = mpl.patches.Rectangle((0.0,1e-3), 10000, 1e-1-1e-3, color='orange')
    ax3.add_patch(rect1)
    ax3.text(10, 1e-2, r'LISA sensitivity band', verticalalignment='center', fontsize=12)    
    ax3.set_xlabel(r'scattering counter')
    ax3.set_ylabel(r'$f_{\rm GW}$')
    ax3.set_xlim(0, 1.05*cc)
    ax3.set_ylim(1e-9, 1e0)
    ax3.set_yscale("log")

    #plot and save fig:
    plt.savefig('BBHpert_test1_fig.pdf', bbox_inches='tight')        
    plt.show()
    
    exit()
#----------------------------------------------------







exit()








#NOTES
#strong encounter zone: looking at the topology paper we know that resint can happen for all
#rp/a < 3. So the limit should be around 3. use fig. 2 in paper and 3b^2/8 = rp/a. (green bands in fig. 2 ill problem: will cont. for ever.)
#I still see a few where sin is bound for rp/a<3.5, so maybe best to choose 3.5.
#CHECK new implemantion of HR and HS in other code. show de as a function of r_p to check.
#new effects:
#BBH can be scattered in and out of high ecc state. Sometimes it leads to merger, sometimes it prevents merger.
#a BBH can evolve through ... 

#MAKE MORE DETAILED GUIDE BELOW!!!!!
#run cluster:
#move all files to new folder and compile
#change settings in BBH_incl_distpert.py and output folder and filename
#change parallel.cmd + folder path in parallel.cmd + ... etc.
#....
#Submit (run)!




#print info:
#t_orb   = (2.*np.pi)*np.sqrt((a_bin_dt**3.)/(G_new_SI*(2.*m_SI)))/sec_year
#t_GR    = (2.*np.pi/3.)*((c_SI**2.)*(1.-e_bin_dt**2.))*(a_bin_dt**(5./2.))/((G_new_SI*(2.*m_SI))**(3./2.))/sec_year
#t_enc   = (rpENC/vdis_SI)/sec_year
#print t_orb, t_GR, t_enc





