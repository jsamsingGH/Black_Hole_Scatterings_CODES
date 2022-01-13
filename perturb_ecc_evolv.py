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
import csv

#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
def func_return_binsin_posvel(inputfunc_binsin_posvel):  
#------------------------------------------------------------------------------------------

    #define:
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
b2_mass      = 20.0
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
b3_mass      = 20.0
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
#------------------------------------------------------------------------------------------
#Define:
#------------------------------------------------------------------------------------------
#mass: (all below are in correct code units (Msun))
mass_tot    = b1_mass + b2_mass + b3_mass
mass_bin    = b1_mass + b2_mass
mass_sin    = b3_mass
mass_red    = mass_bin*mass_sin/mass_tot 
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#SIM: simexp1
#------------------------------------------------------------------------------------------
#Set:
sim_YN      = 0
rperi_abin_min      = 3.0
rperi_abin_max      = 20.0    #rememeber to also adjust r_simsurf_abin: do ALWAYS choose r_simsurf_abin >> rperi_abin_max
nr_rperi            = 50
nr_rnd_fbin         = 1
ecc_ini_arr         = [0.95, 0.99, 0.999]#[0.95, 0.99, 0.999]

#define:
nr_ecc_ini          = len(ecc_ini_arr)
rperi_abin_arr      = 10.**(np.linspace(np.log10(rperi_abin_min), np.log10(rperi_abin_max), nr_rperi))
simexp1_data_arr    = np.zeros((nr_ecc_ini, nr_rperi, nr_rnd_fbin, 10), dtype=np.float64)
simexp1_params_arr  = np.zeros((nr_ecc_ini, nr_rperi, nr_rnd_fbin, 10), dtype=np.float64)
#------------------------------------------------------------------------------------------


#Run settings:


#NOTES!!!!!!
#do such that different simexp can be run.
#mention in code that several varible names are being reused.
#sometimes I see points in the plot that does not go through lines. why?
#the slope is still slightly off.
#include 1PN, 2.5PN...
#FUNCTION: the functions we have added rely fully on global variables! everything has to be
#correctly defined before calling the function. 


#------------------------------------------------------------------------------------------
if (sim_YN == 1):
#------------------------------------------------------------------------------------------
#simulate YES:
#------------------------------------------------------------------------------------------  

    save_filename = raw_input('input save_filename: ')
    #---------------------------------------------------------------------
    #START sim:
    #---------------------------------------------------------------------
    for sc_e, sc_r, sc_f in itertools.product(range(0,nr_ecc_ini), range(0,nr_rperi), range(0,nr_rnd_fbin)):
    #---------------------------------------------------------------------

        #-----------------------------------------------------------------
        #binary-single parameters (SMA, vinf, ...):
        #-----------------------------------------------------------------
        print sc_e, sc_r, sc_f
        #---------------------------------------
        #set initial values:
        #---------------------------------------
        #---------------------
        #obj 1,2:
        #---------------------
        SMAbin_AU           = 0.5                   #input in AU
        ecc_ini             = ecc_ini_arr[sc_e]
        #---------------------
        #obj 3:
        #---------------------
        vinf_kmsec          = 0.001                 #input in km/sec
        rperi_abin          = rperi_abin_arr[sc_r]  #unit of SMA
        r_simsurf_abin      = 250.0  #20.0          #unit of SMA
        #---------------------
        #orbital orientations:
        #---------------------
        fbin_01     = 0.0                           #set from [0:1] or CHOOSEOR CHOOSE RANDOM:  = np.random.random_sample()
        w_rot       = 0.0
        i_rot       = np.pi/2.
        Omega_rot   = -np.pi/4.
        #---------------------
        #---------------------------------------
        #---------------------------------------
        #change units/define etc:
        #---------------------------------------
        SMA_bin     = SMAbin_AU*AU_U                                                #in correct code units
        vinf_sin    = vinf_kmsec*kmsec_U                                            #in correct code units
        r_peri      = rperi_abin*SMA_bin                                            #in correct code units
        b_imp       = r_peri*np.sqrt(1.0 + (2.*mass_tot)/(r_peri*(vinf_sin**2.)))   #in correct code units
        r_simsurf   = r_simsurf_abin*SMA_bin                                        #in correct code units
        Torb_inibin = 2.*np.pi*np.sqrt((SMA_bin**3.)/(mass_bin))                    #orbital time of ini target bin in code units.
        #---------------------------------------
        #-----------------------------------------------------------------        
        
        #-----------------------------------------------------------------
        #calc etc.:
        #-----------------------------------------------------------------
        m1_SI   = b1_mass*M_sun_SI
        m2_SI   = b2_mass*M_sun_SI
        m3_SI   = b3_mass*M_sun_SI
        m12_SI  = m1_SI+m2_SI
        m123_SI = m1_SI+m2_SI+m3_SI
        a_SI    = SMA_bin*R_sun_SI
        rp_SI   = r_peri*R_sun_SI
        #analytical de (parabolic limit. eq.8 in heggie&rasio)
        de_analytical   = -(15.*np.pi/16.)*(((2.*(m3_SI**2.)*(a_SI**3.))/(m123_SI*m12_SI*(rp_SI**3.)))**(1./2.))*(ecc_ini*np.sqrt(1.0 - ecc_ini**2.))*np.sin(2.*Omega_rot)*(np.sin(i_rot)**2.)
        epsSA           = (((m3_SI**2.)/(m12_SI*(m12_SI + m3_SI)))*((a_SI/rp_SI)**3.)*((1.+1.)**(-3.)))**(1./2.) 
        de_HS19         = epsSA*(15.*np.pi/4.)*ecc_ini*np.sqrt(1.-ecc_ini**2.)*np.sin(2.*(-Omega_rot))*(np.sin(i_rot)**2.) + (epsSA**2.)*(3.*np.pi*ecc_ini/512.)*(100.*(1.-ecc_ini**2.)*np.sin(2.*(-Omega_rot))*((5.*np.cos(i_rot)+3.*np.cos(3.*i_rot))*np.cos(2.*(-w_rot)) + 6.*np.sin(i_rot)*np.sin(2.*i_rot)) + 4.*np.cos(2.*i_rot)*(3.*np.pi*(81.*(ecc_ini**2.) - 56.) + 200.*(1.-ecc_ini**2.)*np.cos(2.*(-Omega_rot))*np.sin(2.*(-w_rot))) + 3.*np.pi*(200.*(ecc_ini**2.)*(np.sin(i_rot)**4.)*np.cos(4.*(-Omega_rot)) + 8.*(16.*(ecc_ini**2.) + 9.)*(np.sin(2.*i_rot)**2.)*np.cos(2.*(-Omega_rot)) + (39.*(ecc_ini**2.) + 36.)*np.cos(4.*i_rot) - 299.*(ecc_ini**2.) + 124.))
        #-----------------------------------------------------------------
                
        #-----------------------------------------------------------------
        #Code Input/stopping criteria etc.:
        #-----------------------------------------------------------------
        #simulation time:
        pfac    = 2.*r_peri
        Dfac    = np.tan((np.arccos((pfac/r_simsurf) - 1.0))/2.)
        Tpara_U     = np.sqrt((pfac**3.)/(mass_tot))*(Dfac + (1./3.)*(Dfac**3.))
        nbody_params_arr_2_REAL[1]  = Tpara_U
        nbody_params_arr_2_REAL[4]  = 100000000
        #inspiral:
        nbody_params_arr_2_REAL[8]  = insp_SMA_abin*SMA_bin        
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        #ICs: 3-body pos, vel:
        #-----------------------------------------------------------------
        inputfunc_binsin_posvel = [b1_mass, b2_mass, b3_mass, SMA_bin, ecc_ini, fbin_01, i_rot, w_rot, Omega_rot, vinf_sin, b_imp, r_simsurf]
        RF_binsin_posvel        = func_return_binsin_posvel(inputfunc_binsin_posvel)   #function returns: [b1_posxyz_CM, b2_posxyz_CM, b3_posxyz_CM, b1_velxyz_CM, b2_velxyz_CM, b3_velxyz_CM]        
        #from output define:
        [b1_posxyz_CM, b2_posxyz_CM, b3_posxyz_CM, b1_velxyz_CM, b2_velxyz_CM, b3_velxyz_CM] = RF_binsin_posvel[:]
        #-----------------------------------------------------------------
        

        #-----------------------------------------------------------------
        #Run/Analyze:
        #-----------------------------------------------------------------

        #-----------------------------------------------------------------
        #OPTION 1: TEST/PLOT SECTION!
        #-----------------------------------------------------------------
        if (outputinfo_screenfiles == 1):









            #CLEN UP BELOW. EVERYTHING SHOULD BE OK. JUST TRIM SECTION!
            #NOT TESSTT ANYMORE. OK PART.
            #-----------------------------------------------------------------
            #TESSSSSSSTTTTTTTT:
            #-----------------------------------------------------------------
            #---------------------------------------
            #run N-body code:
            #---------------------------------------
            func_runcode()  #function writes IC .txt file and run code.
            #---------------------------------------
            #Read in data from Nbody solver:
            #---------------------------------------
            #open data:
            tf = open('NbodyTides_dataout_pos.dat',  "r")
            NbodyTides_dataout_pos          = np.loadtxt(tf, dtype=float)
            tf.close()
            b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
            b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
            b3_posxyz   = NbodyTides_dataout_pos[:,6:9]
            #open data:
            tf = open('NbodyTides_dataout_vel.dat',  "r")
            NbodyTides_dataout_vel          = np.loadtxt(tf, dtype=float)
            tf.close()
            b1_velxyz   = NbodyTides_dataout_vel[:,0:3]
            b2_velxyz   = NbodyTides_dataout_vel[:,3:6]
            b3_velxyz   = NbodyTides_dataout_vel[:,6:9]
            #open data:
            tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
            NbodyTides_dataout_a1a2a3          = np.loadtxt(tf, dtype=float)
            tf.close()
            time = NbodyTides_dataout_a1a2a3[:,0]
            nr_sim_steps    = len(time)
            #---------------------------------------
            #Define:
            #---------------------------------------
            #object 1:
            xp1 = b1_posxyz[:,0]
            yp1 = b1_posxyz[:,1]
            zp1 = b1_posxyz[:,2]
            #object 2:
            xp2 = b2_posxyz[:,0]
            yp2 = b2_posxyz[:,1]
            zp2 = b2_posxyz[:,2]
            #object 3:
            xp3 = b3_posxyz[:,0]
            yp3 = b3_posxyz[:,1]
            zp3 = b3_posxyz[:,2]
            #other things..:
            sma_a0      = SMA_bin
            M12         = b1_mass + b2_mass
            mu12        = (b1_mass*b2_mass)/M12
            #---------------------------------------
            #---------------------------------------
            #Analyze:
            #---------------------------------------
            #define:
            plot_arr    = np.zeros((nr_sim_steps, 5), dtype=np.float64)
            totc        = 0
            #calc e: (obj 1,2)
            flip_YN     = 0 #0=NO, 1=YES
            flip_pos    = 0
            flip_time   = 0
            for i in range(0,nr_sim_steps,2):
                r_vec   = np.array([b1_posxyz[i,0]-b2_posxyz[i,0], b1_posxyz[i,1]-b2_posxyz[i,1], b1_posxyz[i,2]-b2_posxyz[i,2]])
                v_vec   = np.array([b1_velxyz[i,0]-b2_velxyz[i,0], b1_velxyz[i,1]-b2_velxyz[i,1], b1_velxyz[i,2]-b2_velxyz[i,2]])
                L_vec   = mu12*(np.cross(r_vec,v_vec))
                r_len   = np.sqrt(sum(r_vec**2.))
                v_len   = np.sqrt(sum(v_vec**2.))
                L_len   = np.sqrt(sum(L_vec**2.))
                E_orb   = (1./2.)*mu12*(v_len**2.) - M12*mu12/r_len
                a_orb   = - M12*mu12/(2.*E_orb)
                ecc_orb = np.sqrt(1. + 2.*E_orb*(L_len**2.)/(mu12*((M12*mu12)**2.)))
                rp_orb  = a_orb*(1.0-ecc_orb)    
                de_orb  = ecc_orb-ecc_ini
                #perform check:
                Lz_sign = L_vec[2]/abs(L_vec[2])
                if (flip_YN == 0 and Lz_sign < 0.0):
                    flip_pos    = i    
                    flip_YN     = 1
                    flip_time   = time[i]
                #save:
                plot_arr[totc,0]    = time[i]
                plot_arr[totc,1]    = r_len
                plot_arr[totc,2]    = a_orb
                plot_arr[totc,3]    = ecc_orb
                plot_arr[totc,4]    = rp_orb
                #update counter:
                totc = totc+1
            #trim:
            plot_arr    = plot_arr[0:totc,:]            
            #---------------------------------------
            #PLOT 1:
            #---------------------------------------
            fig1, ax1 = plt.subplots(figsize=(5, 2))
            fig2, ax2 = plt.subplots(figsize=(5, 2))
            
            if (flip_YN == 1):
                posfN   = np.where(plot_arr[:,0] <  flip_time)[0]
                posfY   = np.where(plot_arr[:,0] >= flip_time)[0]
                #flipN:
                xplot   = (plot_arr[posfN,0]-Tpara_U/2.)/Torb_inibin
                #ax1:
                yplot   = plot_arr[posfN,3] 
                ax1.plot(xplot, yplot, color = 'blue', linewidth = 2.0, alpha = 0.75, linestyle = '-')
                #ax2:
                yplot   = np.log10(plot_arr[posfN,1]/sma_a0)
                ax2.plot(xplot, yplot, color = 'blue', linewidth = 2.0, alpha = 0.75, linestyle = '-')
                yplot   = np.log10(plot_arr[posfN,4]/sma_a0)
                ax2.plot(xplot, yplot, color = 'blue', linewidth = 2.0, alpha = 0.75, linestyle = ':')
                yplot   = np.log10(plot_arr[posfN,2]/sma_a0)
                ax2.plot(xplot, yplot, color = 'blue', linewidth = 2.0, alpha = 0.75, linestyle = '--')

                #flipY:
                xplot   = (plot_arr[posfY,0]-Tpara_U/2.)/Torb_inibin
                #ax1:
                yplot   = plot_arr[posfY,3] 
                ax1.plot(xplot, yplot, color = 'orange', linewidth = 2.0, alpha = 0.75, linestyle = '-')
                #ax2:
                yplot   = np.log10(plot_arr[posfY,1]/sma_a0)
                ax2.plot(xplot, yplot, color = 'orange', linewidth = 2.0, alpha = 0.75, linestyle = '-')
                yplot   = np.log10(plot_arr[posfY,4]/sma_a0)
                ax2.plot(xplot, yplot, color = 'orange', linewidth = 2.0, alpha = 0.75, linestyle = ':')
                yplot   = np.log10(plot_arr[posfY,2]/sma_a0)
                ax2.plot(xplot, yplot, color = 'orange', linewidth = 2.0, alpha = 0.75, linestyle = '--')

                
            #axis settings:
            ax1.set_xlim([-8,8])
            ax1.set_xlabel(r'time [$T_{0}$]')
            ax1.set_ylabel(r'$e$')

            #dummy plots for legend:
            ax2.plot([1e10,1e10], [1e10,1e10], color = 'black', linewidth = 2.0, alpha = 0.75, linestyle = '-',     label = r'$|\bf{r}_1 - \bf{r}_2|$')
            ax2.plot([1e10,1e10], [1e10,1e10], color = 'black', linewidth = 2.0, alpha = 0.75, linestyle = ':',     label = r'$a(1-e)$')
            ax2.plot([1e10,1e10], [1e10,1e10], color = 'black', linewidth = 2.0, alpha = 0.75, linestyle = '--',    label = r'$a$')
            ax2.set_xlim([-8,8])
            ax2.set_ylim([-2.5,1.00])
            ax2.set_xlabel(r'time [$T_{0}$]')
            ax2.set_ylabel(r'log $r$ $[a_0]$')
            ax2.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 3, frameon = False)
            
            fig1.savefig(save_filename+'time_ecc.pdf',    bbox_inches='tight')        
            fig2.savefig(save_filename+'time_rdist.pdf',  bbox_inches='tight')        
            
            plt.show()
            #---------------------------------------
            #PLOT 2:
            #---------------------------------------
            fig = plt.figure(figsize=(10, 6))
            ax  = fig.add_subplot(111, projection='3d')
            
            cm_pos12    = (b1_mass*b1_posxyz + b2_mass*b2_posxyz)/(b1_mass+b2_mass)
            cmx = cm_pos12[:,0]
            cmy = cm_pos12[:,1]
            cmz = cm_pos12[:,2]
            if (flip_YN == 1):
                #obj 1:
                xyz = np.array([xp1-cmx, yp1-cmy, zp1-cmz])/sma_a0
                ax.plot(xyz[0][0:flip_pos],xyz[1][0:flip_pos],xyz[2][0:flip_pos],     linewidth = 1.0, color='blue')
                ax.plot(xyz[0][flip_pos::],xyz[1][flip_pos::],xyz[2][flip_pos::],     linewidth = 1.0, color='orange')
                #obj 2:
                xyz = np.array([xp2-cmx, yp2-cmy, zp2-cmz])/sma_a0
                ax.plot(xyz[0][0:flip_pos],xyz[1][0:flip_pos],xyz[2][0:flip_pos],     linewidth = 1.0, color='blue')
                ax.plot(xyz[0][flip_pos::],xyz[1][flip_pos::],xyz[2][flip_pos::],     linewidth = 1.0, color='orange')
                #obj 3:
                xyz = np.array([xp3-cmx, yp3-cmy, zp3-cmz])/sma_a0
                ax.plot(xyz[0][0:flip_pos],xyz[1][0:flip_pos],xyz[2][0:flip_pos],     linewidth = 1.0, color='blue')
                ax.plot(xyz[0][flip_pos::],xyz[1][flip_pos::],xyz[2][flip_pos::],     linewidth = 1.0, color='orange')

            if (flip_YN == 0):
                #obj 1:
                xyz = np.array([xp1-cmx, yp1-cmy, zp1-cmz])/sma_a0
                ax.plot(xyz[0][:],xyz[1][:],xyz[2][:],     linewidth = 1.0, color='black')
                #obj 2:
                xyz = np.array([xp2-cmx, yp2-cmy, zp2-cmz])/sma_a0
                ax.plot(xyz[0][:],xyz[1][:],xyz[2][:],     linewidth = 1.0, color='black')
                #obj 3:
                xyz = np.array([xp3-cmx, yp3-cmy, zp3-cmz])/sma_a0
                ax.plot(xyz[0][:],xyz[1][:],xyz[2][:],     linewidth = 1.0, color='black')

            #axis settings:
            #ax.view_init(azim=-110, elev=65)
            ax.view_init(azim=-105, elev=29)            
            ax.set_xlabel(r'$x/a_0$')
            ax.set_ylabel(r'$y/a_0$')
            ax.set_zlabel(r'$z/a_0$')
            ax.set_xlim(-1, 1)
            ax.set_ylim(-0.1, 0.1)
            ax.set_zlim(-0.1, 0.1)

            plt.savefig(save_filename+'3bodyint_ex1.pdf', bbox_inches='tight')        
            plt.show()
            
            print de_orb
            print de_analytical
            print de_HS19
            exit()
            #---------------------------------------           
            #-----------------------------------------------------------------    
            #END TESSSSSTTTT!!!!
            #-----------------------------------------------------------------






            #make fig:
            fig = plt.figure(figsize=(6, 10))
            ax1 = fig.add_subplot(311)
            ax2 = fig.add_subplot(312)
            ax3 = fig.add_subplot(313)
            
            #below we loop over these PN settings:
            run_PN_param  = [use_25PN]#[0,1]   # = [use_25PN]
            for pnc in range(0,len(run_PN_param)):
                
                #-------------------------------
                #incl. PN setting param:
                #-------------------------------
                PN_param_val    = run_PN_param[pnc]
                nbody_params_arr_1_INT[1]   = PN_param_val  #nbody_params_arr_1_INT = np.array([0, use_25PN, 1, 1000000000, 10, outputinfo_screenfiles, 3, de_analysis,0,0], dtype='i')
                #-------------------------------
                
                #-------------------------------
                #Write param file to Nbody code and run:
                #-------------------------------
                func_runcode()  #function writes IC .txt file and run code.
                #-------------------------------
                
                #-------------------------------
                #open data:
                #-------------------------------
                #open data:            
                tf = open('Nbody_endsim_info_data.txt',  "r")
                fline_split = tf.readline().split()
                endstate_ID = int(fline_split[0])
                tf.close()
                #open data:
                tf = open('NbodyTides_dataout_pos.dat',  "r")
                NbodyTides_dataout_pos          = np.loadtxt(tf, dtype=np.float64)
                tf.close()
                b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
                b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
                b3_posxyz   = NbodyTides_dataout_pos[:,6:9]
                #open data:
                tf = open('NbodyTides_dataout_vel.dat',  "r")
                NbodyTides_dataout_vel          = np.loadtxt(tf, dtype=np.float64)
                tf.close()
                b1_velxyz   = NbodyTides_dataout_vel[:,0:3]
                b2_velxyz   = NbodyTides_dataout_vel[:,3:6]
                b3_velxyz   = NbodyTides_dataout_vel[:,6:9]
                #open data:
                tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
                NbodyTides_dataout_a1a2a3          = np.loadtxt(tf, dtype=np.float64)
                tf.close()
                sim_time    = NbodyTides_dataout_a1a2a3[:,0]
                nr_sim_steps    = len(sim_time)
                print 'stop reading data'
                print 'endstate_ID', endstate_ID
                #-------------------------------
            
                #-------------------------------
                #calc...:
                #-------------------------------
                it      = nr_sim_steps-1    #index
                m12     = b1_mass + b2_mass
                mu12    = (b1_mass*b2_mass)/m12
                r_vec   = np.array([b1_posxyz[it,0]-b2_posxyz[it,0], b1_posxyz[it,1]-b2_posxyz[it,1], b1_posxyz[it,2]-b2_posxyz[it,2]])
                v_vec   = np.array([b1_velxyz[it,0]-b2_velxyz[it,0], b1_velxyz[it,1]-b2_velxyz[it,1], b1_velxyz[it,2]-b2_velxyz[it,2]])
                L_vec   = mu12*(np.cross(r_vec,v_vec))
                r_len   = np.sqrt(sum(r_vec**2.))
                v_len   = np.sqrt(sum(v_vec**2.))
                L_len   = np.sqrt(sum(L_vec**2.))
                E_orb   = (1./2.)*mu12*(v_len**2.) - m12*mu12/r_len
                a_orb   = - m12*mu12/(2.*E_orb)
                ecc_orb = np.sqrt(1. + 2.*E_orb*(L_len**2.)/(mu12*((m12*mu12)**2.)))
                #define:
                #de_simulation   = ecc_orb - ecc_ini
                #Lz_sign         = L_vec[2]/abs(L_vec[2])
                #a_fin           = a_orb
                #e_fin           = ecc_orb            
                #-------------------------------
                
                #-------------------------------
                #MAKE FIGURE:
                #-------------------------------
                #select plot range:
                plot_dTbin  = 10.0
                Tpara_Tbin  = Tpara_U/Torb_inibin
                time_Tbin   = sim_time/Torb_inibin
                pos_plot_GT = np.where(time_Tbin[:] > (Tpara_Tbin/2. - plot_dTbin))[0] 
                pos_plot_LT = np.where(time_Tbin[:] < (Tpara_Tbin/2. + plot_dTbin))[0]                 
                pos_plot    = list(set(pos_plot_GT).intersection(pos_plot_LT))
                #prepare plot data:
                nrppoints   = 1000
                plotevery   = max([1,int((max(pos_plot)-min(pos_plot))/nrppoints)])
                tot_pps     = len(np.arange(min(pos_plot),max(pos_plot),plotevery))
                plot_data_arr   = np.zeros((tot_pps, 10), dtype=np.float64)
                ic          = 0
                #loop and generate/save plot data:
                for sc in range(min(pos_plot),max(pos_plot),plotevery):
                    #calc/define:
                    sim_t   = sim_time[sc]
                    m12     = b1_mass + b2_mass
                    mu12    = (b1_mass*b2_mass)/m12
                    r_vec   = np.array([b1_posxyz[sc,0]-b2_posxyz[sc,0], b1_posxyz[sc,1]-b2_posxyz[sc,1], b1_posxyz[sc,2]-b2_posxyz[sc,2]])
                    v_vec   = np.array([b1_velxyz[sc,0]-b2_velxyz[sc,0], b1_velxyz[sc,1]-b2_velxyz[sc,1], b1_velxyz[sc,2]-b2_velxyz[sc,2]])
                    L_vec   = mu12*(np.cross(r_vec,v_vec))
                    r_len   = np.sqrt(sum(r_vec**2.))
                    v_len   = np.sqrt(sum(v_vec**2.))
                    L_len   = np.sqrt(sum(L_vec**2.))
                    E_orb   = (1./2.)*mu12*(v_len**2.) - m12*mu12/r_len
                    a_orb   = - m12*mu12/(2.*E_orb)
                    ecc_orb = np.sqrt(1. + 2.*E_orb*(L_len**2.)/(mu12*((m12*mu12)**2.)))
                    Lz_sign = L_vec[2]/abs(L_vec[2])
                    #calc T25PN, T1PN, Torb3:
                    m1_SI   = b1_mass*M_sun_SI
                    m2_SI   = b2_mass*M_sun_SI
                    m3_SI   = b3_mass*M_sun_SI
                    m12_SI  = m1_SI+m2_SI
                    m123_SI = m1_SI+m2_SI+m3_SI
                    a_SI    = a_orb*R_sun_SI
                    rp3_SI  = r_peri*R_sun_SI
                    T25PN   = (768./425.)*((a_SI**4.)/(4.*(64./5.)*(G_new_SI**3.)*m1_SI*m2_SI*m12_SI/(c_SI**5.)))*((1. - ecc_orb**2.)**(7./2.))
                    T1PN    = (2.*np.pi/3.)*((c_SI**2.)*(1.-ecc_orb**2.))*(a_SI**(5./2.))/((G_new_SI*m12_SI)**(3./2.))  #CHECK M=m1+m2 MASS TERM!!!!!
                    Torb3   = np.pi*np.sqrt((rp3_SI**3.)/(G_new_SI*m123_SI))                                            #ONLY APPROXIMATE!!!
                    #save data in array:
                    plot_data_arr[ic,0] = sim_t/Torb_inibin - Tpara_Tbin/2.
                    plot_data_arr[ic,1] = ecc_orb
                    plot_data_arr[ic,2] = (a_orb/SMA_bin)
                    plot_data_arr[ic,3] = (a_orb*(1. - ecc_orb))/(SMA_bin*(1. - ecc_ini))
                    plot_data_arr[ic,4] = (T25PN/Torb3)
                    plot_data_arr[ic,5] = (T1PN/Torb3)
                    plot_data_arr[ic,6] = Lz_sign
                    ic = ic+1
                #-------------------------------
                #plot:
                #-------------------------------
                #plot settings:
                Lz_sign_arr = plot_data_arr[:,6]
                posLsGT0    = np.where(Lz_sign_arr[:] >= 0.0)[0]
                posLsLT0    = np.where(Lz_sign_arr[:] <= 0.0)[0]
                posLsLT0    = np.insert(posLsLT0,0,max(posLsGT0))   #add element to connect GT/LT Lsign curves.
                
                #PN: NO
                if (PN_param_val == 0):
                    lstyle      = '-'
                    psymsize    = 0
                #PN: YES    
                if (PN_param_val == 1):
                    lstyle      = '--'
                    psymsize    = 0
                
                #fig 1:
                pos     = posLsGT0
                pcolor  = 'black'
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,1]
                ax1.plot(px, py,   marker = 'o', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                pos     = posLsLT0
                pcolor  = 'red'
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,1]
                ax1.plot(px, py,   marker = 'o', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                #axis settings, etc.:
                ax1.plot([-100, 100], [1,1], marker='', linestyle=':', color='black')
                xdeltax = 5#plot_dTbin
                ax1.set_xlim(-xdeltax, xdeltax)
                ax1.set_ylim(min(plot_data_arr[:,1]) - 0.001, 1.0 + 0.001)
                ax1.set_xlabel(r'time $(t-t_{\rm p})/T_{\rm ini}$')
                ax1.set_ylabel(r'eccentricity $e$')
                
                #fig 2:
                pos     = posLsGT0
                pcolor  = 'black'
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,2]
                ax2.plot(px, py,   marker = 'o', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                pos     = posLsLT0
                pcolor  = 'red'
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,2]
                ax2.plot(px, py,   marker = 'o', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                #axis settings, etc.:
                ax2.plot([-100, 100], [1,1], marker='', linestyle=':', color='black')
                xdeltax = 5#plot_dTbin
                ax2.set_xlim(-xdeltax, xdeltax)
                ax2.set_ylim(0.1,1.1)
                ax2.set_yscale("log")
                ax2.set_xlabel(r'time $(t-t_{\rm p})/T_{\rm ini}$')
                ax2.set_ylabel(r'semi-major axis $a$')
                
                #fig 3:
                pos     = posLsGT0
                pcolor  = 'black'
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,4]
                ax3.plot(px, py,   marker = '*', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,5]
                ax3.plot(px, py,   marker = 'o', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                pos     = posLsLT0
                pcolor  = 'red'
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,4]
                ax3.plot(px, py,   marker = '*', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                px  = plot_data_arr[pos,0]
                py  = plot_data_arr[pos,5]
                ax3.plot(px, py,   marker = 'o', linestyle = lstyle, markersize=psymsize, markeredgewidth=0.0, color=pcolor)
                #axis settings, etc.:
                ax3.text(-4, plot_data_arr[0,4]*(1.+0.2), r'$2.5PN$',   fontsize=12)
                ax3.text(-4, plot_data_arr[0,5]*(1.+0.2), r'$1PN$',     fontsize=12)
                ax3.plot([-100, 100], [1,1], marker='', linestyle=':', color='black')
                xdeltax = 5#plot_dTbin
                ax3.set_xlim(-xdeltax, xdeltax)
                ax3.set_ylim(0.5,1e8)
                ax3.set_yscale("log")
                ax3.set_xlabel(r'time $(t-t_{\rm p})/T_{\rm ini}$')
                ax3.set_ylabel(r'$\tau/T_{\rm 3}$')
                
                #save plot:
                plt.savefig('sim_ill_1.eps', bbox_inches='tight')
                #-------------------------------    
            
            plt.show()
            exit()   
        #-----------------------------------------------------------------    
            
            
            
        #-----------------------------------------------------------------
        #OPTION 2:
        #-----------------------------------------------------------------
        if (outputinfo_screenfiles == 2):
            
            #-------------------------------
            #Write param file to Nbody code and run:
            #-------------------------------
            func_runcode()  #function writes IC .txt file and run code.
            #-------------------------------

            #open data:            
            tf = open('Nbody_endsim_info_data.txt',  "r")
            fline_split = tf.readline().split()
            endstate_ID = int(fline_split[0])
            tf.close()
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
            #-----------------------------------
            #calc de:
            #-----------------------------------
            m12     = b1_mass + b2_mass
            mu12    = (b1_mass*b2_mass)/m12
            r_vec   = np.array([b1_posxyz[0]-b2_posxyz[0], b1_posxyz[1]-b2_posxyz[1], b1_posxyz[2]-b2_posxyz[2]])
            v_vec   = np.array([b1_velxyz[0]-b2_velxyz[0], b1_velxyz[1]-b2_velxyz[1], b1_velxyz[2]-b2_velxyz[2]])
            L_vec   = mu12*(np.cross(r_vec,v_vec))
            r_len   = np.sqrt(sum(r_vec**2.))
            v_len   = np.sqrt(sum(v_vec**2.))
            L_len   = np.sqrt(sum(L_vec**2.))
            E_orb   = (1./2.)*mu12*(v_len**2.) - m12*mu12/r_len
            a_orb   = - m12*mu12/(2.*E_orb)
            ecc_orb = np.sqrt(1. + 2.*E_orb*(L_len**2.)/(mu12*((m12*mu12)**2.)))
            #define:
            de_simulation   = ecc_orb - ecc_ini
            Lz_sign         = L_vec[2]/abs(L_vec[2])
            a_fin           = a_orb
            e_fin           = ecc_orb
            print de_simulation, Lz_sign, endstate_ID
            #-----------------------------------
        #-----------------------------------------------------------------
        
        
        
        
        
        #-----------------------------------------------------------------
        #save data:
        #-----------------------------------------------------------------
        #simexp1_data_arr    = np.zeros((nr_ecc_ini, nr_rperi, nr_rnd_fbin, 10), dtype=np.float64)
        #simexp1_params_arr  = np.zeros((nr_ecc_ini, nr_rperi, nr_rnd_fbin, 10), dtype=np.float64)
        #sim data:
        simexp1_data_arr[sc_e, sc_r, sc_f, 0]   = de_analytical 
        simexp1_data_arr[sc_e, sc_r, sc_f, 1]   = de_simulation
        simexp1_data_arr[sc_e, sc_r, sc_f, 2]   = Lz_sign
        simexp1_data_arr[sc_e, sc_r, sc_f, 3]   = a_fin
        simexp1_data_arr[sc_e, sc_r, sc_f, 4]   = e_fin
        simexp1_data_arr[sc_e, sc_r, sc_f, 5]   = endstate_ID
        simexp1_data_arr[sc_e, sc_r, sc_f, 6]   = de_HS19
        #sim params:
        simexp1_params_arr[sc_e, sc_r, sc_f, 0] = rperi_abin 
        simexp1_params_arr[sc_e, sc_r, sc_f, 1] = ecc_ini
        simexp1_params_arr[sc_e, sc_r, sc_f, 2] = SMA_bin
        simexp1_params_arr[sc_e, sc_r, sc_f, 3] = b1_mass
        simexp1_params_arr[sc_e, sc_r, sc_f, 4] = b2_mass
        simexp1_params_arr[sc_e, sc_r, sc_f, 5] = b3_mass
        #-----------------------------------------------------------------
    #---------------------------------------------------------------------
    #END SIM.
    #---------------------------------------------------------------------
    
    #---------------------------------------------------------------------
    #save data for later use:
    #---------------------------------------------------------------------
    np.save(save_filename + '_' + 'simexp1_data_arr',     simexp1_data_arr,   allow_pickle=True, fix_imports=True)
    np.save(save_filename + '_' + 'simexp1_params_arr',   simexp1_params_arr, allow_pickle=True, fix_imports=True)
    #---------------------------------------------------------------------

#------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------
#PLOT:
#------------------------------------------------------------------------------------------
#THIS sections gets all of its information from the input data files!!!
plot_filename_arr   = ['deDATA_NR2_PNN', 'deDATA_NR2_PNY']#'deDATA_1_25PNY']#['deDATA_1_PNN', 'deDATA_1_25PNY']#['t3N25', 't3Y25']#'t3Y1225']
nr_datasets         = len(plot_filename_arr)

#NOTES:
#check and change ploting settings below for different datasets.
#CHECK AGAIN INT DATA TYPE!!!!! see below
#finish f1ax2 plot.
#MAKE tGW plot.
#make sure we filter out bad endstates. e.g. =13
#CHECK EVERYTHING AGAIN!!!

#fig 1:
fig1 = plt.figure(figsize=(5, 10))
f1ax1 = fig1.add_subplot(311)
f1ax2 = fig1.add_subplot(312)
f1ax3 = fig1.add_subplot(313)


#loop over datasets:
#---------------------------------------------------------------------
for pc in range(0,nr_datasets):
#---------------------------------------------------------------------
    
    plot_filename       = plot_filename_arr[pc]
    simexp1_data_arr    = np.load(plot_filename + '_' + 'simexp1_data_arr.npy')
    simexp1_params_arr  = np.load(plot_filename + '_' + 'simexp1_params_arr.npy')

    #simexp1_data_arr    = np.zeros((nr_ecc_ini, nr_rperi, nr_rnd_fbin, 10), dtype=np.float64)
    #simexp1_params_arr  = np.zeros((nr_ecc_ini, nr_rperi, nr_rnd_fbin, 10), dtype=np.float64)

    #sim data:
    #simexp1_data_arr[sc_e, sc_r, sc_f, 0]   = de_analytical 
    #simexp1_data_arr[sc_e, sc_r, sc_f, 1]   = de_simulation
    #simexp1_data_arr[sc_e, sc_r, sc_f, 2]   = Lz_sign
    #simexp1_data_arr[sc_e, sc_r, sc_f, 3]   = a_fin
    #simexp1_data_arr[sc_e, sc_r, sc_f, 4]   = e_fin
    #simexp1_data_arr[sc_e, sc_r, sc_f, 5]   = endstate_ID
    #simexp1_data_arr[sc_e, sc_r, sc_f, 6]   = de_HS19
    
    ##sim params:
    #simexp1_params_arr[sc_e, sc_r, sc_f, 0] = rperi_abin 
    #simexp1_params_arr[sc_e, sc_r, sc_f, 1] = ecc_ini
    #simexp1_params_arr[sc_e, sc_r, sc_f, 2] = SMA_bin
    #simexp1_params_arr[sc_e, sc_r, sc_f, 3] = b1_mass
    #simexp1_params_arr[sc_e, sc_r, sc_f, 4] = b2_mass
    #simexp1_params_arr[sc_e, sc_r, sc_f, 5] = b3_mass

    #-----------------------------------------------------------------
    #config data and define:
    #-----------------------------------------------------------------
    nr_e            = len(simexp1_params_arr[:,0,0,0])    
    nr_r            = len(simexp1_params_arr[0,:,0,0])    
    endstates_arr   = [5,10]
    nr_endstates    = len(endstates_arr)
    #ENDSTATE CODES: 1 = TIDES, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 5 = Inspiral, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...

    analyze_simexp1_data_arr            = np.zeros((nr_e, nr_r, nr_endstates, 6), dtype=np.float64)
    analyze_simexp1_data_arr[:,:,:,:]   = np.nan    #initialize

    #in this simexp these are constants:
    SMA_ini = simexp1_params_arr[0, 0, 0, 2]
    b1_mass = simexp1_params_arr[0, 0, 0, 3]
    b2_mass = simexp1_params_arr[0, 0, 0, 4]

    #define for calc:
    m1_SI   = b1_mass*M_sun_SI
    m2_SI   = b2_mass*M_sun_SI
    m12_SI  = m1_SI+m2_SI

    for sc_e, sc_r, in itertools.product(range(0,nr_e), range(0,nr_r)):
        #--------------------
        #initial/analytical results:
        #--------------------
        eini    = simexp1_params_arr[sc_e,  0, 0, 1]
        de_ana  = simexp1_data_arr[sc_e, sc_r, 0, 0]
        de_HS19 = simexp1_data_arr[sc_e, sc_r, 0, 6]
        tGW_ini = (768./425.)*(((SMA_ini*R_sun_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m1_SI*m2_SI*m12_SI/(c_SI**5.)))*((1. - eini**2.)**(7./2.))
        #save in new array:  
        analyze_simexp1_data_arr[sc_e,sc_r,:,0] = de_ana
        analyze_simexp1_data_arr[sc_e,sc_r,:,3] = tGW_ini
        analyze_simexp1_data_arr[sc_e,sc_r,:,5] = de_HS19
        
        #--------------------
        #simulation results:
        #--------------------
        esid_er = simexp1_data_arr[sc_e, sc_r, :, 5]    #CHECK AGAIN INT DATA TYPE!!!!!
        for sc_ids in range(0,nr_endstates):
            esid    = endstates_arr[sc_ids]
            pos_es  = np.where(esid_er[:] == esid)[0] 
            if (len(pos_es) > 0):
                #calc, etc.:
                afin_arr    = simexp1_data_arr[sc_e, sc_r, pos_es, 3]
                efin_arr    = simexp1_data_arr[sc_e, sc_r, pos_es, 4]
                tGW_fin_arr = (768./425.)*(((afin_arr*R_sun_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*m1_SI*m2_SI*m12_SI/(c_SI**5.)))*((1. - efin_arr**2.)**(7./2.))
                #define/calc average quanteties:
                de_sim      = np.mean(simexp1_data_arr[sc_e, sc_r, pos_es, 1])  #average over phase sim (f)  
                Lz_sign     = np.mean(simexp1_data_arr[sc_e, sc_r, pos_es, 2])  #average over phase sim (f)  
                tGW_fin     = np.mean(tGW_fin_arr)
                #save in new array:
                analyze_simexp1_data_arr[sc_e,sc_r,sc_ids,1]   = de_sim
                analyze_simexp1_data_arr[sc_e,sc_r,sc_ids,2]   = Lz_sign
                analyze_simexp1_data_arr[sc_e,sc_r,sc_ids,4]   = tGW_fin
        #--------------------
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    #PLOT:
    #-----------------------------------------------------------------
    #set plot settings for dataset 'plot_filename = plot_filename_arr[pc]':
    if (pc == 0):
        ecc_plot_colors = ['orange', 'firebrick', 'royalblue', 'green']
        pd_psymsize     = 5
        pd_linewidth    = 1.5
        pd_fillstyle    = 'none'
    if (pc == 1):
        ecc_plot_colors = ['black', 'black', 'black', 'black']
        pd_psymsize     = 2
        pd_linewidth    = 0.0
        pd_fillstyle    = 'full'
    if (pc == 2):
        ecc_plot_colors = ['pink', 'pink', 'pink', 'pink']
        pd_psymsize     = 4
        pd_linewidth    = 0.0
        pd_fillstyle    = 'none'

    #loop over initial ecc:    
    for sc_e in range(0,nr_e):

        #--------------------
        #define for plots etc.:
        #--------------------
        eini        = simexp1_params_arr[sc_e,0,0,1]    
        scalefac_y  = 1.0#(eini*np.sqrt(1.-eini**2.))
        rp_parr     = simexp1_params_arr[0,:,0,0]
        #--------------------
        
        #--------------------
        #PLOT 1:
        #--------------------
        #Data:
        #--------------------
        #endstate id:
        esid_index  = 1     
        #analyze_simexp1_data_arr = np.zeros((nr_e, nr_r, nr_endstates, 5), dtype=np.float64)
        #analyze_simexp1_data_arr[sc_e,sc_r,:,0] = de_ana
        #analyze_simexp1_data_arr[sc_e,sc_r,sc_ids,1]   = de_sim
        #analyze_simexp1_data_arr[sc_e,sc_r,sc_ids,2]   = Lz_sign
        #analyze_simexp1_data_arr[sc_e,sc_r,:,3] = tGW_ini
        #analyze_simexp1_data_arr[sc_e,sc_r,sc_ids,4]   = tGW_fin 
        #analyze_simexp1_data_arr[sc_e,sc_r,:,5] = de_HS19       
        de_ana  = analyze_simexp1_data_arr[sc_e,:,0,0]
        de_sim  = analyze_simexp1_data_arr[sc_e,:,esid_index,1]
        Lz_sign = analyze_simexp1_data_arr[sc_e,:,esid_index,2]
        de_HS19 = analyze_simexp1_data_arr[sc_e,:,0,5]
        #pos cut:
        pos_deGT0   = np.where(de_sim[:]    > 0.0)[0] 
        pos_deLT0   = np.where(de_sim[:]    < 0.0)[0]
        pos_LzGT0   = np.where(Lz_sign[:]   > 0.0)[0]
        pos_LzLT0   = np.where(Lz_sign[:]   < 0.0)[0]
        #--------------------
        #Figure 1:
        #--------------------
        symevery = 1
        #Simulation:
        #X>0 (GT0):
        #pos = pos_deGT0
        #f1ax1.plot(rp_parr[pos], abs(de_sim[pos])/scalefac_y, linestyle='-', linewidth=pd_linewidth, color=ecc_plot_colors[sc_e])
        pos = pos_LzGT0
        f1ax1.plot(rp_parr[pos], abs(de_sim[pos])/scalefac_y, marker='^', linestyle='', markersize=pd_psymsize, markeredgewidth=1.0, fillstyle = 'full', markeredgecolor = ecc_plot_colors[sc_e], color=ecc_plot_colors[sc_e],markevery=symevery)
        #X<0 (LT0):
        #pos = pos_deLT0
        #f1ax1.plot(rp_parr[pos], abs(de_sim[pos])/scalefac_y, linestyle='--', linewidth=pd_linewidth, color=ecc_plot_colors[sc_e]) 
        pos = list(set(pos_LzLT0).intersection(pos_deGT0))        
        f1ax1.plot(rp_parr[pos], abs(de_sim[pos])/scalefac_y, marker='v', linestyle='', markersize=pd_psymsize, markeredgewidth=1.0, fillstyle = 'full', markeredgecolor = ecc_plot_colors[sc_e], color=ecc_plot_colors[sc_e],markevery=symevery)

        pos = list(set(pos_LzLT0).intersection(pos_deLT0))        
        f1ax1.plot(rp_parr[pos], abs(de_sim[pos])/scalefac_y, marker='v', linestyle='', markersize=pd_psymsize, markeredgewidth=1.0, fillstyle = 'none', markeredgecolor = ecc_plot_colors[sc_e], color=ecc_plot_colors[sc_e],markevery=symevery)
        #Analytical:
        f1ax1.plot(rp_parr, abs(de_ana)/scalefac_y,     linewidth=pd_linewidth, linestyle='-.', color='grey')
        f1ax1.plot(rp_parr, abs(de_HS19)/scalefac_y,    linewidth=pd_linewidth, linestyle=':',  color='grey')
        #--------------------
        #Figure 2:
        #--------------------
        efin    = eini + de_sim
        f1ax2.plot(rp_parr, efin, linestyle='-', linewidth=pd_linewidth, color=ecc_plot_colors[sc_e])
        f1ax2.plot(rp_parr, efin, marker='o', linestyle='', markersize=pd_psymsize, markeredgewidth=1.0, fillstyle = pd_fillstyle, markeredgecolor = ecc_plot_colors[sc_e], color=ecc_plot_colors[sc_e],markevery=symevery)
        #incl guide lines:
        if (pc == 0): f1ax2.plot([1e-10, 1e10], [eini,eini], marker='', linestyle='--', color=ecc_plot_colors[sc_e])
        #--------------------
        #Figure 3:
        #--------------------
        if (pc == 0): p2_linestyle = '-'
        if (pc == 1): p2_linestyle = ''
        
        esid_index  = 0 #(GW merger (5))
        tGW_ini_arr = analyze_simexp1_data_arr[sc_e,:,0, 3]
        tGW_fin_arr = analyze_simexp1_data_arr[sc_e,:,esid_index, 4]
        f1ax3.plot(rp_parr, tGW_fin_arr/tGW_ini_arr, linestyle=p2_linestyle, linewidth=1.5, color=ecc_plot_colors[sc_e])
        f1ax3.plot(rp_parr, tGW_fin_arr/tGW_ini_arr, marker='x', linestyle='', markersize=pd_psymsize, markeredgewidth=1.0, fillstyle = pd_fillstyle, markeredgecolor = ecc_plot_colors[sc_e], color=ecc_plot_colors[sc_e],markevery=symevery)
        
        esid_index  = 1 #(survive (10))   
        tGW_ini_arr = analyze_simexp1_data_arr[sc_e,:,0, 3]
        tGW_fin_arr = analyze_simexp1_data_arr[sc_e,:,esid_index, 4]
        f1ax3.plot(rp_parr, tGW_fin_arr/tGW_ini_arr, linestyle=p2_linestyle, linewidth=1.5, color=ecc_plot_colors[sc_e])
        f1ax3.plot(rp_parr, tGW_fin_arr/tGW_ini_arr, marker='s', linestyle='', markersize=pd_psymsize, markeredgewidth=1.0, fillstyle = pd_fillstyle, markeredgecolor = ecc_plot_colors[sc_e], color=ecc_plot_colors[sc_e],markevery=symevery)
        #--------------------
        
#---------------------------------------------------------------------        
#---------------------------------------------------------------------
#Include data-file data:
#---------------------------------------------------------------------
#datafile_arr    = ['data_e_0.999_large.csv', 'data_e_0.99_large.csv', 'data_e_0.9_large.csv']
#for fc in range(0,len(datafile_arr)):
#    filereader  = csv.reader(open(datafile_arr[fc], 'rb'))
#    for row in filereader:
#        ts  = str(row)
#        ts  = ts[2:len(ts)-2]
#        ts  = (ts.split())
#        rp_data     = float(ts[0]) 
#        de_data     = float(ts[1]) 
#        scalefac_y  = 1.0
#        #f1ax1.plot(rp_data, abs(de_data)/scalefac_y, marker='*', linestyle='', markersize=5, markeredgewidth=1.0, fillstyle = 'none', markeredgecolor = 'yellow', color='yellow', alpha = 0.5)
#---------------------------------------------------------------------
#axis/legend settings:
#---------------------------------------------------------------------
eini_arr    = simexp1_params_arr[:,0,0,1]
rpini_arr   = simexp1_params_arr[0,:,0,0]   

#f1ax1 plot:  
f1ax1.set_xlim(0.9*min(rpini_arr),1.1*max(rpini_arr))
f1ax1.set_ylim(6e-6,1e-1)
f1ax1.set_xscale("log")
f1ax1.set_yscale("log")
f1ax1.set_xlabel(r'')
f1ax1.set_ylabel(r'$|{\delta}e| = |e_{\rm fin} - e_{0}|$')
#dummy plots for legends:
f1ax1.plot(-1,-1, linestyle='-',  linewidth=5.0, alpha = 1.0, color='orange', label=r'$e_0 = 0.95$')
f1ax1.plot(-1,-1, linestyle='-',  linewidth=5.0, alpha = 1.0, color='firebrick', label=r'$e_0 = 0.99$')
f1ax1.plot(-1,-1, linestyle='-',  linewidth=5.0, alpha = 1.0, color='royalblue', label=r'$e_0 = 0.999$')
f1ax1.plot(-1,-1, linestyle='-.', linewidth=2.0, alpha = 1.0, color='grey', label=r'HR96')
f1ax1.plot(-1,-1, linestyle=':',  linewidth=2.0, alpha = 1.0, color='grey', label=r'HS19')
f1ax1.legend(loc='lower left', numpoints = 1, fontsize = 10.0, ncol = 3, frameon = False)

#f1ax2 plot:
#incl guide line:
f1ax2.plot([1e-10, 1e10], [1,1], marker='', linestyle='-', color='black')
f1ax2.set_xlim(0.9*min(rpini_arr),1.1*max(rpini_arr))
f1ax2.set_ylim(min(eini_arr) - 0.005, 1.0 + 0.005)
f1ax2.set_xscale("log")
f1ax2.set_xlabel(r'')
f1ax2.set_ylabel(r'$e_{\rm fin}$')

#f1ax3 plot:
#incl guide line:
f1ax3.plot([1e-10, 1e10], [1,1], marker='', linestyle='--', color='black')
f1ax3.set_xlim(0.9*min(rpini_arr),1.1*max(rpini_arr))
f1ax3.set_xscale("log")
f1ax3.set_yscale("log")
f1ax3.set_xlabel(r'$r_{\rm p}/a_0$')
f1ax3.set_ylabel(r'$\tau_{\rm fin}/\tau_{0}$')


#save figure:
plt.savefig('deefintfin_etc.eps', bbox_inches='tight')        


plt.show()
exit()
#---------------------------------------------------------------------




#------------------------------------------------------------------------------------------
#---------------------------------------
#Read in data from Nbody solver:
#---------------------------------------
#open data:
tf = open('NbodyTides_dataout_pos.dat',  "r")
NbodyTides_dataout_pos          = np.loadtxt(tf, dtype=float)
tf.close()
b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
b3_posxyz   = NbodyTides_dataout_pos[:,6:9]

#open data:
tf = open('NbodyTides_dataout_vel.dat',  "r")
NbodyTides_dataout_vel          = np.loadtxt(tf, dtype=float)
tf.close()
b1_velxyz   = NbodyTides_dataout_vel[:,0:3]
b2_velxyz   = NbodyTides_dataout_vel[:,3:6]
b3_velxyz   = NbodyTides_dataout_vel[:,6:9]

#open data:
tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
NbodyTides_dataout_a1a2a3          = np.loadtxt(tf, dtype=float)
tf.close()
time = NbodyTides_dataout_a1a2a3[:,0]
nr_sim_steps    = len(time)

#---------------------------------------

#---------------------------------------
#Plot:
#---------------------------------------
sma_a0      = SMA_bin

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#make plot window:
fig = plt.figure(figsize=(6, 5))

#object 1:
xp1 = b1_posxyz[:,0]
yp1 = b1_posxyz[:,1]
zp1 = b1_posxyz[:,2]
x = xp1/sma_a0
y = yp1/sma_a0
fig.add_subplot(111).plot(x, y, linewidth=0.5, linestyle='-', color='black')
#object 2:
xp2 = b2_posxyz[:,0]
yp2 = b2_posxyz[:,1]
zp2 = b2_posxyz[:,2]
x = xp2/sma_a0
y = yp2/sma_a0
fig.add_subplot(111).plot(x, y, linewidth=0.5, linestyle='-', color='black')
fig.add_subplot(111).plot(x[0:1], y[0:1], marker = 'o', markersize=5, color='red')
#object 3:
xp3 = b3_posxyz[:,0]
yp3 = b3_posxyz[:,1]
zp3 = b3_posxyz[:,2]
x = xp3/sma_a0
y = yp3/sma_a0
fig.add_subplot(111).plot(x, y, linewidth=0.5, linestyle='-', color='black')



fig = plt.figure(figsize=(10, 8))
#calc e: (obj 1,2)
flip_YN     = 0 #0=NO, 1=YES
flip_pos    = 0
for i in range(0,nr_sim_steps,10):
    M12     = b1_mass + b2_mass
    mu12    = (b1_mass*b2_mass)/M12
    r_vec   = np.array([b1_posxyz[i,0]-b2_posxyz[i,0], b1_posxyz[i,1]-b2_posxyz[i,1], b1_posxyz[i,2]-b2_posxyz[i,2]])
    v_vec   = np.array([b1_velxyz[i,0]-b2_velxyz[i,0], b1_velxyz[i,1]-b2_velxyz[i,1], b1_velxyz[i,2]-b2_velxyz[i,2]])
    L_vec   = mu12*(np.cross(r_vec,v_vec))
    r_len   = np.sqrt(sum(r_vec**2.))
    v_len   = np.sqrt(sum(v_vec**2.))
    L_len   = np.sqrt(sum(L_vec**2.))
    E_orb   = (1./2.)*mu12*(v_len**2.) - M12*mu12/r_len
    a_orb   = - M12*mu12/(2.*E_orb)
    ecc_orb   = np.sqrt(1. + 2.*E_orb*(L_len**2.)/(mu12*((M12*mu12)**2.)))
    rp_orb  = a_orb*(1.0-ecc_orb)
    #CHECK AGAIN!!!!!
    t1PN    = (2.*np.pi/3.)*((c_SI**2.)*(1.-ecc_orb**2.))*((a_orb*R_sun_SI)**(5./2.))/((G_new_SI*(m1_SI+m2_SI))**(3./2.))    #CHECK M=m1+m2 MASS TERM!!!!!
    tb3     = 2.*np.pi*np.sqrt(((r_peri*R_sun_SI)**3.)/(G_new_SI*m123_SI))   #ONLY APPROXIMATE!!!
    
    #print t1PN/tb3
    de_orb  = ecc_orb-ecc_ini
    #perform check:
    Lz_sign = L_vec[2]/abs(L_vec[2])
    if (flip_YN == 0 and Lz_sign < 0.0):
        flip_pos    = i    
        flip_YN     = 1
    #plot:
    t_Torb  = time[i]/Torb_inibin
    if (t_Torb > 24.0 and t_Torb < 26.0):
        fig.add_subplot(221).plot(t_Torb, ecc_orb,                          marker = 'o', markersize=2, markeredgewidth=0.0, color='black')
        fig.add_subplot(222).plot(t_Torb, L_vec[0]/L_len,                   marker = 'o', markersize=2, markeredgewidth=0.0, color='red')
        fig.add_subplot(222).plot(t_Torb, L_vec[1]/L_len,                   marker = 'o', markersize=2, markeredgewidth=0.0, color='blue')
        fig.add_subplot(222).plot(t_Torb, L_vec[2]/L_len,                   marker = 'o', markersize=2, markeredgewidth=0.0, color='green')
        fig.add_subplot(223).plot(t_Torb, rp_orb/(SMA_bin*(1.0-ecc_ini)),   marker = 'o', markersize=2, markeredgewidth=0.0, color='black')
        


print ecc_orb, ecc_ini, L_vec
print 'numerical:', de_orb, ecc_orb
print 'analytical:', de_analytical, ecc_ini+de_analytical, 1.-((ecc_ini+de_analytical)-1.)

plt.show()


#---------------------------------------
#axis settings/labels/ranges, etc:
#---------------------------------------
ax = fig.add_subplot(111)
ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
ax.set_xlim(-5, 5)
ax.set_ylim(-5, 5)

#show:    
plt.show()


fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')
ax.plot(xp1/sma_a0, yp1/sma_a0, zp1/sma_a0,     linewidth = 1.0, color='black')
ax.plot(xp2/sma_a0, yp2/sma_a0, zp2/sma_a0,     linewidth = 1.0, color='black')
ax.plot(xp3/sma_a0, yp3/sma_a0, zp3/sma_a0,     linewidth = 1.0, color='black')
ax.plot(xp3[0:1]/sma_a0, yp3[0:1]/sma_a0, zp3[0:1]/sma_a0, marker = 'o', markersize=5, color='red')

ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
ax.set_zlabel(r'$z/a$')

plt.show()




fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

cm_pos12    = (b1_mass*b1_posxyz + b2_mass*b2_posxyz)/(b1_mass+b2_mass)
cmx = cm_pos12[:,0]
cmy = cm_pos12[:,1]
cmz = cm_pos12[:,2]
if (flip_YN == 1):
    xyz = np.array([xp1-cmx, yp1-cmy, zp1-cmz])/sma_a0
    ax.plot(xyz[0][0:flip_pos],xyz[1][0:flip_pos],xyz[2][0:flip_pos],     linewidth = 1.0, color='black')
    ax.plot(xyz[0][flip_pos::],xyz[1][flip_pos::],xyz[2][flip_pos::],     linewidth = 1.0, color='blue')

    xyz = np.array([xp2-cmx, yp2-cmy, zp2-cmz])/sma_a0
    ax.plot(xyz[0][0:flip_pos],xyz[1][0:flip_pos],xyz[2][0:flip_pos],     linewidth = 1.0, color='black')
    ax.plot(xyz[0][flip_pos::],xyz[1][flip_pos::],xyz[2][flip_pos::],     linewidth = 1.0, color='blue')

    xyz = np.array([xp3-cmx, yp3-cmy, zp3-cmz])/sma_a0
    ax.plot(xyz[0][0:flip_pos],xyz[1][0:flip_pos],xyz[2][0:flip_pos],     linewidth = 1.0, color='black')
    ax.plot(xyz[0][flip_pos::],xyz[1][flip_pos::],xyz[2][flip_pos::],     linewidth = 1.0, color='blue')

if (flip_YN == 0):
    xyz = np.array([xp1-cmx, yp1-cmy, zp1-cmz])/sma_a0
    ax.plot(xyz[0][:],xyz[1][:],xyz[2][:],     linewidth = 1.0, color='black')

    xyz = np.array([xp2-cmx, yp2-cmy, zp2-cmz])/sma_a0
    ax.plot(xyz[0][:],xyz[1][:],xyz[2][:],     linewidth = 1.0, color='black')

    xyz = np.array([xp3-cmx, yp3-cmy, zp3-cmz])/sma_a0
    ax.plot(xyz[0][:],xyz[1][:],xyz[2][:],     linewidth = 1.0, color='black')


ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
ax.set_zlabel(r'$z/a$')

ax.set_xlim(-1, 1)
ax.set_ylim(-0.1, 0.1)
ax.set_zlim(-0.1, 0.1)

plt.show()
#-----------------------------------------------------------------    

 















#print func_return_binsin_posvel(5)
#exit()


#-----------------------------------------------------------------
#Set up 2-body (1,2) system:
#-----------------------------------------------------------------
#find rnd phase f: (see: http://articles.adsabs.harvard.edu/cgi-bin/nph-iarticle_query?1983ApJ...268..319H&defaultprint=YES&filetype=.pdf)
rnd02pi             = (2.*np.pi)*fbin_01
func                = lambda eccE : rnd02pi - (eccE - ecc_ini*np.sin(eccE))    
eccE_initial_guess  = 1.0
sol_eccE            = fsolve(func, eccE_initial_guess)
bin_f               = 2.*np.arctan((((1.+ecc_ini)/(1.-ecc_ini))**(1./2.))*np.tan(sol_eccE/2.))[0]
#calc/define:
Mb12_SI     = mass_bin*M_sun_SI
pfac_SI     = ((SMA_bin*R_sun_SI)*(1.-ecc_ini**2.))
hfac_SI     = (np.sqrt(G_new_SI*Mb12_SI*pfac_SI))
rdist_SI    = (pfac_SI/(1.+ecc_ini*np.cos(bin_f)))
v_ang_SI    = (hfac_SI/rdist_SI)
v_rad_SI    = ((ecc_ini*np.sin(bin_f)/pfac_SI)*hfac_SI)
#calc pos/vel of reduced-mass object:
posx_U      = - ((rdist_SI*np.cos(bin_f))/R_sun_SI)
posy_U      =   ((rdist_SI*np.sin(bin_f))/R_sun_SI)
velx_U      = - ((v_ang_SI*np.sin(bin_f) - v_rad_SI*np.cos(bin_f))*(kmsec_U/1000.))
vely_U      = - ((v_rad_SI*np.sin(bin_f) + v_ang_SI*np.cos(bin_f))*(kmsec_U/1000.))
#define pos/vel vecs:
posvec_U    = np.array([-posx_U, -posy_U, 0.0])    #the '-'signs here change the bin direction such that the 'pericenter vector' points in the pos-x-axis direction.
velvec_U    = np.array([-velx_U, -vely_U, 0.0])    #the '-'signs here change the bin direction such that the 'pericenter vector' points in the pos-x-axis direction.
#change to binary COM:
b1_posxyz_binCM     = - (b2_mass/mass_bin)*posvec_U
b1_velxyz_binCM     = - (b2_mass/mass_bin)*velvec_U
b2_posxyz_binCM     =   (b1_mass/mass_bin)*posvec_U
b2_velxyz_binCM     =   (b1_mass/mass_bin)*velvec_U
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
E_orb               = (1./2.)*mass_red*(vinf_sin**2.)   # = E_tot
L_orb               = abs(mass_red*vinf_sin*b_imp)      # = L tot
ecc_orb             = np.sqrt(1. + 2.*E_orb*(L_orb**2.)/(mass_red*(mass_tot*mass_red)**(2.)))
a_orb               = -mass_tot*mass_red/(2.*E_orb)
theta_orb_A         = np.arccos(-1./ecc_orb)            #theta (asymptotic)    
#--------------------
#set up system where v_inf is parallel to the x-axis:
#--------------------
r_orb               = r_simsurf                         #constant. Not varied.
theta_orb_r         = np.arccos((a_orb*(1.-ecc_orb**2.) - r_orb)/(ecc_orb*r_orb))
phi_orb_r           = theta_orb_A - theta_orb_r    
rx_orb              = r_orb*np.cos(phi_orb_r)
ry_orb              = r_orb*np.sin(phi_orb_r)
v_alpha             = np.arctan(ecc_orb*np.sin(theta_orb_r)/(1.+ecc_orb*np.cos(theta_orb_r)))
v_beta              = phi_orb_r + v_alpha - np.pi/2.
v_orb_r             = np.sqrt(2.*E_orb/mass_red + 2.*mass_tot/r_orb)
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
pos_CM  = (b1_mass*b1_posxyz_binCM + b2_mass*b2_posxyz_binCM + b3_mass*b3_posxyz_binCM)/mass_tot
vel_CM  = (b1_mass*b1_velxyz_binCM + b2_mass*b2_velxyz_binCM + b3_mass*b3_velxyz_binCM)/mass_tot
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


























