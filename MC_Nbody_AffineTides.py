from matplotlib import pyplot as plt
import numpy as np
import sys
from itertools import combinations
import itertools
from scipy import linalg
from subprocess import call
import subprocess
from sympy import *
from sympy.solvers import solve
from sympy import Symbol
from sympy.abc import x
import matplotlib as mpl
import matplotlib.gridspec as gridspec

#NOTES:
#Mqp fac should be as a function of n_gas.
#what do we do about the 1,2 PN terms?





#------------------------------------------------------------------------------------------
#User input:
#------------------------------------------------------------------------------------------
#Allow usr input: (if run on laptop: ok(=1), if run on e.g. cluster not ok(=0))
allow_terminal_usr_input = 1
#output file names:
if (allow_terminal_usr_input == 1): add_filename = raw_input('add this filename to all MC output files: ')
if (allow_terminal_usr_input == 0): add_filename = 'TigerTest3'
#resonance analysis yes/no:
if (allow_terminal_usr_input == 1): resonance_analysis_yesno = int(raw_input('do resonance analysis yes(1)/no(0): '))
if (allow_terminal_usr_input == 0): resonance_analysis_yesno = 0
#ouput folder:
MCoutput_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/TEST_MC_OUT/' #'/scratch/network/jsamsing/Stellar_Tides/'
#------------------------------------------------------------------------------------------




#polytropic table (From P. Diener): http://adsabs.harvard.edu/abs/1995MNRAS.275..498D
#n      M/(mR^2)
#0.0    0.2
#0.5    0.1629656
#1.0    0.1306910
#1.5    0.1022999
#2.0    0.0774246
#2.5    0.0559025
#3.0    0.0376788
#polytropic values: (WIKI)
# - Neutron stars are well modeled by polytropes with index about in the range between n = 0.5 and n = 1.
# - A polytrope with index n = 1.5 is a good model for fully convective star cores[3][4](like those of red giants), brown dwarfs, giant gaseous planets (like Jupiter), or even for rocky planets.
# - Main sequence stars like our Sun and relativistic degenerate cores like those of white dwarfs are usually modeled by a polytrope with index n = 3, corresponding to the Eddington standard model of stellar structure.




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
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
#------------------------------------------------------------------------------------------
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
b1_mass      = 10.
b1_rad       = Rsch_1Msun_unitRsun*b1_mass #0.013*((1.43/b1_mass)**(1./3.))*((1.-b1_mass/1.43)**(0.447))
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
b2_mass      = 10.
b2_rad       = Rsch_1Msun_unitRsun*b2_mass #Rsch_1Msun_unitRsun*b2_mass     #1.7246e-5
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
b3_mass      = 5.
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
use_25PN                = 1
insp_threshold          = 100.0   # = a*(Ri+Rj).
tidaldisrup_threshold   = 5.0
evolvetides_threshold   = min([0.05, 32.*(b2_mass/b1_mass)*(insp_threshold*(b1_rad+b2_rad)/b1_rad)**(-3.), 32.*(b3_mass/b1_mass)*(insp_threshold*(b1_rad+b3_rad)/b1_rad)**(-3.)]) #works for: TO is only '1' - OR all identical objects. - the 0.05 in front is just to make sure its small.
# [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, ...]
nbody_params_arr_1_INT  = np.array([0, use_25PN, 1, 1000000000, 4, 1, 0,0,0,0], dtype='i')
# [scale_dt, max_sim_time (set later in code!!!), evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, ...]
nbody_params_arr_2_REAL = np.array([0.01, -1.0, evolvetides_threshold, 0.01, -1.0, 0.1, tidaldisrup_threshold, insp_threshold,0,0], dtype='d')
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Make arrays with parameters we vary in this MC (SMA, vinf, ...):
#------------------------------------------------------------------------------------------
#---------------------------------------
#semi-major axis (SMA):
#---------------------------------------
#---------------------
#Sequence of values lin in logspace:
#---------------------
#min_SMA_AU      = 1e-4           #input in AU
#max_SMA_AU      = 1e-3           #input in AU     
#nr_SMA          = 3             #incl min, max vals
#MCvalues_SMA_arr  = AU_U*np.logspace(np.log10(min_SMA_AU),np.log10(max_SMA_AU),nr_SMA)                    #output units are in correct code units.
#---------------------
#A few (>=1) user defined values:
#---------------------
SMA_AU_arr          = [0.1]     #input in AU
MCvalues_SMA_arr    = AU_U*np.array(SMA_AU_arr)
nr_SMA              = len(MCvalues_SMA_arr)
#---------------------
#calc corresponding charac vel (vc) per SMA in arr:
#---------------------
vc_kmsec_perSMA_arr = (1./kmsec_U)*np.sqrt((b1_mass*b2_mass*(b1_mass+b2_mass+b3_mass)/(b3_mass*(b1_mass+b2_mass)))*(1./MCvalues_SMA_arr))
print vc_kmsec_perSMA_arr
#---------------------
#---------------------------------------
#incoming vel of body 3 at infinity (vinf):
#---------------------------------------
#---------------------
#Sequence of values lin in logspace:
#---------------------
#min_vinf_kmsec  = 0.5*vc_kmsec_perSMA_arr[0]  
#max_vinf_kmsec  = 8.0*vc_kmsec_perSMA_arr[0]
##min_vinf_kmsec  = 30.0       #input in km/sec
##max_vinf_kmsec  = 100.0      #input in km/sec
#nr_vinf         = 20         #incl min, max vals 
#MCvalues_vinf_arr = kmsec_U*np.logspace(np.log10(min_vinf_kmsec),np.log10(max_vinf_kmsec),nr_vinf)        #output units are in correct code units.
#---------------------
#A few (>=1) user defined values:
#---------------------
vinf_kmsec_arr      = [0.0001]     #input in km/sec
MCvalues_vinf_arr   = kmsec_U*np.array(vinf_kmsec_arr)
nr_vinf             = len(MCvalues_vinf_arr)
#info: (13.3 AU for 10 km/sec)
#info: (53.3 AU for 5 km/sec) 
#info: (1330 AU for 1 km/sec) 
#print (1./AU_U)*((b1_mass*b2_mass*(b1_mass+b2_mass+b3_mass)/(b3_mass*(b1_mass+b2_mass)))*(1./(MCvalues_vinf_arr**2.)))
#exit()
#---------------------
#---------------------
#opt4: phase-space analysis 1: (see http://adsabs.harvard.edu/abs/1983AJ.....88.1549H)
#---------------------
phase_space_analysis_1_yesno = 0    #ANALYSIS YES/NO.
nr_b_imp            = 50    #set res
nr_f_phase          = 1     #set res
sim_dist_unitSMA    = 20.   #P. Hut: R = 40*r, a=2r, => R = 20*a 
rpmax_unitSMA       = 3.0
plane_gamma         = 0.*np.pi    # = 0 is in same plane. Positive val rotates in positive angular direction. This angle is fixed for a given sim.
bmax_imp_unitSMA    = (np.sqrt(2.*(b1_mass+b2_mass+b3_mass)*(rpmax_unitSMA*MCvalues_SMA_arr[0])/(MCvalues_vinf_arr[0]**2.)))/MCvalues_SMA_arr[0]
#--------------
#b,f full scan:
#--------------
b_imp_unitSMA_arr   = (1./(vinf_kmsec_arr/vc_kmsec_perSMA_arr))*np.linspace(-1.1, -0.9, nr_b_imp)         #impact param b range (set number in front)
f_phase_rad_arr     = np.linspace(2.9, 2.9, nr_f_phase)              #relative phase. use full circle: 2pi.     
#--------------
#b,f zoom scan:
#--------------
#b_imp_unitSMA_arr   = np.linspace(-1.0, -0.5, nr_b_imp)                          #impact param b range (set number in front)
#f_phase_rad_arr     = np.linspace(2.6, 2.6, nr_f_phase)                        #relative phase. use full circle: 2pi.    
#--------------
nr_tot_ph1          = nr_b_imp*nr_f_phase
bf_sampl_arr        = np.array(list(itertools.product(b_imp_unitSMA_arr, f_phase_rad_arr)))
#---------------------
#FIGURE EX:
# 202020,1e-4, 1, 50, -1.1,-0.9,2.9,2.9
#---------------------------------------
#final array with all param combinations: (in correct code units)
#---------------------------------------
MCvalues_allparams_arr    = np.array(list(itertools.product(MCvalues_SMA_arr, MCvalues_vinf_arr)))
nr_comb_allparams         = len(MCvalues_allparams_arr)
#format: MCvalues_allparams_arr[nr combinations(nr = nr_comb_allparams), param values(SMA, vinf, ...)]. The array values are in correct code units.
#---------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#3body scattering settings:
#------------------------------------------------------------------------------------------
#maximum pericenter dist in units of bin SMA:
#print "CHECK RP MAX. HAS NOW BEEN SET FROM TIDAL THRESHOLD"
#rp_max_FtidFbin                 = 0.1   # = Ftid/Fbin (=1 will give max(r_p) = 1.26*a in equal mass case, and = 2.0*a for [0.2,1.4,1.4])
#rp_max_unit_SMAbin              = ((1./rp_max_FtidFbin)*((b1_mass+b2_mass)*b3_mass/(b1_mass*b2_mass)))**(1./3.) 
#print rp_max_unit_SMAbin
rp_max_unit_SMAbin              = b2_mass/(b1_mass + b2_mass)   #0.5
#the rsim surface is where the FtidFbin frac equals this number:
FtidFbin_r_simsurf              = 0.001  #0.1 #For high mass bin set this to lower val.
#nr of rnd samplings around the sphere of the bin per bin-sin conf:
nr_scatterings_per_paramcomb    = 1000
#SPECIAL VALUES FOR PHASE-SPACE ANALYSIS:
if (phase_space_analysis_1_yesno == 1):
    nr_scatterings_per_paramcomb    = nr_tot_ph1
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Time settings:
#------------------------------------------------------------------------------------------
print 'REMEMBER to set time settings in code!!!!'
#Choose which time limit we use (1 = use):
use_time_secs   = 0
use_time_binT0  = 1
#if we 'use_time_secs' we need to specify:
Tsim_max_secs       = 20.
#if we 'use_time_binT0' we need to specify:
Tsim_max_unit_Torb  = 50.0
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


test_pos3_arr = np.zeros((nr_scatterings_per_paramcomb,5), dtype=np.float64)


#------------------------------------------------------------------------------------------
#Define:
#------------------------------------------------------------------------------------------
#mass:
mass_tot    = b1_mass + b2_mass + b3_mass
mass_bin    = b1_mass + b2_mass
mass_sin    = b3_mass
mass_red    = mass_bin*mass_sin/mass_tot 
#counters/numbers:
tot_nr_scatterings  = nr_comb_allparams*nr_scatterings_per_paramcomb
totsc  = 0
#arrays:
MCout_arr_endstate_info_INT   = np.zeros((tot_nr_scatterings,10), dtype=np.int)
MCout_arr_endstate_info_REAL  = np.zeros((tot_nr_scatterings,10), dtype=np.float64)
MCout_arr_scattering_params   = np.zeros((tot_nr_scatterings,10), dtype=np.float64)
#if we do resonance analysis:
if (resonance_analysis_yesno == 1):
    maxnr_binsin_changes    = 100  #max nr of bin changes in a single res interaction (incl. creation, detroy bin, etc..)
    Resonance_data_iniarr   = np.zeros((maxnr_binsin_changes*tot_nr_scatterings,13), dtype=np.float64)
    totnr_res_counter       = 0 #counter
#------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------
#Print info:
#------------------------------------------------------------------------------------------
#nr combinations:
print '======================================================'
print 'MC SETTINGS:'
print '======================================================'
print 'nr_SMA, nr_vinf, nr_comb_allparams, nr_scatterings_per_paramcomb, tot_nr_scatterings'
print nr_SMA, nr_vinf, nr_comb_allparams, nr_scatterings_per_paramcomb, tot_nr_scatterings
#dist from bin to sin at beginning of sim (r_simsurf):
r_simsurf_FtideFbin_unit_SMAbin = ((1./FtidFbin_r_simsurf)*(2.)*((b1_mass+b2_mass)*b3_mass/(b1_mass*b2_mass)))**(1./3.)
r_simsurf_rp_max_unit_SMAbin    = rp_max_unit_SMAbin
r_simsurf_unit_SMAbin           = max([r_simsurf_FtideFbin_unit_SMAbin, r_simsurf_rp_max_unit_SMAbin])
print 'r_simsurf_unit_SMAbin'
print r_simsurf_unit_SMAbin
print 'vc_kmsec_perSMA_arr'
print vc_kmsec_perSMA_arr
print 'we will scan over the following values in [SMA, vinf] (code units):'
print MCvalues_allparams_arr
print 'b1_evoTides_yesno, b2_evoTides_yesno, b3_evoTides_yesno'
print b1_evoTides_yesno, b2_evoTides_yesno, b3_evoTides_yesno
print 'min_SMA, max_SMA'
print min(MCvalues_SMA_arr), max(MCvalues_SMA_arr)
print '======================================================'
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Run MC yes/no:
#------------------------------------------------------------------------------------------
if (allow_terminal_usr_input == 1): runMC_yesno = raw_input('Run MC with these settings? (yes=1, no=0): ')
if (allow_terminal_usr_input == 0): runMC_yesno = 1
if (int(runMC_yesno) == 0): exit()
#------------------------------------------------------------------------------------------






#------------------------------------------------------------------------------------------
#RUN MCMC:
#------------------------------------------------------------------------------------------
for pc in range(0,nr_comb_allparams):

    #-----------------------------------------------------------------
    #Get binary-single parameters (SMA, vinf, ...):
    #-----------------------------------------------------------------
    #Binary SMA
    SMA_bin   = MCvalues_allparams_arr[pc,0]  #already in correct code units
    #vel at infinity of single:
    vinf_sin  = MCvalues_allparams_arr[pc,1]  #already in correct code units
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    #Set simulation time:
    #-----------------------------------------------------------------
    Torb_inibin = 2.*np.pi*np.sqrt((SMA_bin**3.)/(mass_bin))    #orbital time of ini target bin in code units.
    #select and define sim time for this a,v comb scattering:
    if (use_time_secs   == 1):
        nbody_params_arr_2_REAL[1] = 100000000*Torb_inibin          #just put a very high number.
        nbody_params_arr_2_REAL[4] = Tsim_max_secs                  #set real limit                   
    if (use_time_binT0  == 1):
        nbody_params_arr_2_REAL[1] = Tsim_max_unit_Torb*Torb_inibin #set real limit
        nbody_params_arr_2_REAL[4] = 100000000                      #just put a very high number.
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    #Scattering settings:
    #-----------------------------------------------------------------
    #max peri-center dist:
    rp_max      = rp_max_unit_SMAbin*SMA_bin    
    #from rp_max we can now calc max impact param (b):
    b_max       = rp_max*np.sqrt(1.+(2.*mass_tot)/(rp_max*(vinf_sin**2.)))
    #distance to sampling surface:
    r_sampsurf  = 1e20*SMA_bin   #just make sure to put in a very big number.
    #find distance to simulation surface (we analytically propagate sin body3 from sampsurf to simsurf to save sim time):
    r_simsurf_FtideFbin     = ((1./FtidFbin_r_simsurf)*(2.*SMA_bin**3.)*((b1_mass+b2_mass)*b3_mass/(b1_mass*b2_mass)))**(1./3.)
    r_simsurf_rp_max        = rp_max
    r_simsurf               = max([r_simsurf_FtideFbin, r_simsurf_rp_max]) # at rsim the FtidFbin will always be < or = FtidFbin_r_simsurf.
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    #Set body1, body2 up in a binary:
    #-----------------------------------------------------------------
    #calc orbital params in CM of binary:
    v_redmass         = np.sqrt(mass_bin/SMA_bin)
    b1_posxyz_binCM   = np.array([ (b2_mass/mass_bin)*SMA_bin,0,0])
    b2_posxyz_binCM   = np.array([-(b1_mass/mass_bin)*SMA_bin,0,0])
    b1_velxyz_binCM   = np.array([0, (b2_mass/mass_bin)*v_redmass,0])
    b2_velxyz_binCM   = np.array([0,-(b1_mass/mass_bin)*v_redmass,0])
    #-----------------------------------------------------------------

    #-----------------------------------------------------------------
    #Calc and define:
    #-----------------------------------------------------------------
    #Total energy of bin-sin system (where bin is a single par with mass Mbin)
    Etot_binsin = (1./2.)*mass_red*(vinf_sin**2.)
    #make rnd numbers for the rnd incoming scatterings of body 3 below:
    rndnum_arr = np.random.random((nr_scatterings_per_paramcomb,4))
    #-----------------------------------------------------------------
    
    #-----------------------------------------------------------------
    #LOOP over incoming scatterings from rnd angles:
    #-----------------------------------------------------------------
    for sc in range(0, nr_scatterings_per_paramcomb):
        
        #---------------------------------------
        #rnd sampling parameters:
        #---------------------------------------
        #ISOTROPIC:
        b_sampsurf          = b_max*np.sqrt(rndnum_arr[sc,0])	        #P flat in b^2			- rad pos on 0->bmax disc
        psi_sampsurf        = (2.*np.pi)*rndnum_arr[sc,1]			    #P flat in psi			- angle on b disc.
        theta               = np.arccos(2.*rndnum_arr[sc,2]-1.)         #P flat in cos(theta)	- pos on unitsphere
        phi                 = ((2.*np.pi)*rndnum_arr[sc,3])			    #P flat in phi			- pos on unitsphere
        #FIXED ANGLE (or CO-PLANAR) CONFIGURATION: 
        print 'FIXED ANGLE NOW CHOSEN!!!!'
        opening_dtheta_DG   = 5.0
        opening_dtheta_RD   = (opening_dtheta_DG/360.)*(2.*np.pi)       #Fixed angle.
        b_sampsurf          = b_max*(2.*rndnum_arr[sc,0] - 1.0)	        
        psi_sampsurf        = np.pi/2.			   
        theta               = np.pi/2.-opening_dtheta_RD          
        phi                 = ((2.*np.pi)*rndnum_arr[sc,3])                                 
        #---------------------------------------
        
        			    
        #everything rotates according to the "right hand rule".
        
                                  
        #---------------------------------------
        #calc:
        #---------------------------------------
        rb_sampsurf         = np.sqrt(r_sampsurf**2 + b_sampsurf**2)    #dist to point 'b' at sampling surface.
        v_sampsurf          = np.sqrt((2./mass_red)*(Etot_binsin + mass_tot*mass_red/rb_sampsurf))
        #---------------------------------------
                        
        #---------------------------------------
        #(b,psi,v) at simulation surface (r_simsurf):
        #---------------------------------------
        #By conservation of E,L we now propagate v,b onto the rsim surface:
        v_simsurf           = np.sqrt((2./mass_red)*(Etot_binsin + mass_tot*mass_red/r_simsurf))
        b_simsurf           = b_sampsurf*v_sampsurf/v_simsurf
        #angle dist is conserved onto rsim surface:
        psi_simsurf         = psi_sampsurf
        #---------------------------------------
                
        #---------------------------------------
        #Send obj3 in from a pos at the r_sim sphere:
        #---------------------------------------
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
        b3_posxyz_binCM = pos_RyRz              #final pos
        #Final vel: we know the vel is by constr parallel with the final pos simplane unitvec:
        unitvec_RyRz_simplane   = np.array([np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta)])
        b3_velxyz_binCM         = - v_simsurf*unitvec_RyRz_simplane 
        #Final output is: b3_posxyz_binCM, b3_velxyz_binCM
        #---------------------------------------
        
        test_pos3_arr[sc,0] = 360.*(theta/(2.*np.pi))
        test_pos3_arr[sc,1] = 360.*(phi/(2.*np.pi))
        
        #---------------------------------------
        #phase-space analysis 1:
        #---------------------------------------
        #initialize: (we need to define these variables for later saving -- even when the analysis below is not active)
        R_val   = 0.0
        b_val   = 0.0
        f_val   = 0.0
        if (phase_space_analysis_1_yesno == 1):
            #mass_tot    = b1_mass + b2_mass + b3_mass
            #mass_bin    = b1_mass + b2_mass
            #mass_sin    = b3_mass
            #mass_red    = mass_bin*mass_sin/mass_tot         
            #read param vals:
            R_val               = sim_dist_unitSMA*SMA_bin          #constant. Not varied.
            b_val               = bf_sampl_arr[sc,0]*SMA_bin        #b at infinity. correct units
            f_val               = bf_sampl_arr[sc,1]                #binary phase in radians
            #calc: (some eqs taken from this site: http://www.braeunig.us/space/orbmech.htm)
            E_orb               = (1./2.)*mass_red*(vinf_sin**2.)   # = E_tot
            L_orb               = abs(mass_red*vinf_sin*b_val)      # = L tot
            e_orb               = np.sqrt(1. + 2.*E_orb*(L_orb**2.)/(mass_red*(mass_tot*mass_red)**(2.)))
            a_orb               = - mass_tot/(vinf_sin**2.)         # = -Mm_red/(2E_tot)
            theta_orb_A         = np.arccos(-1./e_orb)                                                                  #theta (asymptotic)         
            #--------------------
            #if R_val = r_x (as P. Hut does it)
            #--------------------
            #theta_orb_r         = nsolve(1. + e_orb*cos(x) - (a_orb/R_val)*(1.-e_orb**2.)*cos(theta_orb_A - x), 2.0)    #theta (rad dist r). last val(=2) is just initial guess.
            #phi_orb_r           = theta_orb_A - theta_orb_r    
            #r_orb               = R_val/cos(phi_orb_r)
            #--------------------
            #if R_val = r_orb (much faster since we don't need to call nsolve)
            #--------------------
            theta_orb_r         = np.arccos((a_orb*(1.-e_orb**2.) - R_val)/(e_orb*R_val))
            phi_orb_r           = theta_orb_A - theta_orb_r    
            r_orb               = R_val
            #--------------------
            rx_orb              = r_orb*cos(phi_orb_r)
            ry_orb              = r_orb*sin(phi_orb_r)
            v_alpha             = atan(e_orb*sin(theta_orb_r)/(1.+e_orb*cos(theta_orb_r)))
            v_beta              = phi_orb_r + v_alpha - np.pi/2.
            v_orb_r             = sqrt(2.*E_orb/mass_red + 2.*mass_tot/r_orb)
            vx_orb_r            = abs(v_orb_r*cos(v_beta))  #vx abs val
            vy_orb_r            = abs(v_orb_r*sin(v_beta))  #vy abs val
            #b3 incoming pos,vel:
            if (b_val > 0.0):
                b3_posxyz_binCM     = np.array([rx_orb, ry_orb, 0.0])
                b3_velxyz_binCM     = np.array([-vx_orb_r, -vy_orb_r, 0.0])
            if (b_val == 0.0):
                b3_posxyz_binCM     = np.array([rx_orb, 0.0, 0.0])
                b3_velxyz_binCM     = np.array([-vx_orb_r, 0.0, 0.0])
            if (b_val < 0.0):
                b3_posxyz_binCM     = np.array([rx_orb, -ry_orb, 0.0])
                b3_velxyz_binCM     = np.array([-vx_orb_r, vy_orb_r, 0.0])
            #rotate out of plane (angle gamma):
            Ry	= np.array([[np.cos(-plane_gamma),0,np.sin(-plane_gamma)], [0,1,0], [-np.sin(-plane_gamma),0,np.cos(-plane_gamma)]])
            b3_posxyz_binCM         = np.dot(Ry, b3_posxyz_binCM)   #rotate around y (out of plane)
            b3_velxyz_binCM         = np.dot(Ry, b3_velxyz_binCM)   #rotate around y (out of plane)            
            #b1,b2 binary pos, vel:   
            #initial pos/vel:
            v_redmass           = np.sqrt(mass_bin/SMA_bin)
            b1_posxyz_binCM     = np.array([ (b2_mass/mass_bin)*SMA_bin,0,0])
            b2_posxyz_binCM     = np.array([-(b1_mass/mass_bin)*SMA_bin,0,0])
            b1_velxyz_binCM     = np.array([0, (b2_mass/mass_bin)*v_redmass,0])
            b2_velxyz_binCM     = np.array([0,-(b1_mass/mass_bin)*v_redmass,0])
            #now rotate according to phase f:         
            Rot_M	            = np.array([[np.cos(f_val),-np.sin(f_val),0],[np.sin(f_val),np.cos(f_val),0],[0,0,1]])
            b1_posxyz_binCM     = np.dot(Rot_M, b1_posxyz_binCM)
            b2_posxyz_binCM     = np.dot(Rot_M, b2_posxyz_binCM)
            b1_velxyz_binCM     = np.dot(Rot_M, b1_velxyz_binCM)
            b2_velxyz_binCM     = np.dot(Rot_M, b2_velxyz_binCM)              
            #Final output MUST BE: b3_posxyz_binCM
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
    
        #---------------------------------------
        #Write param file to Nbody code:
        #---------------------------------------
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
        #---------------------------------------
    
        #---------------------------------------
        #Run Nbody_AffineTides_solver.exe with 'MCinput_Nbody.txt' as input:
        #---------------------------------------
        subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
        #---------------------------------------
    
        #exit()
    
        #---------------------------------------
        #Read end-state info from Nbody code and save:
        #---------------------------------------
        #open file:
        tf = open('Nbody_endsim_info_data.txt', "r")
        #From: Nbody_endsim_info_data.txt
        #Nbodyout_endsim_end_state_flag
        #endsim_Return_Info_arr_INT
        #endsim_Return_Info_arr_REAL(:)     = [1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
        #!xtra info:
        #endsim_out_xtra_info_REAL(1:3)		= [rmin_12, rmin_13, rmin_23]
        #endsim_out_xtra_info_REAL(4:5)		= [MRini_IMS_a, MRini_IMS_e]
        #endsim_out_xtra_info_REAL(6:8)		= Return_3body_Info_REAL_XTRAINFO(14:16)	!out_pos_CMij
        #endsim_out_xtra_info_REAL(9)		= endsim_TOUT		
        #!xtra info 2:
        #endsim_out_xtra_2_info_REAL(1:6)	= Return_3body_Info_REAL_XTRAINFO(5:10)
        #endsim_out_xtra_2_info_REAL(7:9)	= Return_3body_Info_REAL_XTRAINFO(11:13)	!Lvec_ij
        #endsim_out_xtra_2_info_REAL(10)	= MRini_IMS_dist_k           
        #!xtra info 3:
        #endsim_out_xtra_3_info_REAL(1:3)	= MRini_IMS_pos_CMij
        #endsim_out_xtra_3_info_REAL(4:6)	= MRini_IMS_vel_CMij
        #endsim_out_xtra_3_info_REAL(7)		= MRini_IMS_time
        #endsim_out_xtra_3_info_REAL(8:10)	= Return_3body_Info_REAL_XTRAINFO(17:19)	!out_vel_CMij
        
        #read: endsim_end_state_flag
        fline_split = tf.readline().split()    
        Nbodyout_endsim_end_state_flag          = int(fline_split[0])
        #ENDSTATE CODES: 1 = TIDES, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...
        
        #read: endsim_Return_Info_arr_INT
        fline_split = tf.readline().split()    
        Nbodyout_endsim_Return_Info_arr_INT     = np.array([int(i)      for i in fline_split[0:10]])    #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, ...]
        
        #read: endsim_Return_Info_arr_REAL
        fline_split = tf.readline().split()    
        Nbodyout_endsim_Return_Info_arr_REAL    = np.array([float(i)    for i in fline_split[0:10]])    #[ [E_kin, E_pot, E_tot, a_bin, e_bin] (bin[i,j]), [E_kin, E_pot, E_tot, a_bin, e_bin] ([bin[i,j],sin[k]]) ]
        
        #read: endsim_out_xtra_2_info_REAL
        fline_split = tf.readline().split()    
        Nbodyout_endsim_out_xtra_2_info_REAL    = np.array([float(i)    for i in fline_split[0:10]])    #[out_pos_CMij_wrt_sink(:), out_vel_CMij_wrt_sink(:), Lvec_ij(:)]
        
        #read: endsim_out_xtra_3_info_REAL
        fline_split = tf.readline().split()    
        Nbodyout_endsim_out_xtra_3_info_REAL    = np.array([float(i)    for i in fline_split[0:10]])    #...
        
        #close file:    
        tf.close()              
        #---------------------------------------
        #Define:
        #---------------------------------------        
        out_bin_i           = Nbodyout_endsim_Return_Info_arr_INT[1]
        out_bin_j           = Nbodyout_endsim_Return_Info_arr_INT[2]
        IMS_rp_counter      = Nbodyout_endsim_Return_Info_arr_INT[5]
        IMS_binsin_counter  = Nbodyout_endsim_Return_Info_arr_INT[6]        
        Lij_vec             = Nbodyout_endsim_out_xtra_2_info_REAL[6:9]
        v_ij_ims_0          = Nbodyout_endsim_out_xtra_3_info_REAL[3:6] 
        v_ij_end_1          = Nbodyout_endsim_out_xtra_3_info_REAL[7:10]
        angle_v01           = np.arccos(np.dot(v_ij_ims_0,v_ij_end_1)/(np.linalg.norm(v_ij_ims_0)*np.linalg.norm(v_ij_end_1)))
        angle_v01_deg       = angle_v01*(360./(2.*np.pi))
        Lij_zang_DG         = 360.*((np.arccos(np.dot(Lij_vec, np.array([0,0,1]))/np.linalg.norm(Lij_vec)))/(2.*np.pi))

        print 'Lij_zang_DG', Lij_zang_DG
        print 'angle_v01_deg', angle_v01_deg
        
        
        

        
        
        #if (Nbodyout_endsim_end_state_flag == 2): test = input('push test')
        
        #test = input('push test')
        
        
        
        
        #print Nbodyout_endsim_Return_Info_arr_INT[:]
        #if (Nbodyout_endsim_end_state_flag == 1 and IMS_binsin_counter > 3 and IMS_rp_counter > 3 and out_bin_i == 1):
        #    exit()
        
        #exit()
        
        #pres1 = raw_input('press 1')

        #if (Nbodyout_endsim_end_state_flag == 3 and out_bin_j == 3):
        #    exit()
        
        #if (Nbodyout_endsim_end_state_flag == 5 and out_bin_i != 1):    # and IMS_binsin_counter > 1 and IMS_rp_counter > 1): # and out_bin_i == 1):
            #exit()

        if (Nbodyout_endsim_end_state_flag == 5):    # and IMS_binsin_counter > 1 and IMS_rp_counter > 1): # and out_bin_i == 1):
            print 'endstate 5'
            print 'angle_v01_deg', angle_v01_deg
            if (angle_v01_deg > 45.0):
                exit()

        #if (Nbodyout_endsim_end_state_flag == 3 and out_bin_j == 2 and Lij_zang_DG > 20.):    # and IMS_binsin_counter > 1 and IMS_rp_counter > 1): # and out_bin_i == 1):
        #    exit()
        
        #if (Nbodyout_endsim_end_state_flag == 5 and IMS_rp_counter > 10 and IMS_binsin_counter > 3 and out_bin_i == 2 and out_bin_j == 3):
        #    exit()

        #if (Nbodyout_endsim_end_state_flag == 5 and IMS_rp_counter > 10 and IMS_binsin_counter > 3): test = input('push test')
        
        #if (Nbodyout_endsim_end_state_flag == 3 and out_bin_j == 3):
        #    exit()
        
        
        
        #---------------------------------------
        #Make 2. merger illustration:
        #---------------------------------------
        make_secondM_ill_yesno = 0
        if (make_secondM_ill_yesno == 1 and Nbodyout_endsim_end_state_flag == 5 and out_bin_i == 2 and out_bin_j == 3):
            #ASSUMES 2,3 MERGES FIRST!!!
            #ASSUMES EQUAL MASS CASE!!!!
            
            
            #------------------------
            #3-body part:
            #------------------------
            #Read data:
            fn_NbodyTides_dataout_pos           = 'NbodyTides_dataout_pos.dat'
            fn_NbodyTides_dataout_vel           = 'NbodyTides_dataout_vel.dat'        
            #pos:
            tf = open(fn_NbodyTides_dataout_pos,  "r")
            NbodyTides_dataout_pos_3b           = np.loadtxt(tf, dtype=float)
            pos_3b                              = NbodyTides_dataout_pos_3b
            tf.close()
            #vel:
            tf = open(fn_NbodyTides_dataout_vel,  "r")
            NbodyTides_dataout_vel_3b           = np.loadtxt(tf, dtype=float)
            vel_3b                              = NbodyTides_dataout_vel_3b
            tf.close()
            #info:
            nr_t_3b                             = len(NbodyTides_dataout_pos_3b[:,0])
            #------------------------
            #------------------------
            
            
            #------------------------
            #2-body part:
            #------------------------
            #set:
            tindex  = nr_t_3b-1
            #set equal mass (em) values:
            em_mass1 = b1_mass
            em_rad1  = b1_rad
            
            #obj1:
            m1_2b       = em_mass1
            r1_2b       = em_rad1
            p1_2b       = pos_3b[tindex,0:3]
            v1_2b       = vel_3b[tindex,0:3]            
            b1_pos_2b   = np.array([p1_2b[0], p1_2b[1], p1_2b[2]])
            b1_vel_2b   = np.array([v1_2b[0], v1_2b[1], v1_2b[2]])
            b1_const_arr_2b = np.array([m1_2b, r1_2b, b1_gas_n, b1_gas_gamma, b1_Mqp, b1_evoTides_yesno, b1_RigidSph_yesno, 0,0,0], dtype='d')
            
            #obj2:
            m2_2b       = (em_mass1 + em_mass1)
            r2_2b       = (em_rad1  + em_rad1)
            p2_2b       = (pos_3b[tindex,3:6] + pos_3b[tindex,6:9])/2.
            v2_2b       = (vel_3b[tindex,3:6] + vel_3b[tindex,6:9])/2.            
            b2_pos_2b   = np.array([p2_2b[0], p2_2b[1], p2_2b[2]])
            b2_vel_2b   = np.array([v2_2b[0], v2_2b[1], v2_2b[2]])
            b2_const_arr_2b = np.array([m2_2b, r2_2b, b1_gas_n, b1_gas_gamma, b1_Mqp, b1_evoTides_yesno, b1_RigidSph_yesno, 0,0,0], dtype='d')
             
            #Make IC file:            
            fn = 'MCinput_Nbody.txt'
            text_file = open(fn, "w")
            #nr particles in the sim   
            text_file.write('2' + '\n')
            #nbody code settings: 
            # [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, ...]
            nbody_params_arr_1_INT_2b  = np.array([0, 1, 0, 1000000000, 10, 1, 0,0,0,0], dtype='i')
            # [scale_dt, max_sim_time (set later in code!!!), evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, ...]
            nbody_params_arr_2_REAL_2b = np.array([0.01, 1e10, evolvetides_threshold, 0.1, 10, 0.01, tidaldisrup_threshold, insp_threshold,0,0], dtype='d')
            #save IC to file:
            np.savetxt(text_file, nbody_params_arr_1_INT_2b[None],   fmt='%10i')
            np.savetxt(text_file, nbody_params_arr_2_REAL_2b[None],  fmt='%10f')
            #body 1:
            np.savetxt(text_file, b1_const_arr_2b[None],  fmt='%10f')
            np.savetxt(text_file, b1_pos_2b[None],  fmt='%10f')
            np.savetxt(text_file, b1_vel_2b[None],  fmt='%10f')
            np.savetxt(text_file, b1_q,  fmt='%10f')
            np.savetxt(text_file, b1_qdot,  fmt='%10f')
            #body 2:
            np.savetxt(text_file, b2_const_arr_2b[None],  fmt='%10f')
            np.savetxt(text_file, b2_pos_2b[None],  fmt='%10f')
            np.savetxt(text_file, b2_vel_2b[None],  fmt='%10f')
            np.savetxt(text_file, b1_q,  fmt='%10f')
            np.savetxt(text_file, b1_qdot,  fmt='%10f')
            #close file:
            text_file.close()    
    
            #Simulate:
            #Run Nbody_AffineTides_solver.exe with 'MCinput_Nbody.txt' as input:
            subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)            
            
            #Read data:
            fn_NbodyTides_dataout_pos           = 'NbodyTides_dataout_pos.dat'
            #pos:
            tf = open(fn_NbodyTides_dataout_pos,  "r")
            NbodyTides_dataout_pos_2b           = np.loadtxt(tf, dtype=float)
            pos_2b                              = NbodyTides_dataout_pos_2b
            tf.close()
            nr_t_2b                             = len(pos_2b[:,0])
            #------------------------
            #------------------------
        
        
            #------------------------
            #plot:
            #------------------------     
            #define:
            a_0 = MCvalues_SMA_arr[0]
            font = {'family' : 'serif'}
            mpl.rc('font', **font)               
            #plt.rcParams['axes.facecolor'] = 'black'            
            f   = plt.figure(figsize=(6,10))
            gs  = gridspec.GridSpec(2, 1,height_ratios=[1,1])
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])
            pos_1GWM    = b2_pos_2b
            #pos_2GWM    = np.array()

            #set zoom box:
            zoom_x12 = [-11.75, -6.5]
            zoom_y12 = [-0.75, 1.5]
            
            
            #------------
            #plot1:
            #------------
            #3-body plot:
            ax1.plot(NbodyTides_dataout_pos_3b[:,0]/a_0, NbodyTides_dataout_pos_3b[:,1]/a_0, color = 'green', linewidth=1, label = 'BH [1]')      #obj 1
            ax1.plot(NbodyTides_dataout_pos_3b[:,3]/a_0, NbodyTides_dataout_pos_3b[:,4]/a_0, color = 'blue',   linewidth=1, label = 'BH [2]')      #obj 2
            ax1.plot(NbodyTides_dataout_pos_3b[:,6]/a_0, NbodyTides_dataout_pos_3b[:,7]/a_0, color = 'red',    linewidth=1, label = 'BH [3]')      #obj 3
            #2-body plot:
            ax1.plot(NbodyTides_dataout_pos_2b[:,0]/a_0, NbodyTides_dataout_pos_2b[:,1]/a_0, color = 'green', linewidth=1)                        #obj 1
            ax1.plot(NbodyTides_dataout_pos_2b[:,3]/a_0, NbodyTides_dataout_pos_2b[:,4]/a_0, color = 'purple', linewidth=1, label = 'BH [2+3]')    #obj 2+3            
            #axis:
            ax1.set_xlim(-13.5, 6.0)
            ax1.set_ylim(-3.0, 2.0)
            #labels:
            ax1.set_title(r'Formation of a Double GW Merger')
            ax1.set_xlabel(r'$x/a_{0}$')
            ax1.set_ylabel(r'$y/a_{0}$')
            #legend:
            #ax1.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = True)
            leg = ax1.legend(loc='lower left', numpoints = 1, fontsize = 10.0, frameon = True)
            leg.get_frame().set_edgecolor('black')
            for t in leg.get_texts():
                t.set_color('black')  
            #zoombox:
            zoombox_xpos    = [zoom_x12[0], zoom_x12[1], zoom_x12[1], zoom_x12[0], zoom_x12[0]]
            zoombox_ypos    = [zoom_y12[0], zoom_y12[0], zoom_y12[1], zoom_y12[1], zoom_y12[0]]
            ax1.plot(zoombox_xpos,zoombox_ypos, color = 'black', linewidth=2, linestyle = ':')
            #------------
            #plot2:
            #------------
            #3-body plot:
            ax2.plot(NbodyTides_dataout_pos_3b[:,0]/a_0, NbodyTides_dataout_pos_3b[:,1]/a_0, color = 'green', linewidth=1, label = 'BH [1]')      #obj 1
            ax2.plot(NbodyTides_dataout_pos_3b[:,3]/a_0, NbodyTides_dataout_pos_3b[:,4]/a_0, color = 'blue',   linewidth=1, label = 'BH [2]')      #obj 2
            ax2.plot(NbodyTides_dataout_pos_3b[:,6]/a_0, NbodyTides_dataout_pos_3b[:,7]/a_0, color = 'red',    linewidth=1, label = 'BH [3]')      #obj 3
            #2-body plot:
            ax2.plot(NbodyTides_dataout_pos_2b[:,0]/a_0, NbodyTides_dataout_pos_2b[:,1]/a_0, color = 'green', linewidth=1)                        #obj 1
            ax2.plot(NbodyTides_dataout_pos_2b[:,3]/a_0, NbodyTides_dataout_pos_2b[:,4]/a_0, color = 'purple', linewidth=1, label = 'BH [2+3]')    #obj 2+3            
            #axis:
            ax2.set_xlim(zoom_x12)
            ax2.set_ylim(zoom_y12)
            #labels:
            ax2.set_xlabel(r'$x/a_{0}$')
            ax2.set_ylabel(r'$y/a_{0}$')
            #add text:
            ax2.text(-7.6, 1.2, r'1. GW Merger', size = 10,
                    horizontalalignment='center',
                    verticalalignment='center')
            ax2.text(-10.9, 0.75, r'2. GW Merger', size = 10,
                    horizontalalignment='center',
                    verticalalignment='center')
            #------------
            #save and show fig:
            #------------
            save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/doubleGWmergers/'
            plt.savefig(save_data_folder + 'doubleGW_ill_1.eps', bbox_inches='tight')
            plt.show()
            exit()
            #------------
            #------------------------
            #------------------------ 
        #---------------------------------------
        
        

        #---------------------------------------
        #Resonance analysis:
        #---------------------------------------
        if (resonance_analysis_yesno == 1):
        #---------------------------------------    
            #open resonance data:
            tf = open('NbodyTides_dataout_resonance.dat',  "r")
            RES_info_arr        = np.loadtxt(tf, dtype=float)
            nr_binsin_changes   = len(open('NbodyTides_dataout_resonance.dat').readlines()) 
            tf.close()
            #if (nr_binsin_changes > 15): exit()
            #save resonance data:
            if (nr_binsin_changes < maxnr_binsin_changes):
                for rc in range(0,nr_binsin_changes):
                    #format: save_arr = [pc,sc,rc,0/1,Return_Info_arr_INT(2:4), Return_Info_arr_REAL(3:5), Return_Info_arr_REAL(8:10)]
                    if (nr_binsin_changes == 1): save_arr = np.concatenate((np.array([pc,sc,rc]),RES_info_arr[:]))
                    if (nr_binsin_changes >  1): save_arr = np.concatenate((np.array([pc,sc,rc]),RES_info_arr[rc,:]))
                    Resonance_data_iniarr[totnr_res_counter,:] = save_arr[:]
                    totnr_res_counter = totnr_res_counter + 1
                #---------------------
                #TEST section:
                #---------------------
                #if (nr_binsin_changes > 1 and Nbodyout_endsim_end_state_flag == 2 and RES_info_arr[nr_binsin_changes-1,0] == 1):
                #    print 'Nbodyout_endsim_Return_Info_arr_INT[5]:  ', Nbodyout_endsim_Return_Info_arr_INT[5]
                #    print 'Nbodyout_endsim_Return_Info_arr_INT[5]:  ', Nbodyout_endsim_Return_Info_arr_INT[5]
                #    print 'Nbodyout_endsim_Return_Info_arr_INT[5]:  ', Nbodyout_endsim_Return_Info_arr_INT[5]
                #    print 'Nbodyout_endsim_Return_Info_arr_INT[5]:  ', Nbodyout_endsim_Return_Info_arr_INT[5]
                #    print 'Nbodyout_endsim_Return_Info_arr_INT[5]:  ', Nbodyout_endsim_Return_Info_arr_INT[5]
                #    if (Nbodyout_endsim_Return_Info_arr_INT[5] > 1): exit()
                #    rperi = RES_info_arr[nr_binsin_changes-1,5]*(1.-RES_info_arr[nr_binsin_changes-1,6])
                #   print 'rperi:   ', rperi
                #    if (rperi > 3.0): exit()
                #if (nr_binsin_changes == 3 and Nbodyout_endsim_end_state_flag == 2 and Nbodyout_endsim_Return_Info_arr_INT[5] == 0): print RES_info_arr
                #if (nr_binsin_changes == 3 and Nbodyout_endsim_end_state_flag == 2 and Nbodyout_endsim_Return_Info_arr_INT[5] == 0): exit()
                #if (nr_binsin_changes >= 2 and Nbodyout_endsim_end_state_flag == 2 and Nbodyout_endsim_Return_Info_arr_INT[5] > 1): print RES_info_arr
                #if (nr_binsin_changes >= 2 and Nbodyout_endsim_end_state_flag == 2 and Nbodyout_endsim_Return_Info_arr_INT[5] > 1): exit()
                #if (nr_binsin_changes > 9 and Nbodyout_endsim_end_state_flag == 3): 
                #    print RES_info_arr
                #    exit()
                #---------------------
            #print info:
            print 'nr_binsin_changes:   ', nr_binsin_changes
        #---------------------------------------
        #---------------------------------------

        #---------------------------------------
        #Save results from MC/Nbody:
        #---------------------------------------
        #scattering param values:
        MCout_arr_scattering_params[totsc,:]  = [SMA_bin,vinf_sin,b_max,0,0,0,0,0,0,0]
        #endstate output from Nbody code:
        MCout_arr_endstate_info_INT[totsc,:]  = Nbodyout_endsim_Return_Info_arr_INT[:]
        MCout_arr_endstate_info_REAL[totsc,:] = Nbodyout_endsim_Return_Info_arr_REAL[:]
        #---------------------------------------

        #---------------------------------------
        #counters:
        #---------------------------------------
        #count tot nr of scatterings:
        totsc = totsc+1
        #---------------------------------------

        #---------------------------------------
        #Print info:
        #---------------------------------------
        print totsc, ' out of :', tot_nr_scatterings
        print 'Output from Nbody:   ', Nbodyout_endsim_end_state_flag, Nbodyout_endsim_Return_Info_arr_INT#, Nbodyout_endsim_Return_Info_arr_REAL
        print 'SMA (AU), vinf (km/sec):   ', SMA_bin/AU_U, vinf_sin/kmsec_U

        #print bf_sampl_arr[sc,0]*(vinf_kmsec_arr/vc_kmsec_perSMA_arr), f_val
        #test = input('push test')        
        #---------------------------------------        
    #-----------------------------------------------------------------    
    #END loop over samplings 'sc'.
    #-----------------------------------------------------------------
#------------------------------------------------------------------------------------------
#END loop over bin-sin params (SMA, vinf, ...) 'pc'.
#------------------------------------------------------------------------------------------



#test_pos3_arr[sc,0] = theta
#test_pos3_arr[sc,1] = phi


#FIG angles:
fig = plt.figure(figsize=(5.0, 4.0))   
ax1 = fig.add_subplot(111)    
#hist settings:
min_H           = 0.0
max_H           = 180.
nr_Hbins        = 100

data_ij         = test_pos3_arr[:,0]
Hdata           = data_ij
Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
ax1.fill_between(Hx, 0, Hy, color = 'grey', step="pre", alpha = 0.35)

data_ij         = test_pos3_arr[:,1]
Hdata           = data_ij
Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
ax1.fill_between(Hx, 0, Hy, color = 'red', step="pre", alpha = 0.35)

plt.show()
exit()


#------------------------------------------------------------------------------------------
#Save results from MC/Nbody run:
#------------------------------------------------------------------------------------------
#---------------------------------------
#Save IC info:
#---------------------------------------
#open file:
tf = open(MCoutput_folder + add_filename+'_'+'MC_IC_all_info.txt',  "w")
#Info body 1,2,3:
np.savetxt(tf, b1_const_arr[None],   fmt='%5f')
np.savetxt(tf, b2_const_arr[None],   fmt='%5f')
np.savetxt(tf, b3_const_arr[None],   fmt='%5f')
#Nbody params:
np.savetxt(tf, nbody_params_arr_1_INT[None],    fmt='%5i')
np.savetxt(tf, nbody_params_arr_2_REAL[None],   fmt='%5f')
#MC params:
MC_setting_params = np.array([nr_SMA, nr_vinf, rp_max_unit_SMAbin, FtidFbin_r_simsurf, nr_scatterings_per_paramcomb, 0,0,0,0,0])
np.savetxt(tf, MC_setting_params[None], fmt='%5f')
#close file:
tf.close()
#---------------------------------------
#---------------------------------------
#Save MC/Nbody output results:
#---------------------------------------
#save: MCout_arr_endstate_info_INT
tf = open(MCoutput_folder + add_filename+'_'+'MCout_arr_endstate_info_INT.txt',     "w")
np.savetxt(tf, MCout_arr_endstate_info_INT,   fmt='%5i')
tf.close()
#save: MCout_arr_endstate_info_REAL
tf = open(MCoutput_folder + add_filename+'_'+'MCout_arr_endstate_info_REAL.txt',    "w")
np.savetxt(tf, MCout_arr_endstate_info_REAL,  fmt='%5f')
tf.close()
#save: MCout_arr_scattering_params
tf = open(MCoutput_folder + add_filename+'_'+'MCout_arr_scattering_params.txt',     "w")
np.savetxt(tf, MCout_arr_scattering_params,   fmt='%5f')
tf.close()
#---------------------------------------
#---------------------------------------
#Save resonance analysis output results:
#---------------------------------------
if (resonance_analysis_yesno == 1):
    #make final array:
    Resonance_data_finalarr = Resonance_data_iniarr[0:totnr_res_counter,:]
    #save final array:
    tf = open(MCoutput_folder + add_filename+'_'+'MCout_arr_resonance_analysis.txt',    "w")
    np.savetxt(tf, Resonance_data_finalarr,   fmt='%5f')
    tf.close()
#---------------------------------------
#------------------------------------------------------------------------------------------


















