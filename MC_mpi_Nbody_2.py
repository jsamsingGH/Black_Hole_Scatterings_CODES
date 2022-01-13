#!/usr/bin/env python
"""
Following MPI template is from:
https://github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py
Demonstrate the task-pull paradigm for high-throughput computing
using mpi4py. Task pull is an efficient way to perform a large number of
independent tasks when there are more tasks than processors, especially
when the run times vary for each task. 
"""

from mpi4py import MPI
import testname
import numpy as np
import sys
from subprocess import call
import subprocess
import time
import itertools
from sympy import *
from sympy.solvers import solve
from sympy import Symbol
from sympy.abc import x



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def func_make_MC_IC_3body_binsin():
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
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
    #Specify Properties for the 3 objects:
    #------------------------------------------------------------------------------------------
    #Order: Body 1,2 will be in a binary and Body 3 is incoming.
    #---------------------------------------
    #IC Obj 1:  (binary 1)
    #---------------------------------------
    b1_mass      = 30.0
    b1_rad       = Rsch_1Msun_unitRsun*b1_mass  #0.013*((1.43/b1_mass)**(1./3.))*((1.-b1_mass/1.43)**(0.447)) #1.7246e-5
    b1_gas_n     = 3.0
    b1_gas_gamma = 5./3.
    b1_Mqp       = 0.0376788*(b1_mass*(b1_rad**2))
    b1_evoTides_yesno = 0
    b1_RigidSph_yesno = 1
    b1_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
    b1_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
    b1_const_arr = np.array([b1_mass, b1_rad, b1_gas_n, b1_gas_gamma, b1_Mqp, b1_evoTides_yesno, b1_RigidSph_yesno, 0,0,0], dtype='d')
    #---------------------------------------
    #IC Obj 2:  (binary 2)
    #---------------------------------------
    b2_mass      = 10.0
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
    b3_mass      = 10.0
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
    sampldist_IC12          = 1     # =1: ISOTROPIC dist. =2 RESTRICTED (CO-PLANAR) dist.
    opening_dtheta_DG       = 5.0   #orbital inclination in degrees (DG) wrt. to coplanar config. For ONLY coplanar set therefore = 0!
    use_12PN                = 0
    use_25PN                = 1
    insp_threshold          = 5.0       #a_min = insp_threshold*(Ri+Rj).
    insp_SMA_abin           = (1./10.)  #units of a_bin
    tidaldisrup_threshold   = 3.0
    evolvetides_threshold   = min([0.05, 32.*(b2_mass/b1_mass)*(insp_threshold*(b1_rad+b2_rad)/b1_rad)**(-3.), 32.*(b3_mass/b1_mass)*(insp_threshold*(b1_rad+b3_rad)/b1_rad)**(-3.)]) #works for: TO is only '1' - OR all identical objects. - the 0.05 in front is just to make sure its small.
    # [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, n_particles=3, de_analysis, ...]
    nbody_params_arr_1_INT  = np.array([use_12PN, use_25PN, 1, 1000000000, 10, 0, 3, 0,0,0], dtype='i')
    # [scale_dt, max_sim_time (set later in code!!!), evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, insp_SMA (set later in code!!!), ...]
    nbody_params_arr_2_REAL = np.array([0.01, -1.0, evolvetides_threshold, 0.01, -1.0, 0.5, tidaldisrup_threshold, insp_threshold, -1.0, 0], dtype='d')    
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #Make arrays with parameters we vary in this MC (SMA, vinf, ...):
    #------------------------------------------------------------------------------------------   
    #---------------------------------------
    #semi-major axis (SMA):
    #---------------------------------------
    #---------------------
    #opt1: Sequence of values lin in logspace:
    #---------------------
    #min_SMA_AU      = 0.0001       #input in AU
    #max_SMA_AU      = 0.001         #input in AU     
    #nr_SMA          = 10           #incl min, max vals
    #MCvalues_SMA_arr  = AU_U*np.logspace(np.log10(min_SMA_AU),np.log10(max_SMA_AU),nr_SMA)                    #output units are in correct code units.
    #---------------------
    #opt2: A few (>=1) user defined values:
    #---------------------
    SMA_AU_arr          = [0.1]     #input in AU
    MCvalues_SMA_arr    = AU_U*np.array(SMA_AU_arr)
    nr_SMA              = len(MCvalues_SMA_arr)
    #---------------------
    #opt3: Time period distribution lin in logspace:
    #---------------------
    #min_Pdays   = 0.01  #in days
    #max_Pdays   = 100.  #in days
    #nr_SMA      = 10    #incl min, max vals
    #Pdays_arr           = np.logspace(np.log10(min_Pdays),np.log10(max_Pdays),nr_SMA)
    #SMA_SI_arr          = ((G_new_SI*((b1_mass+b2_mass)*M_sun_SI)*((86400.*Pdays_arr)**2.))/(4.*np.pi**2.))**(1./3.) 
    #MCvalues_SMA_arr    = SMA_SI_arr/R_sun_SI
    #---------------------
    #calc corresponding charac vel (vc) per SMA in arr:
    #---------------------
    vc_kmsec_perSMA_arr = (1./kmsec_U)*np.sqrt((b1_mass*b2_mass*(b1_mass+b2_mass+b3_mass)/(b3_mass*(b1_mass+b2_mass)))*(1./MCvalues_SMA_arr))
    #---------------------
    #---------------------------------------
    #incoming vel of body 3 at infinity (vinf):
    #---------------------------------------
    #---------------------
    #opt1: Sequence of values lin in logspace:
    #---------------------
    #min_vinf_kmsec  = 0.5*vc_kmsec_perSMA_arr[0]  
    #max_vinf_kmsec  = 8.0*vc_kmsec_perSMA_arr[0]
    #nr_vinf         = 10         #incl min, max vals 
    #MCvalues_vinf_arr = kmsec_U*np.logspace(np.log10(min_vinf_kmsec),np.log10(max_vinf_kmsec),nr_vinf)        #output units are in correct code units.
    #---------------------
    #opt2: A few (>=1) user defined values:
    #---------------------
    vinf_kmsec_arr      = [1.0] #[10.0]     #input in km/sec
    MCvalues_vinf_arr   = kmsec_U*np.array(vinf_kmsec_arr)
    nr_vinf             = len(MCvalues_vinf_arr)
    #info: (13.3 AU for 10 km/sec)
    #info: (53.3 AU for 5 km/sec) 
    #info: (1330 AU for 1 km/sec) 
    #print (1./AU_U)*((b1_mass*b2_mass*(b1_mass+b2_mass+b3_mass)/(b3_mass*(b1_mass+b2_mass)))*(1./(MCvalues_vinf_arr**2.)))
    #---------------------
    #---------------------
    #opt4: phase-space analysis 1: (see http://adsabs.harvard.edu/abs/1983AJ.....88.1549H)
    phase_space_analysis_yesno      = 0    #ANALYSIS YES/NO.
    #---------------------
    #overall settings:    
    #---------------------
    #choose option:
    phase_space_analysis_op1        = 1
    phase_space_analysis_op2        = 0
    res_pa              = 50    #set res (same for b/f/g)
    sim_dist_unitSMA    = 20.   #P. Hut: R = 40*r, a=2r, => R = 20*a (DONT TOUCH!!!)
    nr_tot_pa           = res_pa*res_pa
    #---------------------
    #OPT1: vary: f,b fix: g
    #---------------------
    if (phase_space_analysis_op1 == 1):
        #fix:
        plane_gamma         = 0.0*np.pi    # = 0 is in same plane. Positive val rotates in positive angular direction. This angle is fixed for a given sim.
        #OPT: full scan:
        rpmax_unitSMA       = 3.5
        bmax_imp_unitSMA    = (np.sqrt(2.*(b1_mass+b2_mass+b3_mass)*(rpmax_unitSMA*MCvalues_SMA_arr[0])/(MCvalues_vinf_arr[0]**2.)))/MCvalues_SMA_arr[0]
        b_imp_unitSMA_arr   = bmax_imp_unitSMA*np.linspace(-1.0, 1.0, res_pa)       #impact param b range (set number in front)
        f_phase_rad_arr     = (2.*np.pi)*np.linspace(0.0, 1.0, res_pa)              #relative phase. use full circle: 2pi.     
        #OPT: zoom:
        #b_imp_unitSMA_arr   = np.linspace(1.0, 5.0, res_pa)                        #impact param b range (set number in front)
        #f_phase_rad_arr     = np.linspace(2.0, 4.0, res_pa)                        #relative phase. use full circle: 2pi.    
        #final val array:
        pa_sampl_arr        = np.array(list(itertools.product(b_imp_unitSMA_arr, f_phase_rad_arr)))
    #---------------------
    #OPT2: vary: f,g fix: b
    #---------------------
    if (phase_space_analysis_op2 == 1):
        #fix:
        rpmax_unitSMA       = 1e-10 #do not set exactly to zero as some of the trigonometric funtions later in the code then break down.
        b_imp_unitSMA       = abs((np.sqrt(2.*(b1_mass+b2_mass+b3_mass)*(rpmax_unitSMA*MCvalues_SMA_arr[0])/(MCvalues_vinf_arr[0]**2.)))/MCvalues_SMA_arr[0])
        #OPT: full scan:
        g_phase_rad_arr     = np.linspace(-np.pi/2., np.pi/2., res_pa)          #GAMMA (g) relative phase. use full circle: 2pi.     
        f_phase_rad_arr     = np.linspace(0.0, 2.*np.pi, res_pa)                #PHASE (f) relative phase. use full circle: 2pi.     
        #final val array:
        pa_sampl_arr        = np.array(list(itertools.product(g_phase_rad_arr, f_phase_rad_arr)))
    #---------------------    
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #3body scattering settings:
    #------------------------------------------------------------------------------------------
    #maximum pericenter dist set by a tidal threshold:
    rp_max_FtidFbin                     = 1.0   # = Ftid/Fbin (=1 will give max(r_p) = 1.26*a in equal mass case, and = 2.0*a for [0.2,1.4,1.4])
    #the rsim surface is where the FtidFbin frac equals this number:
    r_simsurf_FtidFbin                  = 0.01 # For high mass bin set this to lower val. (same idea as for rp_max)
    #nr of rnd samplings around the sphere of the bin per bin-sin comb:
    nr_scatterings_per_paramcomb        = 100
    #SPECIAL VALUES FOR PHASE-SPACE ANALYSIS:
    if (phase_space_analysis_yesno == 1):
        nr_scatterings_per_paramcomb    = nr_tot_pa
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #Define:
    #------------------------------------------------------------------------------------------
    #mass:
    mass_tot    = b1_mass + b2_mass + b3_mass
    mass_bin    = b1_mass + b2_mass
    mass_sin    = b3_mass
    mass_red    = mass_bin*mass_sin/mass_tot 
    #final output list:
    MCinfo_ICNbody_list = []
    #numbers:
    nr_tot_scatterings  = nr_SMA*nr_vinf*nr_scatterings_per_paramcomb
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #Time settings:
    #------------------------------------------------------------------------------------------
    #Rough estimate of max 1 bin-sin passage time:
    #By simple scalings we see that the min(a) and min(v) lead to the maximum t/Torb.
    d_ini       = min(MCvalues_SMA_arr)*((1./r_simsurf_FtidFbin)*(mass_bin*b3_mass/(b1_mass*b2_mass)))**(1./3.) #dist from bin to sin at the beginning of the sim (r_simsurf_FtidFbin)
    T_1_passage = 2.*d_ini/min(MCvalues_vinf_arr)                                                               #time it takes to pass 2d assuming constant v (initial at r_simsurf_FtidFbin)
    Torb_ini    = 2.*np.pi*np.sqrt((min(MCvalues_SMA_arr)**3.)/(mass_bin))                                      #orbital time of ini target bin in code units.
    T_1_passage_unit_Torb = T_1_passage/Torb_ini                                                                #1 passage time in units of Torb initial.
    #for m=1, a=1, m3 = 10^6, v=50,  r_simsurf_FtidFbin = 0.001, T_1_passage_unit_Torb = 340
    #for similar mass bin-sin interactions:
    Tsim_max_unit_Torb  = 100.                          #max sim time in nr initial orbital times
    #for bin-SMBH interactions:
    #Tsim_max_unit_Torb  = 2.*T_1_passage_unit_Torb       #max sim time in nr initial orbital times
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------

    
    #------------------------------------------------------------------------------------------
    #make MC IC info list:
    #------------------------------------------------------------------------------------------
    #allocate:
    MC_settings_list_INT    = np.zeros(5,dtype='i')
    MC_settings_list_REAL   = np.zeros(5,dtype='d')
    #insert info:
    MC_settings_list_INT[:]     = [nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, nr_tot_scatterings, 0]    #where nr_tot_scatterings is simply =nr_SMA*nr_vinf*nr_scatterings_per_paramcomb
    MC_settings_list_REAL[:]    = [0.0, 0.0, 0.0, 0.0, 0.0]
    #put into one list:
    MC_settings_list = [MC_settings_list_INT, MC_settings_list_REAL]
    #------------------------------------------------------------------------------------------

    
    #------------------------------------------------------------------------------------------
    #Make MC ICs:
    #------------------------------------------------------------------------------------------
    #initilize IC ID counter:
    icidc = int(0)
    #Loop over all a,v combinations:
    for ac in range(0,nr_SMA):                               # ac-loop over SMA
        for vc in range(0,nr_vinf):                          # vc-loop over vinf
                        
            #-----------------------------------------------------------------
            #Get binary-single parameters (SMA, vinf, ...):
            #-----------------------------------------------------------------
            #Binary SMA
            SMA_bin     = MCvalues_SMA_arr[ac]      #already in correct code units
            #vel at infinity of single:
            vinf_sin    = MCvalues_vinf_arr[vc]     #already in correct code units            
            #-----------------------------------------------------------------

            #-----------------------------------------------------------------
            #Code Input/stopping criteria etc.:
            #-----------------------------------------------------------------
            Torb_inibin = 2.*np.pi*np.sqrt((SMA_bin**3.)/(mass_bin))        #orbital time of ini target bin in code units.
            #select and define sim time for this a,v comb scattering:
            nbody_params_arr_2_REAL[1]  = Tsim_max_unit_Torb*Torb_inibin    #set limit
            nbody_params_arr_2_REAL[4]  = 100000000                         #just put a very high number. Might be set to remaining sim time at Master process. See mpi modules.
            nbody_params_arr_2_REAL[8]  = insp_SMA_abin*SMA_bin             #inspiral criteria: fraction of initial SMA        
            #-----------------------------------------------------------------

            #-----------------------------------------------------------------
            #Scattering settings:
            #-----------------------------------------------------------------
            #The max peri-center dist is set by a tidal threshold:
            rp_max      = SMA_bin*((1./rp_max_FtidFbin)*((b1_mass+b2_mass)*b3_mass/(b1_mass*b2_mass)))**(1./3.)            
            #from rp_max we can now calc max impact param (b):
            b_max       = rp_max*np.sqrt(1.+(2.*mass_tot)/(rp_max*(vinf_sin**2.)))
            #distance to sampling surface:
            r_sampsurf  = 1e20*SMA_bin   #just make sure to put in a very big number.
            #find distance to simulation surface (we analytically propagate sin body3 from sampsurf to simsurf to save sim time):
            r_simsurf_FtideFbin     = SMA_bin*((1./r_simsurf_FtidFbin)*((b1_mass+b2_mass)*b3_mass/(b1_mass*b2_mass)))**(1./3.)
            r_simsurf_rp_max        = rp_max
            r_simsurf               = max([r_simsurf_FtideFbin, r_simsurf_rp_max]) # at rsim the FtidFbin will always be < or = r_simsurf_FtidFbin.
            #-----------------------------------------------------------------

            #-----------------------------------------------------------------
            #Set body1, body2 up in a CIRCULAR binary with COM in (0,0,0):
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
            for sc in range(0,nr_scatterings_per_paramcomb):
                                                
                #---------------------------------------
                #sampling parameters:
                #---------------------------------------
                #ISOTROPIC CONFIGURATION:
                if (sampldist_IC12 == 1):
                    b_sampsurf          = b_max*np.sqrt(rndnum_arr[sc,0])	        #P flat in b^2			- rad pos on 0->bmax disc
                    psi_sampsurf        = (2.*np.pi)*rndnum_arr[sc,1]			    #P flat in psi			- angle on b disc.
                    theta               = np.arccos(2.*rndnum_arr[sc,2]-1.)         #P flat in cos(theta)	- pos on unitsphere
                    phi                 = ((2.*np.pi)*rndnum_arr[sc,3])			    #P flat in phi			- pos on unitsphere
                #FIXED ANGLE (or CO-PLANAR) CONFIGURATION:
                if (sampldist_IC12 == 2):       
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

                #mass_tot    = b1_mass + b2_mass + b3_mass
                #mass_bin    = b1_mass + b2_mass
                #mass_sin    = b3_mass
                #mass_red    = mass_bin*mass_sin/mass_tot          
                
                #---------------------------------------
                #phase-space analysis 1:
                #---------------------------------------
                #initialize: (we need to define these variables for later saving -- even when the analysis below is not active)
                R_val   = 0.0
                b_val   = 0.0
                f_val   = 0.0
                g_val   = 0.0                
                if (phase_space_analysis_yesno == 1):
                    #--------------------
                    #OPT1: vary: f,b fix: g
                    #--------------------
                    if (phase_space_analysis_op1 == 1):
                        #param vals:
                        R_val               = sim_dist_unitSMA*SMA_bin          #constant. Not varied.
                        b_val               = pa_sampl_arr[sc,0]*SMA_bin        #b at infinity. correct units
                        f_val               = pa_sampl_arr[sc,1]                #binary phase in radians
                        g_val               = plane_gamma                       #orbital plane angle
                    #--------------------
                    #OPT2: vary: f,g fix: b
                    #--------------------
                    if (phase_space_analysis_op2 == 1):
                        #param vals:
                        R_val               = sim_dist_unitSMA*SMA_bin          #constant. Not varied.
                        b_val               = b_imp_unitSMA*SMA_bin             #b at infinity. correct units
                        f_val               = pa_sampl_arr[sc,1]                #binary phase in radians
                        g_val               = pa_sampl_arr[sc,0]                #orbital plane angle                    
                    #--------------------
                    #calc: (some eqs taken from this site: http://www.braeunig.us/space/orbmech.htm)
                    #--------------------
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
                    Ry	= np.array([[np.cos(-g_val),0,np.sin(-g_val)], [0,1,0], [-np.sin(-g_val),0,np.cos(-g_val)]])
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
                    #Final output MUST BE: b3_posxyz_binCM, b3_velxyz_binCM
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
                #calc info for later saving:
                #---------------------------------------
                #Total energy:                
                E_tot_binsinsystem = (b1_mass*b2_mass)/(2.*SMA_bin) + (1./2.)*mass_red*(vinf_sin**2.)
                #Total angular momentum: 
                L_tot_binsinsystem = len((b1_mass*np.cross(b1_posxyz_CM, b1_velxyz_CM)) + (b2_mass*np.cross(b2_posxyz_CM, b2_velxyz_CM)) + (b3_mass*np.cross(b3_posxyz_CM, b3_velxyz_CM)))
                #---------------------------------------
                
                #---------------------------------------
                #IC list:
                #---------------------------------------
                IN_IC_code_version                          = 2       #MUST = 2 for parallel version !!!!
                IN_dimlen_IC_nbody                          = 3*9     #(for n_particles=3, = 27)
                #define arrays with a specific format:
                IN_IC_simparams_INT                         = np.zeros(10,dtype='i')
                IN_IC_simparams_REAL                        = np.zeros(10,dtype='d')
                IN_IC_nbody_const_posvel_qqdot_etc_arr      = np.zeros((IN_dimlen_IC_nbody,10),dtype='d')
                #sim parameters:
                IN_IC_simparams_INT[:]                      = nbody_params_arr_1_INT[:]
                IN_IC_simparams_REAL[:]                     = nbody_params_arr_2_REAL[:]
                #obj 1:
                IN_IC_nbody_const_posvel_qqdot_etc_arr[0,:]     = b1_const_arr[:]   #[1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[1,0:3]   = b1_posxyz_CM[:]   #[-0.538861, 93.044669, 64.479373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[2,0:3]   = b1_velxyz_CM[:]   #[0.004368, 0.117192, -0.027313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[3,0:3]   = b1_q[0,:]         #[1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[4,0:3]   = b1_q[1,:]         #[0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[5,0:3]   = b1_q[2,:]         #[0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[6,0:3]   = b1_qdot[0,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[7,0:3]   = b1_qdot[1,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[8,0:3]   = b1_qdot[2,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                #obj 2:
                IN_IC_nbody_const_posvel_qqdot_etc_arr[9,:]     = b2_const_arr[:]   #[1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[10,0:3]  = b2_posxyz_CM[:]   #[-22.038986, 93.044669, 64.479373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[11,0:3]  = b2_velxyz_CM[:]   #[0.004368, -0.187805, -0.027313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[12,0:3]  = b2_q[0,:]         #[1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[13,0:3]  = b2_q[1,:]         #[0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[14,0:3]  = b2_q[2,:]         #[0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[15,0:3]  = b2_qdot[0,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[16,0:3]  = b2_qdot[1,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[17,0:3]  = b2_qdot[2,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                #obj 3:
                IN_IC_nbody_const_posvel_qqdot_etc_arr[18,:]    = b3_const_arr[:]   #[1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[19,0:3]  = b3_posxyz_CM[:]   #[22.577848, -186.089337, -128.958746, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[20,0:3]  = b3_velxyz_CM[:]   #[-0.008735, 0.070613, 0.054626, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[21,0:3]  = b3_q[0,:]         #[1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[22,0:3]  = b3_q[1,:]         #[0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[23,0:3]  = b3_q[2,:]         #[0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[24,0:3]  = b3_qdot[0,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[25,0:3]  = b3_qdot[1,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                IN_IC_nbody_const_posvel_qqdot_etc_arr[26,0:3]  = b3_qdot[2,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
                #Put into one list:
                outlist_IC_Nbody            = [IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL, IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr] 
                #---------------------------------------
                #MC param info list:
                #---------------------------------------
                #allocate:
                MC_info_INT     = np.zeros(10,dtype='i')
                MC_info_REAL    = np.zeros(10,dtype='d')
                #input info:
                MC_info_INT[:]  = [icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
                MC_info_REAL[:] = [SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, R_val, b_val, f_val, g_val]
                outlist_MC_info_INT_REAL    = [MC_info_INT,MC_info_REAL]
                #---------------------------------------
                #append MC,IC list to final output list:
                #---------------------------------------
                MCinfo_ICNbody = [outlist_MC_info_INT_REAL, outlist_IC_Nbody]
                MCinfo_ICNbody_list.append(MCinfo_ICNbody)
                #---------------------------------------
                #update ic id counter
                #---------------------------------------
                icidc = icidc + 1
                #---------------------------------------
            #-----------------------------------------------------------------
                
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #Return final info from function:
    #------------------------------------------------------------------------------------------
    final_MC_output_list = [MC_settings_list, MCinfo_ICNbody_list]
    return final_MC_output_list
    #------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def func_run_Nbody(ICnbody):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    #input format: IC_inputNbody = ICnbody = [IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL, IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr]
    #unpack input list:
    IN_IC_code_version                      = ICnbody[0]
    IN_IC_simparams_INT                     = ICnbody[1]
    IN_IC_simparams_REAL                    = ICnbody[2]
    IN_dimlen_IC_nbody                      = ICnbody[3]
    IN_IC_nbody_const_posvel_qqdot_etc_arr  = ICnbody[4]
    #call Nbody code with these input IC:
    OUT_Nbody = testname.interfacesub(IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL, IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr)
    return OUT_Nbody    #OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL, OUT_out_xtra_3_info_REAL
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def enum(*sequential, **named):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    """
    This function is also from the MPI template example.
    Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main program: make IC list, Run Nbody using MPI and write output to .txt files
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
Basics: This full code will be send to and executed by each processor.
At the time the processor receives the code, the processor gets a unique 
'rank' number. The list of these number goes from 0 to N, where N is the total
number of processors. These numbers (or ids) will not change during the computation.
We use this rank number to pass info between the processors.
Normally you set processor 0 (rank '0') to be the 'Master' which sends and receives
info from the 'Worker' processers (rank >0). What part of the code a specific processor
runs is controlled by a simple 'if rank ==...' statement.
"""


#------------------------------------------------------------------------------------------
#Define and initialize:
#------------------------------------------------------------------------------------------
# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START', 'WAKEUP')

# Initializations and preliminaries
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#WORKER process executes code below
#------------------------------------------------------------------------------------------
if rank != 0:
#------------------------------------------------------------------------------------------       
   
    #---------------------------------------
    #If Worker is active:
    #---------------------------------------
    while True:
        #-------------------------
        #send info: READY to Master (dest=0):
        #-------------------------
        time.sleep(1.0*float(np.random.rand(1)))    #SET TIME BUFFER
        comm.send(None, dest=0, tag=tags.READY)
        #-------------------------
        #Receive info from Master on what to do:
        #-------------------------
        data    = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag     = status.Get_tag()
        #-------------------------
        #If Master says 'Start':
        #-------------------------
        if tag == tags.START:
            # split data file into info and ICs:
            MC_info_INT_REAL    = data[0] 
            IC_Nbody            = data[1]
            #Run Nbody code with 'IC_Nbody' as input:
            output_Nbody    = func_run_Nbody(IC_Nbody)
            result          = [MC_info_INT_REAL, output_Nbody]
            comm.send(result, dest=0, tag=tags.DONE)
        #-------------------------    
        #If Master says 'Exit' break loop:
        #-------------------------
        if tag == tags.EXIT:
            break
    #---------------------------------------
    
    #---------------------------------------
    #If Worker has been stopped (Exit):
    #---------------------------------------        
    #Send info to Master that Worker has now stopped (Exit): 
    comm.send(None, dest=0, tag=tags.EXIT)  
    #---------------------------------------
    
#------------------------------------------------------------------------------------------
#if rank != 0:    
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#MASTER process executes code below
#------------------------------------------------------------------------------------------
if rank == 0:
#------------------------------------------------------------------------------------------ 
    
    #---------------------------------------
    #wake up all processers: (necessary for e.g. SLURM)
    #---------------------------------------
    for nps in range(1,size):   #size = total number of processes
        comm.send(None, dest=nps, tag=tags.WAKEUP)
        print 'wake up processor:    ', nps
    #---------------------------------------
    
    #---------------------------------------
    #Input to Master process:
    #---------------------------------------
    #-------------------------
    #save data:
    #-------------------------
    data_name   = 'TEST4_' #WD10WD10WD10_NTWGR_0001to01AU_1kmsec_50000'#'WD14NSNS_NT_0001to01AU_1kmsec_50000'#'1MS1MS1MS_NTNGR_05to5AU_1kmsec_50000'#'WD14NSNS_WTWGR_0001to1AU_1kmsec_100000'
    #data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'
    data_folder = '/groups/astro/jsamsing/DATA_OUTPUT/'
    #-------------------------
    #max sim time in secs:
    #-------------------------
    #when using a scheduler (SLURM): PUT SAME TIME AS GIVEN TO SLURM!!
    sim_time_maxupper   = 1.0*3600.                             #in secs
    sim_time_buffer     = 0.02*sim_time_maxupper + 2.*60        #in secs
    max_sim_secs        = sim_time_maxupper - sim_time_buffer   
    '''
    NOTES on time: We give each task a physical time limit (N*initial orbital times)
    and a wall clock time limit which is the time we have left of the total sim time.
    If a task does not finish before the total sim time is over, the task will
    end with endstate id 12 (wall clock limit). IMPORTANT: Its important that all tasks have been
    distributed and have start running before the sim time is over. The code still works
    if this is not the case, but its really hard to use the dataset. The main idea
    with the wall clock limit is simply to stop the (often small) fraction of interactions
    that just keeps going forever.    
    '''
    #-------------------------
    #---------------------------------------
    
    #---------------------------------------
    #Generate bin-sin ICs: (all settings are set in the function)
    #---------------------------------------
    out = func_make_MC_IC_3body_binsin()
    #split into different lists:
    MC_settings_list_INT    = out[0][0] # MC_settings_list_INT[:]     = [nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, nr_tot_scatterings, 0]
    MC_settings_list_REAL   = out[0][1] # MC_settings_list_REAL[:]    = [0.0, 0.0, 0.0, 0.0, 0.0]
    MC_output_list          = out[1]    # [outlist_MC_info_INT_REAL, outlist_IC_Nbody], where: outlist_MC_info_INT_REAL = [[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0],[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]] 
    #define:
    nr_tot_scatterings      = MC_settings_list_INT[3]
    all_data                = [] 
    #---------------------------------------    

    #---------------------------------------
    #initialize:
    #---------------------------------------
    taskc               = 0         #task counter - counting the number of tasks (or IC scatterings) that have been send to Workers.
    nr_tasks_completed  = 0         #info counter.
    num_workers         = size - 1
    closed_workers      = 0
    time_start          = time.time()
    #---------------------------------------
    
    #---------------------------------------
    #SIMULATE - if not all workers are closed:
    #---------------------------------------
    while closed_workers < num_workers:        
        #-------------------------
        #receive info from Worker:
        #-------------------------
        data    = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source  = status.Get_source()
        tag     = status.Get_tag()
        #-------------------------
        #if Worker is Ready:
        #-------------------------
        if tag == tags.READY:
            #-------------------------
            #if we have more tasks, send Worker a task:
            #-------------------------
            if taskc < nr_tot_scatterings:
                #calc time left of sim:
                time_now            = time.time() - time_start
                time_left           = max_sim_secs - time_now
                #set max wall clock time for Nbody code:
                data_for_worker             = MC_output_list[taskc]
                data_for_worker[1][2][4]    = time_left + (sim_time_buffer/2.)*float(np.random.rand(1))   #1000.*float(np.random.rand(1)) 
                #send data to Worker:
                data_send = data_for_worker
                comm.send(data_send, dest=source, tag=tags.START)
                #update taskc by +1 because we now succesfully send a task:
                taskc += 1
                #TEST: print info:
                print 'sent:    ', taskc, nr_tot_scatterings
                print time_now, time_left, max_sim_secs
            #-------------------------
            #if we do not have more tasks, tell Worker to stop:
            #-------------------------
            else:
                comm.send(None, dest=source, tag=tags.EXIT)
            #-------------------------    
        #-------------------------
        #if Worker is Done:
        #-------------------------
        if tag == tags.DONE:
            result = data
            #collect data:
            all_data.append(result)
            #TEST: print info:
            nr_tasks_completed +=1  #nr completed tasks
            print 'comp:    ', nr_tasks_completed, nr_tot_scatterings
            print result[1][0]
            print (time.time() - time_start)
        #-------------------------    
        #if Worker is closed down:
        #-------------------------
        if tag == tags.EXIT:
            closed_workers += 1
	    print 'closed_workers:  ', closed_workers 
    #--------------------------------------
    
    #---------------------------------------
    #organize output results:
    #---------------------------------------
    print 'SAVING DATA 1:', (time.time() - time_start)
    MC_output_Nbody     = [None]*nr_tot_scatterings 
    for i in range(0,nr_tot_scatterings):
        data_i              = all_data[i]                       # [outlist_MC_info_INT_REAL, output_Nbody]
        #where: outlist_MC_info_INT_REAL = [[icidc,ac,vc,sc,0],[SMA_bin,vinf_sin,b_max,0.0,0.0]] and output_Nbody: [OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL, OUT_out_xtra_3_info_REAL]
        data_i_ori_id       = data_i[0][0][0]                   # =icidc which is the ori IC ID. This runs from 0->nr_tot_scatterings in the (correct) order we made the ICs in.
        data_i_output_Nbody = data_i[1]                         # output_Nbody for this data_i: [OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL, OUT_out_xtra_3_info_REAL]
        MC_output_Nbody[data_i_ori_id] = data_i_output_Nbody    # Here we use 'icidc' (data_i_ori_id) to re-shuffle our data to the initial order.
    #Now the list: MC_output_Nbody has the ori and same order as the IC was created in.
    #---------------------------------------
    #Unpack Nbody data for saving:
    #---------------------------------------
    print 'SAVING DATA 2:', (time.time() - time_start)    
    output_Nbody_endstate_INT       = [] #1x10 per one MC scattering
    output_Nbody_endstate_REAL      = [] #1x10 per one MC scattering
    output_Nbody_xtra_info_INT      = [] #1x10 per one MC scattering
    output_Nbody_xtra_info_REAL     = [] #1x10 per one MC scattering
    output_Nbody_xtra_2_info_REAL   = [] #1x10 per one MC scattering
    output_Nbody_xtra_3_info_REAL   = [] #1x10 per one MC scattering
    for i in range(0,nr_tot_scatterings):
        #MC_output_Nbody[i][:] = [OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL, OUT_out_xtra_3_info_REAL]
        output_Nbody_endstate_INT.append(MC_output_Nbody[i][0])
        output_Nbody_endstate_REAL.append(MC_output_Nbody[i][1])
        output_Nbody_xtra_info_INT.append(MC_output_Nbody[i][2])
        output_Nbody_xtra_info_REAL.append(MC_output_Nbody[i][3])
        output_Nbody_xtra_2_info_REAL.append(MC_output_Nbody[i][4])
        output_Nbody_xtra_3_info_REAL.append(MC_output_Nbody[i][5])
    #Unpack info from 'MC_output_list':
    #MC_output_list          = out[1]    # [outlist_MC_info_INT_REAL, outlist_IC_Nbody], where: outlist_MC_info_INT_REAL = [[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0],[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]] 
    print 'SAVING DATA 3:', (time.time() - time_start)
    outlist_MC_info_INT     = []
    outlist_MC_info_REAL    = []
    for i in range(0,nr_tot_scatterings):
        outlist_MC_info_INT.append(MC_output_list[i][0][0])
        outlist_MC_info_REAL.append(MC_output_list[i][0][1])
    #---------------------------------------          
    #save to .txt files:
    #---------------------------------------
    print 'SAVING DATA 4:', (time.time() - time_start)
    #write data to files in folder:
    tf = open(data_folder+data_name+'MC_settings_list_INT.txt', "w")
    np.savetxt(tf, MC_settings_list_INT[None],   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'MC_settings_list_REAL.txt', "w")
    np.savetxt(tf, MC_settings_list_REAL[None],   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'outlist_MC_info_INT.txt', "w")
    np.savetxt(tf, outlist_MC_info_INT,   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'outlist_MC_info_REAL.txt', "w")
    np.savetxt(tf, outlist_MC_info_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_endstate_INT.txt', "w")
    np.savetxt(tf, output_Nbody_endstate_INT,   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_endstate_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_endstate_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_info_INT.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_info_INT,   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_info_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_info_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_2_info_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_2_info_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_3_info_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_3_info_REAL,   fmt='%5f')
    tf.close()
    #---------------------------------------
    
#------------------------------------------------------------------------------------------    
#if rank == 0:   
#------------------------------------------------------------------------------------------   


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
