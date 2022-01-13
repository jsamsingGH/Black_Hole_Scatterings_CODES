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
import matplotlib as mpl
from matplotlib import cm
import pylab as pl
from matplotlib.patches import Ellipse
from scipy.optimize import fsolve
import matplotlib.gridspec as gridspec


usrinput_sim_or_analyze = raw_input('simulate(=1), analyze(=2): ')
usrinput_sim_or_analyze = int(usrinput_sim_or_analyze)
add_filename            = raw_input('name of in/out data: ')
usrinput_figure         = raw_input('type nr of figure you will look at: ')
usrinput_figure         = int(usrinput_figure)

#data sets:
#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
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
#Specify Properties for the 2 objects:
#------------------------------------------------------------------------------------------
#---------------------------------------
#IC Obj 1:  (binary 1)  TIDAL OBJECT
#---------------------------------------
b1_mass      = 1.0
b1_rad       = 1.0  #0.013*((1.43/b1_mass)**(1./3.))*((1.-b1_mass/1.43)**(0.447))
b1_gas_n     = 3.0
b1_gas_gamma = 5./3.
b1_Mqp       = 0.0376788*(b1_mass*(b1_rad**2))
b1_evoTides_yesno = 1
b1_RigidSph_yesno = 0
b1_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
b1_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
b1_const_arr = np.array([b1_mass, b1_rad, b1_gas_n, b1_gas_gamma, b1_Mqp, b1_evoTides_yesno, b1_RigidSph_yesno, 0,0,0])
#---------------------------------------
#IC Obj 2:  (binary 2)  TIDAL OBJECT or COMPACT OBJECT
#---------------------------------------
b2_mass      = 1.4
b2_rad       = 1.7246e-5
b2_gas_n     = 3.0
b2_gas_gamma = 5./3.
b2_Mqp       = 0.1022999*(b2_mass*(b2_rad**2))
b2_evoTides_yesno = 0
b2_RigidSph_yesno = 1
b2_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])
b2_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]])
b2_const_arr = np.array([b2_mass, b2_rad, b2_gas_n, b2_gas_gamma, b2_Mqp, b2_evoTides_yesno, b2_RigidSph_yesno, 0,0,0])
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Nbody code settings:
#------------------------------------------------------------------------------------------
nbody_params_arr_1_INT  = np.array([0, 0, 0, 100000000, 10, 1, 0,0,0,0])                            # [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, ...]
nbody_params_arr_2_REAL = np.array([0.01, -1.0, 0.0001, 0.0, 1000000000, 0.5, 3.0, 1.0,0,0])        # [scale_dt, max_sim_time (set later in code!!!), evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec (set later in code!!!), IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, ...]
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Time settings:
#------------------------------------------------------------------------------------------
print 'REMEMBER to set time settings in code!!!!'
#Choose which time limit we use (1 = use):
use_fixed_Torb  = int(1)
use_binsin_a0_T = int(0)
#settings for time limits:
Tsim_unit_Torb  = 20.
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Binary-Single example:
#------------------------------------------------------------------------------------------
binsin_ex_a0_AU = 0.1
binsin_ex_a0    = AU_U*binsin_ex_a0_AU
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Make arrays with (a,e):
#------------------------------------------------------------------------------------------
#---------------------------------------
#semi-major axis (SMA):
#---------------------------------------
#SMA Array:
#min_SMA_AU      = [0.1]#1.1*binsin_ex_a0_AU   #0.025          #input in AU
#max_SMA_AU      = [0.1]#2.0*binsin_ex_a0_AU   #0.1            #input in AU     
#nr_SMA          = 1#60             #incl min, max vals
#scan_SMA_arr    = AU_U*np.linspace(min_SMA_AU, max_SMA_AU, nr_SMA)                                         #output units are in correct code units.
#SMA Single values:
SMA_AU_arr      = [0.1] #input in AU
scan_SMA_arr    = AU_U*np.array(SMA_AU_arr)
nr_SMA          = len(scan_SMA_arr)
#---------------------------------------
#---------------------------------------
#peri-center distance:
#---------------------------------------
r_coll          = (b1_rad+b2_rad)
min_rp          = 1.48*r_coll
max_rp          = 2.0*r_coll
nr_rp           = 70             #incl min, max vals
scan_rp_arr     = np.linspace(min_rp, max_rp, nr_rp)                                                       #output units are in correct code units.
#---------------------------------------
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------
#Make final scan combinations of (a,rp):
#------------------------------------------------------------------------------------------
scan_arp_comb_arr    = np.array(list(itertools.product(scan_SMA_arr, scan_rp_arr)))
nr_arp_comb          = len(scan_arp_comb_arr)
#format:
#scan_arp_comb_arr[nr combinations(nr = nr_arp_comb), param values(a, rp, ...)]
#------------------------------------------------------------------------------------------


#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


#----------------------------------------------------------------------------------------------
#Simulate:
#----------------------------------------------------------------------------------------------
if (usrinput_sim_or_analyze == 1):

    #------------------------------------------------------------------------------------------
    #---------------------------------------
    #open file for final output (DATA and INFO):
    #---------------------------------------
    output_DATA_file            = open(add_filename+'_'+'ae_evol_sim_all_DATA_output.txt',  "w")
    output_INFO_file            = open(add_filename+'_'+'ae_evol_sim_all_INFO_output.txt',  "w")
    output_APO_DATA_file        = open(add_filename+'_'+'ae_evol_APOCENTER_all_DATA_output.txt',  "w")
    output_APO_INFO_file        = open(add_filename+'_'+'ae_evol_APOCENTER_all_INFO_output.txt',  "w")
    output_dE_firstpassage_file = open(add_filename+'_'+'ae_evol_dE_firstpassage_output.txt',  "w")   
    #---------------------------------------
    #Initialize counters:
    #---------------------------------------
    totindex_data_start     = -1
    totindex_data_stop      = -1
    totAPOindex_data_start  = -1
    totAPOindex_data_stop   = -1
    #---------------------------------------
    #------------------------------------------------------------------------------------------
    #RUN:
    #------------------------------------------------------------------------------------------
    for pc in range(0,nr_arp_comb):
        
        #---------------------------------------
        #Get binary parameters:
        #---------------------------------------
        #Binary SMA
        ini_SMA_bin = scan_arp_comb_arr[pc,0]  #already in correct code units
        #Binary peri-center dist:
        ini_rp_bin  = scan_arp_comb_arr[pc,1]  #already in correct code units
        #From this calc eccentricity:
        ini_ecc_bin = 1.-(ini_rp_bin/ini_SMA_bin)
        #---------------------------------------

        #---------------------------------------
        #Set up Binary with input (a,e):
        #---------------------------------------
        #In CM of obj 1:
        Mtot    = (b1_mass+b2_mass)
        r_apo   = ini_SMA_bin*(1.+ini_ecc_bin) 
        v_apo   = np.sqrt((Mtot/ini_SMA_bin)*((1.-ini_ecc_bin)/(1.+ini_ecc_bin)))
        #In CM of binary:
        b1_posxyz_binCM   = np.array([ (b2_mass/Mtot)*r_apo,0,0])
        b2_posxyz_binCM   = np.array([-(b1_mass/Mtot)*r_apo,0,0])
        b1_velxyz_binCM   = np.array([0, (b2_mass/Mtot)*v_apo,0])
        b2_velxyz_binCM   = np.array([0,-(b1_mass/Mtot)*v_apo,0])
        #---------------------------------------

        #---------------------------------------
        #Set simulation time:
        #---------------------------------------
        #calc:
        Torb_bin    = 2.*np.pi*np.sqrt((ini_SMA_bin**3.)/(Mtot)) 
        T_iso       = ((2./(1.-(1./(ini_SMA_bin/binsin_ex_a0))))**(3./2.))*2.*np.pi*np.sqrt((binsin_ex_a0**3.)/(b1_mass+b1_mass+b1_mass))
        #set:
        if (use_fixed_Torb  == 1):
            Tsim        = Tsim_unit_Torb*Torb_bin
        if (use_binsin_a0_T == 1):
            print 'EQUAL MASS ONLY FOR THIS EXAMPLE!'
            Tsim        = T_iso
        #set sim time:
        nbody_params_arr_2_REAL[1] = Tsim
        #---------------------------------------
        
        
        #---------------------------------------
        #Print info:
        #---------------------------------------
        print 'INFO:'
        print 'Tsim/Torb:   ', Tsim/Torb_bin
        print 'nr, out of tot nr    :', pc+1, nr_arp_comb
        print 'a,rp,e: ', ini_SMA_bin, ini_rp_bin, ini_ecc_bin
        print 'T_sim: ', Tsim    
        #---------------------------------------
    
    
        #---------------------------------------
        #Write param file to Nbody code:
        #---------------------------------------
        #open file:
        fn = 'MCinput_Nbody.txt'
        text_file = open(fn, "w")
        #nr particles in the sim   
        text_file.write('2' + '\n')
        #nbody code settings: 
        np.savetxt(text_file, nbody_params_arr_1_INT,   fmt='%10i')
        np.savetxt(text_file, nbody_params_arr_2_REAL,  fmt='%10f')
        #body 1:
        np.savetxt(text_file, b1_const_arr,  fmt='%10f')
        np.savetxt(text_file, b1_posxyz_binCM,  fmt='%10f')
        np.savetxt(text_file, b1_velxyz_binCM,  fmt='%10f')
        np.savetxt(text_file, b1_q,  fmt='%10f')
        np.savetxt(text_file, b1_qdot,  fmt='%10f')
        #body 2:
        np.savetxt(text_file, b2_const_arr,  fmt='%10f')
        np.savetxt(text_file, b2_posxyz_binCM,  fmt='%10f')
        np.savetxt(text_file, b2_velxyz_binCM,  fmt='%10f')
        np.savetxt(text_file, b2_q,  fmt='%10f')
        np.savetxt(text_file, b2_qdot,  fmt='%10f')
        #close file:
        text_file.close()    
        #---------------------------------------

        #---------------------------------------
        #Run Nbody_AffineTides_solver.exe with 'MCinput_Nbody.txt' as input:
        #---------------------------------------
        subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
        #---------------------------------------

        #---------------------------------------
        #Read in data from Nbody solver:
        #---------------------------------------
        #File names:
        fn_NbodyTides_dataout_pos           = 'NbodyTides_dataout_pos.dat'
        fn_NbodyTides_dataout_a1a2a3        = 'NbodyTides_dataout_a1a2a3.dat' 
        fn_NbodyTides_dataout_Eself         = 'NbodyTides_dataout_Eself.dat'
        fn_NbodyTides_dataout_Etot          = 'NbodyTides_dataout_Etot.dat'
        fn_NbodyTides_dataout_binij_info    = 'NbodyTides_dataout_binij_info.dat' 
        fn_Nbody_endsim_info_data           = 'Nbody_endsim_info_data.txt'
        #open data:
        tf = open(fn_NbodyTides_dataout_pos,  "r")
        NbodyTides_dataout_pos          = np.loadtxt(tf, dtype=float)
        tf.close()
        tf = open(fn_NbodyTides_dataout_a1a2a3,  "r")
        NbodyTides_dataout_a1a2a3       = np.loadtxt(tf, dtype=float)
        tf.close()
        tf = open(fn_NbodyTides_dataout_Eself,  "r")
        NbodyTides_dataout_Eself        = np.loadtxt(tf, dtype=float)
        tf.close()
        tf = open(fn_NbodyTides_dataout_Etot,  "r")
        NbodyTides_dataout_Etot         = np.loadtxt(tf, dtype=float)
        tf.close()
        tf = open(fn_NbodyTides_dataout_binij_info,  "r")
        NbodyTides_dataout_binij_info   = np.loadtxt(tf, dtype=float)
        tf.close()
        tf = open(fn_Nbody_endsim_info_data,  "r")
        fline_split = tf.readline().split()
        Nbodyout_endsim_end_state_flag  = int(fline_split[0])
        tf.close()
        #---------------------------------------
        #Define
        #---------------------------------------
        #---------------------
        #Format (for 2 objs):
        #---------------------
    	#write(10,*) pos(1,:),  pos(2,:)
    	#write(11,*) time_t, body_all_a1a2a3(1,:),  body_all_a1a2a3(2,:)
    	#write(12,*) time_t, body_all_Eselftot(1),  body_all_Eselftot(2)
    	#write(13,*) time_t, E_tot_kin, E_tot_pot_pmass, E_tot_pot_tides, &
    	# 			E_tot_internal_terms, E_tot_external_terms, E_tot_system
    	#write(14,*) time_t, bin_i, bin_j, a_bin_ij, e_bin_ij, rperi_ij, T_ij, T_ji
        #---------------------
        #pos:
        #---------------------
        b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
        b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
        #---------------------
        #time:
        #---------------------
        time = NbodyTides_dataout_a1a2a3[:,0]
        nr_sim_steps    = len(time)
        #---------------------
        #stellar axes:
        #---------------------
        b1_a1a2a3   =  NbodyTides_dataout_a1a2a3[:,1:4]
        b2_a1a2a3   =  NbodyTides_dataout_a1a2a3[:,4:7]
        #---------------------
        #tot self energy each star:
        #---------------------
        b1_Eselftot = NbodyTides_dataout_Eself[:,1]
        b2_Eselftot = NbodyTides_dataout_Eself[:,2]
        #---------------------
        #total energies:
        #---------------------
        E_tot_kin               = NbodyTides_dataout_Etot[:,1]
        E_tot_pot_pmass         = NbodyTides_dataout_Etot[:,2]
        E_tot_pot_tides         = NbodyTides_dataout_Etot[:,3]
        E_tot_internal_terms    = NbodyTides_dataout_Etot[:,4]
        E_tot_external_terms    = NbodyTides_dataout_Etot[:,5]
        E_tot_system            = NbodyTides_dataout_Etot[:,6]
        #---------------------
        #Binary info:
        #---------------------
        bin_i       = NbodyTides_dataout_binij_info[:,1]
        bin_j       = NbodyTides_dataout_binij_info[:,2]
        a_bin_ij    = NbodyTides_dataout_binij_info[:,3]
        e_bin_ij    = NbodyTides_dataout_binij_info[:,4]
        rperi_ij    = NbodyTides_dataout_binij_info[:,5]
        T_ij        = NbodyTides_dataout_binij_info[:,6]
        T_ji        = NbodyTides_dataout_binij_info[:,7]
        #---------------------
        #---------------------------------------
        #Calc bin evolution info to all times t:
        #---------------------------------------
        #calc:
        mtot12      = b1_mass+b2_mass
        mred12      = b1_mass*b2_mass/mtot12
        Eorb_arr    = -b1_mass*b2_mass/(2.*a_bin_ij)
        Lorb_arr    = mred12*rperi_ij*np.sqrt(2.*mtot12*(1./(rperi_ij) - 1./(2.*a_bin_ij)))
        dEorb_arr   = Eorb_arr-Eorb_arr[0]
        dLorb_arr   = Lorb_arr-Lorb_arr[0]
        pos12xyz_t  = b1_posxyz - b2_posxyz                                                     #rel pos of b1 and b2 at all t.
        rdist12_t   = np.sqrt(pos12xyz_t[:,0]**2 + pos12xyz_t[:,1]**2 + pos12xyz_t[:,2]**2)     #r dist between b1 and b2 at all t.
        fractid_t   = abs((b2_mass/b1_mass)*(((rdist12_t/b1_rad)-1.)**(-2.) - ((rdist12_t/b1_rad)+1.)**(-2.)))         
        max_a_t     = np.array([max([b1_a1a2a3[i,0],b1_a1a2a3[i,1],b1_a1a2a3[i,2]]) for i in range(0,nr_sim_steps)])
        #---------------------------------------
        #Linear estimates of dE and dL:
        #---------------------------------------
        #calc dE/dL: (ASSUMING IDENTICAL OBJECTS)
        io_gas_n   = b1_gas_n
        io_mass    = b1_mass
        io_radius  = b1_rad 
        Epot_sph_star       = - (3./(5.-io_gas_n))*(io_mass*io_mass/io_radius)
        MInertia_sph_star   = b1_Mqp
        dEodL_lin_fac    = (abs(Epot_sph_star)/np.sqrt(15.))*(1./np.sqrt(MInertia_sph_star*abs(Epot_sph_star)))
        #calc PT dE (orbit) from initial rp: (ASSUMES obj1 IS THE ONLY TIDAL OBJECT - but works for any mass ratio)
        rp_ini  = rperi_ij[0]
        eta_imp = ((b1_mass/(b1_mass+b2_mass))**(1./2.))*((rp_ini/b1_rad)**(3./2.)) #(1./np.sqrt(2.))*((rp_ini/io_radius)**(3./2.))
        log10eta    = np.log10(eta_imp)
        #For n=1.5:
        if (b1_gas_n == 1.5):
            print 'PT: n=', b1_gas_n
            fA = -0.397
            fB = 1.678
            fC = 1.277
            fD = -12.42
            fE = 9.446
            fF = -5.550
        #For n=3.0:
        if (b1_gas_n == 3.):
            print 'PT: n=', b1_gas_n
            fA = -1.124
            fB = 0.877
            fC = -13.37
            fD = 21.55
            fE = -16.48
            fF = 4.124
        log10T2     = fA + fB*(log10eta**(1.)) + fC*(log10eta**(2.)) + fD*(log10eta**(3.)) + fE*(log10eta**(4.)) + fF*(log10eta**(5))
        T2_PT          = 10.**(log10T2)
        dEorb_PT_tides  = ((b2_mass**2.)/b1_rad)*((b1_rad/rp_ini)**(6.))*T2_PT
        #---------------------------------------        
        #---------------------------------------
        #APO-CENTER analysis:
        #---------------------------------------
        #initialize:
        Ncounter_apo    = 0
        apo_center_bin_info = np.array([Ncounter_apo, time[0], a_bin_ij[0], rperi_ij[0], e_bin_ij[0]])  #we start sim at apo-center.
        #loop over all t to find each apo-center:
        for tc in range(1,nr_sim_steps-1):
            rdist_tm1   = rdist12_t[tc-1]
            rdist_t0    = rdist12_t[tc+0]
            rdist_tp1   = rdist12_t[tc+1]
            #apo-center:
            if (rdist_t0 > rdist_tm1 and rdist_t0 > rdist_tp1):
                #get info:
                Ncounter_apo    = Ncounter_apo+1
                t_apo           = time[tc]
                SMA_apo         = a_bin_ij[tc]
                rperi_apo       = rperi_ij[tc]
                ecc_apo         = e_bin_ij[tc]
                #save:
                apo_center_save_arr = np.array([Ncounter_apo, t_apo, SMA_apo, rperi_apo, ecc_apo])
                apo_center_bin_info = np.vstack([apo_center_bin_info,apo_center_save_arr])
        #apo_center_bin_info = apo_center_bin_info[1::,:]
        #define:
        nr_apo_centers      = Ncounter_apo+1
        Ncounter_apo_arr    = apo_center_bin_info[:,0]
        time_apo_arr        = apo_center_bin_info[:,1]
        SMA_apo_arr         = apo_center_bin_info[:,2]
        rperi_apo_arr       = apo_center_bin_info[:,3]
        ecc_apo_arr         = apo_center_bin_info[:,4]
        #calc:
        mtot12  = b1_mass+b2_mass
        mred12  = b1_mass*b2_mass/mtot12
        Eorb_apo_arr        = -b1_mass*b2_mass/(2.*SMA_apo_arr)
        Lorb_apo_arr        = mred12*rperi_apo_arr*np.sqrt(2.*mtot12*(1./(rperi_apo_arr) - 1./(2.*SMA_apo_arr)))
        dEorb_apo_arr       = Eorb_apo_arr-Eorb_apo_arr[0]
        dLorb_apo_arr       = Lorb_apo_arr-Lorb_apo_arr[0]
        pos_time_1pass      = np.where(time < time_apo_arr[1])[0]
        max_fractid_1pass   = max(fractid_t[pos_time_1pass])
        max_max_a_1pass     = max(max_a_t[pos_time_1pass])
        #save info:
        dE_firstpassage = dEorb_apo_arr[1] - dEorb_apo_arr[0]
        T2_SIM  = abs(dE_firstpassage)/(((b2_mass**2.)/b1_rad)*(b1_rad/ini_rp_bin)**(6.))
        dE_firstpassage_save_arr = np.array([ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides, max_fractid_1pass, max_max_a_1pass, T2_PT, T2_SIM, eta_imp])
        #---------------------------------------
        #---------------------------------------    
        
        
        #---------------------------------------
        #Energy tests:
        #---------------------------------------
        Escale1 = abs(dE_firstpassage)*(b1_rad/(b2_mass**2.))*(ini_rp_bin/b1_rad)**(6.)
        print 'Escale1:  ', Escale1
        xf  = 3.5
        Mf  = b2_mass*((b1_mass+b2_mass)/b1_mass)**(xf/4.)
        Rf  = b1_rad
         
        Escale2 = abs(dE_firstpassage)*(((Mf**2.)/Rf)*(Rf/ini_rp_bin)**(6.+3.*xf/2.))**(-1.)
        print 'Escale2:  ', Escale2

        print 'Ecorrect:    ', abs(dE_firstpassage)/((((b2_mass**2.)/b1_rad)*(b1_rad/ini_rp_bin)**(6.))*T2_PT)
        print 'dE_firstpassage: ', dE_firstpassage
        #---------------------------------------
        
        
        #-------------------------------------------
        #-------------------------------------------
        #Figure 3:
        #-------------------------------------------
        #-------------------------------------------
        if (usrinput_figure == 3):

            #---------------------------------------
            #PAPER figure:
            #---------------------------------------
            f = plt.figure(figsize=(5,6.5))
            gs = gridspec.GridSpec(2, 1,height_ratios=[3,1])
            ax1 = plt.subplot(gs[0])
            ax2 = plt.subplot(gs[1])
            
            color_vals = [0,100,150]
        
            time_apo_unit_iniTorb = time_apo_arr/Torb_bin
        
            #subplot 1:
            ax1.plot(time_apo_unit_iniTorb, SMA_apo_arr/SMA_apo_arr[0], drawstyle='steps-mid', linewidth=1.5, color='black', linestyle='-', label=r'$a/a_{0}$')
            ax1.plot(time_apo_unit_iniTorb, rperi_apo_arr/rperi_apo_arr[0], drawstyle='steps-mid', linewidth=1.5, color='black', linestyle='--', label=r'$r_{p}/r_{p,0}$')
            ax1.plot(time_apo_unit_iniTorb, ecc_apo_arr, drawstyle='steps-mid', linewidth=1.5, color='black', linestyle=':', label=r'$e$')      
            ax1.plot(time_apo_unit_iniTorb, dEorb_apo_arr, drawstyle='steps-mid', linewidth=1.5, color='purple', linestyle=':', label=r'dE')      
          
            #ax1.plot(time_apo_unit_iniTorb, SMA_apo_arr/SMA_apo_arr[0], drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color=plt.cm.gnuplot(color_vals[0]), label=r'semi-major axis $a(t)/a(t=0)$')
            #ax1.plot(time_apo_unit_iniTorb, rperi_apo_arr/rperi_apo_arr[0], drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color=plt.cm.gnuplot(color_vals[1]), label=r'peri-center $r_{p}(t)/r_{p}(t=0)$')
            #ax1.plot(time_apo_unit_iniTorb, ecc_apo_arr, drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color=plt.cm.gnuplot(color_vals[2]), label=r'eccentricity $e(t)$')      
            #TEST lines:
            #ax1.plot(time/Torb_bin, b1_a1a2a3[:,0], drawstyle='steps-mid', alpha=0.25, linewidth=1.5, color='red')      
            #ax1.plot(time/Torb_bin, b1_a1a2a3[:,1], drawstyle='steps-mid', alpha=0.25, linewidth=1.5, color='orange')      
            #ax1.plot(time/Torb_bin, b1_a1a2a3[:,2], drawstyle='steps-mid', alpha=0.25, linewidth=1.5, color='blue')      
            #ax1.plot(time/Torb_bin, rdist12_t/(b1_a1a2a3[:,2]+b2_a1a2a3[:,2]), drawstyle='steps-mid', alpha=0.25, linewidth=1.5, color='orange')
            #ax1.plot(time/Torb_bin, Lorb_arr/Lorb_arr[0], drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color='green')      
            #ax1.plot([0,100], [0.75,0.75], color='green', linestyle='--', linewidth=1.5)
            #ax1.plot(time/Torb_bin, rdist12_t/a_bin_ij[0], drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color=plt.cm.gnuplot(color_vals[0]))
            #ax1.plot(time/Torb_bin, e_bin_ij, drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color=plt.cm.gnuplot(color_vals[2]))
            #ax1.plot(time/Torb_bin, 1.+dEorb_arr/2., drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color='red')
            #ax1.plot(time/Torb_bin, a_bin_ij/a_bin_ij[0], drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color='black')
            #ax1.plot(time/Torb_bin, a_bin_ij/(b1_a1a2a3[:,2]+b2_a1a2a3[:,2]), drawstyle='steps-mid', alpha=0.5, linewidth=1.5, color='yellow')
            #ax1.plot(time/Torb_bin, T_ij[:], drawstyle='steps-mid', alpha=0.25, linewidth=1.5, color='black')      
            
            #upper lower boarder lines on SMA a:
            #upper:
            ax1.plot([0,Tsim_unit_Torb], [1,1], color='black', linestyle='-', linewidth=0.25)#, label=r'analytical min,max $a(t)/a(t=0)$')
            #lower:        
            psolv_p0 = 1.0
            psolv_p1 = -np.sqrt(2.)*(np.sqrt(ini_rp_bin/io_radius))
            psolv_p2 = 0.0
            psolv_p3 = (15.*MInertia_sph_star*io_mass/(2.*abs(Epot_sph_star)*(io_radius**3.)))**(1./2.)
            psolv_coeff = [psolv_p0, psolv_p1, psolv_p2, psolv_p3]
            r_circular  = ((max(np.roots(psolv_coeff)))**2.)*io_radius
            ax1.plot([0,Tsim_unit_Torb], [r_circular,r_circular]/SMA_apo_arr[0], color='black', linestyle='-', linewidth=0.25)
            print 'min a:', r_circular
           
            #analytical estimate for SMA a(t):    
            Mtot12  = b1_mass + b2_mass
            mu12    = b1_mass*b2_mass/Mtot12
            gamma   = (1./2.)*np.sqrt(Mtot12*mu12)*(dEorb_PT_tides/Mtot12)*(1./np.pi)*(mu12**(-3./2.))   #(dEorb_PT_tides/(io_mass**(3./2.)))*(1./(np.pi*np.sqrt(2.)))
            t_arr   = np.arange(0,1.1*np.amax(time_apo_arr),(np.amax(time_apo_arr)/10000.)) 
            SMA_a_t_analesti = (np.sqrt(SMA_apo_arr[0])-gamma*t_arr)**2.        
            posmin  = 0
            posmax  = min([np.amin(np.append(np.where(SMA_a_t_analesti < r_circular),1e10)), np.amax(np.where(SMA_a_t_analesti > r_circular))])+1     
            plot_t_arr_unit_iniTorb = t_arr[posmin:posmax]/Torb_bin
            plot_SMA_a_t_analesti   = SMA_a_t_analesti[posmin:posmax]
            nrpos                   = len(plot_t_arr_unit_iniTorb)
            ax1.plot(t_arr/Torb_bin, SMA_a_t_analesti/SMA_apo_arr[0], alpha=0.5, linewidth=1.0, color=plt.cm.gnuplot(color_vals[0]), linestyle='--', label=r'analytical estimate $a(t)/a(t=0)$')
            plot_e_t_analesti       = 1. - (ini_rp_bin/plot_SMA_a_t_analesti)   #in this ana esti we assume rp (dE) to be constant (eq to ini val) during the evol.
            #ax1.plot(plot_t_arr_unit_iniTorb, plot_e_t_analesti, alpha=0.5, linewidth=1.0, color=plt.cm.gnuplot(color_vals[2]), linestyle='--', label=r'analytical estimate $e(t)$')
            #ax1.plot(plot_t_arr_unit_iniTorb[nrpos-1], plot_SMA_a_t_analesti[nrpos-1]/SMA_apo_arr[0], marker='x', color=plt.cm.gnuplot(color_vals[0]))
      
            #isolation time: T_iso/T_orb:
            T_iso_unit_iniTorb = T_iso/Torb_bin
            #ax1.plot([T_iso_unit_iniTorb,T_iso_unit_iniTorb], [-10,10], color='black', linestyle='-.', linewidth=0.5, label=r'isolation time $T_{iso}$'+' ('+'$a/a_{0}$ = '+str(round(ini_SMA_bin/binsin_ex_a0,3))+')')
            #print (1. - abs(dEorb_PT_tides/Eorb_arr[0])*np.sqrt(4./(3.*(ini_SMA_bin/binsin_ex_a0 - 1.)**3)))**2.
        
            #insert figure showing orbits:
            subfig = plt.axes([.26, .625, 0.4, .25], axisbg='white')
            subfig.plot(pos12xyz_t[::10,0],pos12xyz_t[::10,1], color='black', linewidth=0.5)
            #subfig.plot(b1_posxyz[::10,0],b1_posxyz[::10,1], color='black', linewidth=0.5)
            #subfig.plot(b2_posxyz[::10,0],b2_posxyz[::10,1], color='red', linewidth=0.5)
            subfig.add_patch(Ellipse((0, 0), b1_rad, b1_rad, color='black'))
            subfig.add_patch(Ellipse((pos12xyz_t[0,0], pos12xyz_t[0,1]), b2_rad, b2_rad, color='black'))
            subfig.axis('equal')
            subfig.set_xlabel(r'x pos $[R_{\odot}]$')
            subfig.set_ylabel(r'y pos $[R_{\odot}]$')
            
            #plt.setp(subfig, xticks=[], yticks=[])
                
            #set x,y limits:
            ax1.set_ylim(0,3.0)
            ax1.set_xlim(0,np.max(1.1*time/Torb_bin))
            #ax1.set_xlim(0,np.amax(np.append(1.1*time/Torb_bin,1.1*T_iso_unit_iniTorb)))
        
            #set labels:
            #ax1.set_xlabel(r'time $[T_{orb}(t=0)]$')
            ax1.set_title(r'Evolution of a tidal capture')
            ax1.set_ylabel(r'orbital parameters $a$, $r_{p}$, $e$')
            ax1.legend(loc='upper right', numpoints = 1, fontsize = 12.0, frameon = False)
            
            pt = 1#50 #plot thinning (plot every ...)
            ax2.plot(time[::pt]/Torb_bin, b1_a1a2a3[::pt,0], drawstyle='steps-mid', alpha=0.50, linewidth=0.5, color='red')      
            ax2.plot(time[::pt]/Torb_bin, b1_a1a2a3[::pt,1], drawstyle='steps-mid', alpha=0.50, linewidth=0.5, color='orange')      
            ax2.plot(time[::pt]/Torb_bin, b1_a1a2a3[::pt,2], drawstyle='steps-mid', alpha=0.50, linewidth=0.5, color='blue')      
            ax2.set_xlabel(r'time $[T_{orb,0}]$')
            ax2.set_ylabel(r'$a_{1}$, $a_{2}$, $a_{3}$ $[R_{star}]$')
            ax2.set_xlim(0,np.max(1.1*time/Torb_bin))
            #ax2.set_xlim(0,1.2)
            ax2.set_ylim(0.5,2.0)
            
            #show and save:
            plt.savefig('2body_evolution_example.eps', bbox_inches='tight')
            
            plt.show()
            
            #exit()
            #save_yes_no = raw_input('if save as pdf and quit press 1: ')
            #if (int(save_yes_no) == 1):
            #    figname = raw_input('figure name: ')
            #    plt.show()
            #    plt.savefig(figname+'.eps', bbox_inches='tight')
            #---------------------------------------
        
        
        
        
        
        
        
        
        
            #---------------------------------------
            #TEST FIGURES:
            #---------------------------------------
            print 'CALC dE TIDES:'
            print dEorb_PT_tides
            print abs(Eorb_apo_arr[0]-Eorb_apo_arr[1])
        
            fig = plt.figure(figsize=(8.0, 3.0))

            pos12xyz_t =  b1_posxyz - b2_posxyz
            fig.add_subplot(331).plot(pos12xyz_t[:,0],pos12xyz_t[:,1])
        
            fig.add_subplot(332).plot(Ncounter_apo_arr, SMA_apo_arr/SMA_apo_arr[0], drawstyle='steps-mid', color='blue')
            fig.add_subplot(332).plot(Ncounter_apo_arr, rperi_apo_arr/rperi_apo_arr[0], drawstyle='steps-mid', color='red')
            fig.add_subplot(332).plot(Ncounter_apo_arr, ecc_apo_arr, drawstyle='steps-mid', color='black')
            fig.add_subplot(332).plot(Ncounter_apo_arr, Lorb_apo_arr/Lorb_apo_arr[0], drawstyle='steps-mid', color='green')
            #IDENTICAL OBJECTS CASE:
            dEodL        = (dEorb_apo_arr/dLorb_apo_arr)
            fig.add_subplot(332).plot(Ncounter_apo_arr, abs(dEodL/dEodL_lin_fac), drawstyle='steps-mid', color='orange')
            fig.add_subplot(332).set_ylim(0,2)
            a2 = (1./2.)*(rperi_apo_arr**2./(rperi_apo_arr - rperi_apo_arr[0]**2*(1./rperi_apo_arr[0] - 1./(2.*SMA_apo_arr[0]))))
            fig.add_subplot(332).plot(Ncounter_apo_arr, a2/SMA_apo_arr[0], drawstyle='steps-mid', color='yellow')
                 
        
            print Lorb_apo_arr[0]**2, rperi_apo_arr[0]
        
            fig.add_subplot(332).plot(Ncounter_apo_arr, (rperi_apo_arr[0] + abs(2*Lorb_apo_arr[0]*dEorb_apo_arr/dEodL_lin_fac))/rperi_apo_arr[0], drawstyle='steps-mid', color='gray')
        
            fig.add_subplot(333).plot(time_apo_arr, SMA_apo_arr/SMA_apo_arr[0], drawstyle='steps-mid', color='blue')
            fig.add_subplot(333).plot(time_apo_arr, rperi_apo_arr/rperi_apo_arr[0], drawstyle='steps-mid', color='red')
            fig.add_subplot(333).plot(time_apo_arr, ecc_apo_arr, drawstyle='steps-mid', color='black')
            fig.add_subplot(333).plot(time_apo_arr, Lorb_apo_arr/Lorb_apo_arr[0], drawstyle='steps-mid', color='green')
            #IDENTICAL OBJECTS CASE:
            dEodL        = (dEorb_apo_arr/dLorb_apo_arr)
            fig.add_subplot(333).plot(time_apo_arr, abs(dEodL/dEodL_lin_fac), drawstyle='steps-mid', color='orange')
        
            gamma = (dEorb_PT_tides/(io_mass**(3./2.)))*(1./(np.pi*np.sqrt(2.)))
            at = (np.sqrt(SMA_apo_arr[0])-gamma*time_apo_arr)**2.
            fig.add_subplot(333).plot(time_apo_arr, at/SMA_apo_arr[0], drawstyle='steps-mid', color='purple')
    
            fig.add_subplot(333).set_ylim(0,2)
        
            fig.add_subplot(334).plot(time, a_bin_ij/a_bin_ij[0], color='blue')
            fig.add_subplot(334).plot(time, rperi_ij/rperi_ij[0], color='red')
            fig.add_subplot(334).plot(time, e_bin_ij, color='black')
            fig.add_subplot(334).plot(time, Lorb_arr/Lorb_arr[0], color='green')
            fig.add_subplot(334).set_ylim(0,2)
        
            fig.add_subplot(335).plot(time, dEorb_arr, color='blue')
            fig.add_subplot(335).plot(time, dLorb_arr, color='black')
        
            fig.add_subplot(336).plot(time, rperi_ij, color='blue')
            fig.add_subplot(336).plot(time, rdist12_t, color='red')
            fig.add_subplot(336).plot(time, a_bin_ij, color='black')
            fig.add_subplot(336).plot(time_apo_arr, SMA_apo_arr, drawstyle='steps-mid', color='black')
        
            rc_arr  = np.arange(2,10,0.01)
            rp_arr  = (1./2.)*((np.sqrt(rc_arr) + 0.95/rc_arr)**2.)
            fig.add_subplot(337).plot(rc_arr, rp_arr, color='black')
            fig.add_subplot(337).set_xlabel(r'final circulation radius')
            fig.add_subplot(337).set_ylabel(r'initial peri-center distance')
            
           
            #fig.add_subplot(338).plot(time, rdist12_t/(b1_a1a2a3[:,0]+b2_a1a2a3[:,0]), color='black')
            #fig.add_subplot(338).plot(time, rdist12_t/(b1_a1a2a3[:,1]+b2_a1a2a3[:,1]), color='red')
            #fig.add_subplot(338).plot(time, rdist12_t/(b1_a1a2a3[:,2]+b2_a1a2a3[:,2]), color='blue')
            
            plt.show()
            print Nbodyout_endsim_end_state_flag
            test = raw_input('test')
            
            
            fig = plt.figure(figsize=(5, 5))
            
            #fig.add_subplot(339).plot(time, E_tot_kin, color='black')
            #fig.add_subplot(339).plot(time, E_tot_pot_pmass, color='red')
            #fig.add_subplot(339).plot(time, E_tot_pot_tides, color='blue')
            #fig.add_subplot(339).plot(time, E_tot_internal_terms, color='orange')
            #fig.add_subplot(339).plot(time, E_tot_external_terms, color='green')
            #fig.add_subplot(339).plot(time, E_tot_system, color='gray')
              
            #fig.add_subplot(339).plot(time, b1_a1a2a3[:,0], color='black')
            #fig.add_subplot(339).plot(time, b1_a1a2a3[:,1], color='red')
            #fig.add_subplot(339).plot(time, b1_a1a2a3[:,2], color='yellow')
               
            phi_tot = Lorb_arr**2/(rdist12_t**2) + E_tot_pot_pmass + E_tot_pot_tides
            fig.add_subplot(111).plot(rdist12_t, phi_tot, color='black')
            fig.add_subplot(111).plot(rdist12_t, E_tot_pot_pmass, color='red')
            fig.add_subplot(111).plot(rdist12_t, E_tot_pot_tides, color='blue')
            #fig.add_subplot(111).plot(rdist12_t, Lorb_arr**2/(rdist12_t**2) + E_tot_pot_pmass, color='gray')
        
            dEtides_arr = np.arange(0.0,0.2,0.0001)
            dLtides_arr = dEtides_arr/dEodL_lin_fac
            E_arr       = Eorb_arr[0] - dEtides_arr
            L_arr       = Lorb_arr[0] - dLtides_arr
            e_arr       = np.sqrt(1+4.*E_arr*(L_arr**2))
            a_arr       = -1./(2*E_arr)
            print e_arr
            rp_arr = a_arr*(1-e_arr)
            ra_arr = a_arr*(1+e_arr)
            fig.add_subplot(111).plot(rp_arr, E_arr, color='red')
            fig.add_subplot(111).plot(ra_arr, E_arr, color='red')
              
            r_arr = np.arange(1,20,0.01)
            rc_in = 3.5
            Veff_arr_1 = (1./2.)*(rc_in/(r_arr**2))
            Veff_arr_2 = - 1./r_arr
            Veff_arr_3 = - 2.*0.51*(1.0**5)/(r_arr**6.)
            Veff_arr = Veff_arr_1 + Veff_arr_2 + Veff_arr_3
            #fig.add_subplot(339).plot(r_arr, Veff_arr, color='orange')
            #fig.add_subplot(339).plot(r_arr, Veff_arr_1, color='red')
            #fig.add_subplot(339).plot(r_arr, Veff_arr_2, color='blue')
            #fig.add_subplot(339).plot(r_arr, Veff_arr_3, color='gray')
            #fig.add_subplot(339).plot(r_arr, Veff_arr_1+Veff_arr_2, color='purple')
            
            fig.add_subplot(111).set_xlim(0,30)
            fig.add_subplot(111).set_ylim(-0.20, 0.05)
            
            plt.show()
            save_yes_no = raw_input('if save as pdf and quit press 1: ')
            save_yes_no = int(save_yes_no)
            if (save_yes_no == 1):
                figname = raw_input('figure name: ')
                fig.savefig(figname+'.eps')
                   
            #---------------------------------------
            #---------------------------------------   
    
    
    
    
    
        #---------------------------------------
        #Save output to file:
        #---------------------------------------
        #---------------------
        #save DATA:
        #---------------------
        save_data_arr   = np.zeros((nr_sim_steps,5), dtype=np.float64)
        save_data_arr[:,0] = time[:]
        save_data_arr[:,1] = a_bin_ij[:]
        save_data_arr[:,2] = e_bin_ij[:]
        np.savetxt(output_DATA_file, save_data_arr, fmt='%10f')
        #---------------------
        #save INFO:
        #---------------------
        totindex_data_start = totindex_data_stop + 1
        totindex_data_stop  = totindex_data_start + nr_sim_steps - 1
        initial_sma         = a_bin_ij[0]
        final_sma           = a_bin_ij[nr_sim_steps-1] 
        save_info_arr = np.array([ini_SMA_bin, ini_rp_bin, ini_ecc_bin, totindex_data_start, totindex_data_stop, Nbodyout_endsim_end_state_flag, initial_sma, final_sma])
        np.savetxt(output_INFO_file, save_info_arr[None], fmt='%5f')
        print 'HERE:', save_info_arr, totindex_data_start, totindex_data_stop
        print 'a_ini/a_fin: ', initial_sma/final_sma
        #---------------------
        #---------------------------------------
        
        #---------------------------------------
        #Save APOCENTER DATA to file:
        #---------------------------------------
        #---------------------
        #save APOCENTER-DATA:
        #---------------------
        save_data_arr   = np.zeros((nr_apo_centers,5), dtype=np.float64)
        save_data_arr[:,0] = Ncounter_apo_arr[:]
        save_data_arr[:,1] = time_apo_arr[:]
        save_data_arr[:,2] = SMA_apo_arr[:]
        save_data_arr[:,3] = rperi_apo_arr[:]
        save_data_arr[:,4] = ecc_apo_arr[:]
        np.savetxt(output_APO_DATA_file, save_data_arr, fmt='%10f')
        #---------------------
        #save APOCENTER-INFO:
        #---------------------
        totAPOindex_data_start = totAPOindex_data_stop + 1
        totAPOindex_data_stop  = totAPOindex_data_start + nr_apo_centers - 1
        save_info_arr = np.array([ini_SMA_bin, ini_rp_bin, ini_ecc_bin, totAPOindex_data_start, totAPOindex_data_stop, Nbodyout_endsim_end_state_flag])
        np.savetxt(output_APO_INFO_file, save_info_arr[None], fmt='%5f')
        #---------------------
        #first passage dE:
        #---------------------
        #output_dE_firstpassage_file = open(add_filename+'_'+'ae_evol_dE_firstpassage_output.txt',  "w")   
        #dE_firstpassage_save_arr = np.array([ini_SMA_bin, ini_rp_bin, dE_firstpassage])
        np.savetxt(output_dE_firstpassage_file, dE_firstpassage_save_arr[None], fmt='%5f') 
        #---------------------
        #---------------------------------------


    #------------------------------------------------------------------------------------------
    #---------------------------------------
    #close output file:
    #---------------------------------------
    output_DATA_file.close()    #data
    output_INFO_file.close()    #info
    output_APO_DATA_file.close()
    output_APO_INFO_file.close()
    output_dE_firstpassage_file.close()
    #---------------------------------------
    #------------------------------------------------------------------------------------------
#----------------------------------------------------------------------------------------------




#exit()

#------------------------------------------------------------------------------------------
#Analyze data:
#------------------------------------------------------------------------------------------
#---------------------------------------
#Set general plot settings incl. font size etc.:
#---------------------------------------
font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 10.5}
mpl.rc('font', **font)
#---------------------------------------

#---------------------------------------
#read data and info files:
#---------------------------------------
#----------------
#FULL SIMULATION:
#----------------
#DATA:	
input_DATA_file         = open(add_filename+'_'+'ae_evol_sim_all_DATA_output.txt',  "r")
if (usrinput_figure == 1): ae_evol_sim_all_DATA    = np.loadtxt(input_DATA_file, dtype=float)  #[time,a,e]
if (usrinput_figure == 2): print 'we dont read full data for FIG: ', usrinput_figure
input_DATA_file.close()
#INFO:
input_INFO_file         = open(add_filename+'_'+'ae_evol_sim_all_INFO_output.txt',  "r")
ae_evol_sim_all_INFO    = np.loadtxt(input_INFO_file, dtype=float)  #[ini_SMA_bin, ini_rp_bin, ini_ecc_bin, totindex_data_start, totindex_data_stop, Nbodyout_endsim_end_state_flag]
input_INFO_file.close()
#----------------
#APOCENTER ANALYSIS:
#----------------
if (usrinput_figure == 4):
    #DATA:
    input_APO_DATA_file         = open(add_filename+'_'+'ae_evol_APOCENTER_all_DATA_output.txt',  "r")
    ae_evol_APOCENTER_all_DATA  = np.loadtxt(input_APO_DATA_file, dtype=float)  #[Ncounter_apo_arr,time_apo_arr,SMA_apo_arr,rperi_apo_arr,ecc_apo_arr]
    input_APO_DATA_file.close()
    #INFO:
    input_APO_INFO_file         = open(add_filename+'_'+'ae_evol_APOCENTER_all_INFO_output.txt',  "r")
    ae_evol_APOCENTER_all_INFO  = np.loadtxt(input_APO_INFO_file, dtype=float)  #[ini_SMA_bin, ini_rp_bin, ini_ecc_bin, totindex_data_start, totindex_data_stop, Nbodyout_endsim_end_state_flag]
    input_APO_INFO_file.close()
#----------------    
#---------------------------------------
#---------------------------------------
#Define:
#---------------------------------------
nr_tot_arp_comb = len(ae_evol_sim_all_INFO)
#from info file restore a,rp scan arrays etc.:
scan_SMA_arr    = np.unique(ae_evol_sim_all_INFO[:,0])
scan_rp_arr     = np.unique(ae_evol_sim_all_INFO[:,1])
nr_SMA          = len(scan_SMA_arr)
nr_rp           = len(scan_rp_arr)
if (nr_SMA  > 1): delta_SMA = scan_SMA_arr[1]-scan_SMA_arr[0]   #assumes constant step size in a.  Change if necesseary.
if (nr_rp   > 1): delta_rp  = scan_rp_arr[1]-scan_rp_arr[0]     #assumes constant step size in rp. Change if necesseary.
print 'scan_SMA_arr, scan_rp_arr, nr_SMA, nr_rp :'
print scan_SMA_arr, scan_rp_arr, nr_SMA, nr_rp
#---------------------------------------






#-------------------------------------------
#plot (dE,rp)
#-------------------------------------------

#DATA:
input_dE_firstpassage_file  = open(add_filename+'_'+'ae_evol_dE_firstpassage_output.txt',  "r")
ae_evol_dE_firstpassage_arr = np.loadtxt(input_dE_firstpassage_file, dtype=float)  #[ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides/2.]
input_dE_firstpassage_file.close()
#dE_firstpassage = dEorb_apo_arr[1] - dEorb_apo_arr[0]
#dE_firstpassage_save_arr = np.array([ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides, max_fractid_1pass, max_max_a_1pass, T2_PT, T2_SIM, eta_imp])

#NOTE: FOR THIS PLOT CHOOSE:
# CHOOSE: obj 1: MS(sun), obj 2: CO(M_sun, R=0).
# - just one value for a (does not really play any role)
# - several values for rp
# - same mass objects is best (standard 1R,1M sun).
# - one of the objects to be a compact obj to probe below 2R!!! ASSUME THIS BELOW!!

#define from data file: ae_evol_dE_firstpassage_arr [nr (a,rp)comb, [ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides, max_fractid_1pass, max_max_a_1pass, T2_PT, T2_SIM, eta_imp]]
rp_arr          = ae_evol_dE_firstpassage_arr[:,1]
dE1pass_arr     = abs(ae_evol_dE_firstpassage_arr[:,2])
dE1pass_PT_arr  = abs(ae_evol_dE_firstpassage_arr[:,3])
T2_PT_arr       = abs(ae_evol_dE_firstpassage_arr[:,6])
T2_SIM_arr      = abs(ae_evol_dE_firstpassage_arr[:,7])
eta_arr         = abs(ae_evol_dE_firstpassage_arr[:,8])

#define:
#dE1pass_beta_model  = 0.1*(1./rp_arr)**(6.)

fig = plt.figure(figsize=(8, 10))

#plot 1: (dE, rp)
#affine sim:
fig.add_subplot(211).plot(rp_arr, dE1pass_arr, marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, color='black', label=r'TEST')
#PT analytical:
fig.add_subplot(211).plot(rp_arr, dE1pass_PT_arr, linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75, color='black', label=r'PT')
#beta model:
#fig.add_subplot(211).plot(rp_arr, dE1pass_beta_model, marker='o', linestyle='--', markersize=5.0, linewidth=1.5, alpha=0.75, color='black')
#axis settings:
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)
plt.yscale('log')
plt.xlabel(r'peri-center $r_{p} [R_{sun}]$')
plt.ylabel(r'$\Delta{E}$ single passage')

#plot 2:
fig.add_subplot(212).plot(np.log10(eta_arr), np.log10(T2_PT_arr), marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, label= 'T2_PT')
fig.add_subplot(212).plot(np.log10(eta_arr), np.log10(T2_SIM_arr), marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, label= 'T2_SIM')

plt.show()
exit()

##plot 2:
#a0bin_AU_arr    = np.array([0.1, 1.0, 10.0])
#a0bin_arr       = AU_U*a0bin_AU_arr
#E0bin_arr       = b1_mass*b2_mass/(2.*a0bin_arr)
#nr_a0           = len(a0bin_AU_arr)
#for ac in range(0,nr_a0):
#    E0bin       = E0bin_arr[ac]
#    Efrac_dE_E0 = dE1pass_arr/E0bin
#    Vfrac_dE_E0 = np.sqrt(Efrac_dE_E0)
#    fig.add_subplot(212).plot(rp_arr, Efrac_dE_E0, marker='o', linestyle='-', markersize=5.0, linewidth=1.5, alpha=0.75, label= str("{0:.1f}".format(a0bin_AU_arr[ac])))
##axis settings:
#plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)
#plt.xlabel(r'peri-center $r_{p} [R_{sun}]$')
#plt.ylabel(r'$dE/E0$')
#plt.yscale('log')

##plt.savefig(add_filename + '_dE_rp_plots.eps')
#plt.show()
##exit()
#-------------------------------------------




#-------------------------------------------
#plot (dE,rp): OTHER VERSION (maybe for paper)
#-------------------------------------------
#DATA:
input_dE_firstpassage_file  = open(add_filename+'_'+'ae_evol_dE_firstpassage_output.txt',  "r")
ae_evol_dE_firstpassage_arr = np.loadtxt(input_dE_firstpassage_file, dtype=float)  #[ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides/2.]
input_dE_firstpassage_file.close()
#dE_firstpassage = dEorb_apo_arr[1] - dEorb_apo_arr[0]
#dE_firstpassage_save_arr = np.array([ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides/2., max_fractid_1pass, max_max_a_1pass])

#NOTE: FOR THIS PLOT CHOOSE:
# CHOOSE: obj 1: MS(sun), obj 2: CO(M_sun, R=0).
# - just one value for a (does not really play any role)
# - several values for rp
# - same mass objects is best (standard 1R,1M sun).
# - one of the objects to be a compact obj to probe below 2R!!! ASSUME THIS BELOW!!

#define from data file: ae_evol_dE_firstpassage_arr [nr (a,rp)comb, [ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides/2.]]
rp_arr                  = ae_evol_dE_firstpassage_arr[:,1]
dE1pass_arr             = abs(ae_evol_dE_firstpassage_arr[:,2])
max_fractid_1pass_arr   = ae_evol_dE_firstpassage_arr[:,4]
max_max_a_1pass_arr     = ae_evol_dE_firstpassage_arr[:,5]

#make figure:
fig = plt.figure(figsize=(5, 6))

#PLOT 1:
#compare dE with Eorb:
a0bin_AU_arr    = np.array([0.001, 0.01, 0.1, 1.0])
a0bin_arr       = AU_U*a0bin_AU_arr
E0bin_arr       = b1_mass*b2_mass/(2.*a0bin_arr)
nr_a0           = len(a0bin_AU_arr)
color_vals      = np.linspace(50, 200, num=nr_a0, endpoint=True).astype(int)

print color_vals
for ac in range(0,nr_a0):
    E0bin       = E0bin_arr[ac]
    Efrac_dE_E0 = dE1pass_arr/E0bin
    fig.add_subplot(211).plot(rp_arr/b1_rad, np.log10(Efrac_dE_E0), marker='o', linestyle='-', markersize=2.0, linewidth=1.5, alpha=0.75, color=plt.cm.cubehelix(color_vals[ac]), label= 'a [AU] = ' + str("{0:.3f}".format(a0bin_AU_arr[ac])))
#region fill:
fig.add_subplot(211).fill_between([2,4], np.log10(0.1),np.log10(1.0), facecolor='yellow', alpha=0.25)
#axis settings:
plt.legend(loc='upper right', numpoints = 1, fontsize = 8.0, frameon = False, ncol=2)
plt.xlabel(r'pericenter $r_{p}\ [R_{WD}]$')
plt.ylabel(r'log $\Delta{E}_{tid}(r_{p})/E_{orb}(a)$')
fig.add_subplot(211).set_xlim(2.0,4.0)
#plt.yscale('log')

#PLOT 2:
#dE_firstpassage_save_arr = np.array([ini_SMA_bin, ini_rp_bin, dE_firstpassage, dEorb_PT_tides/2., max_fractid_1pass, max_max_a_1pass])
E_TO = b1_mass*b1_mass/b1_rad
fig.add_subplot(212).plot(rp_arr/b1_rad, max_fractid_1pass_arr, marker='o', linestyle='-', markersize=2.0, linewidth=1.5, alpha=0.75, color='black', label= r'$F_{tides}/F_{WD}$')
fig.add_subplot(212).plot(rp_arr/b1_rad, max_max_a_1pass_arr, marker='o', linestyle='--', markersize=2.0, linewidth=1.5, alpha=0.75, color='black', label= r'$max(a_1,a_2,a_3)$')
fig.add_subplot(212).plot(rp_arr/b1_rad, dE1pass_arr/E_TO, marker='o', linestyle=':', markersize=2.0, linewidth=1.5, alpha=0.75, color='black', label= r'$\Delta{E}_{tid}/E_{WD}$')
#axis settings:
plt.legend(loc='upper right', numpoints = 1, fontsize = 8.0, frameon = False)
plt.xlabel(r'pericenter $r_{p}\ [R_{WD}]$')
fig.add_subplot(212).set_xlim(2.0,4.0)
fig.add_subplot(212).set_ylim(-1,3.0)

#plt.savefig(add_filename + '_dE_rp_plots.eps')
plt.show()
exit()
#-------------------------------------------




#-------------------------------------------
#-------------------------------------------
#Figure 4:
#-------------------------------------------
#-------------------------------------------
if (usrinput_figure == 4):
    
    #a is const in this ex:
    ini_SMA_bin = ae_evol_APOCENTER_all_INFO[0,0]
    ini_Torb    = 2.*np.pi*np.sqrt((ini_SMA_bin**3.)/(b1_mass+b2_mass))        

    #figure settings:
    pl.figure(figsize=(5.0, 7.5))
    ax1 = pl.subplot(211) 
    ax2 = pl.subplot(212)    
       
    #loop over different r_p (a is const in this ex):
    for j in range(0,nr_rp):              
        #---------------------
        #read in info/data for sim a(i), rp(j):
        #---------------------
        #initial values for rp sim:
        ini_rp_bin  = ae_evol_APOCENTER_all_INFO[j,1]
        #sim data indices:
        index_min   = ae_evol_APOCENTER_all_INFO[j,3]
        index_max   = ae_evol_APOCENTER_all_INFO[j,4]
        #read in APOCENTER data:
        apo_data_N  = ae_evol_APOCENTER_all_DATA[index_min:index_max+1,0]
        apo_data_t  = ae_evol_APOCENTER_all_DATA[index_min:index_max+1,1]
        apo_data_a  = ae_evol_APOCENTER_all_DATA[index_min:index_max+1,2]
        #---------------------
        #Calc:
        #---------------------
        nr_apo_center   = len(apo_data_N)
        #---------------------
        #plot:
        #---------------------
        for k in range(0,nr_apo_center-1):
            
            #PLOT as a function of TIME:
            Torb_a  = 2.*np.pi*np.sqrt((apo_data_a[k]**3.)/(b1_mass+b2_mass))
            #rectange: width, hight
            recw = Torb_a/ini_Torb + 0.1    #make a little wider to remove small gaps in plotting.
            rech = 1.1*delta_rp             #make a little wider to remove small gaps in plotting.
            #rectange: left corner pos x,y
            recx = apo_data_t[k]/ini_Torb   - 0.5*recw
            recy = scan_rp_arr[j]           - 0.5*rech
            #color settings:
            cnorm   = 1.0
            cval    = int(((1-(apo_data_a[k]/ini_SMA_bin))/cnorm)*250)  
            #plot rectangle:
            rect    = mpl.patches.Rectangle((recx,recy), recw,rech, color=plt.cm.BrBG(cval), linewidth=0)
            ax1.add_patch(rect)

            #PLOT as a function of NR APO:
            recw = 1.2          #make a little wider to remove small gaps in plotting.
            rech = 1.1*delta_rp #make a little wider to remove small gaps in plotting.
            #rectange: left corner pos x,y
            recx = apo_data_N[k]    - 0.5*recw
            recy = scan_rp_arr[j]   - 0.5*rech
            #color settings:
            cnorm   = 1.0
            cval    = int(((1-(apo_data_a[k]/ini_SMA_bin))/cnorm)*250)  
            #plot rectangle:
            rect    = mpl.patches.Rectangle((recx,recy), recw,rech, color=plt.cm.BrBG(cval), linewidth=0)
            ax2.add_patch(rect)
            
            #print info:
            print j, nr_rp       
        #---------------------
        
    
    
    #---------------------
    #oplot analytical guide lines:
    #---------------------
    
    #PLOT time at some input a:
    #calc PT dE (orbit) from initial rp (ASSUMING n=1.5 polytrope (other n, change coefficients in fitting formula)!!!!. 'io' stand for Identical Objects.)
    io_mass     = b1_mass
    io_radius   = b1_rad
    rp_arr      = np.arange(2,10,0.001)
    eta_imp     = (1./np.sqrt(2.))*((rp_arr/io_radius)**(3./2.))
    log10eta    = np.log10(eta_imp)
    log10T2     = (-0.397) + (1.678)*(log10eta**(1.)) + (1.277)*(log10eta**(2.)) + (-12.42)*(log10eta**(3.)) + (9.446)*(log10eta**(4.)) + (-5.550)*(log10eta**(5))
    T2          = 10**(log10T2)
    dEorb_PT_tides_Arr  = (1./2.)*(io_mass**2/io_radius)*(T2/(eta_imp**4))    #the energy the orbit looses per peri-center passage (=2*dE_star), assuming equal obj stars.
    #calc gamma scale fac:
    t_evolv_gamma   = (dEorb_PT_tides_Arr/(io_mass**(3./2.)))*(1./(np.pi*np.sqrt(2.)))
    #calc t time arr:
    input_frac_a0_0 = 0.90
    input_frac_a0_1 = 0.75
    input_frac_a0_2 = 0.50
    input_frac_a0_3 = 2.*rp_arr/ini_SMA_bin
    time_input_frac_a0_0_arr    =  (np.sqrt(ini_SMA_bin) - np.sqrt(input_frac_a0_0*ini_SMA_bin))/t_evolv_gamma
    time_input_frac_a0_1_arr    =  (np.sqrt(ini_SMA_bin) - np.sqrt(input_frac_a0_1*ini_SMA_bin))/t_evolv_gamma
    time_input_frac_a0_2_arr    =  (np.sqrt(ini_SMA_bin) - np.sqrt(input_frac_a0_2*ini_SMA_bin))/t_evolv_gamma
    time_input_frac_a0_3_arr    =  (np.sqrt(ini_SMA_bin) - np.sqrt(input_frac_a0_3*ini_SMA_bin))/t_evolv_gamma
    #plot:
    ax1.plot(time_input_frac_a0_0_arr/ini_Torb, rp_arr, color='black', linewidth=1.0, linestyle='--')
    ax1.plot(time_input_frac_a0_1_arr/ini_Torb, rp_arr, color='black', linewidth=1.0, linestyle='--')
    ax1.plot(time_input_frac_a0_2_arr/ini_Torb, rp_arr, color='black', linewidth=1.0, linestyle='--')
    ax1.plot(time_input_frac_a0_3_arr/ini_Torb, rp_arr, color='black', linewidth=1.5)
    #PLOT isolation time:
    T_iso       = ((2./(1.-(1./(ini_SMA_bin/binsin_ex_a0))))**(3./2.))*2.*np.pi*np.sqrt((binsin_ex_a0**3.)/(b1_mass+b1_mass+b1_mass))
    ax1.plot([T_iso/ini_Torb, T_iso/ini_Torb], [-10,10], color='black', linewidth=1.5, linestyle='-.')
    #---------------------
    
    #---------------------
    #plot settings: 
    #---------------------
    max_apo_data_t  = max(ae_evol_APOCENTER_all_DATA[:,1])
    max_N_apo       = max(ae_evol_APOCENTER_all_DATA[:,0])
    
    ax1.patch.set_facecolor('salmon')    
    ax1.set_xlim(0,(max_apo_data_t/ini_Torb)-5)
    ax1.set_ylim(2.0,4.0)
    ax1.set_title(r'Orbital evolution of semi-major axis $a/a_{0}$')
    ax1.set_xlabel(r'time $t/T_{orbit}$')
    ax1.set_ylabel(r'initial peri-center distance $r_{p}$')
   
    ax2.patch.set_facecolor('salmon')    
    ax2.set_xlim(0,max_N_apo)
    ax2.set_ylim(2.0,4.0)
    #ax2.set_title(r'Orbital evolution of semi-major axis $a/a_{0}$')
    ax2.set_xlabel(r'nr apo-center passages')
    ax2.set_ylabel(r'initial peri-center distance $r_{p}$')
    #---------------------

    #---------------------
    #show and save fig:
    #---------------------
    save_yes_no = raw_input('if save as pdf and quit press 1: ')
    save_yes_no = int(save_yes_no)
    if (save_yes_no == 1):
        figname = raw_input('figure name: ')
        plt.savefig(figname+'.eps')
    plt.show()
    #---------------------
    
#-------------------------------------------



#-------------------------------------------







#-------------------------------------------
#-------------------------------------------
#Figure 1:
#-------------------------------------------
#-------------------------------------------
if (usrinput_figure == 1):

    #NOTES:
    # We use mass, radii etc. below from the input in the beginning. THAT IS NOT GOOD
    # AS THESE NUMBERS MIGHT CHANGE. Put it into and info file.
    # Should we plot the trajectories also?
    # ADD: colorbar, legends til line examples, maybe T/T_star extra y axis, ...
    # Experiment with res and find good line exmaples.
    # For high res we still see small gaps between rectangles.
    # Should we also make a similar study where we vary r_peri and keep a fixed?

    #---------------------------------------
    #Define plot windows:
    #---------------------------------------
    ax_1 = pl.subplot(221)
    ax_2 = pl.subplot(222)
    ax_3 = pl.subplot(223)
    ax_4 = pl.subplot(224)
    #---------------------------------------
    pc = 0
    for i in range(0,nr_SMA):
        for j in range(0,nr_rp):              
            #---------------------
            #read in info/data for sim a(i), rp(j):
            #---------------------
            #initial values for a,rp sim:
            bin_sma         = ae_evol_sim_all_INFO[pc,0]
            bin_rp          = ae_evol_sim_all_INFO[pc,1]
            bin_ecc         = ae_evol_sim_all_INFO[pc,2]
            #sim data indices:
            sim_index_min   = ae_evol_sim_all_INFO[pc,3]
            sim_index_max   = ae_evol_sim_all_INFO[pc,4]
            #read in sim data:
            sim_data_t  = ae_evol_sim_all_DATA[sim_index_min:sim_index_max+1,0]
            sim_data_a  = ae_evol_sim_all_DATA[sim_index_min:sim_index_max+1,1]
            sim_data_e  = ae_evol_sim_all_DATA[sim_index_min:sim_index_max+1,2]
            #---------------------
            #Calc:
            #---------------------
            #Initial orbital time:
            Torb_bin_all_SMA    = 2.*np.pi*np.sqrt((scan_SMA_arr**3.)/(b1_mass+b2_mass))
            Torb_bin            = Torb_bin_all_SMA[i]
            #Energy:
            sim_data_dE = -(b1_mass*b2_mass)*(1./(2.*sim_data_a[:]) - 1./(2.*sim_data_a[0]))    # dE = E(t) - E(t=0)
            E_star      = -min([(b1_mass*b1_mass/b1_rad), (b2_mass*b2_mass/b2_rad)])    #This is the E for the star with the smallest E.
            T_star      = max([np.sqrt((b1_rad**3.)/b1_mass), np.sqrt((b2_rad**3.)/b2_mass)])
            #define:
            sim_data_dE_unit_E_star = sim_data_dE/E_star
            sim_data_a_unit_a0      = sim_data_a/sim_data_a[0]
            nr_tot_sim_steps    = int(len(sim_data_t))
            print 'nr_tot_sim_steps:    ', nr_tot_sim_steps
            #---------------------
        
            #---------------------
            #Thin dataset for plotting:
            #---------------------
            nr_plot_res = 100       #input
            plotevery           = int(1.*nr_tot_sim_steps/(1.*nr_plot_res))
            t_index_ori         = np.arange(0,nr_tot_sim_steps)
            t_index_thin_out    = t_index_ori[0::plotevery]
            nr_t_index_thin_out = len(t_index_thin_out)
            #---------------------
        
            plot_index_rp       = 0
            plot_indices_sma    = [0,4,6,11,29,49]
            color_vals_sma      = [25,50,75,100,125,150]
                
            #---------------------
            #PLOT CHANGE IN ENERGY:
            #---------------------
            if (j == plot_index_rp): #fixed rp
                #Plot (ax_1) contour plot:
                for k in range(0,nr_t_index_thin_out-1):
                    t   = t_index_thin_out[k]   #t(k)
                    t1  = t_index_thin_out[k+1] #t(k+1)
                    #rectange: width, hight
                    recw = 1.2*(sim_data_t[t1]-sim_data_t[t])/Torb_bin  #make a little wider to remove small gaps in plotting.
                    rech = delta_SMA
                    #rectange: left corner pos x,y
                    recx = sim_data_t[t]/Torb_bin
                    recy = scan_SMA_arr[i]-0.5*rech
                    #color settings:
                    cnorm   = 0.015 
                    cval    = int((sim_data_dE_unit_E_star[t]/cnorm)*250)  
                    #plot rectangle:
                    rect    = mpl.patches.Rectangle((recx,recy), recw,rech, color=plt.cm.YlGnBu(cval), linewidth=0)
                    ax_1.add_patch(rect)
                    #plot guide lines:
                    ax_1.plot([0,Tsim_unit_Torb],[recy+0*delta_SMA,recy+0*delta_SMA], color='lightgrey', lw=0.2)
                    ax_1.plot([0,Tsim_unit_Torb],[recy+1*delta_SMA,recy+1*delta_SMA], color='lightgrey', lw=0.2)                  
                #add color bar:
                #how do we do that...?    
                #Plot (ax_2) a few examples (lines):
                if (i == plot_indices_sma[0]):     ax_2.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_dE_unit_E_star[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[0]), lw=1.0)
                if (i == plot_indices_sma[1]):     ax_2.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_dE_unit_E_star[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[1]), lw=1.0)
                if (i == plot_indices_sma[2]):     ax_2.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_dE_unit_E_star[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[2]), lw=1.0)
                if (i == plot_indices_sma[3]):     ax_2.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_dE_unit_E_star[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[3]), lw=1.0)
                if (i == plot_indices_sma[4]):     ax_2.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_dE_unit_E_star[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[4]), lw=1.0)
                if (i == plot_indices_sma[5]):     ax_2.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_dE_unit_E_star[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[5]), lw=1.0)
            #--------------------- 
        
            #---------------------
            #PLOT CHANGE IN SMA:
            #---------------------
            if (j == plot_index_rp): #fixed rp
                #Plot (ax_3) contour plot:
                for k in range(0,nr_t_index_thin_out-1):
                    t   = t_index_thin_out[k]   #t(k)
                    t1  = t_index_thin_out[k+1] #t(k+1)
                    #rectange: width, hight
                    recw = 1.2*(sim_data_t[t1]-sim_data_t[t])/Torb_bin  #make a little wider to remove small gaps in plotting.
                    rech = delta_SMA
                    #rectange: left corner pos x,y
                    recx = sim_data_t[t]/Torb_bin
                    recy = scan_SMA_arr[i]-0.5*rech
                    #color settings:
                    cnorm   = 0.5 
                    cval    = int(((1-sim_data_a_unit_a0[t])/cnorm)*250)  
                    #plot rectangle:
                    rect    = mpl.patches.Rectangle((recx,recy), recw,rech, color=plt.cm.YlGnBu(cval), linewidth=0)
                    ax_3.add_patch(rect)
                    #plot guide lines:
                    ax_3.plot([0,Tsim_unit_Torb],[recy+0*delta_SMA,recy+0*delta_SMA], color='lightgrey', lw=0.2)
                    ax_3.plot([0,Tsim_unit_Torb],[recy+1*delta_SMA,recy+1*delta_SMA], color='lightgrey', lw=0.2)                  
                #add color bar:
                #how do we do that...?    
                #Plot (ax_4) a few examples (lines):
                if (i == plot_indices_sma[0]):     ax_4.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_a_unit_a0[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[0]), lw=1.0)
                if (i == plot_indices_sma[1]):     ax_4.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_a_unit_a0[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[1]), lw=1.0)
                if (i == plot_indices_sma[2]):     ax_4.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_a_unit_a0[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[2]), lw=1.0)
                if (i == plot_indices_sma[3]):     ax_4.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_a_unit_a0[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[3]), lw=1.0)
                if (i == plot_indices_sma[4]):     ax_4.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_a_unit_a0[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[4]), lw=1.0)
                if (i == plot_indices_sma[5]):     ax_4.plot(sim_data_t[0::plotevery]/Torb_bin, sim_data_a_unit_a0[0::plotevery], alpha=0.75, color=plt.cm.Dark2(color_vals_sma[5]), lw=1.0)
            #--------------------- 
       
            #---------------------
            #update counter:
            #---------------------
            pc = pc+1
            #---------------------        
    #---------------------------------------
    #axis settings, save and show figure:
    #---------------------------------------
    ax_1.set_xlim(0,Tsim_unit_Torb)
    ax_1.set_ylim(min(scan_SMA_arr[:])-delta_SMA,max(scan_SMA_arr[:])+delta_SMA)
    ax_1.set_title(r'ENERGY ($r_{p} = 3R_{star}$)')
    ax_1.set_xlabel(r'time $t/T_{orbit}$')
    ax_1.set_ylabel(r'initial SMA $a$')

    ax_2.set_xlim(0,Tsim_unit_Torb)
    ax_2.set_ylim(0,0.015)
    ax_2.set_title(r'ENERGY ($r_{p} = 3R_{star}$)')
    ax_2.set_xlabel(r'time $t/T_{orbit}$')
    ax_2.set_ylabel(r'$(E_{orbit}(t)-E_{orbit}(t=0))/E_{star}$')

    ax_3.set_xlim(0,Tsim_unit_Torb)
    ax_3.set_ylim(min(scan_SMA_arr[:])-delta_SMA,max(scan_SMA_arr[:])+delta_SMA)
    ax_3.set_title(r'Semi-Major Axis ($r_{p} = 3R_{star}$)')
    ax_3.set_xlabel(r'time $t/T_{orbit}$')
    ax_3.set_ylabel(r'initial SMA $a$')

    ax_4.set_xlim(0,Tsim_unit_Torb)
    ax_4.set_ylim(0,1.1)
    ax_4.set_title(r'Semi-Major Axis ($r_{p} = 3R_{star}$)')
    ax_4.set_xlabel(r'time $t/T_{orbit}$')
    ax_4.set_ylabel(r'SMA $a(t)/a(t=0)$')

    plt.subplots_adjust(wspace=0.3)
    plt.subplots_adjust(hspace=0.3)

    plt.savefig('test_binary_evol.eps', bbox_inches='tight')

    plt.show()
    exit()    
    #---------------------------------------
    #---------------------------------------
#-------------------------------------------










#-------------------------------------------
#-------------------------------------------
#Figure 2:
#-------------------------------------------
#-------------------------------------------
if (usrinput_figure == 2):
    
    #NOTES:
    #We use a0 below. But we still use the val from input at the top. Not good.
    # - MUST put into some kind of input file!!!!
    # for this fig we DONT READ ALL DATA. Change that option in some way.
    # save initial, final e in the sim loop.
    # LATER: we need to track a,e IMS. Put that into the NbodyTides code.
    # on this fig overplot a few r_p lines.
    # Check all codes again.
    # Make figure as a split figure: left zoom out, right zoom in.

    fig = plt.figure(figsize=(3.5, 3.0)) 

    pc = 0
    for i in range(0,nr_SMA):
        for j in range(0,nr_rp):              
            #---------------------
            #read in info/data for sim a(i), rp(j):
            #---------------------
            #initial values for a,rp sim:
            bin_sma         = ae_evol_sim_all_INFO[pc,0]
            bin_rp          = ae_evol_sim_all_INFO[pc,1]
            bin_ecc         = ae_evol_sim_all_INFO[pc,2]
            #sim data indices:
            #sim_index_min   = ae_evol_sim_all_INFO[pc,3]
            #sim_index_max   = ae_evol_sim_all_INFO[pc,4]
            #outcome type:
            outcome_flag    = int(ae_evol_sim_all_INFO[pc,5])
            #initial initial/final values:
            initial_ecc     = bin_ecc
            initial_sma     = bin_sma
            final_sma       = ae_evol_sim_all_INFO[pc,7]
            initial_rp      = initial_sma*(1.-initial_ecc)
            #read in sim data:
            #sim_data_t  = ae_evol_sim_all_DATA[sim_index_min:sim_index_max+1,0]
            #sim_data_a  = ae_evol_sim_all_DATA[sim_index_min:sim_index_max+1,1]
            #sim_data_e  = ae_evol_sim_all_DATA[sim_index_min:sim_index_max+1,2]
            #---------------------
            #Define:
            #---------------------
            #nr_tot_sim_steps    = int(len(sim_data_t))
            #final_sma           = sim_data_a[nr_tot_sim_steps-1]
            #initial_sma         = sim_data_a[0]
            #initial_ecc         = sim_data_e[0]
            #---------------------
            #Calc:
            #---------------------
            final_aeff          = 1./((1./binsin_ex_a0) + (1./final_sma) - (1./initial_sma)) 
            final_aeff_unit_a0  = final_aeff/binsin_ex_a0
            #---------------------              
            print bin_sma, bin_ecc
            print 'final_aeff_unit_a0:  ', final_aeff_unit_a0       
            #---------------------
            #PLOT: 1
            #---------------------
            #If IMS binary survives:
            if (outcome_flag == 10):
                cnorm = 0.5
                colorval = int(((1-final_aeff_unit_a0)/cnorm)*250)
                fig.add_subplot(111).scatter(initial_sma/binsin_ex_a0, initial_ecc, s=3.00, alpha=1.0, edgecolor='black', lw=0.2, c=plt.cm.YlGnBu(colorval), marker='o')
            #If IMS collide in less than T_iso:
            if (outcome_flag == 1 or outcome_flag == 2):
                fig.add_subplot(111).scatter(initial_sma/binsin_ex_a0, initial_ecc, s=3.00, alpha=1.0, edgecolor='black', lw=0.2, c=plt.cm.Reds(150), marker='o')
            #---------------------       
            ##---------------------
            ##PLOT: 2
            ##---------------------
            ##If IMS binary survives:
            #if (outcome_flag == 10):
            #    cnorm = 0.5
            #    colorval = int(((1-final_aeff_unit_a0)/cnorm)*250)
            #    fig.add_subplot(122).scatter(initial_sma/binsin_ex_a0, initial_rp, s=10.00, alpha=1.0, edgecolor='black', lw=0.2, c=plt.cm.YlGnBu(colorval), marker='o')
            ##If IMS collide in less than T_iso:
            #if (outcome_flag == 1 or outcome_flag == 2):
            #    fig.add_subplot(122).scatter(initial_sma/binsin_ex_a0, initial_rp, s=10.00, alpha=1.0, edgecolor='black', lw=0.2, c=plt.cm.Reds(150), marker='o')
            ##---------------------
            #---------------------
            #update counter:
            #---------------------
            pc = pc+1
            #---------------------        
    #---------------------------------------
    #oplot shaded areas:
    #---------------------------------------
    coll_a_arr = np.arange(0,2*max(scan_SMA_arr),0.001)
    coll_e_arr = 1. - r_coll/coll_a_arr
    fig.add_subplot(111).fill_between(coll_a_arr/binsin_ex_a0, coll_e_arr, 1, alpha=0.75, color=plt.cm.Reds(150))
    #pa_arr      = np.arange(1,2,0.001)
    #pe_arr      = 1. - (r_coll/binsin_ex_a0)/pa_arr
    #fig.add_subplot(111).fill_between(pa_arr, 0, pe_arr, alpha=0.5, color='gray')
    #fig.add_subplot(111).plot([1.0,1.0], [0,1], alpha=0.25, linewidth=0.5, linestyle='-', color=plt.cm.Greys(250))
    #---------------------------------------
    
    #---------------------------------------
    #oplot guide lines:
    #---------------------------------------
    #---------------------------------------
    #object data:
    #---------------------------------------
    inp_mass    = b1_const_arr[0]
    inp_radius  = b1_const_arr[1]
    inp_a0      = binsin_ex_a0
    #---------------------------------------
    #---------------------------------------
    #collision boundary:
    #---------------------------------------
    r_coll = 2.*inp_radius
    a_pa = np.arange((r_coll/inp_a0),2,0.01)
    e_pa = 1.0 - (r_coll/inp_a0)/a_pa
    fig.add_subplot(111).plot(a_pa,e_pa, color='black', alpha=1.0, linewidth=1.0, linestyle='-')
    #fig.add_subplot(111).fill_between(a_pa, e_pa, 1.0, alpha=0.25, color=plt.cm.Reds(150))
    #---------------------------------------
    #---------------------------------------
    #orbital instability boundary:
    #---------------------------------------
    r_orbital_inst_coll = 2.4
    a_pa = np.arange((r_orbital_inst_coll/inp_a0),2,0.01)
    e_pa = 1.0 - (r_orbital_inst_coll/inp_a0)/a_pa
    fig.add_subplot(111).plot(a_pa,e_pa, color='black', alpha=1.0, linewidth=1.0, linestyle='-')
    #---------------------------------------
    #---------------------------------------
    #a/a_0 = 1 vertical boundary:
    #---------------------------------------
    fig.add_subplot(111).plot([1,1],[0,1], color='black', alpha=1.0, linewidth=1.0, linestyle='-')
    #---------------------------------------
    #---------------------------------------
    #Analytical limits:
    #---------------------------------------
    #make input rp array:
    rp_arr          = np.arange(2.0,20.0,0.01)
    nr_rp           = len(rp_arr)
    #make arrays for different solutions:
    #afoai:
    sol_afoai_aprime_arr    = np.zeros((nr_rp), dtype=np.float64)
    sol_afoai_ecc_arr       = np.zeros((nr_rp), dtype=np.float64)
    sol_afoai_aprime_arr[:]     = -1    #set to -1 for later checking
    #afeq2rp:
    sol_afeq2rp_aprime_arr  = np.zeros((nr_rp), dtype=np.float64)
    sol_afeq2rp_ecc_arr     = np.zeros((nr_rp), dtype=np.float64)
    sol_afeq2rp_aprime_arr[:]   = -1    #set to -1 for later checking
    #aeffoa0:
    sol_aeffoa0_aprime_arr  = np.zeros((nr_rp), dtype=np.float64)
    sol_aeffoa0_ecc_arr     = np.zeros((nr_rp), dtype=np.float64)
    sol_aeffoa0_aprime_arr[:]   = -1    #set to -1 for later checking
    #input frac limits:
    f_afoai     = 0.5
    f_aeffoa0   = 0.95
    #---------------------------------------
    for rpc in range(0,nr_rp):
    #---------------------------------------    
        #---------------------
        rp      = rp_arr[rpc]   #peri-center in code units.
        #---------------------
        #calc PT dE (orbit) from rp (ASSUMING n=1.5 polytrope (other n, change coefficients in fitting formula)!!!! see paper: )
        #---------------------
        eta_imp         = (1./np.sqrt(2.))*((rp/inp_radius)**(3./2.))
        log10eta        = np.log10(eta_imp)
        T2              = 10**((-0.397) + (1.678)*(log10eta**(1.)) + (1.277)*(log10eta**(2.)) + (-12.42)*(log10eta**(3.)) + (9.446)*(log10eta**(4.)) + (-5.550)*(log10eta**(5)))
        dEorb_PT_tides  = (1./2.)*((inp_mass**2.)/inp_radius)*(T2/(eta_imp**4))    #the energy the orbit looses per peri-center passage (=2*dE_star), assuming equal obj stars.
        #---------------------
        #define:
        #---------------------
        dE      = dEorb_PT_tides
        E0      = (inp_mass**2.)/(2.*inp_a0)
        F_rp    = (2./np.sqrt(3.))*(dE/E0)
        #---------------------
        #solve for: a_f/a_i
        #---------------------
        kfac2   = ((1./F_rp)*(1.-np.sqrt(f_afoai)))**2.
        coeff   = [(kfac2), -(1+3.*kfac2), (3.*kfac2), (-kfac2)]
        roots   = np.roots(coeff)
        aprime  = roots[0]
        ecc     = 1. - (rp/inp_a0)/aprime
        sol_afoai_aprime_arr[rpc]       = aprime
        sol_afoai_ecc_arr[rpc]          = ecc
        #---------------------
        #solve for: a_f = 2*rp
        #---------------------
        func = lambda x : x*(1. - F_rp*(x/((x-1.)**(3./2.))))**(2.) - 2.*(rp/inp_a0)
        x_iniguess  = 1.001
        x_solution  = fsolve(func, x_iniguess)
        aprime      = x_solution
        ecc         = 1. - (rp/inp_a0)/aprime
        if (func(aprime) < 1e-10 and ecc > 0.75):
            sol_afeq2rp_aprime_arr[rpc]     = aprime
            sol_afeq2rp_ecc_arr[rpc]        = ecc
        #---------------------
        #solve for: a_eff/a_0
        #---------------------
        func        = lambda x : 1. + 1./(x*(1. - F_rp*(x/((x-1.)**(3./2.))))**(2.)) - 1./x - 1./f_aeffoa0
        func_check  = lambda x : 1. - F_rp*(x/((x-1.)**(3./2.)))    #for checking solution: this MUST be > 0 by contruction.
        x_iniguess  = 1.001
        x_solution  = fsolve(func, x_iniguess)
        aprime      = x_solution
        ecc         = 1. - (rp/inp_a0)/aprime
        if (func(aprime) < 1e-10 and ecc > 0.75 and func_check(aprime) > 0.0):
            sol_aeffoa0_aprime_arr[rpc]     = aprime
            sol_aeffoa0_ecc_arr[rpc]        = ecc
    #---------------------------------------        
    #finalize solution arrays:
    #---------------------------------------
    pos = np.where(sol_afeq2rp_aprime_arr[:] != -1)
    finsol_afeq2rp_aprime_arr   = sol_afeq2rp_aprime_arr[pos[0]]
    finsol_afeq2rp_ecc_arr      = sol_afeq2rp_ecc_arr[pos[0]]

    pos = np.where(sol_aeffoa0_aprime_arr[:] != -1)
    finsol_aeffoa0_aprime_arr   = sol_aeffoa0_aprime_arr[pos[0]]
    finsol_aeffoa0_ecc_arr      = sol_aeffoa0_ecc_arr[pos[0]]
    #---------------------------------------
    #---------------------------------------    
    #plot solutions:
    #---------------------------------------
    #fig.add_subplot(111).plot(sol_afoai_aprime_arr,sol_afoai_ecc_arr)
    fig.add_subplot(111).plot(finsol_afeq2rp_aprime_arr,finsol_afeq2rp_ecc_arr, alpha=1.0, color='black', linewidth=1.0, linestyle='--')
    fig.add_subplot(111).plot(finsol_aeffoa0_aprime_arr,finsol_aeffoa0_ecc_arr, alpha=1.0, color='black', linewidth=1.0, linestyle='-.')
    #---------------------------------------
    
    #---------------------------------------
    #axis settings, save and show figure:
    #---------------------------------------
    fig.add_subplot(111).set_xlim(0.5,2.01)
    fig.add_subplot(111).set_ylim(0.8,1.0)
    fig.add_subplot(111).set_title(r'$(a_{eff}/a_{0}$) for IMS binaries with $a_{0} = 0.1$ au')
    fig.add_subplot(111).set_xlabel(r'semi-major axis $a/a_{0}$')
    fig.add_subplot(111).set_ylabel(r'eccentricity $e$')
    
    #fig.add_subplot(122).set_xlim(1.0,2.01)
    #fig.add_subplot(122).set_ylim(1.5,4.5)
    #fig.add_subplot(122).set_title(r'Tidal Energy loss ($a_{eff}/a_{0}$) of IMS binaries ($a_{0} = 0.1$ au)')
    #fig.add_subplot(122).set_xlabel(r'SMA $a_{IMS}/a_{0}$')
    #fig.add_subplot(122).set_ylabel(r'peri-center distance $r_{p}$')

    plt.savefig(add_filename+'_'+'aeffa0IMSbins.eps', bbox_inches='tight')

    plt.show()
    exit()
    #---------------------------------------
    #---------------------------------------
#-------------------------------------------






























































