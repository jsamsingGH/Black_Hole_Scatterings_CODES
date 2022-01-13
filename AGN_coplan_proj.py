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
from scipy.integrate import solve_ivp
from matplotlib import rcParams
from scipy.interpolate import make_interp_spline, BSpline


#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
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
#------------------------------------------------------------------------------------------


N   = 20.
k   = (5.*(c_SI**5.)/(512.*(G_new_SI**3.)))
m   = 20.*M_sun_SI
a   = 1.0*AU_SI
tau = (10**5.)*(sec_year)
dcp = 2./3.
fHz = 10.
ef  = 0.1
tcs = k*(a**4./m**3.)
print tcs*(1-0.999**2)**(7./2.)
#exit()

T0s = np.sqrt(a**3./(G_new_SI*(2.*m)))
P2  = (tau**(1./7.))*(k**(-1./7.))*m**(3./7.)*a**(-4./7.)
P3  = N*((k*(G_new_SI**(1./2.)))**(-1./7.))*m**(5./14.)*a**(-5./14.)
a23 = (N**(-14./3.))*(tau**(2./3.))*(G_new_SI**(1./3.))*(m**(1./3.))
a1  = ((5.*(1.-dcp))/(14.*N))**(-14./5.)*((k*(G_new_SI**(1./2.)))**(-2./5.))*m
at  = (tau*(k**(-1.))*(m**3.))**(1./4.)
Fef = ((ef**(12./19.))/(1.+ef))*(1. + (121./304.)*(ef**2.))**(870./2299.)
rf  = ((2.*G_new_SI*m)/((fHz**2.)*(np.pi**2.)))**(1./3.) 
rc  = rf*(1./(2.*Fef))*(425./304.)**(870./2299.) 
Pef = N*((2.*rc)/(a))**(1./2.)
Rsc = (2.*G_new_SI*m)/(c_SI**2.)
print P2
print P3
print a23/AU_SI
print a1/AU_SI
print at/AU_SI
print Pef
print a/Rsc
print np.sqrt(G_new_SI*(50.*M_sun_SI)/(1.0*AU_SI)) 

print np.sqrt(G_new_SI*(M_sun_SI)/(R_sun_SI))


print (tcs/T0s)**(1./7.)
print (tcs/tau)**(1./7.)
print (tau/T0s)**(1./7.)

print (((((1./0.025)**(-8./3.))/(1000.*0.1))**(3./5.))*(c_SI**3.)/G_new_SI)/M_sun_SI



#exit()






mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

analyze_data_YN = 0
plotetc_data_YN = 1

m1_SI       = 20.*M_sun_SI#40.*M_sun_SI
m2_SI       = 20.*M_sun_SI#20.*M_sun_SI
m3_SI       = 20.*M_sun_SI#20.*M_sun_SI
m123_SI_arr = np.array([m1_SI, m2_SI, m3_SI])

tint_SI     = (10.**5.)*sec_year

EMcase_YN   = 1 #equal mass (EM) case, YES/NO

obj121323_arr       =  np.zeros((3, 2), dtype=np.int)
obj121323_arr[0,:]  = [1,2]
obj121323_arr[1,:]  = [1,3]
obj121323_arr[2,:]  = [2,3]


#--------------------------------------------------------------
if (analyze_data_YN == 1):
#--------------------------------------------------------------
    #----------------------------------------------------------
    #INPUT:
    #----------------------------------------------------------
    dataset_arr = ['RVA_1AU_FA05_1000000_402020']#['RVA2_1AU_FA5_1000000_402020']#['RVA_1AU_FA05_1000000_402020']
    #['RVA_1AU_FA0_50000_202020', 'RVA_1AU_FA003_50000_202020', 'RVA_1AU_FA01_50000_202020', 'RVA_1AU_FA03_50000_202020', 'RVA_1AU_FA1_50000_202020', 'RVA_1AU_FA2_50000_202020', 'RVA_1AU_FA3_50000_202020', 'RVA_1AU_FA4_50000_202020', 'RVA_1AU_FA5_50000_202020', 'RVA_1AU_FA10_50000_202020', 'RVA_1AU_FA15_50000_202020', 'RVA_1AU_FA25_50000_202020', 'RVA_1AU_FA45_50000_202020', 'RVA_1AU_FA75_50000_202020', 'RVA_1AU_FA90_50000_202020', 'RVA_1AU_ISO_50000_202020']
    #['RVA_1AU_FA5_1000000_402020']
    #['T1_1e21e2AU10_OA0_10000_202020', 'T1_1e21e2AU10_ISO_10000_202020']    
    #['RVA_1AU_FA0_50000_202020', 'RVA_1AU_FA01_50000_202020', 'RVA_1AU_FA03_50000_202020', 'RVA_1AU_FA1_50000_202020', 'RVA_1AU_FA2_50000_202020', 'RVA_1AU_FA3_50000_202020', 'RVA_1AU_FA4_50000_202020', 'RVA_1AU_FA5_50000_202020', 'RVA_1AU_FA10_50000_202020', 'RVA_1AU_FA15_50000_202020', 'RVA_1AU_FA25_50000_202020', 'RVA_1AU_FA45_50000_202020', 'RVA_1AU_FA75_50000_202020', 'RVA_1AU_FA90_50000_202020', 'RVA_1AU_ISO_50000_202020']
   
    #['T1_1e21e2AU10_OA0_10000_202020', 'T1_1e21e2AU10_ISO_10000_202020']#['T1_1AUOA5_100000_402020'] # ['T1_OA1_202020'] #['T1_OA5_202020'] #['T1_OA10_202020']#['T5_BH202020_1e21e2_10_10000_ISO', 'T5_BH202020_1e21e2_10_10000_CP'] #['T5_BH202020_1e21e2_10_10000_CP', 'T5_BH202020_1e21e2_10_10000_CP']#['T4_BH202020_1e21e2_10_10000_CP', 'T4_BH202020_1e21e2_10_10000_CP']#['T3_BH202020_1e250_10_10000_CP', 'T3_BH202020_1e250_10_10000_CP']#['T3_BH202020_1e150_10_10000_CP', 'T3_BH202020_1e150_10_10000_CP']#['T3_BH402020_1AU_50000_CP']
    #['T3_BH202020_1e21e2_10_10000_ISO', 'T3_BH202020_1e21e2_10_10000_CP']#['T1_BH202020_1e21e2_CP', 'T1_BH202020_1e21e2_CP'] #['S1ISO_', 'S1CP_']#['S422ISO_', 'S422CP_'] #['S1ISO_', 'S1CP_']
    data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'
    #----------------------------------------------------------
    #Define:
    #----------------------------------------------------------
    nrdataset   = len(dataset_arr)
    #Get info from 1st dataset:
    tf = open(data_folder+dataset_arr[0]+'MC_settings_list_INT.txt', "r")
    MC_settings_list_INT    = np.loadtxt(tf, dtype=int)         #[nr_SMA, nr_vinf, nr_sim, 0, 0]
    nr_SMA                  = MC_settings_list_INT[0]
    nr_vinf                 = 1
    nr_sim                  = MC_settings_list_INT[2]
    print nr_SMA, nr_vinf, nr_sim, 'nr_SMA, nr_vinf, nr_sim'
    #define saving array:
    info_per_DATA_SMA_arr       = np.zeros((nrdataset, nr_SMA, 5), dtype=np.float64)
    info_per_DATA_SMA_OC_arr    = np.zeros((nrdataset, nr_SMA, 3, 5), dtype=np.float64)
    HIST_data_arr               = np.zeros((nrdataset, nr_SMA, 3, 5, 5, 2, 50), dtype=np.float64)   
    savetestarr                 = np.zeros((nrdataset, nr_SMA, nr_sim, 10), dtype=np.float64)
    #----------------------------------------------------------

    #----------------------------------------------------------
    #LOOP over datasets:
    #----------------------------------------------------------
    save_data_filename = raw_input('Save data filename: ')
    for dc in range(0,nrdataset):
    
        #get data_name:
        data_name = dataset_arr[dc]
        print 'ANALYZING dataset: ', dc+1

        #------------------------------------------------------
        #Read in data from files:
        #------------------------------------------------------
        #read data:
        tf = open(data_folder+data_name+'MC_settings_list_INT.txt', "r")
        MC_settings_list_INT        = np.loadtxt(tf, dtype=int)         #[nr_SMA, nr_vinf, nr_sim, 0, 0]
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
        #------------------------------------------------------
        #------------------------------------------------------
        #Define:
        #------------------------------------------------------
        #MC_settings_list_INT[:]     = [nr_SMA, nr_vinf, nr_sim, nr_tot_scatterings, 0]    #where nr_tot_scatterings is simply =nr_SMA*nr_vinf*nr_sim
        #nr_sim      = MC_settings_list_INT[2]
        #------------------------------------------------------
        #reshape input arrays for easy analysis:
        #------------------------------------------------------
        RS_avs_outlist_MC_info_INT              = outlist_MC_info_INT.reshape(nr_SMA, nr_vinf, nr_sim, 10)              #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
        RS_avs_outlist_MC_info_REAL             = outlist_MC_info_REAL.reshape(nr_SMA, nr_vinf, nr_sim, 10)             #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]
        RS_avs_output_Nbody_endstate_INT        = output_Nbody_endstate_INT.reshape(nr_SMA, nr_vinf, nr_sim, 10)        #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
        RS_avs_output_Nbody_endstate_REAL       = output_Nbody_endstate_REAL.reshape(nr_SMA, nr_vinf, nr_sim, 10)       #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
        RS_avs_output_Nbody_xtra_info_REAL      = output_Nbody_xtra_info_REAL.reshape(nr_SMA, nr_vinf, nr_sim, 10)      #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e,...]
        RS_avs_output_Nbody_xtra_2_info_REAL    = output_Nbody_xtra_2_info_REAL.reshape(nr_SMA, nr_vinf, nr_sim, 10)    #[out_pos_CMij_wrt_sink(:), out_vel_CMij_wrt_sink(:), Lvec_ij(:)]
        #------------------------------------------------------
        #Define:
        #------------------------------------------------------
        SMA_arr                                 = RS_avs_outlist_MC_info_REAL[:,0,0,0]
        #------------------------------------------------------
        #Analyze:
        #------------------------------------------------------
        #ENDSTATE CODES: 1 = TDE, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 5 = Inspiral, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...
        #IF ENDSTATE code is < 10 then the code reached a physical OK endstate - if >=10 then the code did not converge towards a physical solution.
        #Return_3body_Info_REAL(1:10)		= [E_kin(ij), E_pot(ij), E_tot(ij), a_bin(ij), e_bin(ij), E_kin(ijk), E_pot(ijk), E_tot(ijk), a_bin(ijk), e_bin(ijk)]
        #RS_avs_output_Nbody_endstate_INT format: [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
                
        for ac in range(0,nr_SMA):

            print dc, ac, 'nrdataset, nr_SMA'
            
            #define:
            id_arr      = RS_avs_output_Nbody_endstate_INT[ac,0,:,0]
            
            io_arr      = RS_avs_output_Nbody_endstate_INT[ac,0,:,1]
            jo_arr      = RS_avs_output_Nbody_endstate_INT[ac,0,:,2]
            ko_arr      = 6 - (io_arr + jo_arr)
            
            mi_SI_arr   = m123_SI_arr[io_arr-1]
            mj_SI_arr   = m123_SI_arr[jo_arr-1]
            mk_SI_arr   = m123_SI_arr[ko_arr-1]
                        
            aij_SI_arr  = RS_avs_output_Nbody_endstate_REAL[ac,0,:,3]*R_sun_SI
            eij_SI_arr  = RS_avs_output_Nbody_endstate_REAL[ac,0,:,4]
            
            Lij_ang_arr = np.zeros(nr_sim, dtype=np.float64)
            z_axis_vec  = np.array([0,0,1])
            for sc in range(0,nr_sim):
                Lij_vec         = RS_avs_output_Nbody_xtra_2_info_REAL[ac,0,sc,6:9]
                Lij_ang         = np.arccos(np.dot(Lij_vec, z_axis_vec)/(np.linalg.norm(Lij_vec)*np.linalg.norm(z_axis_vec)))
                Lij_ang_DG      = (Lij_ang/(2.*np.pi))*360.           
                Lij_ang_arr[sc] = Lij_ang_DG
            
            #calc
            tGWij_SI_arr    = ((768./425.)*((aij_SI_arr**4.)/(4.*(64./5.)*(G_new_SI**3.)*mi_SI_arr*mj_SI_arr*(mi_SI_arr+mj_SI_arr)/(c_SI**5.)))*((1.-eij_SI_arr**2.)**(7./2.)))            
            rpij_SI_arr     = aij_SI_arr*(1.-eij_SI_arr)
            fGWij_SI_arr    = (1./np.pi)*np.sqrt(G_new_SI*(mi_SI_arr + mj_SI_arr)/(rpij_SI_arr**3.))
            
            #define:
            pos_id3         = np.where(id_arr[:] == 3)[0]               #binary-single ejection
            pos_id5         = np.where(id_arr[:] == 5)[0]               #3-body merger      
            pos_tlim        = np.where(tGWij_SI_arr < tint_SI)[0]       #time limit
            pos_id3_tlim    = list(set(pos_id3).intersection(pos_tlim))
            pos_id5_tlim    = list(set(pos_id5).intersection(pos_tlim))

            nr_id3          = len(pos_id3)            
            nr_id5          = len(pos_id5)

            #----------------------------------------
            #get info per dataset, SMA, scattering, ...
            #----------------------------------------
            #save info:
            savetestarr[dc, ac, :, 0]   = id_arr[:]
            
            #find and save obj. endstate combination (12, 13, 23):
            for oc in range(0,3):
                obj_i           = obj121323_arr[oc,0]
                obj_j           = obj121323_arr[oc,1]
                pos_i           = np.where(io_arr[:] == obj_i)[0]
                pos_j           = np.where(jo_arr[:] == obj_j)[0]
                pos_ij          = list(set(pos_i).intersection(pos_j))      
                savetestarr[dc, ac, pos_ij, 1]   = (oc+1)   #comb is here referred to as: `1,2,3'        
          
            #save info:
            savetestarr[dc, ac, :, 2]           = tGWij_SI_arr[:]
            savetestarr[dc, ac, pos_tlim, 3]    = 1
            savetestarr[dc, ac, :, 4]           = fGWij_SI_arr[:]
            savetestarr[dc, ac, :, 9]           = Lij_ang_arr[:]
            
            #---------------------------
            #calc ecc at x fGW:
            #---------------------------
            for sc in range(0,nr_sim):
                                
                #input fGW limit:
                fGW_limit_A = 0.01
                fGW_limit_B = 10.
                
                #define:
                rp0     = rpij_SI_arr[sc]
                e0      = eij_SI_arr[sc]
                fGW0    = fGWij_SI_arr[sc]
                mi      = mi_SI_arr[sc]
                mj      = mj_SI_arr[sc]
                
                #CASE A: calc ecc at fGW:
                fGW_limit   = fGW_limit_A
                ecc_A       = 0.0
                if (fGW0 < fGW_limit):
                    #propagating ini a,e to a,e(fGW_limit)
                    c0      = rp0*(1.+e0)*(e0**(-12./19.))*((1.+(121./304.)*e0**2.)**(-870./2299.))
                    func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mi+mj)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
                    ecc_initial_guess   = 1e-8
                    ecc_A               = fsolve(func, ecc_initial_guess)
                if (fGW0 > fGW_limit):  ecc_A = -1
                
                #CASE B: calc ecc at fGW:
                fGW_limit   = fGW_limit_B
                ecc_B       = 0.0
                if (fGW0 < fGW_limit):
                    #propagating ini a,e to a,e(fGW_limit)
                    c0      = rp0*(1.+e0)*(e0**(-12./19.))*((1.+(121./304.)*e0**2.)**(-870./2299.))
                    func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mi+mj)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
                    ecc_initial_guess   = 1e-8
                    ecc_B               = fsolve(func, ecc_initial_guess)
                    #print tGWij_SI_arr[sc]/tint_SI, id_arr[sc], fGW0, fGW_limit, e0, ecc_B
                if (fGW0 > fGW_limit):  ecc_B = -1
                
                #save info:
                savetestarr[dc, ac, sc, 7]           = ecc_A
                savetestarr[dc, ac, sc, 8]           = ecc_B
                
                
                #NOW TEST AND INCL ALSO ECC AT FGW FOR 2nd MERGER BELOW! GET RID OF HIST SECTION?
            #---------------------------
            
            
            #---------------------------
            #calc 2nd GW merger info:
            #---------------------------
            for id5c in range(0,nr_id5):
                #get pos:
                pij     = pos_id5[id5c]
                #mass m1 (mi + mj), m2 (mk):
                m2_ij   = mi_SI_arr[pij] + mj_SI_arr[pij]   #SI
                m1_k    = mk_SI_arr[pij]                    #SI
                #kick vel: DO THIS LATER!
                vel_vec_kick_SI = 0.0*np.array([1,1,1])     
                #define:
                m_1     = m1_k
                m_2     = m2_ij
                m_12    = m_1+m_2
                mu_12   = m_1*m_2/m_12
                #before kick (no kick NK):
                NK_rel_pos_vec_SI   = RS_avs_output_Nbody_xtra_2_info_REAL[ac,0,pij,0:3]*R_sun_SI
                NK_rel_vel_vec_SI   = RS_avs_output_Nbody_xtra_2_info_REAL[ac,0,pij,3:6]*(1000./kmsec_U)
                #after kick:
                rel_pos_vec_SI = NK_rel_pos_vec_SI
                rel_vel_vec_SI = NK_rel_vel_vec_SI + vel_vec_kick_SI
                #calc orb params of post-merger system:
                rel_L_vec_SI    = mu_12*np.cross(rel_pos_vec_SI, rel_vel_vec_SI)
                rel_r_SI        = np.sqrt(sum(rel_pos_vec_SI**2.))
                rel_v_SI        = np.sqrt(sum(rel_vel_vec_SI**2.))
                rel_L_SI        = np.sqrt(sum(rel_L_vec_SI**2.))
                rel_E_SI        = (1./2.)*mu_12*(rel_v_SI**2.) - G_new_SI*m_12*mu_12/rel_r_SI      
                a_SI            = - G_new_SI*m_12*mu_12/(2.*rel_E_SI)
                e_SI            = np.sqrt(1. + (2.*rel_E_SI*(rel_L_SI**2.))/((G_new_SI**2.)*(m_12**2.)*(mu_12**3.)))
                #calc fGW:
                rp_SI           = a_SI*(1.-e_SI)
                fGW12_SI        = (1./np.pi)*np.sqrt(G_new_SI*m_12/(rp_SI**3.))
                #flag:
                if (a_SI > 0.0):  flagval = 1
                if (a_SI < 0.0):  flagval = -1
                savetestarr[dc, ac, pij, 5] = flagval
                savetestarr[dc, ac, pij, 6] = fGW12_SI
            #---------------------------
            #----------------------------------------    
            #----------------------------------------
            #INFO:
            #savetestarr[dc, ac, :, 0]   = id_arr[:]
            #savetestarr[dc, ac, :, 1]   = 1,2,3#comb is here referred to as: `1,2,3'        
            #savetestarr[dc, ac, :, 2]   = tGWij_SI_arr[:]
            #savetestarr[dc, ac, :, 3]   = 1/0 #=1 IF tGW<tlim
            #savetestarr[dc, ac, :, 4]   = fGWij_SI_arr #1st GW merger
            #savetestarr[dc, ac, :, 5]   = 1/0/-1       #=1 if 2nd binary is OK and bound. only use = 1 cases here.
            #savetestarr[dc, ac, :, 6]   = fGW12_SI     #2nd GW merger   

            #pos_id3         = np.where(id_arr[:] == 3)[0]               #binary-single ejection
            #pos_id5         = np.where(id_arr[:] == 5)[0]               #3-body merger      
            #pos_tlim        = np.where(tGWij_SI_arr < tint_SI)[0]       #time limit
            #pos_id3_tlim    = list(set(pos_id3).intersection(pos_tlim))
            #pos_id5_tlim    = list(set(pos_id5).intersection(pos_tlim))

            #define:
            pos_2ndOK       = np.where(savetestarr[dc, ac, :, 5] == 1)[0]
            fGW12_SI_arr    = savetestarr[dc, ac, :, 6]
            
            #general info:
            pos_NAN         = np.where(id_arr[:] >= 10)[0]              #sim did not finish
            nr_NAN          = len(pos_NAN)
            nr_simOK        = nr_sim - nr_NAN
            info_per_DATA_SMA_arr[dc, ac, 0]   = 1.*nr_NAN              #sim did not finish
            info_per_DATA_SMA_arr[dc, ac, 1]   = 1.*nr_simOK  
            info_per_DATA_SMA_arr[dc, ac, 2]   = SMA_arr[ac]
            
            for oc in range(0,3):
                #define:
                pos_oc123                   = np.where(savetestarr[dc, ac, :, 1] == (oc+1))[0]
                pos_id3_tlim_oc123          = list(set(pos_id3_tlim).intersection(pos_oc123))
                pos_id5_tlim_oc123          = list(set(pos_id5_tlim).intersection(pos_oc123))
                pos_id5_tlim_oc123_2ndOK    = list(set(pos_id5_tlim_oc123).intersection(pos_2ndOK)) 
                
                #cross section info:
                nr_id3_tlim_oc123       = len(pos_id3_tlim_oc123)
                nr_id5_tlim_oc123       = len(pos_id5_tlim_oc123)
                nr_id5_tlim_oc123_2ndOK = len(pos_id5_tlim_oc123_2ndOK)
                info_per_DATA_SMA_OC_arr[dc, ac, oc, 0]   = 1.*nr_id3_tlim_oc123        #binary-single ejection
                info_per_DATA_SMA_OC_arr[dc, ac, oc, 1]   = 1.*nr_id5_tlim_oc123        #3-body merger   
                info_per_DATA_SMA_OC_arr[dc, ac, oc, 2]   = 1.*nr_id5_tlim_oc123_2ndOK  #2nd binary   
                
                #HISTOGRAMS:
                #define:
                pos2bd  = pos_id3_tlim_oc123        #short notation
                pos3bd  = pos_id5_tlim_oc123        #short notation
                pos2nd  = pos_id5_tlim_oc123_2ndOK  #short notation
                #fGW histograms:
                Htype_index = 0 # =0 for fGW, =1 for ...
                min_H       = -5
                max_H       = 3 
                nr_Hbins    = 50
                #Hist 1: fGW-2body
                outc_index      = 0                                                  
                data_ij         = fGWij_SI_arr[pos2bd]
                Hdata           = np.log10(data_ij)
                hist, bin_edges = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
                bin_centers     = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
                HIST_data_arr[dc,ac,oc, outc_index, Htype_index, 0,:]   = bin_centers
                HIST_data_arr[dc,ac,oc, outc_index, Htype_index, 1,:]   = hist
                #Hist 2: fGW-3body  
                outc_index      = 1                                                                                                  
                data_ij         = fGWij_SI_arr[pos3bd]
                Hdata           = np.log10(data_ij)
                hist, bin_edges = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
                bin_centers     = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
                HIST_data_arr[dc,ac,oc, outc_index, Htype_index, 0,:]   = bin_centers
                HIST_data_arr[dc,ac,oc, outc_index, Htype_index, 1,:]   = hist
                #Hist 3: fGW-2ndBIN   
                outc_index      = 2                                                                                                                                                 
                data_ij         = fGW12_SI_arr[pos2nd]
                Hdata           = np.log10(data_ij)
                hist, bin_edges = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
                bin_centers     = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
                HIST_data_arr[dc,ac,oc, outc_index, Htype_index, 0,:]   = bin_centers
                HIST_data_arr[dc,ac,oc, outc_index, Htype_index, 1,:]   = hist
            
            
            #impossible to find out this system...:
            #HIST_data_arr = np.zeros((nrdataset, nr_SMA, 3, 5, 5, 2, 50), dtype=np.float64)        [nr dataset, nr SMA, objc (12,13,23), outc id (2-body, 3-body,..), hist type (fGW, ecc,..), x/y axis, hist vals]         
                        
        #------------------------------------------------------
    #----------------------------------------------------------
    #write to file:
    #----------------------------------------------------------
    saveoutput_arr  = info_per_DATA_SMA_arr  
    np.savez(save_data_filename + 'info_per_DATA_SMA_arr', saveoutput_arr)

    saveoutput_arr  = info_per_DATA_SMA_OC_arr  
    np.savez(save_data_filename + 'info_per_DATA_SMA_OC_arr', saveoutput_arr)
    
    saveoutput_arr  = HIST_data_arr  
    np.savez(save_data_filename + 'HIST_data_arr', saveoutput_arr)

    saveoutput_arr  = savetestarr  
    np.savez(save_data_filename + 'savetestarr', saveoutput_arr)
    #----------------------------------------------------------
#--------------------------------------------------------------
#--------------------------------------------------------------

#NOTES:
#cut on ejection velocity, or?
#clean up and CHECK, prepare ecc dist for 12 13 23,   
#Hills radius?
#calc fpeak for 2nd merger.
#use MR2,a for calc fpeak for 3-body merger.
#separate retro- vs. prograde interactions.
#start sim with adn model. from that generate correct P(e).
#funtional fits.
#extand sim to larger SMA?
#different mass ratios
#the outcome will be a 1m+2m binary that will keep interacting OR merge.
#make most of the letter on the equal mass case. Finish with arguing that 1+2 systems are naturally build up and interacts. Therefore end with a 1+2 - 1 case and show the results are still interesting.
#how is the hardening working?
#MYABE USE: MRIMS (ini_a_SI    = output_Nbody_xtra_info_REAL[c,3]*R_sun_SI...) to make sure we catch the IMS inspirals at the best fGW... ????


#--------------------------------------------------------------
if (plotetc_data_YN == 1):
#--------------------------------------------------------------
    #input data file for plotting etc.:
    open_data_filename          = 'R_0_ISO'#'RVA_1AU_FA05_1000000_402020'
    #'RVA_1AU_FA5_1000000_402020'
    #'RVA_0_003_01_03_1_2_3_4_5_10_15_25_45_75_90_ISO'
    #'RVA_1AU_FA05_1000000_402020'
    #RVA_1AU_FA5_1000000_402020
    #R_0_ISO
    #'RVA_0_003_01_03_1_2_3_4_5_10_15_25_45_75_90_ISO' #'T1_1AUOA5_100000_402020' #'T1_OA5_202020' #'T1_OA10_202020' #T3_BH402020_1AU_100000_CP'#T5_BH202020_1e21e2_10_10000' # 'T5_BH202020_1e21e2_10_10000_CP'
    #T4_BH202020_1e21e2_10_10000_CP'#T3_BH202020_1e250_10_10000_CP'#T3_BH202020_1e150_10_10000_CP'#T3_BH402020_1AU_50000_CP'#T3_BH202020_1e21e2_10_10000'#'nt1'#'ST4'#raw_input('open filename: ')

    #datasets used for paper:
    #'RVA_1AU_FA5_1000000_402020'
    #R_0_ISO
    #RVA_0_01_03_1_2_3_4_5_10_15_25_45_75_90_ISO
    
    output_arr                  = np.load(open_data_filename + 'info_per_DATA_SMA_arr' + '.npz')['arr_0']
    info_per_DATA_SMA_arr       = output_arr
    
    output_arr                  = np.load(open_data_filename + 'info_per_DATA_SMA_OC_arr' + '.npz')['arr_0']
    info_per_DATA_SMA_OC_arr    = output_arr
    
    output_arr                  = np.load(open_data_filename + 'HIST_data_arr' + '.npz')['arr_0']
    HIST_data_arr               = output_arr

    output_arr                  = np.load(open_data_filename + 'savetestarr' + '.npz')['arr_0']
    savetestarr                 = output_arr
    
    #DEFINE:
    SMA_arr     = info_per_DATA_SMA_arr[0,:,2]
    nr_SMA      = len(SMA_arr)
    nrdataset   = len(info_per_DATA_SMA_arr[:,0,0])
    SMA_arr_SI  = SMA_arr*R_sun_SI
    SMA_arr_AU  = SMA_arr_SI/AU_SI
        
    #info_per_DATA_SMA_arr[dc, ac, 0]   = 1.*nr_NAN                     #sim did not finish
    #info_per_DATA_SMA_arr[dc, ac, 1]   = 1.*nr_simOK  
    #info_per_DATA_SMA_arr[dc, ac, 2]   = SMA_arr[ac]
    
    #info_per_DATA_SMA_OC_arr[dc, ac, oc, 0]   = 1.*nr_id3_tlim_ij      #binary-single ejection
    #info_per_DATA_SMA_OC_arr[dc, ac, oc, 1]   = 1.*nr_id5_tlim_ij      #3-body merger
        
    #HIST_data_arr = np.zeros((nrdataset, nr_SMA, 3, 5, 5, 2, 50), dtype=np.float64)        [nr dataset, nr SMA, objc (12,13,23),       outc id (2-body, 3-body,..), hist type (fGW, ecc,..), x/y axis, hist vals]         
                
    #savetestarr[dc, ac, :, 0]   = id_arr[:]
    #savetestarr[dc, ac, :, 1]   = 1,2,3#comb is here referred to as: `1,2,3'        
    #savetestarr[dc, ac, :, 2]   = tGWij_SI_arr[:]
    #savetestarr[dc, ac, :, 3]   = 1/0 #=1 IF tGW<tlim
    #savetestarr[dc, ac, :, 4]   = fGWij_SI_arr #1st GW merger
    #savetestarr[dc, ac, :, 5]   = 1/0/-1       #=1 if 2nd binary is OK and bound. only use = 1 cases here.
    #savetestarr[dc, ac, :, 6]   = fGW12_SI     #2nd GW merger
    #savetestarr[dc, ac, :, 7]           = ecc_A (if > fGW_limit then = -1)   
    #savetestarr[dc, ac, :, 8]           = ecc_B (if > fGW_limit then = -1)
    #savetestarr[dc, ac, :, 9]           = Lij_ang_arr[:]
       
        
    #----------------------------------------------------------
    #PLOT:
    #----------------------------------------------------------
    
    
    
    #illustration script (for students, etc.):
#    open_data_filename          = 'RVA_0_003_01_03_1_2_3_4_5_10_15_25_45_75_90_ISO'
#    data_all_array              = np.load(open_data_filename + 'savetestarr' + '.npz')['arr_0']
#    #define:
#    sec_year        = 31536000.
#    dc  = 1         #dataset index.
#    id_arr          = data_all_array[dc, 0, :, 0]           
#    obij_arr        = data_all_array[dc, 0, :, 1]
#    tGWij_SI_arr    = data_all_array[dc, 0, :, 2] 
#    fGWij_SI_arr    = data_all_array[dc, 0, :, 4] 
#   ecc_A_arr       = data_all_array[dc, 0, :, 7]   #(if > fGW_limit then = -1)
#    ecc_B_arr       = data_all_array[dc, 0, :, 8]   #(if > fGW_limit then = -1) 
#    Lij_ang_arr     = data_all_array[dc, 0, :, 9]   #in degrees
#    pos_LTtlim        = np.where(tGWij_SI_arr < (10**5)*sec_year) 
#    pos_GTtlim        = np.where(tGWij_SI_arr > (10**5)*sec_year)     
#    #plot:
#    fig = plt.figure(figsize=(8.0, 8.0))
#    ax1 = fig.add_subplot(111)
#    ax1.plot(fGWij_SI_arr[pos_LTtlim], Lij_ang_arr[pos_LTtlim], marker = '.', markersize = 2.0, linewidth = 0, color = 'red')
#    ax1.plot(fGWij_SI_arr[pos_GTtlim], Lij_ang_arr[pos_GTtlim], marker = '.', markersize = 0.5, linewidth = 0, color = 'black')
#    ax1.set_xlim(1e-8,1e3)
#    ax1.set_xscale('log')
#    #show:
#    plt.show()    
#    exit()    
    
    
    
    
    
    
    
    
    


    #Fig P(a): (EQUAL MASS ONLY)
    mBH_SI      = m1_SI
    RsBH_SI     = 2.*G_new_SI*mBH_SI/c_SI**2.
    Nims        = 20.

    EMC_nrid_arr    = np.add(np.add(info_per_DATA_SMA_OC_arr[:,:,0,:], info_per_DATA_SMA_OC_arr[:,:,1,:]), info_per_DATA_SMA_OC_arr[:,:,2,:])             

    #Figure 1:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    #DATA: (loop over dataset)
    dc_plotYN  = [1,1]
    dc_ls_arr  = ['', '']          
    dc_lw_arr  = [0.5, 0.5]           
    dc_al_arr  = [0.5, 0.5]            
    dc_ms_arr  = [8, 6]        
    dc_mt_arr  = ['P', 'o']
    legnd_arr  = [r' (2D)', r' (3D)']
    
    P3b_NUM_arr_DS12    = np.zeros((2, len(SMA_arr_AU)), dtype=np.float64)
    P2b_NUM_arr_DS12    = np.zeros((2, len(SMA_arr_AU)), dtype=np.float64)
    
    for dc in range(0,nrdataset):
        if (dc_plotYN[dc] == 1):    
            sma_l = 0
            sma_u = len(SMA_arr_AU)
            #3-body mergers:
            nr_3b       = 1.*EMC_nrid_arr[dc,:,1]
            nr_simOK    = 1.*info_per_DATA_SMA_arr[dc,:,1]
            P3b_arr     = nr_3b/nr_simOK
            P3b_arr_err = np.sqrt(nr_3b)/nr_simOK
            legnd       = legnd_arr[dc]
            ax1.errorbar(   SMA_arr_AU, P3b_arr, yerr=P3b_arr_err, label = r'3-body merg.'+ legnd, color = 'deepskyblue', markeredgecolor = 'black', marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
            P3b_NUM_arr_DS12[dc,:]  = P3b_arr[:]
            #2-body mergers:
            nr_2b       = 1.*EMC_nrid_arr[dc,:,0]
            nr_simOK    = 1.*info_per_DATA_SMA_arr[dc,:,1]
            P2b_arr     = nr_2b/nr_simOK
            P2b_arr_err = np.sqrt(nr_2b)/nr_simOK
            legnd       = legnd_arr[dc]
            ax1.errorbar(   SMA_arr_AU, P2b_arr, yerr=P2b_arr_err, label = r'2-body merg.'+ legnd, color = 'deeppink', markeredgecolor = 'black',    marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
            P2b_NUM_arr_DS12[dc,:]  = P2b_arr[:]
            #ecc mergers:
            P_eccB_arr      = np.zeros(nr_SMA, dtype=np.float64)
            P_eccB_arr_err  = np.zeros(nr_SMA, dtype=np.float64)            
            for ac in range(0,nr_SMA):
                nr_simOK        = 1.*info_per_DATA_SMA_arr[dc, ac, 1]
                ecc_B_arr       = savetestarr[dc, ac, :, 8]
                nr_ecc_B        = 1.*len(np.append(np.where(ecc_B_arr[:] > 0.1)[0], np.where(ecc_B_arr[:] == -1.0)[0]))
                P_eccB_arr[ac]      = nr_ecc_B/nr_simOK
                P_eccB_arr_err[ac]  = np.sqrt(nr_ecc_B)/nr_simOK
            ax1.errorbar(   SMA_arr_AU, P_eccB_arr, yerr=P_eccB_arr_err, label = r'$e>0.1$ at $>10\ Hz$'+ legnd, color = 'black', markeredgecolor = 'black', marker = dc_mt_arr[dc], markersize = 0.75*dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
            ax1.plot(       SMA_arr_AU, P_eccB_arr, color = 'black', linestyle = ':', alpha = 0.5)            
             

    #ANALYTICAL calc:
    SMA_AU_plotarr  = 10.**(np.linspace(-4., 5.0, 1000))
    SMA_SI_plotarr  = SMA_AU_plotarr*AU_SI
    tGW_SI          = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((SMA_SI_plotarr**4.)/(mBH_SI**3.))
    Torb_a          = (2.*np.pi*np.sqrt(SMA_SI_plotarr**3./(2.*G_new_SI*mBH_SI)))

    #3-body merger:
    #ISO:
    P3b_arr         = 1.2*Nims*(Torb_a/tGW_SI)**(2./7.)
    pos             = np.where(P3b_arr > 1.0)[0]
    P3b_arr[pos]    = 1.0 
    P3b_arr_ISO     = P3b_arr
    ax1.plot(SMA_AU_plotarr, P3b_arr_ISO,       color = 'deepskyblue',      linestyle = '-', linewidth = 1.0, alpha = 0.75)
    #CP:
    P3b_arr         = 1.0*Nims*(Torb_a/tGW_SI)**(1./7.) #(1.05)*np.sqrt(2.)*Nims*(RsBH_SI/SMA_SI_plotarr)**(5./14.) 
    pos             = np.where(P3b_arr > 1.0)[0]
    P3b_arr[pos]    = 1.0 
    P3b_arr_CP      = P3b_arr
    ax1.plot(SMA_AU_plotarr, P3b_arr_CP,       color = 'deepskyblue',      linestyle = '--', linewidth = 1.0, alpha = 0.75)    
    #3b area:
    #ax1.fill_between(SMA_AU_plotarr, P3b_arr_ISO, P3b_arr_CP, color = 'deepskyblue', alpha = 0.15)
    ax1.fill_between(SMA_arr_AU, P3b_NUM_arr_DS12[0,:], P3b_NUM_arr_DS12[1,:], color = 'deepskyblue', alpha = 0.15)#, step='mid')


    #2-body merger:
    #ISO:
    P2b_arr         = 1.3*(tint_SI/tGW_SI)**(2./7.)
    pos             = np.where(P2b_arr > 1.0)[0]
    P2b_arr[pos]    = 1.0 
    P2b_arr_ISO     = P2b_arr
    ax1.plot(SMA_AU_plotarr, P2b_arr_ISO,       color = 'deeppink',      linestyle = '-', linewidth = 1.0, alpha = 0.75)
    #CP:
    P2b_arr         = 1.2*(tint_SI/tGW_SI)**(1./7.)
    pos             = np.where(P2b_arr > 1.0)[0]
    P2b_arr[pos]    = 1.0 
    P2b_arr_CP      = P2b_arr
    ax1.plot(SMA_AU_plotarr, P2b_arr_CP,       color = 'deeppink',      linestyle = '--', linewidth = 1.0, alpha = 0.75)
    #2b area:
    #ax1.fill_between(SMA_AU_plotarr, P2b_arr_ISO, P2b_arr_CP, color = 'deeppink', alpha = 0.15)
    ax1.fill_between(SMA_arr_AU, P2b_NUM_arr_DS12[0,:], P2b_NUM_arr_DS12[1,:], color = 'deeppink', alpha = 0.15)#, step='mid')


    #axis/plot settings:
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = True, ncol=1, framealpha = 1.0)
    ax1.text(0.1, 0.95, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11, color = 'black')
    ax1.text(0.1, 0.88, r'$t_{\rm GW} < 10^{5}\ yrs.$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11, color = 'black')
    ax1.set_xlabel(r'semi-major axis $a$ [AU]')
    ax1.set_ylabel(r'Probability')
    ax1.set_xlim(0.02, 2*(1e2))
    ax1.set_ylim(1e-3, 2.0)
    ax1.set_facecolor('xkcd:white')

    plt.savefig(open_data_filename + '_Pfig_ISOCPL.pdf', bbox_inches='tight')
    plt.show()

    exit()
 
    
    
    
    
    
    
    
    
    #MANY THINGS IN THE SCRIPT HAVE TO CHANGE WHEN SWITHCING BETWEEN
    #DIFFERET PLOTS/DATASET. MAKE IT MORE EASY TO AVIOD MISTAKES!      
    
    #SETTINGS FOR FIG 1,2:
    #define:
    color123_arr        = ['blue', 'forestgreen', 'red']
    label123_arr        = ['12', '13', '23']
    label123_2b_arr     = ['[1,2]-2b', '[1,3]-2b', '[2,3]-2b']
    label123_3b_arr     = ['[1,2]-3b', '[1,3]-3b', '[2,3]-3b']
    label123_all_arr    = ['all-[1,2]', 'all-[1,3]', 'all-[2,3]']

    #choose dataset (index):
    dc  = 0
    ac  = 0
    #define:
    pos_tlim            = np.where(savetestarr[dc, ac, :, 3] == 1)[0] 
    pos_2b              = np.where(savetestarr[dc, ac, :, 0] == 3)[0]
    pos_3b              = np.where(savetestarr[dc, ac, :, 0] == 5)[0]
    pos_2ndOK           = np.where(savetestarr[dc, ac, :, 5] == 1)[0]
    pos_2b_tlim         = list(set(pos_2b).intersection(pos_tlim)) 
    pos_3b_tlim         = list(set(pos_3b).intersection(pos_tlim)) 
    pos_3b_tlim_2ndOK   = list(set(pos_3b_tlim).intersection(pos_2ndOK)) 
    nr_simOK            = info_per_DATA_SMA_arr[dc, ac, 1]  
    ecc_A_arr           = savetestarr[dc, ac, :, 7]
    ecc_B_arr           = savetestarr[dc, ac, :, 8]
    Lij_ang_arr         = savetestarr[dc, ac, :, 9]
    pos_eccB            = np.append(np.where(ecc_B_arr[:] > 0.1)[0], np.where(ecc_B_arr[:] == -1.0)[0])
    fGW_1bin_arr        = savetestarr[dc, ac, :, 4]
    
 
 
    fig = plt.figure(figsize=(10.0, 4.0))        
    ax1 = fig.add_subplot(121)    
    ax2 = fig.add_subplot(122)    
    
    #FIG A:
    #hist settings:
    min_H           = -3.
    max_H           = 3.
    nr_Hbins        = 50
    #define data:   
    fGWij_SI_arr    = savetestarr[dc, ac, :, 4] #1st GW merger   
    #all:
    data_ij         = fGWij_SI_arr[pos_tlim]
    Hdata           = np.log10(data_ij)
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    #ax1.step(Hx, Hy, label = 'all', color = 'grey', alpha=0.75, linewidth = 1)
    ax1.fill_between(Hx, 0, Hy, label = r'all merg.', color = 'grey', step="pre", alpha = 0.35)
    maxH_all        = max(Hy)
    #all-3body
    data_ij         = fGWij_SI_arr[pos_3b_tlim]
    Hdata           = np.log10(data_ij)
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    ax1.fill_between(Hx, 0, Hy, label = r'$3$-body merg.', color = 'deepskyblue', step="pre", alpha = 0.75, hatch = '\\\\\\')
    ax1.step(Hx, Hy, color = 'deepskyblue', alpha=0.25, linewidth = 2)    

    #all highly-eccentric
    data_ij         = fGWij_SI_arr[pos_eccB]
    Hdata           = np.log10(data_ij)
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    #ax1.fill_between(Hx, 0, Hy, color = 'black', step="pre", alpha = 0.35, hatch = '\\\\\\')
    ax1.step(Hx, Hy, color = 'black', alpha=0.5, linestyle = ':', linewidth = 1.5, label = r'$e>0.1$ at $>10\ Hz$')    

    #all-2body
    data_ij         = fGWij_SI_arr[pos_2b_tlim]
    Hdata           = np.log10(data_ij)
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    ax1.fill_between(Hx, 0, Hy, label = r'$2$-body merg.', color = 'deeppink', step="pre", alpha = 0.75, hatch = '///') 
    ax1.step(Hx, Hy, color = 'deeppink', alpha=0.25, linewidth = 2)
    #Loop over 12, 13, 23:
    for oc in range(0,3):
        pos_oc          = np.where(savetestarr[dc, ac, :, 1] == (oc+1))[0]
        pos_2b_tlim_oc  = list(set(pos_2b_tlim).intersection(pos_oc))  
        pos_3b_tlim_oc  = list(set(pos_3b_tlim).intersection(pos_oc))  
        #plot 2body:
        data_ij         = fGWij_SI_arr[pos_2b_tlim_oc]
        Hdata           = np.log10(data_ij)
        Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
        Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
        #ax1.step(Hx, Hy, label = label123_2b_arr[oc], color = color123_arr[oc])
        #plot 3body:
        data_ij         = fGWij_SI_arr[pos_3b_tlim_oc]
        Hdata           = np.log10(data_ij)
        Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
        Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
        #ax1.step(Hx, Hy, label = label123_3b_arr[oc], color = color123_arr[oc], alpha= 1, linestyle = '--')
    #labels, etc..:
    print SMA_arr_AU[ac]
    print len(pos_tlim)
    ax1.text(0.05, 0.95, r'[$40M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.88, r'$a = 1$ AU, $e = 0$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.81, r'$t_{\rm GW} < 10^{5}\ yrs.$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.74, r'$\psi = 0.5$ deg.', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    
    ax1.set_xlim(min_H, max_H)
    ax1.set_ylim(0.0, 1.5*maxH_all)  
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 9, frameon = True, ncol=1, framealpha = 1.0)
    ax1.set_xlabel(r'log $f_{GW}$ [Hz]')
    ax1.set_ylabel(r'Number')
    #ax1.set_title(r'GW peak frequency at formation')
    #plt.savefig(open_data_filename + 'fGW_dist_a.pdf', bbox_inches='tight')
    #plt.show()

    #FIG B:
    #hist settings:
    min_H           = -4
    max_H           = 0 
    nr_Hbins        = 50
    normfac         = 1.*len(pos_tlim)
    #define data:
    ecc_fGW_arr         = savetestarr[dc, ac, :, 8] #at 10 Hz
    pos                 = np.where(ecc_fGW_arr[:] < 0.0)[0]
    ecc_fGW_arr[pos]    = 1.0    
    #all:
    Hdata           = np.log10(ecc_fGW_arr[pos_tlim])
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    CHy             = np.cumsum(Hy[::-1])[::-1]/normfac
    ax2.fill_between(Hx, 0, CHy, label = r'all merg.', color = 'grey', step="pre", alpha = 0.35)    
    #3body:
    Hdata           = np.log10(ecc_fGW_arr[pos_3b_tlim])
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    CHy             = np.cumsum(Hy[::-1])[::-1]/normfac
    ax2.fill_between(Hx, 0, CHy, label = r'$3$-body merg.', color = 'deepskyblue', step="pre", alpha = 0.75, hatch = '\\\\\\')
    ax2.step(Hx, CHy, color = 'deepskyblue', alpha=0.25, linewidth = 2)        
    #2body:
    Hdata           = np.log10(ecc_fGW_arr[pos_2b_tlim])
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    CHy             = np.cumsum(Hy[::-1])[::-1]/normfac
    ax2.fill_between(Hx, 0, CHy, label = r'$2$-body merg.', color = 'deeppink', step="pre", alpha = 0.75, hatch = '///')
    ax2.step(Hx, CHy, color = 'deeppink', alpha=0.25, linewidth = 2)        
    #Loop over 12, 13, 23:
    for oc in range(0,3):
        pos_oc          = np.where(savetestarr[dc, ac, :, 1] == (oc+1))[0]
        pos_tlim_oc     = list(set(pos_tlim).intersection(pos_oc))  
        #plot:
        Hdata           = np.log10(ecc_fGW_arr[pos_tlim_oc])
        Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
        Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
        CHy             = np.cumsum(Hy[::-1])[::-1]/normfac
        #ax2.step(Hx, CHy, label = label123_all_arr[oc], color = color123_arr[oc], alpha=1, linestyle = ':')
    #plot settings:
    ax2.set_xlim(min_H, max_H)
    ax2.set_ylim(0.0, 1.3)  
    ax2.text(0.05, 0.95, r'[$40M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes, fontsize = 11)
    ax2.text(0.05, 0.88, r'$a = 1$ AU, $e = 0$', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes, fontsize = 11)
    ax2.text(0.05, 0.81, r'$t_{\rm GW} < 10^{5}\ yrs.$', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes, fontsize = 11)
    ax2.text(0.05, 0.74, r'$\psi = 0.5$ deg.', horizontalalignment='left', verticalalignment='center', transform=ax2.transAxes, fontsize = 11)
    ax2.set_xlim(min_H, max_H)
    ax2.legend(loc='upper right', numpoints = 1, fontsize = 9, frameon = True, ncol=1, framealpha = 1.0)
    ax2.set_xlabel(r'log $e$')
    ax2.set_ylabel(r'Probability $e(\geq 10Hz) > e$')
    #ax2.set_title(r'orbital eccentricity at $10\ Hz$')
    plt.savefig(open_data_filename + 'ecc_fGW_Cdist.pdf', bbox_inches='tight')
    plt.show()
    
    exit()
     
    
    
    
    
    
    
    
    
    
    
    


    #FIG angles:
    fig = plt.figure(figsize=(5.0, 4.0))   
    ax1 = fig.add_subplot(111)
        
    #hist settings:
    min_H           = -5.0
    max_H           = 185.
    nr_Hbins        = 100

    #THIS WORKS ONLY FOR ecc e0 at 10 hz!
    pos_mask        = np.where(ecc_B_arr[:] == -1.0)[0]
    pos_e_GT_e0_Xhz = np.append(np.where(ecc_B_arr[:] > 0.1)[0], pos_mask)
    pos_e_LT_e0_Xhz = list(set(np.where(ecc_B_arr[:]  < 0.1)[0]).difference(pos_mask))
    #add merger time limit:
    pos_e_GT_e0_Xhz_tlim    = list(set(pos_e_GT_e0_Xhz).intersection(pos_tlim)) 
    pos_e_LT_e0_Xhz_tlim    = list(set(pos_e_LT_e0_Xhz).intersection(pos_tlim)) 
    
    Nfactor = 100.0    
    #plot:
    pos             = pos_e_LT_e0_Xhz_tlim
    data_ij         = Lij_ang_arr[pos]
    Hdata           = data_ij
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=True)
    Hy              = Hy*Nfactor
    Hy_e_LT_e0_Xhz_tlim = Hy
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    ax1.fill_between(Hx, 0, Hy, color = 'orange', step="pre", alpha = 0.5, label = r'$e<0.1$ at $10\ Hz$')#, hatch = '///')
    ax1.step(Hx, Hy, color = 'black', alpha=1, linewidth = 1)
    
    #plot:
    pos             = pos_e_GT_e0_Xhz_tlim
    data_ij         = Lij_ang_arr[pos]
    Hdata           = data_ij
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=True)
    Hy              = Hy*Nfactor
    Hy_e_GT_e0_Xhz_tlim = Hy
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    ax1.fill_between(Hx, 0, Hy, color = 'black', step="pre", alpha = 0.5, label = r'$e>0.1$ at $10\ Hz$')#, hatch = '\\\\\\')
    #ax1.step(Hx, Hy, color = 'black', alpha=0.25, linewidth = 1)        

    #plot:
    iso_dist        = np.sin(2.*np.pi*(Hx/360.))
    iso_dist_norm   = (iso_dist/(sum(iso_dist)*(Hx[1] - Hx[0])))*Nfactor
    ax1.plot(Hx, iso_dist_norm, color = 'black', linestyle = ':', linewidth = 1, label = r'isotropic')
    #maxISO          = max(iso_dist_norm)

    ax1.text(0.05, 0.95, r'[$40M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.88, r'$a = 1$ AU, $e = 0$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.81, r'$t_{\rm GW} < 10^{5}\ yrs.$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.74, r'$\psi = 0.5$ deg.', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    
    ax1.set_xlabel(r'$\theta$ [degrees]')
    ax1.set_ylabel(r'Number [normalized]')
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = True, ncol=1, framealpha = 1.0)
    
    ax1.set_xlim(-5, 180)  
    ax1.set_ylim(0.*max(Hy_e_LT_e0_Xhz_tlim), 1.35*max(Hy_e_LT_e0_Xhz_tlim))  
    
    plt.savefig(open_data_filename + 'angle_z_dist.pdf', bbox_inches='tight')
    plt.show()

    exit()




    
    
    #Fig P(theta): (EQUAL MASS ONLY)
    #RVA_0_003_01_03_1_2_3_4_5_10_15_25_45_75_90_ISO
    
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    
    P_eccB_arr      = np.zeros(nrdataset, dtype=np.float64)
    P_eccB_arr_err  = np.zeros(nrdataset, dtype=np.float64)
    P_Mtot_arr      = np.zeros(nrdataset, dtype=np.float64)                
    P_cir_arr       = np.zeros(nrdataset, dtype=np.float64)                
                    
    
    #DATA: (loop over dataset)
    ac = 0
    for dc in range(0, nrdataset):
        #ecc mergers:
        nr_simOK        = 1.*info_per_DATA_SMA_arr[dc, ac, 1]
        ecc_B_arr       = savetestarr[dc, ac, :, 8]
        nr_ecc_B        = 1.*len(np.append(np.where(ecc_B_arr[:] > 0.1)[0], np.where(ecc_B_arr[:] == -1.0)[0]))
        P_eccB_arr[dc]      = nr_ecc_B/nr_simOK
        P_eccB_arr_err[dc]  = np.sqrt(nr_ecc_B)/nr_simOK
        #all mergers:
        merg_y1no0_arr  = savetestarr[dc, ac, :, 3] 
        nr_Mtot         = 1.*len(np.where(merg_y1no0_arr[:] == 1)[0])
        P_Mtot_arr[dc]  = nr_Mtot/nr_simOK
        
    #NOTE: LAST should be ISO for COMPARISON and FIRST should be CO-PLANAR!
    
    plotmask_arr    = np.array([1,1,1,1,1,1,0,0,1,1,0,1,0,0,1])
    lastpos         = nrdataset-1    
    theta_arr       = np.array([0., 0.03, 0.1, 0.3, 1., 2., 3., 4., 5., 10., 15., 25., 45., 75., 90.]) #nr. here has to be = nrdataset - 1 (last MUST BE ISO)    
    x_values        = theta_arr[0:lastpos]

    #P eccentric mergers:
    y_values        = P_eccB_arr[0:lastpos]
    yerr_values     = P_eccB_arr_err[0:lastpos]
    #prepare for masking arrays - 'conventional' arrays won't do it
    y_values        = np.ma.array(y_values)
    #mask based on: plotmask_arr
    y_values        = np.ma.masked_where(plotmask_arr[:] != 1, y_values)
    y_P_ecc_arr     = y_values
    ax1.errorbar(   x_values, y_values, yerr=yerr_values, marker = 'X', markersize = 8, color = 'lime', markeredgecolor = 'black', linestyle = '', alpha = 0.5, label = r'$N_{m}(e>0.1:\geq 10Hz)/N_{sc}$')
    #ax1.plot(       theta_arr[0:lastpos], P_eccB_arr[0:lastpos], marker = 'X', markersize = 8, color = 'black',  markeredgecolor = 'black', linestyle = '', alpha = 0.5)#, linestyle = ':', alpha = 0.5)            
    #ax1.plot(       theta_arr[0:lastpos], 0.0125*(theta_arr[0:lastpos])**(-0.7), linestyle = ':')            

    #CPL for comparison:
    #ax1.plot([-10,100], [P_eccB_arr[0], P_eccB_arr[0]], color = 'black', linestyle = '--', alpha = 0.5, linewidth = 2, label = 'Co-planar (2D)')            
    #ax1.fill_between([-10,100], (P_eccB_arr[0]+P_eccB_arr_err[0]), (P_eccB_arr[0]-P_eccB_arr_err[0]), color = 'black', alpha = 0.25)
    #ISO for comparison:
    #ax1.plot([-10,100], [P_eccB_arr[lastpos], P_eccB_arr[lastpos]], color = 'black', linestyle = '-', alpha = 0.5, linewidth = 1, label = 'Isotropic (3D)')            
    #ax1.fill_between([-10,100], (P_eccB_arr[lastpos]+P_eccB_arr_err[lastpos]), (P_eccB_arr[lastpos]-P_eccB_arr_err[lastpos]), color = 'lime', alpha = 0.25)
    
    y_values        = P_Mtot_arr[0:lastpos]
    y_values        = np.ma.array(y_values)
    y_values        = np.ma.masked_where(plotmask_arr[:] != 1, y_values)
    y_P_all_arr     = y_values
    ax1.plot(x_values, y_P_ecc_arr/y_P_all_arr, marker = 'X', markersize = 8, color = 'black', markeredgecolor = 'black', linestyle = '', alpha = 0.5, label = r'$N_{m}(e>0.1:\geq 10Hz)/N_{m}(t_{GW} < 10^{5}\ yrs.)$')
    
    #settings etc:
    ax1.set_xlabel(r'$\psi$ [degrees]')
    ax1.set_ylabel(r'Probability')
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlim(0.02, 110)
    ax1.set_ylim(1e-3, 5.0)
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = True, ncol=1, framealpha = 1.0)
    ax1.text(0.05, 0.95-0.8, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    ax1.text(0.05, 0.88-0.8, r'$a = 1$ AU, $e = 0$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11)
    
    plt.savefig(open_data_filename + '_PeccM_thetaz.pdf', bbox_inches='tight')
    plt.show()
    
    exit()
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    






 
 
 
 
 
 
 
 
 
 
 
 


 

























    #Fig P(a): (EQUAL MASS ONLY)
    mBH_SI      = m1_SI
    RsBH_SI     = 2.*G_new_SI*mBH_SI/c_SI**2.
    Nims        = 20.

    EMC_nrid_arr    = np.add(np.add(info_per_DATA_SMA_OC_arr[:,:,0,:], info_per_DATA_SMA_OC_arr[:,:,1,:]), info_per_DATA_SMA_OC_arr[:,:,2,:])             

    #Figure 1:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    #DATA: (loop over dataset)
    dc_plotYN  = [1,1]
    dc_ls_arr  = ['', '']          
    dc_lw_arr  = [0.5, 0.5]           
    dc_al_arr  = [0.5, 0.5]            
    dc_ms_arr  = [8, 6]        
    dc_mt_arr  = ['P', 'o']
    legnd_arr  = [r' (2D)', r' (3D)']
    
    for dc in range(0,nrdataset):
        if (dc_plotYN[dc] == 1):    
            sma_l = 0
            sma_u = len(SMA_arr_AU)
            #3-body mergers:
            nr_3b       = 1.*EMC_nrid_arr[dc,:,1]
            nr_simOK    = 1.*info_per_DATA_SMA_arr[dc,:,1]
            P3b_arr     = nr_3b/nr_simOK
            P3b_arr_err = np.sqrt(nr_3b)/nr_simOK
            legnd       = legnd_arr[dc]
            ax1.errorbar(   SMA_arr_AU, P3b_arr, yerr=P3b_arr_err, label = r'3-body merg.'+ legnd, color = 'deepskyblue', markeredgecolor = 'black', marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
            #2-body mergers:
            nr_2b       = 1.*EMC_nrid_arr[dc,:,0]
            nr_simOK    = 1.*info_per_DATA_SMA_arr[dc,:,1]
            P2b_arr     = nr_2b/nr_simOK
            P2b_arr_err = np.sqrt(nr_2b)/nr_simOK
            legnd       = legnd_arr[dc]
            ax1.errorbar(   SMA_arr_AU, P2b_arr, yerr=P2b_arr_err, label = r'2-body merg.'+ legnd, color = 'deeppink', markeredgecolor = 'black',    marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
            #ecc mergers:
            P_eccB_arr      = np.zeros(nr_SMA, dtype=np.float64)
            P_eccB_arr_err  = np.zeros(nr_SMA, dtype=np.float64)            
            for ac in range(0,nr_SMA):
                nr_simOK        = 1.*info_per_DATA_SMA_arr[dc, ac, 1]
                ecc_B_arr       = savetestarr[dc, ac, :, 8]
                nr_ecc_B        = 1.*len(np.append(np.where(ecc_B_arr[:] > 0.1)[0], np.where(ecc_B_arr[:] == -1.0)[0]))
                P_eccB_arr[ac]      = nr_ecc_B/nr_simOK
                P_eccB_arr_err[ac]  = np.sqrt(nr_ecc_B)/nr_simOK
            ax1.errorbar(   SMA_arr_AU, P_eccB_arr, yerr=P_eccB_arr_err, label = r'$e>0.1$ at $>10\ Hz$'+ legnd, color = 'black', markeredgecolor = 'black', marker = dc_mt_arr[dc], markersize = 0.75*dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
            ax1.plot(       SMA_arr_AU, P_eccB_arr, color = 'black', linestyle = ':', alpha = 0.5)            
             

    #ANALYTICAL calc:
    SMA_AU_plotarr  = 10.**(np.linspace(-4., 5.0, 1000))
    SMA_SI_plotarr  = SMA_AU_plotarr*AU_SI
    tGW_SI          = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((SMA_SI_plotarr**4.)/(mBH_SI**3.))
    Torb_a          = (2.*np.pi*np.sqrt(SMA_SI_plotarr**3./(2.*G_new_SI*mBH_SI)))

    #3-body merger:
    #ISO:
    P3b_arr         = 1.2*Nims*(Torb_a/tGW_SI)**(2./7.)
    pos             = np.where(P3b_arr > 1.0)[0]
    P3b_arr[pos]    = 1.0 
    P3b_arr_ISO     = P3b_arr
    ax1.plot(SMA_AU_plotarr, P3b_arr_ISO,       color = 'deepskyblue',      linestyle = '-', linewidth = 2.0, alpha = 0.75)
    #CP:
    P3b_arr         = 1.0*Nims*(Torb_a/tGW_SI)**(1./7.) #(1.05)*np.sqrt(2.)*Nims*(RsBH_SI/SMA_SI_plotarr)**(5./14.) 
    pos             = np.where(P3b_arr > 1.0)[0]
    P3b_arr[pos]    = 1.0 
    P3b_arr_CP      = P3b_arr
    ax1.plot(SMA_AU_plotarr, P3b_arr_CP,       color = 'deepskyblue',      linestyle = '--', linewidth = 2.0, alpha = 0.75)    
    #3b area:
    ax1.fill_between(SMA_AU_plotarr, P3b_arr_ISO, P3b_arr_CP, color = 'deepskyblue', alpha = 0.15)

    #2-body merger:
    #ISO:
    P2b_arr         = 1.3*(tint_SI/tGW_SI)**(2./7.)
    pos             = np.where(P2b_arr > 1.0)[0]
    P2b_arr[pos]    = 1.0 
    P2b_arr_ISO     = P2b_arr
    ax1.plot(SMA_AU_plotarr, P2b_arr_ISO,       color = 'deeppink',      linestyle = '-', linewidth = 2.0, alpha = 0.75)
    #CP:
    P2b_arr         = 1.2*(tint_SI/tGW_SI)**(1./7.)
    pos             = np.where(P2b_arr > 1.0)[0]
    P2b_arr[pos]    = 1.0 
    P2b_arr_CP      = P2b_arr
    ax1.plot(SMA_AU_plotarr, P2b_arr_CP,       color = 'deeppink',      linestyle = '--', linewidth = 2.0, alpha = 0.75)
    #2b area:
    ax1.fill_between(SMA_AU_plotarr, P2b_arr_ISO, P2b_arr_CP, color = 'deeppink', alpha = 0.15)

    #ax1.plot([0,0], [0,0],           label = r'analytical solutions',       color = 'black',      linestyle = ':', linewidth = 1.0, alpha = 1.0)

    #axis/plot settings:
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = True, ncol=1, framealpha = 1.0)
    ax1.text(0.1, 0.95, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11, color = 'black')
    ax1.text(0.1, 0.88, r'$t_{\rm GW} < 10^{5}\ yrs.$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 11, color = 'black')
    ax1.set_xlabel(r'semi-major axis $a$ [AU]')
    ax1.set_ylabel(r'Probability')
    ax1.set_xlim((0.75*(1e-2)), 2*(1e2))
    ax1.set_ylim(1e-3, 2.0)
    ax1.set_facecolor('xkcd:white')

    plt.savefig(open_data_filename + '_Pfig_ISOCPL.pdf', bbox_inches='tight')
    plt.show()

    exit()
 
 
 
 
 
 
 
 
 
 

 
    



 
 
 
    
    
    
    
    
    
    
    
    
    
    
    #plot 2:
    P_eccB_arr      = np.zeros(nr_SMA, dtype=np.float64)
    P_eccB_arr_err  = np.zeros(nr_SMA, dtype=np.float64)            
    testP_arr = np.zeros(nr_SMA, dtype=np.float64) 
    for ac in range(0,nr_SMA):
     
        
        nr_tlim         = len(np.where(savetestarr[dc, ac, :, 3] == 1)[0]) 
        nr_simOK        = info_per_DATA_SMA_arr[dc, ac, 1]
        ecc_B_arr       = savetestarr[dc, ac, :, 8]
        nr_ecc_B        = len(np.append(np.where(ecc_B_arr[:] > 0.1)[0], np.where(ecc_B_arr[:] == -1.0)[0]))
        P_eccB_arr[ac]      = (1.*nr_ecc_B)/(1.*nr_simOK)
        P_eccB_arr_err[ac]  = np.sqrt(1.*nr_ecc_B)/(1.*nr_simOK)
        #pos_3b          = np.where(savetestarr[dc, ac, :, 0] == 5)[0]
        #pos_tlim        = np.where(savetestarr[dc, ac, :, 3] == 1)[0] 
        #pos_3b_tlim     = list(set(pos_3b).intersection(pos_tlim)) 
        #nr_3b_tlim      = len(pos_3b_tlim)
        #testP_arr[ac]   = (1.*nr_3b_tlim)/(1.*nr_simOK)
        
    ax2.errorbar(SMA_arr_AU, P_eccB_arr,  yerr=P_eccB_arr_err, color = 'black')
    
    
    
   
    
    
    #            ax2.errorbar(SMA_arr_AU, P_eccB_arr,  yerr=P_eccB_arr_err, color = 'black')

    
    
    
    
    
    
    
    
        
    
         
    
    
    
    #FIG P(2/3-body mergers):
    if (EMcase_YN == 1):

         mBH_SI      = m1_SI
         RsBH_SI     = 2.*G_new_SI*mBH_SI/c_SI**2.
         Nims        = 20.

         EMC_nrid_arr    = np.add(np.add(info_per_DATA_SMA_OC_arr[:,:,0,:], info_per_DATA_SMA_OC_arr[:,:,1,:]), info_per_DATA_SMA_OC_arr[:,:,2,:])             

         #Figure 1:
         fig = plt.figure(figsize=(5.0, 4.0))
         ax1 = fig.add_subplot(111)
         ax1.set_xscale('log')
         ax1.set_yscale('log')

         #DATA: (loop over dataset)
         dc_ls_arr   = ['', '']   #ISO, CPL
         dc_lw_arr   = [1,2]         #ISO, CPL
         dc_al_arr   = [0.75,0.85]     #ISO, CPL        
         dc_ms_arr   = [5,20]        #ISO, CPL
         dc_mt_arr   = ['x','.']     #ISO, CPL
         legnd_arr   = [r'3D (cluster)', r'2D (disk)']
         for dc in range(0,nrdataset):
             if (dc == 0):
                 sma_l = 0
                 sma_u = len(SMA_arr_AU) #5
             if (dc == 1):
                 sma_l = 0
                 sma_u = len(SMA_arr_AU)    
             #3-body mergers:
             nr_id5      = EMC_nrid_arr[dc,:,1]
             nr_simOK    = info_per_DATA_SMA_arr[dc,:,1]
             P3b_arr     = nr_id5/nr_simOK
             P3b_arr_err = np.sqrt(nr_id5)/nr_simOK
             legnd       = legnd_arr[dc]
             #ax1.errorbar(SMA_arr_AU, P3b_arr,  yerr=P3b_arr_err, label = '3-body '+ legnd, color = 'deepskyblue',  marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
             ax1.plot(SMA_arr_AU[sma_l:sma_u], P3b_arr[sma_l:sma_u], label = '3-body merg. '+ legnd, color = 'deepskyblue',  marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
             #2-body mergers:
             nr_id3      = EMC_nrid_arr[dc,:,0]
             nr_simOK    = info_per_DATA_SMA_arr[dc,:,1]
             P2b_arr     = nr_id3/nr_simOK
             P2b_arr_err = np.sqrt(nr_id3)/nr_simOK
             legnd       = legnd_arr[dc]
             #ax1.errorbar(SMA_arr_AU, P2b_arr,  yerr=P2b_arr_err, label = '2-body '+ legnd, color = 'deeppink',     marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])
             ax1.plot(SMA_arr_AU[sma_l:sma_u], P2b_arr[sma_l:sma_u], label = '2-body merg. '+ legnd, color = 'deeppink',     marker = dc_mt_arr[dc], markersize = dc_ms_arr[dc], linestyle = dc_ls_arr[dc], linewidth = dc_lw_arr[dc], alpha = dc_al_arr[dc])            
    
         #ANALYTICAL calc:
         SMA_AU_plotarr  = 10.**(np.linspace(-4., 5.0, 1000))
         SMA_SI_plotarr  = SMA_AU_plotarr*AU_SI
         tGW_SI          = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((SMA_SI_plotarr**4.)/(mBH_SI**3.))
         Torb_a          = (2.*np.pi*np.sqrt(SMA_SI_plotarr**3./(2.*G_new_SI*mBH_SI)))

         #2-body merger:
         #ISO:
         P2b_arr         = 1.3*(tint_SI/tGW_SI)**(2./7.)
         pos             = np.where(P2b_arr > 1.0)[0]
         P2b_arr[pos]    = 1.0 
         ax1.plot(SMA_AU_plotarr, P2b_arr,       color = 'deeppink',      linestyle = ':', linewidth = 0.5, alpha = 0.75)
         #CP:
         P2b_arr         = 1.2*(tint_SI/tGW_SI)**(1./7.)
         pos             = np.where(P2b_arr > 1.0)[0]
         P2b_arr[pos]    = 1.0 
         ax1.plot(SMA_AU_plotarr, P2b_arr,       color = 'deeppink',      linestyle = ':', linewidth = 1.5, alpha = 0.75)

         #3-body merger:
         #ISO:
         P3b_arr         = 1.1*Nims*(Torb_a/tGW_SI)**(2./7.)
         pos             = np.where(P3b_arr > 1.0)[0]
         P3b_arr[pos]    = 1.0 
         ax1.plot(SMA_AU_plotarr, P3b_arr,       color = 'deepskyblue',      linestyle = ':', linewidth = 0.5, alpha = 0.75)
         #CP:
         P3b_arr         = 0.9*Nims*(Torb_a/tGW_SI)**(1./7.) #(1.05)*np.sqrt(2.)*Nims*(RsBH_SI/SMA_SI_plotarr)**(5./14.) 
         pos             = np.where(P3b_arr > 1.0)[0]
         P3b_arr[pos]    = 1.0 
         ax1.plot(SMA_AU_plotarr, P3b_arr,       color = 'deepskyblue',      linestyle = ':', linewidth = 1.5, alpha = 0.75)    

         #TEST: INTEGRATED PROB
         #ISO:
         d_bs                = 7./9.
         P2b_arr             = (tint_SI/tGW_SI)**(2./7.)
         P2b_INT_arr         = (1./(1.-d_bs))*(7./8.)*(P2b_arr - 0.0)
         pos                 = np.where(P2b_INT_arr > 1.0)[0]
         P2b_INT_arr[pos]    = 1.0
         SMA_AU_2bP1_ISO     = SMA_AU_plotarr[max(pos)]
         #ax1.plot(SMA_AU_plotarr, P2b_INT_arr,           label = r'INT 2b-ISO (CALC)',       color = dc_cl_arr[0],      linestyle = '-.', linewidth = 1.0, alpha = 0.75)
         #CP:
         d_bs                = 2./3.
         P2b_arr             = (tint_SI/tGW_SI)**(1./7.)
         P2b_INT_arr         = (1./(1.-d_bs))*(7./4.)*(P2b_arr - 0.0)
         pos                 = np.where(P2b_INT_arr > 1.0)[0]
         P2b_INT_arr[pos]    = 1.0
         SMA_AU_2bP1_CP      = SMA_AU_plotarr[max(pos)]
         #ax1.plot(SMA_AU_plotarr, P2b_INT_arr,           label = r'INT 2b-CPL (CALC)',       color = dc_cl_arr[1],      linestyle = '-.', linewidth = 1.0, alpha = 0.75)

         ax1.plot([0,0], [0,0],           label = r'analytical solutions',       color = 'black',      linestyle = ':', linewidth = 1.0, alpha = 1.0)

         #GUIDE LINES:
         #ax1.plot([1e-5,1e2], [1,1],       color = 'black',     linestyle = '-', linewidth = 3.0)
         #ax1.plot([0.0075,0.0075], [1e-3, 1.0],       color = 'grey',     linestyle = '-', linewidth = 1.0)
         #ax1.axvspan(1e-6, 0.0075, alpha=0.25, color='grey')
         #ax1.plot([SMA_AU_2bP1_ISO,SMA_AU_2bP1_ISO], [1e-3, 1.0],        color = dc_cl_arr[0],     linestyle = '-', linewidth = 1.0)
         #ax1.axvspan(SMA_AU_2bP1_ISO, SMA_AU_2bP1_CP, alpha=0.25, color=dc_cl_arr[0])
         #ax1.plot([SMA_AU_2bP1_CP,SMA_AU_2bP1_CP],   [1e-3, 1.0],        color = dc_cl_arr[1],     linestyle = '-', linewidth = 1.0)
         #ax1.axvspan(SMA_AU_2bP1_CP, 100, alpha=0.25, color=dc_cl_arr[1])

         #axis/plot settings:
         ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = True, ncol=1, framealpha = 1.0)
         ax1.text(0.1, 0.95, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 10, color = 'black')
         ax1.text(0.2, 0.89, r'$t_{\rm GW} < 10^{5}\ yrs.$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 10, color = 'black')
         ax1.set_xlabel(r'semi-major axis $a$ [AU]')
         ax1.set_ylabel(r'Probability')
         #ax1.set_xlim(0.75*min(SMA_arr_AU), 2.0*max(SMA_arr_AU))
         ax1.set_xlim((0.75*(1e-2)), 2.0*(1e3))
         ax1.set_ylim(1e-3, 1.2)

         ax1.set_facecolor('xkcd:white')

         plt.savefig(open_data_filename + '_Pfig_ISOCPL.pdf', bbox_inches='tight')
         plt.show()

         exit()
     
     
     
     
    
    
    
    #FIG 3:
    #ASSUMES EQUAL MASS CASE!
    mBH_SI      = m1_SI
    RsBH_SI     = 2.*G_new_SI*mBH_SI/c_SI**2.
    Nims        = 20.    
  
    #info_per_DATA_SMA_OC_arr[dc, ac, oc, 0]   = 1.*nr_id3_tlim_ij      #binary-single ejection
    #info_per_DATA_SMA_OC_arr[dc, ac, oc, 1]   = 1.*nr_id5_tlim_ij      #3-body merger
 
    EMC_nrid_arr    = np.add(np.add(info_per_DATA_SMA_OC_arr[:,:,0,:], info_per_DATA_SMA_OC_arr[:,:,1,:]), info_per_DATA_SMA_OC_arr[:,:,2,:])             
    
    #Figure 1:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    ax1.set_xscale('log')
    ax1.set_yscale('log')

    ax1.axvspan(1e-6, 0.01, alpha=0.25, color='grey', hatch="XX")
    ax1.axvspan(0.01, 1e6 , alpha=0.25, color='limegreen')
    ax1.axhspan(1, 1e6 , alpha=1, color='white')

    #DATA: (loop over dataset)
    dc_cl_arr   = ['salmon', 'black']  #ISO, CPL
    dc_LW_arr   = [1,2]
    dc_ALP_arr  = [0.5, 1.0]
    legnd3b_arr = [r'(D) 3b-IS', r'(D) 3b-CP']
    legnd2b_arr = [r'(D) 2b-IS', r'(D) 2b-CP']
    for dc in range(0,nrdataset):
        #3body:
        nr_3b       = EMC_nrid_arr[dc,:,1]
        nr_simOK    = info_per_DATA_SMA_arr[dc,:,1]
        P3b_arr     = nr_3b/nr_simOK
        P3b_arr_err = np.sqrt(nr_3b)/nr_simOK
        ax1.errorbar(SMA_arr_AU, P3b_arr,  yerr=P3b_arr_err, label = legnd3b_arr[dc], color = dc_cl_arr[dc],     marker = '.', markersize = 10, linestyle = '--', linewidth = dc_LW_arr[dc])
        #2body:
        nr_2b       = EMC_nrid_arr[dc,:,0]
        nr_simOK    = info_per_DATA_SMA_arr[dc,:,1]
        P2b_arr     = nr_2b/nr_simOK
        P2b_arr_err = np.sqrt(nr_2b)/nr_simOK
        ax1.errorbar(SMA_arr_AU, P2b_arr,  yerr=P2b_arr_err, label = legnd2b_arr[dc], color = dc_cl_arr[dc],     marker = '.', markersize = 10, linestyle = ':', linewidth = dc_LW_arr[dc])            
   
    #ANALYTICAL calc:
    SMA_AU_plotarr  = 10.**(np.linspace(-4.0, 4.0, 1000))
    SMA_SI_plotarr  = SMA_AU_plotarr*AU_SI
    tGW_SI          = (768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((SMA_SI_plotarr**4.)/(mBH_SI**3.))
    #2-body merger:
    #ISO:
    P2b_arr         = 1.3*(tint_SI/tGW_SI)**(2./7.)
    pos             = np.where(P2b_arr > 1.0)[0]
    P2b_arr[pos]    = 1.0 
    ax1.plot(SMA_AU_plotarr, P2b_arr,           label = r'(C) 2b-IS',       color = dc_cl_arr[0],      linestyle = ':', linewidth = 1.0, alpha = 0.5)
    #CP:
    P2b_arr         = 1.1*(tint_SI/tGW_SI)**(1./7.)
    pos             = np.where(P2b_arr > 1.0)[0]
    P2b_arr[pos]    = 1.0 
    ax1.plot(SMA_AU_plotarr, P2b_arr,           label = r'(C) 2b-CP',       color = dc_cl_arr[1],      linestyle = ':', linewidth = 1.0, alpha = 0.5)
    #3-body merger:
    #ISO:
    P3b_arr  = (1.7)*2.*Nims*(RsBH_SI/SMA_SI_plotarr)**(5./7.) 
    pos             = np.where(P3b_arr > 1.0)[0]
    P3b_arr[pos]    = 1.0 
    ax1.plot(SMA_AU_plotarr, P3b_arr,           label = r'(C) 3b-IS',       color = dc_cl_arr[0],      linestyle = '--', linewidth = 1.0, alpha = 0.5)
    #CP:
    P3b_arr  = (0.9)*np.sqrt(2.)*Nims*(RsBH_SI/SMA_SI_plotarr)**(5./14.) 
    pos             = np.where(P3b_arr > 1.0)[0]
    P3b_arr[pos]    = 1.0 
    ax1.plot(SMA_AU_plotarr, P3b_arr,           label = r'(C) 3b-CP',       color = dc_cl_arr[1],      linestyle = '--', linewidth = 1.0, alpha = 0.5)    
    
    #axis/plot settings:
    ax1.text(0.35, 0.95, r'[$20M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 9)
    ax1.text(0.35, 0.89, r'$e_0 = 0$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 9)
    ax1.plot([1e-5,1e4], [1,1],       color = 'black',     linestyle = '-', linewidth = 1)
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = True, ncol=1, framealpha = 1.0)
    ax1.set_xlabel(r'semi-major axis [AU]')
    ax1.set_ylabel(r'outcome probability')
    ax1.set_title(r'GW merger probability')
    #ax1.set_xlim(0.75*min(SMA_arr_AU), 2.0*max(SMA_arr_AU))
    ax1.set_xlim(0.75*(1e-3), 2.0*(1e3))
    ax1.set_ylim(1e-3, 1)

    plt.savefig(open_data_filename + '_Pfig_ISOCPL.pdf', bbox_inches='tight')
    plt.show()

    exit()
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    #FIG 4: NOTE: maybe make as contour, or mix of scatter and contour.
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    ax1.set_xscale('log')
    ax1.set_yscale('log')    
    #define data:
    ecc_fGW_arr         = savetestarr[dc, ac, :, 8] #ecc at 10 Hz
    pos                 = np.where(ecc_fGW_arr[:] < 0.0)[0]
    ecc_fGW_arr[pos]    = 1.0
    pos_ecc01_10Hz      = np.where(ecc_fGW_arr[:] > 0.1)[0]
    pos_ecc01_10Hz_tlim = list(set(pos_ecc01_10Hz).intersection(pos_tlim))
    nr_ecc01_10Hz_tlim  = 1.*len(pos_ecc01_10Hz_tlim)
    #relevant 1st/2nd merger dataset:
    pos_ecc01_10Hz_3b_tlim_2ndOK    = list(set(pos_ecc01_10Hz).intersection(pos_3b_tlim_2ndOK))
    fGW_1st = savetestarr[dc, ac, pos_ecc01_10Hz_3b_tlim_2ndOK, 4]
    fGW_2nd = savetestarr[dc, ac, pos_ecc01_10Hz_3b_tlim_2ndOK, 6]
    #boundaries and nrs:
    LISA_l  = 1e-3
    LISA_u  = 1e-1
    LIGO_l  = 10
    pos_2ndLISA = np.where((fGW_2nd > LISA_l) & (fGW_2nd < LISA_u))[0]
    pos_2ndLIGO = np.where((fGW_2nd > LIGO_l))[0] 
    print 1.*len(pos_2ndLISA)/nr_ecc01_10Hz_tlim
    print 1.*len(pos_2ndLIGO)/nr_ecc01_10Hz_tlim
    #plot:
    xp = fGW_1st
    yp = fGW_2nd
    ax1.plot(xp, yp, marker = '.', linestyle = '', color = 'teal')
    #guide lines:
    ax1.axhspan(LISA_l, LISA_u ,    alpha=0.25, color='red')
    ax1.axhspan(LIGO_l, 1e10,       alpha=0.25, color='grey')
    ax1.text(100, 1e-2,  r'2nd LISA $(\sim 10\%)$', horizontalalignment='left', verticalalignment='center', transform=ax1.transData, fontsize = 12, color = 'black')
    ax1.text(100, 1e2,   r'2nd LIGO $(\sim 1\%)$',  horizontalalignment='left', verticalalignment='center', transform=ax1.transData, fontsize = 12, color = 'black')
    ax1.text(0.05, 0.95, r'[$40M_{\odot}$, $20M_{\odot}$] $\leftarrow$ $20M_{\odot}$', horizontalalignment='left', verticalalignment='center', transform=ax1.transAxes, fontsize = 9)
    #limits/labels etc.
    #ax1.set_xlim(...)
    ax1.set_ylim(1e-5, 1e3)  
    ax1.set_xlabel(r'$f_{GW}$ ($e_{10} > 0.1$) [$1$st merger]')
    ax1.set_ylabel(r'$f_{GW}$ [$2$nd merger]')
    ax1.set_title(r'$f_{GW}$ at formation -- $1$st and $2$nd GW merger')
    #plot and save:
    plt.savefig(open_data_filename + 'fGW_1st2nd.pdf', bbox_inches='tight')
    plt.show()
    exit()
    
    
    
    
    
    
    
    
    
    
    
    
    
    #fGW plot 2:
    normfac         = 1.*len(pos_tlim)
    
    data_ij         = fGWij_SI_arr[pos_tlim]
    Hdata           = np.log10(data_ij)
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    CHy             = np.cumsum(Hy[::-1])[::-1]/normfac #backwards cumsum
    ax2.step(Hx, CHy, label = 'all', color = 'black', alpha=1)
    
    data_ij         = fGW12_SI_arr[pos_3b_tlim_2ndOK]
    Hdata           = np.log10(data_ij)
    Hy, bin_edges   = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx              = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    CHy             = np.cumsum(Hy[::-1])[::-1]/normfac #backwards cumsum
    ax2.step(Hx, CHy, label = 'all', color = 'black', alpha=1)
    

    #EVERYTHING SEEMS OK! BUT BELOW: find another way of changing to ISO dataset: we don't want mistakes
    #and mitchach in dc = ...
    #Now think about figure. 2nd merger: normalized to nr ecc in LIGO? looks better?
    #WE NEED: incl. real GW kick vel. AND: maybe calc ecc at some fgw.
    #make '3 sections' in paper: populaitons (cross section: 3body dominates, etc...), eccentricity (), double GW mergers ().
    
    
    #oplot dc = 0 case (ISO):
    dc_ISO              = 0
    fGWij_SI_arr_ISO    = savetestarr[dc_ISO, ac, :, 4] #1st GW merger
    pos_tlim_ISO        = np.where(savetestarr[dc_ISO, ac, :, 3] == 1)[0] 
    data_ij             = fGWij_SI_arr_ISO[pos_tlim_ISO]
    Hdata               = np.log10(data_ij)
    Hy, bin_edges       = np.histogram(Hdata, bins = nr_Hbins, range=[min_H, max_H], density=False)
    Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    CHy                 = 1-np.cumsum(Hy)/(1.*max(np.cumsum(Hy)))
    ax2.step(Hx, CHy, label = 'all', color = 'red', alpha=1)
    
    
    
    
    plt.savefig(open_data_filename + 'fGW_dist_a.pdf', bbox_inches='tight')
    plt.show()

    exit()






    fig = plt.figure(figsize=(7.0, 7.0))
    ax1 = fig.add_subplot(111)
    ac = 3
    pos = np.where(savetestarr[dc, ac, :,   5] == 1)[0]
    xp = np.log10(savetestarr[ dc, ac, pos, 4])
    yp = np.log10(savetestarr[ dc, ac, pos, 6])
    ax1.plot(xp, yp, marker = '.', linestyle = '')
    plt.show()





    #ecc dist:
    
    #Fig. 1:
    #choose SMA (index):
    ac  = 4
    #define fig:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    #Loop over 12, 13, 23:
    for oc in range(0,3): 
        H_X_data_ij         = HIST_data_arr[dc,ac,oc, 0,1,0,:]
        H_Y_data_ij         = HIST_data_arr[dc,ac,oc, 0,1,1,:]
        ax1.step(H_X_data_ij, H_Y_data_ij, label = label123_arr[oc], color = color123_arr[oc])

    ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = False, ncol=1, framealpha = 1.0)
    
    ax1.set_xlabel(r'ecc')
    ax1.set_ylabel(r'Number')
    
    plt.savefig(open_data_filename + 'ecc.pdf', bbox_inches='tight')
    plt.show()




    #Fig. 1:
    #fGW distribution at a given SMA:
    #choose SMA (index):
    ac  = 4
    #define fig:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    #Loop over 12, 13, 23:
    for oc in range(0,3):
        H_X_data_ij         = HIST_data_arr[dc,ac,oc, 0,0,0,:]
        H_Y_data_ij         = np.add(HIST_data_arr[dc,ac,oc, 0,0,1,:], HIST_data_arr[dc,ac,oc, 1,0,1,:])
        ax1.step(H_X_data_ij, H_Y_data_ij, label = label123_arr[oc], color = color123_arr[oc])

    ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = False, ncol=1, framealpha = 1.0)
    
    ax1.set_xlabel(r'fGW [Hz]')
    ax1.set_ylabel(r'Number')
    
    plt.savefig(open_data_filename + 'fGW_dist_a.pdf', bbox_inches='tight')
    plt.show()
    
    
    
    
    
    #Fig. 2:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    #define:
    nr_simOK    = info_per_DATA_SMA_arr[dc,:,1]
    #Loop over 12, 13, 23 (individual prob)
    for oc in range(0,3):
    
        nr_2b       = info_per_DATA_SMA_OC_arr[dc,:,oc,0]
        nr_3b       = info_per_DATA_SMA_OC_arr[dc,:,oc,1]

        P2b_arr     = nr_2b/nr_simOK
        P2b_arr_err = np.sqrt(nr_2b)/nr_simOK
        ax1.errorbar(SMA_arr_AU, P2b_arr,  yerr=P2b_arr_err, linestyle = '--', label = label123_2b_arr[oc], color = color123_arr[oc])

        P3b_arr     = nr_3b/nr_simOK
        P3b_arr_err = np.sqrt(nr_3b)/nr_simOK
        ax1.errorbar(SMA_arr_AU, P3b_arr,  yerr=P3b_arr_err, linestyle = '-', label = label123_3b_arr[oc], color = color123_arr[oc])        
    
    #sum over 12, 13, 23 (tot prob)
    TOT123_nrid_arr = np.add(np.add(info_per_DATA_SMA_OC_arr[:,:,0,:], info_per_DATA_SMA_OC_arr[:,:,1,:]), info_per_DATA_SMA_OC_arr[:,:,2,:])             
    TOT_nrid_arr    = np.add(TOT123_nrid_arr[:,:,0], TOT123_nrid_arr[:,:,1])
    
    ax1.plot(SMA_arr_AU, TOT123_nrid_arr[dc,:,0]/nr_simOK[:],   color = 'black', linestyle = '-',   label = '2b', linewidth = 2.0, alpha = 0.75)    
    ax1.plot(SMA_arr_AU, TOT123_nrid_arr[dc,:,1]/nr_simOK[:],   color = 'black', linestyle = '-',   label = '3b', linewidth = 2.0, alpha = 0.75)    
    ax1.plot(SMA_arr_AU, TOT_nrid_arr[dc,:]/nr_simOK[:],        color = 'black', linestyle = '--',  label = 'all', linewidth = 2.0, alpha = 0.75)    

    ax1.legend(loc='lower left', numpoints = 1, fontsize = 8, frameon = False, ncol=3, framealpha = 1.0)
            
    ax1.set_xscale('log')
    ax1.set_yscale('log')
    ax1.set_xlabel(r'semi-major axis [AU]')
    ax1.set_ylabel(r'outcome probability')

    plt.savefig(open_data_filename + 'Prob123_a.pdf', bbox_inches='tight')
    plt.show()
    
    exit()





    
    
    
    
    
    
    
    
    
    

#--------------------------------------------------------------
#--------------------------------------------------------------










#2nd BBH:
#NOTES: make vk an array and show several results.
testarrfGW  = np.zeros((nr_id5, 5), dtype=np.float64)
vkick_kmsec = 10.0
vkick_SI    = vkick_kmsec*1000.
rndnum_arr  = np.random.random((nr_id5,2))
rnd_phi     = ((2.*np.pi)*rndnum_arr[:,0])			    #P flat in phi			- pos on unitsphere
rnd_theta   = np.arccos(2.*rndnum_arr[:,1]-1.)          #P flat in cos(theta)	- pos on unitsphere
m_i         = 1.*mBH_SI
m_j         = 2.*mBH_SI
m_ij        = m_i+m_j
mu_ij       = m_i*m_j/m_ij
for tb5 in range(0,nr_id5):
    p5  = pos_id5[tb5]
    #before kick (no kick NK):
    NK_rel_pos_vec_SI   = RS_avs_output_Nbody_xtra_2_info_REAL[ac,0,p5,0:3]*R_sun_SI
    NK_rel_vel_vec_SI   = RS_avs_output_Nbody_xtra_2_info_REAL[ac,0,p5,3:6]*(1000./kmsec_U)
    #kick velocity vec (isotropic vkick dist):       
    phi                 = rnd_phi[tb5]
    theta               = rnd_theta[tb5]
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
    fGW_SI  = (1./np.pi)*np.sqrt(2.*G_new_SI*mBH_SI/(rp_SI**3.))
    testarrfGW[tb5, 0]  = fGW_SI



print 'dc, ac, nr_id3, nr_id5, nr_NAN: ', dc, ac, nr_id3, nr_id5, nr_NAN

#PLOT TEST: ECC DIST AT A GIVEN SMA(ac = ...)
if (ac == 3):
    #define fig window:
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    #outcome: 2-body merger
    pos_id      = pos_id3
    a_SI_arr    = RS_avs_output_Nbody_endstate_REAL[ac,0,pos_id,3]*R_sun_SI
    ecc_arr     = RS_avs_output_Nbody_endstate_REAL[ac,0,pos_id,4]
    rp_SI_arr   = a_SI_arr*(1.-ecc_arr)
    fGW_arr     = (1./np.pi)*np.sqrt(2.*G_new_SI*mBH_SI/(rp_SI_arr**3.))    #in Hz
    tGW_SI_arr  = ((768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((a_SI_arr**4.)/(mBH_SI**3.)))*((1.-(ecc_arr**2.))**(7./2.))
    tint_SI     = (0.1*(10.**6.))*sec_year
    pos_tGW     = np.where(tGW_SI_arr < tint_SI)[0]
    fGW_arr_C_tGW = fGW_arr[pos_tGW]   
    ax1.hist(np.log10(fGW_arr_C_tGW), bins = 50, range=[-3, 3], alpha = 0.5, label = '2-body merger')
    #outcome: 3-body merger
    pos_id      = pos_id5
    a_SI_arr    = RS_avs_output_Nbody_endstate_REAL[ac,0,pos_id,3]*R_sun_SI
    ecc_arr     = RS_avs_output_Nbody_endstate_REAL[ac,0,pos_id,4]
    rp_SI_arr   = a_SI_arr*(1.-ecc_arr)
    fGW_arr     = (1./np.pi)*np.sqrt(2.*G_new_SI*mBH_SI/(rp_SI_arr**3.))    #in Hz
    tGW_SI_arr  = ((768./425.)*((5.*(c_SI**5.))/(512.*(G_new_SI**3.)))*((a_SI_arr**4.)/(mBH_SI**3.)))*((1.-(ecc_arr**2.))**(7./2.))
    tint_SI     = (0.1*(10.**6.))*sec_year
    pos_tGW     = np.where(tGW_SI_arr < tint_SI)[0]
    fGW_arr_C_tGW = fGW_arr[pos_tGW]   
    ax1.hist(np.log10(fGW_arr_C_tGW), bins = 50, range=[-3, 3], alpha = 0.5, label = '3-body merger')


    #TEST:
    ax1.hist(np.log10(testarrfGW[:, 0]), bins = 50, range=[-6, 3], alpha = 0.5, label = '2nd merger')


    #plot settings, etc.:
    ax1.legend(loc='upper right', numpoints = 1, fontsize = 8, frameon = False, ncol=1, framealpha = 1.0)
    ax1.set_xlim(-6,3)
    ax1.set_xlabel(r'log $f_{GW}$ [Hz]')
    ax1.set_ylabel(r'Counts')
    #save/show:
    plt.savefig(data_name + '_fGW_dist.pdf', bbox_inches='tight')
    plt.show()            






  