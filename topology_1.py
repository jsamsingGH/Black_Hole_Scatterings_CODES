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
yr_sec     = 31536000.0    
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif',
        'size'   : 35}
mpl.rc('font', **font)


#data name:
#data_name = 'TOP_WGR_101010Msun_a1e2_v001vc_fbfull_gamma00pi_fbres1000'
#data_name = 'TOP_WGR_101010Msun_a1e4_v001vc_fbfull_gamma00pi_fbres500'
data_name = 'test_NO3'


M1 = 10.
M2 = 10.
M3 = 10.

#low res data: SOME OF THESE ARE NOT LABELED CORRECTLY!!

#'TOP_WGR_101010Msun_a1e4_v01vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e4_v01vc_fbfull_gamma14pi_fbres50'
#'TOP_WGR_101010Msun_a1e4_v01vc_fbfull_gamma12pi_fbres50'

#'TOP_WGR_101010Msun_a1e2_v01vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v01vc_fbfull_gamma14pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v01vc_fbfull_gamma12pi_fbres50'

#'TOP_WGR_101010Msun_a1e2_v001vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v01vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v025vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v05vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v075vc_fbfull_gamma00pi_fbres50'
#'TOP_WGR_101010Msun_a1e2_v10vc_fbfull_gamma00pi_fbres50'



#high res data:

#TOP_WGR_101010Msun_a1e2_v001vc_fgfull_bimpL0_fbres1000
#TOP_WGR_101010Msun_a1e2_v001vc_fgfull_bimp0_fbres1000

#'TOP_WGR_101010Msun_a1e0_v001vc_fbfull_gamma00pi_fbres500'
#'TOP_WGR_101010Msun_a1e0_v01vc_fbfull_gamma00pi_fbres500'
#'TOP_WGR_101010Msun_a1e0_v05vc_fbfull_gamma00pi_fbres500'

#'TOP_WGR_101010Msun_a1e4_v001vc_fbfull_gamma00pi_fbres500'
#'TOP_WGR_101010Msun_a1e4_v001vc_fbfull_gamma14pi_fbres500'
#'TOP_WGR_101010Msun_a1e4_v001vc_fbfull_gamma12pi_fbres500'

#'TOP_NGR_101010Msun_a1e4_v001vc_zoomR1_gamma00pi_fbres500'
#'TOP_WGR_101010Msun_a1e4_v001vc_zoomR1_gamma00pi_fbres500'

#'TOP_WGR_101010Msun_a1e4_v001vc_zoomf0pib350350_gamma00pi_fbres500'
#'TOP_NGR_101010Msun_a1e4_v001vc_zoomf0pib350350_gamma00pi_fbres500'

#'TOP_WGR_101010Msun_a1e0_v001vc_zoomR2_gamma00pi_fbres500'
#'TOP_WGR_101010Msun_a1e0_v001vc_zoomR4_gamma00pi_fbres500'



#'TOP_WGR_101010Msun_a1e2_v001vc_fbfull_gamma00pi_fbres1000'

#'TOP_test_laptop_WGR_pg14'#'TOP_noGR_TEST1zoom'#'TOP_test2_laptop_noGR'#'TOP_test_laptop_WGR_pg12'#'TOP_test2_laptop_noGR'#'TOP_test_laptop_WGR_pg14'#'TOP_test3_laptop_noGR'#'TOP_test2_laptop_noGR'

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
nr_tot  = MC_settings_list_INT[2]
res_pa  = int(np.sqrt(nr_tot))
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
xy_outlist_MC_info_REAL         = outlist_MC_info_REAL.reshape(res_pa, res_pa, 10)                  #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, r_val, b_val, f_val, 0.0]   
xy_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(res_pa, res_pa, 10)             #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
xy_output_Nbody_endstate_REAL   = output_Nbody_endstate_REAL.reshape(res_pa, res_pa, 10)            #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
xy_output_Nbody_xtra_info_REAL  = output_Nbody_xtra_info_REAL.reshape(res_pa, res_pa, 10)           #[rmin_12, rmin_13, rmin_23]
#----------------------------------------------------------


#----------------------------------------------------------
#Analysis:
#----------------------------------------------------------

#---------------------------
#analyze: 
#---------------------------
#allocate:
binssin_endstates_xy_arr            = np.zeros((res_pa, res_pa), dtype=int)
binssin_GW_COLL_endstates_xy_arr    = np.zeros((res_pa, res_pa), dtype=int)
binssin_info_xy_arr                 = np.zeros((res_pa, res_pa, 10), dtype=float)
#test_arr                            = np.zeros((nr_tot,2), dtype=float)

#counter:
ac = 0

#loop over imp b and phase f:
for xc in range(0, res_pa):                               
    for yc in range(0, res_pa):       
        
        #define:
        endstate_id         = xy_output_Nbody_endstate_INT[xc,yc,0]
        sin_k               = xy_output_Nbody_endstate_INT[xc,yc,3]     #dont work for inspirals. Works fine for other endstates.
        endbin_a            = xy_output_Nbody_endstate_REAL[xc,yc,3]
        endbin_e            = xy_output_Nbody_endstate_REAL[xc,yc,4]
        rmin_121323         = min([xy_output_Nbody_xtra_info_REAL[xc,yc,0], xy_output_Nbody_xtra_info_REAL[xc,yc,1], xy_output_Nbody_xtra_info_REAL[xc,yc,2]])
        nr_resonanses       = xy_output_Nbody_endstate_INT[xc,yc,6]
        sin_wrt_endbin_E    = xy_output_Nbody_endstate_REAL[xc,yc,7]
        sin_wrt_endbin_a    = xy_output_Nbody_endstate_REAL[xc,yc,8]
        sin_wrt_endbin_ecc  = xy_output_Nbody_endstate_REAL[xc,yc,9]
        sin_wrt_endbin_rp   = sin_wrt_endbin_a*(1.-sin_wrt_endbin_ecc)
        #calc second-merger insp time (EQUAL MASS CASE):
        m_i     = M_sun_SI*(M1)
        m_j     = M_sun_SI*(M1+M1)
        m_ij    = m_i+m_j
        mu_ij   = m_i*m_j/m_ij
        Rs      = (2.*G_new_SI*m_ij/(c_SI**2.))
        Ms      = mu_ij
        eps     = 85.*np.pi/96.
        beta    = 7./2.
        a_SI    = R_sun_SI*sin_wrt_endbin_a
        a0_SI   = R_sun_SI*sma_a0
        rp_SI   = a_SI*(1.-sin_wrt_endbin_ecc)
        Torb_a_secs         = 2.*np.pi*np.sqrt((a_SI**3.)/(G_new_SI*(m_i+m_j)))        
        Torb_a0_secs        = 2.*np.pi*np.sqrt((a0_SI**3.)/(G_new_SI*(m_i+m_i)))
        tinsp_ini_bin_secs  = (a0_SI**(4.))/(4.*(64./5.)*(G_new_SI**3.)*m_i*m_i*(m_i+m_i)/(c_SI**5.))       
        tinsp_secondM_secs  = (2.*np.pi*(rp_SI**(beta))*np.sqrt(a_SI)*(m_ij*mu_ij*(Rs**(1.-beta))/(eps*(Ms**(2.))*np.sqrt(G_new_SI*m_ij)))) + 0.5*Torb_a_secs 
        tinsp_secondM_years = tinsp_secondM_secs/yr_sec
        tinsp_secondM_unit_inibint  = tinsp_secondM_secs/tinsp_ini_bin_secs
        tinsp_secondM_unit_Torb0    = tinsp_secondM_secs/Torb_a0_secs
        print tinsp_secondM_unit_inibint, a_SI/a0_SI, tinsp_secondM_unit_Torb0
        ##TEST:
        #vel_vc          = np.sqrt((3./2.)*M1/sma_a0)            #EQUAL MASS LIMIT!
        #vel_v           = xy_outlist_MC_info_REAL[0,0,1]/vel_vc #EQUAL MASS LIMIT!
        #if (endstate_id == 5 or endstate_id == 2): test_arr[ac,0]  = xy_outlist_MC_info_REAL[xc,yc,7]/(sma_a0/vel_v)
        #if (endstate_id == 5 or endstate_id == 2): test_arr[ac,1] = (rp_SI/a0_SI)/((3./16.)*(np.sin((1./4.)*np.pi))**(2.))
        
        
        #endstate ids:
        if (endstate_id == 3 and sin_k == 1):
            binssin_endstates_xy_arr[xc,yc] = 1 #endstate: binsin   (sin k = 1)
        if (endstate_id == 3 and sin_k == 2):   
            binssin_endstates_xy_arr[xc,yc] = 2 #endstate: binsin   (sin k = 2)
        if (endstate_id == 3 and sin_k == 3):   
            binssin_endstates_xy_arr[xc,yc] = 3 #endstate: binsin   (sin k = 3)
        if (endstate_id == 2):                 
            binssin_endstates_xy_arr[xc,yc] = 4 #collision
        if (endstate_id == 5):                 
            binssin_endstates_xy_arr[xc,yc] = 5 #inspiral
        if (endstate_id >= 10):                
            binssin_endstates_xy_arr[xc,yc] = 6 #noid

        #endstate ids SPLITTING UP GW INSPIRALS AND COLLISIONS:
        out_bin_i         = xy_output_Nbody_endstate_INT[xc,yc,1]
        out_bin_j         = xy_output_Nbody_endstate_INT[xc,yc,2]        
        #GW inspirals  
        if (endstate_id == 5 and out_bin_i == 2 and out_bin_j == 3):
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 1 #endstate: GW   (sin k = 1)
        if (endstate_id == 5 and out_bin_i == 1 and out_bin_j == 3):
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 2 #endstate: GW   (sin k = 2)
        if (endstate_id == 5 and out_bin_i == 1 and out_bin_j == 2):
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 3 #endstate: GW   (sin k = 3)
        #GW inspirals        
        if (endstate_id == 2 and sin_k == 1):   
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 4 #endstate: COLL (sin k = 1)
        if (endstate_id == 2 and sin_k == 2):   
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 5 #endstate: COLL (sin k = 2)
        if (endstate_id == 2 and sin_k == 3):   
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 6 #endstate: COLL (sin k = 3)
        #BS:
        if (endstate_id == 3 or endstate_id >= 10):
            binssin_GW_COLL_endstates_xy_arr[xc,yc] = 0
            
        
        
        #IN THE HB LIMIT: IS A IONIZATION IS SEEN ITS LIKELY TO BE A COLLISION!!! NOT A PERFECT SOLUTION BUT OK FOR NOW.      
        if (endstate_id == 4):
            binssin_endstates_xy_arr[xc,yc] = 4
            
                             
              
        #param values a,e, etc:      
        if (endstate_id == 3):
            binssin_info_xy_arr[xc,yc,0]    = endbin_a
            binssin_info_xy_arr[xc,yc,1]    = endbin_e
        
        binssin_info_xy_arr[xc,yc,2]        = 10.**10.    #initilize
        if (endstate_id == 3 or endstate_id == 2 or endstate_id == 5):
            binssin_info_xy_arr[xc,yc,2]    = rmin_121323
            
    
        binssin_info_xy_arr[xc,yc,3]    = nr_resonanses

        #rp for single wrt coll/insp bin:
        binssin_info_xy_arr[xc,yc,4]        = 100.*sma_a0       #initialize value in regions outside coll/insp. rp in these regions are not defined (sin k in not bound): just put a large value.
        binssin_info_xy_arr[xc,yc,6]        = 1e10
        if (endstate_id == 2 or endstate_id == 5):              #replace initialize val with correct rp when endstate is coll/insp.
            if (a_SI > 0.0): binssin_info_xy_arr[xc,yc,4]    = sin_wrt_endbin_rp
            if (a_SI > 0.0): binssin_info_xy_arr[xc,yc,6]    = tinsp_secondM_unit_Torb0
            
        binssin_info_xy_arr[xc,yc,5]    = sin_wrt_endbin_E/(M1*M2/(2.*sma_a0))  #endstate bin-sin energy E_ijk in units of initial orbital energy.
        
        
        
        #TEST FOR ART-OF-SCIENCE!!!!!:
        if (endstate_id == 2 or endstate_id == 5):
            binssin_info_xy_arr[xc,yc,7]    = 1e-10
        if (endstate_id == 3):
            binssin_info_xy_arr[xc,yc,7]    = ((endbin_a/sma_a0)**(4.))*(1.-endbin_e**2)**(7./2.)
        if (endstate_id >= 10):
            binssin_info_xy_arr[xc,yc,7]    = 1.0

        
        
        #print and update counter:
        ac = ac+1
        print ac






##TEST:
#fig = plt.figure(figsize=(10, 8))
#ax  = plt.subplot(111)
##N-body res:
#ax.plot(test_arr[:,0], test_arr[:,1], marker='o',  linestyle='none', color='black')
#plt.show()
#exit()




#---------------------------
#PLOT -- OPT1: vary: f,b fix: g
#---------------------------
figs_seperate_yesno = 1

#define:
vel_vc              = np.sqrt((3./2.)*M1/sma_a0)            #EQUAL MASS LIMIT!
vel_v               = xy_outlist_MC_info_REAL[0,0,1]/vel_vc #EQUAL MASS LIMIT!
imp_b_unita0v_arr   = xy_outlist_MC_info_REAL[:,0,7]/(sma_a0/vel_v)
phase_f_unitrad_arr = xy_outlist_MC_info_REAL[0,:,8]

#min max f,b;
minf    = min(phase_f_unitrad_arr)
maxf    = max(phase_f_unitrad_arr)
minb    = min(imp_b_unita0v_arr)
maxb    = max(imp_b_unita0v_arr)
#make x,y grid:
X, Y    = np.meshgrid(phase_f_unitrad_arr, imp_b_unita0v_arr)

save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/TOPOLOGY_STUDY/'

#this fig is used if we plot all figs together:
if (figs_seperate_yesno == 0):  fig = plt.figure(figsize=(16, 10))

#binssin_info_xy_arr[xc,yc,0]    = endbin_a
#binssin_info_xy_arr[xc,yc,1]    = endbin_e
#binssin_info_xy_arr[xc,yc,2]    = rmin_121323
#binssin_info_xy_arr[xc,yc,3]    = nr_resonanses









#endstates:
#binssin_endstates_xy_arr[xc,yc] = 1 #endstate: binsin   (sin k = 1)
#binssin_endstates_xy_arr[xc,yc] = 2 #endstate: binsin   (sin k = 2)
#binssin_endstates_xy_arr[xc,yc] = 3 #endstate: binsin   (sin k = 3)
#binssin_endstates_xy_arr[xc,yc] = 4 #collision
#binssin_endstates_xy_arr[xc,yc] = 5 #inspiral
#binssin_endstates_xy_arr[xc,yc] = 6 #noid
if (figs_seperate_yesno == 0):  ax  = plt.subplot(236)
#if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(30, 25)) #use for large plot
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
cMap    = c.ListedColormap(['#ff7f00', '#377eb8', '#ffff99', 'black', 'crimson', '#4daf4a' ])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
Z       = binssin_endstates_xy_arr[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=1.0, vmax=6.0)
#color bar:
bounds = np.array([1,2,3,4,5,6,7])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .5)
cbar.ax.set_yticklabels([r'BS[23]', r'BS[13]', r'BS[12]', r'Coll', r'GW', r'VLint'], rotation=270, fontsize = 25)
#axis settings:
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Endstates')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'endstate_ids.eps')#, bbox_inches = 'tight')
    plt.show()

exit()

















#endstates:
#binssin_endstates_xy_arr[xc,yc] = 1 #endstate: binsin   (sin k = 1)
#binssin_endstates_xy_arr[xc,yc] = 2 #endstate: binsin   (sin k = 2)
#binssin_endstates_xy_arr[xc,yc] = 3 #endstate: binsin   (sin k = 3)
#binssin_endstates_xy_arr[xc,yc] = 4 #collision
#binssin_endstates_xy_arr[xc,yc] = 5 #inspiral
#binssin_endstates_xy_arr[xc,yc] = 6 #noid
if (figs_seperate_yesno == 0):  ax  = plt.subplot(236)
#if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(30, 25)) #use for large plot
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
cMap    = c.ListedColormap(['#ff7f00', '#377eb8', '#ffff99', 'black', 'crimson', '#4daf4a' ])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
Z       = binssin_endstates_xy_arr[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=1.0, vmax=6.0)
#color bar:
bounds = np.array([1,2,3,4,5,6,7])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .5)
cbar.ax.set_yticklabels([r'BS[23]', r'BS[13]', r'BS[12]', r'Coll', r'GW', r'VLint'], rotation=270, fontsize = 25)
#axis settings:
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Endstates')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'endstate_ids.eps')#, bbox_inches = 'tight')
    plt.show()

exit()





#nr IMS:
if (figs_seperate_yesno == 0):  ax  = plt.subplot(234)
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
maxnrIMS    = 6
cmap        = cm.get_cmap('Greys', maxnrIMS)
Z       = binssin_info_xy_arr[:,:,3]
quad    = plt.pcolormesh(X, Y, Z, cmap=cmap, linewidth=0, rasterized=True, vmin=0, vmax=maxnrIMS-1)
#color bar:
bounds = np.array([0,1,2,3,4,5,6])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cmap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .5)
cbar.ax.set_yticklabels([r'0', r'1', r'2', r'3', r'4', r'$\geq$5'])
#axis settings:
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Number of intermediate states')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'nrRES.eps')#, bbox_inches='tight')
    plt.show()
    
exit()










#endstates:
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 0 #endstate: BS     
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 1 #endstate: GW     (sin k = 1)
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 2 #endstate: GW     (sin k = 2)
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 3 #endstate: GW     (sin k = 3)
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 4 #endstate: COLL   (sin k = 1)
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 5 #endstate: COLL   (sin k = 2)
#binssin_GW_COLL_endstates_xy_arr[xc,yc] = 6 #endstate: COLL   (sin k = 3)
if (figs_seperate_yesno == 0):  ax  = plt.subplot(236)
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
cMap    = c.ListedColormap(["white", "#CC1821", '#FC6B4B', '#FCBBA1', "#525353", "#979793", "#D9D9DA"])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
Z       = binssin_GW_COLL_endstates_xy_arr[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=0.0, vmax=6.0)
#color bar:
bounds = np.array([0,1,2,3,4,5,6,7])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals)
cbar.set_ticks(vals + .5)
cbar.ax.set_yticklabels([r'BS[ij]', r'GW[23]', r'GW[13]', r'GW[12]', r'Coll[23]', r'Coll[13]', r'Coll[12]'], rotation=270, fontsize = 25)
#axis settings:
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Collisions and GW Inspirals')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'GW_COLL_endstates.eps')#, bbox_inches='tight')
    plt.show()

exit()










if (figs_seperate_yesno == 0):  ax  = plt.subplot(231)
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
Z       = binssin_info_xy_arr[:,:,5]
quad    = plt.pcolormesh(X, Y, Z, cmap='RdBu', linewidth=0, rasterized=True, vmin=-1.0, vmax=1.0)
plt.colorbar()
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Binary-single energy')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'endEijk.eps')#, bbox_inches='tight')
    plt.show()

exit()





if (figs_seperate_yesno == 0):  ax  = plt.subplot(232)
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
Z       = binssin_info_xy_arr[:,:,1]
quad    = plt.pcolormesh(X, Y, Z, cmap='viridis', linewidth=0, rasterized=True, vmin=0.0, vmax=1.0)
plt.colorbar()
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'eccentricity')
if (figs_seperate_yesno == 1):
    #save and show fig:
    plt.savefig(save_data_folder+data_name+'_' + 'endecc.eps', bbox_inches='tight')
    plt.show()





if (figs_seperate_yesno == 0):  ax  = plt.subplot(235)
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
Z       = np.log10(binssin_info_xy_arr[:,:,6])
quad    = plt.pcolormesh(X, Y, Z, cmap='magma_r', linewidth=0, rasterized=True, vmin=0.0, vmax=5.0)
plt.colorbar()
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Double GW Mergers')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'secondmerger.eps')#, bbox_inches='tight')
    plt.show()

exit()






if (figs_seperate_yesno == 0):  ax  = plt.subplot(233)
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(15, 12.5))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
Z       = np.log10(binssin_info_xy_arr[:,:,2]/sma_a0)
quad    = plt.pcolormesh(X, Y, Z, cmap='viridis', linewidth=0, rasterized=True, vmin=-5.0, vmax=0.0)
plt.colorbar()
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
ax.set_xlabel(r'binary phase $f$')
ax.set_ylabel(r'impact parameter $b$')
ax.set_title(r'Minimum distance')
if (figs_seperate_yesno == 1):
    #save and show fig:
    fig.tight_layout()
    plt.savefig(save_data_folder+data_name+'_' + 'rmin121323.eps')#, bbox_inches='tight')
    plt.show()

exit()






plt.show()

exit()
#---------------------------







#FOR ART-OF-SCIENCE!!!
if (figs_seperate_yesno == 1):  fig = plt.figure(figsize=(50,50))
if (figs_seperate_yesno == 1):  ax  = plt.subplot(111)
Z       = - np.log10(binssin_info_xy_arr[:,:,7])
print np.log10(binssin_info_xy_arr[:,:,7])
quad    = plt.pcolormesh(X, Y, Z, cmap='magma', linewidth=0, rasterized=True, vmin=0.0, vmax=10.0)
#plt.colorbar()
ax.set_xlim(minf,maxf)
ax.set_ylim(minb,maxb)
#ax.set_xlabel(r'$f$ [radians]')
#ax.set_ylabel(r'$b$ [$a_{0}/v$]')
#ax.set_title(r'Binary-Single Energy')
#final axis settings:
ax.tick_params(
    axis='both',           # changes apply to the x,y-axis
    which='both',       # both major and minor ticks are affected
    bottom='on',        # ticks along the bottom edge are off
    top='off',           # ticks along the top edge are off
    labelbottom='off',  # labels along the bottom edge are off
    right='off',
    left='off',
    labelleft='off')
if (figs_seperate_yesno == 1):
    #save and show fig:
    plt.savefig(save_data_folder+data_name+'_' + 'test.eps', bbox_inches='tight')
    plt.show()

exit()














#---------------------------
#PLOT -- OPT2: vary: f,g fix: b
#---------------------------
#define:
phase_f_unitrad_arr = xy_outlist_MC_info_REAL[0,:,8] - np.pi
phase_g_unitrad_arr = xy_outlist_MC_info_REAL[:,0,9]

#make x,y grid:
X, Y    = np.meshgrid(phase_f_unitrad_arr, phase_g_unitrad_arr)

save_data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/TOPOLOGY_STUDY/'

#binssin_info_xy_arr[xc,yc,0]    = endbin_a
#binssin_info_xy_arr[xc,yc,1]    = endbin_e
#binssin_info_xy_arr[xc,yc,2]    = rmin_121323
#binssin_info_xy_arr[xc,yc,3]    = nr_resonanses

#endstates:
#binssin_endstates_xy_arr[xc,yc] = 1 #endstate: binsin   (sin k = 1)
#binssin_endstates_xy_arr[xc,yc] = 2 #endstate: binsin   (sin k = 2)
#binssin_endstates_xy_arr[xc,yc] = 3 #endstate: binsin   (sin k = 3)
#binssin_endstates_xy_arr[xc,yc] = 4 #collision
#binssin_endstates_xy_arr[xc,yc] = 5 #inspiral
#binssin_endstates_xy_arr[xc,yc] = 6 #noid
fig = plt.figure(figsize=(30, 25))
ax  = plt.subplot(111, projection = 'mollweide')
cMap    = c.ListedColormap(['#ff7f00', '#377eb8', '#ffff99', 'black', 'crimson', '#4daf4a' ])  #give the endstates colors (I dont know exactly how this works, but it does. I used an exmaple from the web.)
Z       = binssin_endstates_xy_arr[:,:]
quad    = plt.pcolormesh(X, Y, Z, cmap=cMap, linewidth=0, rasterized=True, vmin=1.0, vmax=6.0)

#color bar:
bounds = np.array([1,2,3,4,5,6,7])
vals = bounds[:-1]
cbar = fig.colorbar(quad, cmap=cMap, boundaries=bounds, values=vals, orientation='horizontal')
cbar.set_ticks(vals + .5)
cbar.ax.set_xticklabels([r'BS[23]', r'BS[13]', r'BS[12]', r'Coll', r'GWinsp', r'VLint'], fontsize = 25)#, rotation=270, fontsize = 25)

#axis settings:
ax.set_xticklabels([])
ax.set_yticklabels([])
ax.set_title(r'Endstates')
#save and show fig:
plt.savefig(save_data_folder+data_name+'_' + 'endstate_ids.eps', bbox_inches='tight')
plt.show()

exit()
#---------------------------










