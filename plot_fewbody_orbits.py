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
from matplotlib.patches import Rectangle
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition, mark_inset)




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
tU_tSI     = np.sqrt(R_sun_SI**3./(G_new_SI*M_sun_SI))    
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'

#---------------------------------------
#input names of Nbody IC files:
#---------------------------------------
#file_name    = 'fig_paper_eccBH.txt'
#file_name    = 'fig_paper_eccBH_NOGR.txt'
#file_name    = 'AGN10.txt'
#file_name    = 'testfGW10.txt'  #50, 20, 20, 1-3
#file_name    = 'testfGW12.txt'  #60, 60, 20, 1-2

#file_name    = 'testfGW16.txt'   #30, 30, 8 OK 1e-2 (insp 2,3)
#file_name    = 'testfGW20.txt'   #30, 30, 8 OK 1e-1 (insp )


#BEST:
#file_name    = 'testfGW16.txt'#'test.txt'   #30, 30, 8 
#file_name    = 'testfGW_20.txt'#'test.txt'   #30, 30, 8 1e-3, insp 1,3 
#file_name    = 'testfGW_21.txt'#'test.txt'   #30, 30, 8 1e-1, insp 1,2 
#file_name    = 'testfGW_22.txt'#'test.txt'   #30, 30, 8 1e-1, insp 2,3 
#file_name    = 'testfGW_23.txt'#'test.txt'   #30, 30, 8 1e-1, insp 2,3 

#file_name    = 'testfGW_25.txt'#'test.txt'   #10, 10, 3 1e-1, insp 1,3 
file_name     = 'testfGW_25_NOGR.txt'#'test.txt'   #10, 10, 3 1e-1, insp 1,3, NOGR!


#file_name    = 'testfGW_26.txt'#'test.txt'   #20, 20, 20 1e-1, insp 1,2 
#file_name    = 'testfGW_27.txt'#'test.txt'   #20, 20, 10 1e-1, insp 2,3 
#file_name    = 'testfGW_29.txt'#'test.txt'   #20, 20, 8 1e-1, insp 2,3 
#file_name    = 'testfGW_30.txt'#'test.txt'   #15, 15, 5 1e-1, insp 1,3 
#file_name    = 'testfGW_31.txt'#'test.txt'   #10, 10, 5 1e-1, insp ??

#2 GW sources (first burst in one band, then in a second...) LOOK AT IT!!!
#file_name    = 'testfGW_28.txt'#'test.txt'   #40, 40, 40 1e-1, insp 2,3 


#testfGW10010011e1.txt quite amazing example of 100+100+1 in near coplanar interaction.



#file_name    = 'MCinput_Nbody.txt'


#NOTES:
#difficult to have perfect 2-body macth template. better fit for c0?
#change observer position to optimize romer.
#can we plot or do better than accumulated romer?
#can we plot angle instead of dt?
#


#---------------------------------------
QUICK_look_YN       = 1
ANALYZE_etc_YN      = 0
SIM_2Body_YN        = 0
#---------------------------------------
m1  = 10.#30.
m2  = 10.#30.
m3  = 3.#8.

obj_i   = 1
obj_j   = 3

obj_k   = 6-(obj_i+obj_j)
m123_arr    = [m1, m2, m3]
mi          = m123_arr[obj_i-1]
mj          = m123_arr[obj_j-1]
mk          = m123_arr[obj_k-1]
#---------------------------------------


#---------------------------------------------------------------------
if (QUICK_look_YN == 1):
#---------------------------------------------------------------------
    #simulate if necessary:
    subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + file_name, shell=True)
    #open data:
    tf = open('NbodyTides_dataout_pos.dat',  "r")
    NbodyTides_dataout_pos          = np.loadtxt(tf, dtype=float)
    tf.close()
    b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
    b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
    b3_posxyz   = NbodyTides_dataout_pos[:,6:9]
    
    #PLOT: Orbits xyz:
    penr    = 10
    fig = plt.figure(figsize=(6, 4))
    ax = fig.add_subplot(111)

    #3-body:
    xp = b1_posxyz[:,0]
    yp = b1_posxyz[:,1]
    ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.5, color = 'pink', rasterized=False)

    xp = b2_posxyz[:,0]
    yp = b2_posxyz[:,1]
    ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'grey', rasterized=False)

    xp = b3_posxyz[:,0]
    yp = b3_posxyz[:,1]
    ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'cadetblue', rasterized=False)

    ax.set_xlim(-100, 100)
    ax.set_ylim(-200, 200)

    ax.set_xlabel(r'x-pos')
    ax.set_ylabel(r'y-pos')

    plt.show()
    exit()
#---------------------------------------------------------------------
    







def func_abin_ebin(pos_1, pos_2, vel_1, vel_2, mass_1, mass_2):
    pos			= pos_1 - pos_2
    vel			= vel_1 - vel_2 		
    Mtot		= mass_1+mass_2
    mred		= mass_1*mass_2/Mtot
    Lvec		= mred*np.cross(pos,vel)
    len_L		= np.linalg.norm(Lvec)			
    len_r 		= np.linalg.norm(pos)
    len_v		= np.linalg.norm(vel)
    z_axis_vec  = np.array([0,0,1])		
    #calc orbital params wrt bin CM:
    E_kin	= (1./2.)*mred*(len_v**2.)
    E_pot	= - Mtot*mred/len_r
    E_tot	= E_kin + E_pot
    a_bin	= - Mtot*mred/(2.*E_tot)
    e_bin	= np.sqrt(1. + (2.*E_tot*(len_L**2.))/(mred*((Mtot*mred)**2.)))
    Lz_ang  = np.arccos(np.dot(Lvec, z_axis_vec)/(np.linalg.norm(Lvec)*np.linalg.norm(z_axis_vec)))
    return [a_bin,e_bin,Lz_ang]; 
    

def GW_aet(v, t, mi, mj):
    fm1_SI  = mi*M_sun_SI
    fm2_SI  = mj*M_sun_SI
    fa      = v[0]
    fe      = v[1] 
    dadt    = ((-64./5.)*(G_new_SI**3.)*fm1_SI*fm2_SI*(fm1_SI+fm2_SI)/((c_SI**5.)*(fa**3.)*((1.-fe**2.)**(7./2.))))*(1.+(73./24.)*fe**2. + (37./96.)*fe**4.)
    dedt    = ((-304./15.)*fe*(G_new_SI**3.)*fm1_SI*fm2_SI*(fm1_SI+fm2_SI)/((c_SI**5.)*(fa**4.)*((1.-fe**2.)**(5./2.))))*(1.+(121./304.)*fe**2.)    
    return [dadt, dedt]
    
    
#---------------------------------------------------------------------
if (ANALYZE_etc_YN == 1):
#---------------------------------------------------------------------
    #---------------------------------------
    #run Nbody code with that input file:
    #---------------------------------------
    subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + file_name, shell=True)
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
    #calc info:
    sma_a0      = np.sqrt(sum((b1_posxyz[0,:]-b2_posxyz[0,:])**2.))
    print sma_a0/AU_U
    #open time info:
    tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
    NbodyTides_dataout_a1a2a3       = np.loadtxt(tf, dtype=float)
    tf.close()
    time = NbodyTides_dataout_a1a2a3[:,0]
    max_time        = 1.0*max(time)
    #give files unique name:
    subprocess.call('cp ' + 'Nbody_endsim_info_data.txt ' + file_name + '_Nbody_endsim_info_data.txt', shell=True)
    #---------------------------------------

    #-----------------------------------------------------------------
    #a,e vs Peters
    #-----------------------------------------------------------------
    nrs         = len(NbodyTides_dataout_pos[:,0])
    output_arr  = np.zeros((nrs, 15), dtype=np.float64) 

    for nc in range(0,nrs):
        
        print nc, nrs
        #define:
        pxyz_i  = NbodyTides_dataout_pos[nc, 0+3*(obj_i-1):3+3*(obj_i-1)]
        vxyz_i  = NbodyTides_dataout_vel[nc, 0+3*(obj_i-1):3+3*(obj_i-1)]
        pxyz_j  = NbodyTides_dataout_pos[nc, 0+3*(obj_j-1):3+3*(obj_j-1)]
        vxyz_j  = NbodyTides_dataout_vel[nc, 0+3*(obj_j-1):3+3*(obj_j-1)]
        time_t  = time[nc]
        #calc:
        out_func    = func_abin_ebin(pxyz_i, pxyz_j, vxyz_i, vxyz_j, mi, mj)    #return [a_bin,e_bin,Lz_ang];
        aij_eij     = out_func[0:2]
        Lzangij     = out_func[2]
        rij_SI      = np.linalg.norm(pxyz_i-pxyz_j)*R_sun_SI
        fGWij       = (1./np.pi)*np.sqrt((G_new_SI*(mi+mj)*M_sun_SI)/(rij_SI**3.))
        #save:
        output_arr[nc,0]    = aij_eij[0]
        output_arr[nc,1]    = aij_eij[1]
        output_arr[nc,2]    = fGWij
        output_arr[nc,3]    = time_t
        output_arr[nc,4]    = np.sqrt(sum((pxyz_i-pxyz_j)**2.))
        output_arr[nc,5]    = Lzangij
        #angular mumentum test:
        L1_vec  = m1*np.cross(b1_posxyz[nc,:],b1_velxyz[nc,:])
        L2_vec  = m2*np.cross(b2_posxyz[nc,:],b2_velxyz[nc,:])
        L3_vec  = m3*np.cross(b3_posxyz[nc,:],b3_velxyz[nc,:])
        Lx_vec  = L1_vec[0] + L2_vec[0] + L3_vec[0]
        Ly_vec  = L1_vec[1] + L2_vec[1] + L3_vec[1]
        Lz_vec  = L1_vec[2] + L2_vec[2] + L3_vec[2]
        Lxy_vec = np.linalg.norm(np.array([Lx_vec, Ly_vec]))        
        Lt_vec  = np.linalg.norm(np.array([Lx_vec, Ly_vec, Lz_vec])) 
        Lij_vec = np.linalg.norm((mi*mj/(mi+mj))*np.cross(pxyz_i-pxyz_j,vxyz_i-vxyz_j))
        output_arr[nc,6]    = Lx_vec        
        output_arr[nc,7]    = Ly_vec        
        output_arr[nc,8]    = Lz_vec
        output_arr[nc,9]    = Lxy_vec                
        output_arr[nc,10]   = Lt_vec        
        output_arr[nc,11]   = Lij_vec        
        
    #save output file:
    data_save_posvel    = [b1_posxyz,b2_posxyz,b3_posxyz, b1_velxyz,b2_velxyz,b3_velxyz]
    data_save_time      = [time]
    data_save_outputarr = [output_arr]
    
    np.savez(file_name + '_data_save_posvel',       data_save_posvel)
    np.savez(file_name + '_data_save_time',         data_save_time)
    np.savez(file_name + '_data_save_outputarr',    data_save_outputarr)
#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#Open and define data:
#---------------------------------------------------------------------
#open:
data_save_posvel        = np.load(file_name + '_data_save_posvel' + '.npz')['arr_0']
data_save_time          = np.load(file_name + '_data_save_time' + '.npz')['arr_0']
data_save_outputarr     = np.load(file_name + '_data_save_outputarr' + '.npz')['arr_0']
Nbody_endsim_info_data  = open(file_name + '_Nbody_endsim_info_data.txt',  "r")
#define:
[b1_posxyz,b2_posxyz,b3_posxyz, b1_velxyz,b2_velxyz,b3_velxyz]  = data_save_posvel
[time]                                                          = data_save_time
[output_arr]                                                    = data_save_outputarr
sma_a0      = np.sqrt(sum((b1_posxyz[0,:]-b2_posxyz[0,:])**2.))
max_time    = 1.0*max(time)
nrs         = len(b1_posxyz[:,0])
Lines = Nbody_endsim_info_data.readlines()
endsim_Return_Info_arr_REAL = [float(x) for x in Lines[2].split()]
endsim_out_xtra_2_info_REAL = [float(x) for x in Lines[3].split()]
endsim_out_xtra_3_info_REAL = [float(x) for x in Lines[4].split()]
endsim_out_xtra_info_REAL   = [float(x) for x in Lines[5].split()]
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#Define:
#---------------------------------------------------------------------
#from sim data:
data_a_arr  = output_arr[:,0] 
data_e_arr  = output_arr[:,1] 
data_f_arr  = output_arr[:,2] 
data_t_arr  = output_arr[:,3] 
data_r_arr  = output_arr[:,4] 
data_Lz_arr = output_arr[:,5] 
data_Lx     = output_arr[:,6] 
data_Ly     = output_arr[:,7] 
data_Lz     = output_arr[:,8] 
data_Lxy    = output_arr[:,9] 
data_Lt     = output_arr[:,10] 
data_Lij    = output_arr[:,11] 

posxyz_i    = data_save_posvel[obj_i-1]
posxyz_j    = data_save_posvel[obj_j-1]
posxyz_k    = data_save_posvel[obj_k-1]

velxyz_i    = data_save_posvel[obj_i+3-1]
velxyz_j    = data_save_posvel[obj_j+3-1]

#From: Nbody_endsim_info_data.txt
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
a_ims_0     = endsim_out_xtra_info_REAL[3] 
e_ims_0     = endsim_out_xtra_info_REAL[4] 
r_ims_k_0   = endsim_out_xtra_2_info_REAL[9] 
a_ims_k     = endsim_Return_Info_arr_REAL[8]
rp_ims_0    = a_ims_0*(1.-e_ims_0)
rp_ims_0_SI = rp_ims_0*R_sun_SI
mi_SI       = mi*M_sun_SI
mj_SI       = mj*M_sun_SI
mk_SI       = mk*M_sun_SI
mij_SI      = mi_SI + mj_SI
mijk_SI     = mi_SI + mj_SI + mk_SI
a_ims_k_SI  = a_ims_k*R_sun_SI
time_ims_0  = endsim_out_xtra_3_info_REAL[6] 
time_end_1  = endsim_out_xtra_info_REAL[8] 
v_ij_ims_0  = endsim_out_xtra_3_info_REAL[3:6] 
v_ij_end_1  = endsim_out_xtra_3_info_REAL[7:10]
p_ij_ims_0  = endsim_out_xtra_3_info_REAL[0:3]
p_ij_end_1  = endsim_out_xtra_info_REAL[5:8]
pos_ims_0   = np.where(time[:] == time_ims_0)[0]
dt_ims0_end1    = time_end_1-time_ims_0
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#Reference Properties for 3B merger:
#---------------------------------------------------------------------
a_URij_lim  = 50.0 #This is our reference distance between ij in units of their total Schwarzschild Radius Rij.  

Rij_SI      = 2.*G_new_SI*(mi+mj)*M_sun_SI/(c_SI**2.)
a_arr_SI    = data_a_arr*R_sun_SI
pos_lim     = int(max(np.where(a_arr_SI/Rij_SI > a_URij_lim)[0])) # a reference distance at which we take the "3-body merger" to be well approximated by an isolated 2-body merger.

#peak finder:
pos_3B_r_peaks, _   = find_peaks(data_r_arr)                                                     
t_3B_peakvals       = data_t_arr[pos_3B_r_peaks]                                                
r_3B_peakvals       = data_r_arr[pos_3B_r_peaks]
a_3B_peakvals       = data_a_arr[pos_3B_r_peaks]
e_3B_peakvals       = data_e_arr[pos_3B_r_peaks]

#define pos_ref:
pos_3B_p0           = int(min(np.where(t_3B_peakvals > data_t_arr[pos_lim])[0])) #index _p0 refers to the pos of the first peak following pos_lim. This pos is RELATIVE TO pos_3B_r_peaks
pos_ref             = pos_3B_r_peaks[pos_3B_p0]
#---------------------------------------------------------------------
#Analysis:
#---------------------------------------------------------------------
#time merger vs time orbit:
t_ims_GW    = (768./425.)*(5./256.)*(2.**(7./2.))*(((6.*np.sqrt(2.))/(85.*np.pi))**(1./2.))*(c_SI**(15./2.)/(G_new_SI**(17./4.)))*((rp_ims_0_SI**(21./4.))/((mi_SI**(3./2.))*(mj_SI**(3./2.))*(mij_SI**(5./4.))))
T_orb_k     = 2.*np.pi*np.sqrt((a_ims_k_SI**3.)/(G_new_SI*(mijk_SI)))
print t_ims_GW/T_orb_k, 't_ims_GW/T_orb_k'

#curvature angle:
angle_v01   = np.arccos(np.dot(v_ij_ims_0,v_ij_end_1)/(np.linalg.norm(v_ij_ims_0)*np.linalg.norm(v_ij_end_1)))
print angle_v01*(360./(2.*np.pi)), 'angle_v01 in degree'

#calc dt_ROMER (max val):
dist_lin_ims0end1   = np.linalg.norm(v_ij_end_1)*dt_ims0_end1
pxyz_lin_ims0       = p_ij_end_1 - dist_lin_ims0end1*(v_ij_end_1/np.linalg.norm(v_ij_end_1))
pxyz_3bd_ims0       = p_ij_ims_0
dist_ims0_lin_3bd   = np.linalg.norm(pxyz_lin_ims0 - pxyz_3bd_ims0)
dt_ims0_lin_3bd_SI  = (dist_ims0_lin_3bd*R_sun_SI)/c_SI
print dt_ims0_lin_3bd_SI, 'dt_ims0_lin_3bd_SI'

#find time between 1 and 2 burst after IMS:
pos_ims_1b          = int(min(np.where(t_3B_peakvals[:] > time_ims_0)[0]))
pos_ims_2b          = pos_ims_1b+1
dt_ims_1b2b         = np.abs(t_3B_peakvals[pos_ims_1b] - t_3B_peakvals[pos_ims_2b])
dt_ims_1b2b_SI      = dt_ims_1b2b*tU_tSI
print dt_ims_1b2b_SI, 'dt_ims_1b2b_SI'

#find effect from tides:
#[b1_posxyz,b2_posxyz,b3_posxyz, b1_velxyz,b2_velxyz,b3_velxyz]  = data_save_posvel
posxyz_CMij = (posxyz_i*mi + posxyz_j*mj)/(mi+mj)
dist_CMij_k = np.array([np.linalg.norm(posxyz_CMij[x,:] - posxyz_k[x,:]) for x in range(0,nrs)])
F_tid_ij_k  = (2.*((mi+mj)*mk)/(dist_CMij_k**3.))*(data_a_arr*(1.+data_e_arr)) 
F_bin_ij    = mi*mj/((data_a_arr*(1.+data_e_arr))**2.)
F_tid_bin   = F_tid_ij_k/F_bin_ij

#angle between vel vectors:
velxyz_CMij = (velxyz_i*mi + velxyz_j*mj)/(mi+mj)
angle_vij_vijend1 = np.array([np.arccos(np.dot(velxyz_CMij[x,:], v_ij_end_1)/(np.linalg.norm(velxyz_CMij[x,:])*np.linalg.norm(v_ij_end_1))) for x in range(0,nrs)])
#---------------------------------------------------------------------


#---------------------------------------------------------------------
#2-body (single-single capture) path:
#---------------------------------------------------------------------
if (SIM_2Body_YN == 1):
#---------------------------------------------------------------------
    #Simulate backwards from pos_ref:
    ini_posxyz_i_2B =  posxyz_i[pos_ref,:]
    ini_posxyz_j_2B =  posxyz_j[pos_ref,:]
    ini_velxyz_i_2B = -velxyz_i[pos_ref,:]
    ini_velxyz_j_2B = -velxyz_j[pos_ref,:]
    sim_time_2B     =  dt_ims0_end1
    #---------------------------------------
    #Write param file to Nbody code:
    #---------------------------------------
    #IC Obj 1:  (binary 1 = i)
    #---------------------------------------
    b1_mass      = mi
    b1_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
    b1_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
    b1_const_arr = np.array([b1_mass, Rsch_1Msun_unitRsun*b1_mass, 1.5, 5./3., 1, 0, 1, 0,0,0], dtype='d')
    #---------------------------------------
    #IC Obj 2:  (binary 2 = j)
    #---------------------------------------
    b2_mass      = mj
    b2_q         = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
    b2_qdot      = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
    b2_const_arr = np.array([b2_mass, Rsch_1Msun_unitRsun*b2_mass, 1.5, 5./3., 1, 0, 1, 0,0,0], dtype='d')
    #---------------------------------------    
    #Nbody code settings:
    #---------------------------------------
    #important parameters:
    insp_threshold          = 10.0  # = a*(Ri+Rj).
    PN25sign_param          = -1    # switch sign on 2.5PN term for backwards evolution!!!! (-1 = switch (backwards), = 0 for normal (forward))
    # [use_12PN, use_25PN, Identify_3Body_endstate, max_sim_nrsteps, nfac_DF, outputinfo_screenfiles, nr_nbody_particles, de_analysis, PN25sign_param, ...]
    nbody_params_arr_1_INT  = np.array([0, 1, 0, 1000000000, 4, 1, 0, 0, PN25sign_param,0], dtype='i')
    # [scale_dt, max_sim_time, evolvetides_threshold, ENDbinsingle_threshold, max_simtime_sec, IMSbinsingle_threshold, tidaldisrup_threshold, insp_threshold, insp_SMA, ...]
    nbody_params_arr_2_REAL = np.array([0.01, sim_time_2B, -1, 0.01, 1e10, 0.1, 5., insp_threshold,0,0], dtype='d')
    #---------------------------------------
    #define:
    b1_posxyz_CM = ini_posxyz_i_2B
    b1_velxyz_CM = ini_velxyz_i_2B
    b2_posxyz_CM = ini_posxyz_j_2B
    b2_velxyz_CM = ini_velxyz_j_2B
    #open file:
    fn = 'MCinput_Nbody.txt'
    text_file = open(fn, "w")
    #nr particles in the sim   
    text_file.write('2' + '\n')
    #nbody code settings: 
    np.savetxt(text_file, nbody_params_arr_1_INT[None],   fmt='%10i')
    np.savetxt(text_file, nbody_params_arr_2_REAL[None],  fmt='%10.10f')
    #body 1:
    np.savetxt(text_file, b1_const_arr[None],  fmt='%10.10f')
    np.savetxt(text_file, b1_posxyz_CM[None],  fmt='%10.10f')
    np.savetxt(text_file, b1_velxyz_CM[None],  fmt='%10.10f')
    np.savetxt(text_file, b1_q,  fmt='%10f')
    np.savetxt(text_file, b1_qdot,  fmt='%10f')
    #body 2:
    np.savetxt(text_file, b2_const_arr[None],  fmt='%10.10f')
    np.savetxt(text_file, b2_posxyz_CM[None],  fmt='%10.10f')
    np.savetxt(text_file, b2_velxyz_CM[None],  fmt='%10.10f')
    np.savetxt(text_file, b2_q,  fmt='%10.10f')
    np.savetxt(text_file, b2_qdot,  fmt='%10.10f')
    #close file:
    text_file.close()    
    #---------------------------------------
    #Run Nbody_AffineTides_solver.exe with 'MCinput_Nbody.txt' as input:
    #---------------------------------------
    subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + fn, shell=True)
    #---------------------------------------
    #Read in data from Nbody solver:
    #---------------------------------------
    #open data:
    tf = open('NbodyTides_dataout_pos.dat',  "r")
    NbodyTides_dataout_pos          = np.loadtxt(tf, dtype=float)
    tf.close()
    posxyz_i_2B = NbodyTides_dataout_pos[:,0:3]
    posxyz_j_2B = NbodyTides_dataout_pos[:,3:6]
    #open data:
    tf = open('NbodyTides_dataout_vel.dat',  "r")
    NbodyTides_dataout_vel          = np.loadtxt(tf, dtype=float)
    tf.close()
    velxyz_i_2B = NbodyTides_dataout_vel[:,0:3]
    velxyz_j_2B = NbodyTides_dataout_vel[:,3:6]
    #open time info:
    tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
    NbodyTides_dataout_a1a2a3       = np.loadtxt(tf, dtype=float)
    tf.close()
    time_2B     = NbodyTides_dataout_a1a2a3[:,0]
    
    posxyz_i_2B = posxyz_i_2B[::-1,:]
    posxyz_j_2B = posxyz_j_2B[::-1,:]
    velxyz_i_2B = velxyz_i_2B[::-1,:]
    velxyz_j_2B = velxyz_j_2B[::-1,:]
    time_2B     = time_2B[::-1]
    #---------------------------------------------------------------------
    #---------------------------------------
    #define:
    #---------------------------------------
    nrs             = len(time_2B[:])
    output_arr_2B   = np.zeros((nrs, 5), dtype=np.float64) 
    for nc in range(0,nrs):
        print nc, nrs
        pxyz_i  = posxyz_i_2B[nc,:]
        pxyz_j  = posxyz_j_2B[nc,:]
        vxyz_i  = velxyz_i_2B[nc,:]
        vxyz_j  = velxyz_j_2B[nc,:]
        out_func    = func_abin_ebin(pxyz_i, pxyz_j, vxyz_i, vxyz_j, mi, mj)
        aij_eij     = out_func[0:2] 
        rij         = np.linalg.norm(pxyz_i-pxyz_j)
        rij_SI      = rij*R_sun_SI
        fGWij       = (1./np.pi)*np.sqrt((G_new_SI*(mi+mj)*M_sun_SI)/(rij_SI**3.))
        #save info:
        output_arr_2B[nc,0] = rij
        output_arr_2B[nc,1] = aij_eij[0]    #aij
        output_arr_2B[nc,2] = aij_eij[1]    #eij
        output_arr_2B[nc,3] = fGWij
    #---------------------------------------
    #save output data:
    #---------------------------------------
    data_save_posvel_2B = [posxyz_i_2B,posxyz_j_2B, velxyz_i_2B,velxyz_j_2B]
    np.savez(file_name + '_data_save_2B_posvel',    data_save_posvel_2B)
    np.savez(file_name + '_data_save_2B_time',      time_2B)
    np.savez(file_name + '_data_save_2B_outputarr', output_arr_2B)
    #---------------------------------------
#---------------------------------------------------------------------
#load data and define:
#---------------------------------------------------------------------
#open:
data_save_posvel_2B = np.load(file_name + '_data_save_2B_posvel' + '.npz')['arr_0']
time_2B             = np.load(file_name + '_data_save_2B_time' + '.npz')['arr_0']
output_arr_2B       = np.load(file_name + '_data_save_2B_outputarr' + '.npz')['arr_0'] 
#define:
[posxyz_i_2B,posxyz_j_2B, velxyz_i_2B,velxyz_j_2B]  = data_save_posvel_2B
data_2B_t_arr   = time_2B
data_2B_r_arr   = output_arr_2B[:,0]
data_2B_a_arr   = output_arr_2B[:,1]
data_2B_e_arr   = output_arr_2B[:,2]
data_2B_f_arr   = output_arr_2B[:,3]
#---------------------------------------------------------------------
#---------------------------------------------------------------------


#---------------------------------------------------------------------
#Compare 2-body and 3-body:
#---------------------------------------------------------------------
#time relative to pos_ref (t = 0 at pos_ref where the 2-body ss runs from)
t_RF_3B_R0      = data_t_arr - data_t_arr[pos_ref]
t_RF_2B_R0      = -data_2B_t_arr
t_SI_RF_3B_R0   = t_RF_3B_R0*tU_tSI
t_SI_RF_2B_R0   = t_RF_2B_R0*tU_tSI
#---------------------------------------------------------------------
#ROMER DELAY (RD):
#---------------------------------------------------------------------
#simulation is in 3B COM,
#---------------------------------------
#3-body
#---------------------------------------
#pos, vel of ij COM:
pos_ijCOM           = (posxyz_i*mi + posxyz_j*mj)/(mi+mj) 
vel_ijCOM           = (velxyz_i*mi + velxyz_j*mj)/(mi+mj)
#Position of observer (sometimes we indicate frame with `3C' (short for 3bodyCOM) and sometimes not. Anyway, in this section we always assume 3C frame if nothing else is stated.):
obs_norm_vec        = pos_ijCOM[pos_ref,:]/np.linalg.norm(pos_ijCOM[pos_ref,:])
obs_dist_Ua0        = 10000000.
obs_pos_vec         = (obs_dist_Ua0*sma_a0)*obs_norm_vec    #position vector of observer.
#dist from ijCOM to obs:
dist_ij_obs_arr     = np.zeros(nrs, dtype=np.float64)       #dist from ij to obs for all steps 
for nc in range(0,nrs):
    dist_obs_ij     = np.linalg.norm(pos_ijCOM[nc,:]-obs_pos_vec[:])
    dist_ij_obs_arr[nc] = dist_obs_ij
#convert dist to time:
dist_ij_obs_arr_SI  = dist_ij_obs_arr*R_sun_SI
time_ij_obs_arr_SI  = dist_ij_obs_arr_SI/(c_SI) #time (secs) it takes for light to travel from ijCOM to observer at each step 'nc'.
#---------------------------------------
#---------------------------------------
#2-body
#---------------------------------------
#pos, vel of 2body ij COM:
pos_ijCOM_2B        = (posxyz_i_2B*mi + posxyz_j_2B*mj)/(mi+mj) 
vel_ijCOM_2B        = (velxyz_i_2B*mi + velxyz_j_2B*mj)/(mi+mj)
#dist from ijCOM_2B to obs:
nrs_2B              = len(t_RF_2B_R0[:])
dist_2B_ij_obs_arr  = np.zeros(nrs_2B, dtype=np.float64)       #dist from ij to obs for all steps 
for nc in range(0,nrs_2B):
    dist_2B_obs_ij          = np.linalg.norm(pos_ijCOM_2B[nc,:]-obs_pos_vec[:])
    dist_2B_ij_obs_arr[nc]  = dist_2B_obs_ij
#convert dist to time:
dist_2B_ij_obs_arr_SI  = dist_2B_ij_obs_arr*R_sun_SI
time_2B_ij_obs_arr_SI  = dist_2B_ij_obs_arr_SI/(c_SI) #time (secs) it takes for light to travel from ijCOM to observer at each step 'nc'.
#---------------------------------------
#time shifts:
#---------------------------------------
t_SI_obs_3B     = t_SI_RF_3B_R0 + time_ij_obs_arr_SI
t_SI_obs_2B     = t_SI_RF_2B_R0 + time_2B_ij_obs_arr_SI
t_SI_obs_R0     = t_SI_obs_3B[pos_ref]
t_SI_obs_3B_R0  = t_SI_obs_3B - t_SI_obs_R0
t_SI_obs_2B_R0  = t_SI_obs_2B - t_SI_obs_R0
#---------------------------------------------------------------------
#f_GW peak timings:
#---------------------------------------------------------------------
pos_2B_f_peaks, _       = find_peaks(data_2B_f_arr)
pos_3B_f_peaks, _       = find_peaks(data_f_arr)
nr_peaks_PA             = len(pos_2B_f_peaks)
#finding peak analysis set (PAset):
pos_f0_3B               = max(np.where(pos_3B_f_peaks < pos_ref)[0])
pos_2B_f_peaks_PAset    = pos_2B_f_peaks                                                          
pos_3B_f_peaks_PAset    = pos_3B_f_peaks[pos_f0_3B+1-nr_peaks_PA:pos_f0_3B+1]                                                          
print 'nr_peaks_PA', nr_peaks_PA
#define:
t_SI_RF_2B_R0_PAset     = t_SI_RF_2B_R0[pos_2B_f_peaks_PAset]
t_SI_obs_2B_R0_PAset    = t_SI_obs_2B_R0[pos_2B_f_peaks_PAset]
t_SI_RF_3B_R0_PAset     = t_SI_RF_3B_R0[pos_3B_f_peaks_PAset]
t_SI_obs_3B_R0_PAset    = t_SI_obs_3B_R0[pos_3B_f_peaks_PAset]
#---------------------------------------------------------------------
#2B and 3B Time difference:
#---------------------------------------------------------------------
#Define linear Doppler factor (line of sight relative to obs) at ref point:
vel_ijCOM_ref   = vel_ijCOM[pos_ref]
v_LOS_SI        = np.dot(obs_norm_vec, vel_ijCOM_ref)*np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)
Dfac            = (1.-v_LOS_SI/c_SI)
#define:
dt_SI_obs_3B2B_PAset    = t_SI_obs_3B_R0_PAset      - t_SI_obs_2B_R0_PAset      #total difference
dt_SI_RD_3B2B_PAset     = t_SI_obs_3B_R0_PAset      - Dfac*t_SI_RF_3B_R0_PAset  #contribution from Romer Delay (RD)
dt_SI_tides_3B2B_PAset  = Dfac*t_SI_RF_3B_R0_PAset  - t_SI_obs_2B_R0_PAset      #contribution from tidal foces (tides)

#peak analysis relative time (PART):
deltaT_burst_SI_obs_3B_R0_PARTset   = t_SI_obs_3B_R0_PAset[1:nr_peaks_PA] - t_SI_obs_3B_R0_PAset[0:nr_peaks_PA-1]
t_SI_obs_3B_R0_PARTset              = t_SI_obs_3B_R0_PAset[0:nr_peaks_PA-1] 
frac_dt_obs_3B2B_PARTset            = dt_SI_obs_3B2B_PAset[0:nr_peaks_PA-1]/deltaT_burst_SI_obs_3B_R0_PARTset
#---------------------------------------------------------------------
#---------------------------------------------------------------------




#COMMENTS:
#LOOK AT THE x,y figs: the plane of the binary changes very much in the 3D case, WHY??? Loook at the 3D plot case
#and figure it out...
#calculate effect from tides.
#... Lorenz Swick for tides?
#check that the backwards sim is ok. Can we prove it? Test it at least with different initial starting points...




#---------------------------------------------------------------------
#PLOTS:
#---------------------------------------------------------------------



#PLOT: Orbits xyz:
penr    = 2
fig = plt.figure(figsize=(6, 4))
ax = fig.add_subplot(111)

posf    = np.where(data_f_arr[:] > 10)[0]

#3-body:
xp = posxyz_k[:,0]
yp = posxyz_k[:,1]
ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.5, color = 'pink', rasterized=False)
ax.plot(xp[pos_ims_0[0]:], yp[pos_ims_0[0]:], linewidth = 2.0, linestyle='-', marker = '', alpha = 1, color = 'pink', rasterized=False)
ax.plot(xp[pos_ims_0[0]],  yp[pos_ims_0[0]], marker = 'o', markersize = 8, color = 'pink', rasterized=False)

xp = posxyz_i[:,0]
yp = posxyz_i[:,1]
ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'grey', rasterized=False)

xp = posxyz_j[:,0]
yp = posxyz_j[:,1]
ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'cadetblue', rasterized=False)

xp = pos_ijCOM[:,0]
yp = pos_ijCOM[:,1]
ax.plot(xp[pos_ims_0[0]:], yp[pos_ims_0[0]:], linewidth = 2.0, linestyle='--', marker = '', alpha = 1, color = 'black', rasterized=False)
ax.plot(xp[pos_ims_0[0]],  yp[pos_ims_0[0]], marker = 'o', markersize = 8, color = 'black', rasterized=False)

#2-body:
xp = posxyz_i_2B[:,0]
yp = posxyz_i_2B[:,1]
ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.5, color = 'coral', rasterized=False)

xp = posxyz_j_2B[:,0]
yp = posxyz_j_2B[:,1]
ax.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.5, color = 'coral', rasterized=False)

xp = pos_ijCOM_2B[:,0]
yp = pos_ijCOM_2B[:,1]
ax.plot(xp[0:], yp[0:], linewidth = 2.0, linestyle='--', marker = '', alpha = 1, color = 'coral', rasterized=False)
ax.plot(xp[0],  yp[0], marker = 'o', markersize = 8, color = 'coral', rasterized=False)

#ax.set_xlim(-10*sma_a0, 10*sma_a0)
#ax.set_ylim(-10*sma_a0, 10*sma_a0)
ax.set_xlim(-100, 100)
ax.set_ylim(-200, 200)

ax.set_xlabel(r'x-pos')
ax.set_ylabel(r'y-pos')

plt.savefig(file_name + '_ijk_orbits_xyz.pdf', bbox_inches='tight')
plt.show()



#PLOT:
ptmin       = -150000/1000.
ptmax       = 1000/1000.

fig = plt.figure(figsize=(4, 8))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)

#Fig 1:
xp = t_SI_obs_2B_R0/1000.
yp = data_2B_f_arr
ax1.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'coral', rasterized=False)
ax1.plot(xp[pos_2B_f_peaks_PAset], yp[pos_2B_f_peaks_PAset], linewidth = 1.0, linestyle='', marker = 'o', markersize = 5, alpha = 1, color = 'coral', rasterized=False)

xp = t_SI_obs_3B_R0/1000.
yp = data_f_arr
ax1.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.75, color = 'black', rasterized=False)
ax1.plot(xp[pos_3B_f_peaks_PAset], yp[pos_3B_f_peaks_PAset], linewidth = 1.0, linestyle='', marker = 'o', markersize = 2, alpha = 0.75, color = 'black', rasterized=False)

#labels, limits.. etc:
ax1.legend(loc='upper left', numpoints = 1, fontsize = 10, frameon = False, ncol=1, framealpha = 1.0)
#ax1.set_xlabel(r'time [sec]')
ax1.set_ylabel(r'$f_{GW}$ [Hz]')
ax1.set_xlim(ptmin, ptmax)
ax1.set_ylim(-1., 15.)


#Fig 2:
#total difference between obs 3B and 2B:
xp  = t_SI_obs_3B_R0_PAset/1000.
yp  = dt_SI_obs_3B2B_PAset
ax2.plot(xp, yp, marker = 'o', markersize = 8, linestyle='', alpha = 0.75, color='c', label = r'$\delta t$ Obs.')

#contribution from 'Doppler':
xp  = t_SI_obs_3B_R0_PAset/1000.
yp  = dt_SI_RD_3B2B_PAset
ax2.plot(xp, yp, marker = 'o', markersize = 3, linestyle='', alpha = 0.75, color='black', label = r'$\delta t$ Doppler')

#contribution from 'tides':
xp  = t_SI_obs_3B_R0_PAset/1000.
yp  = dt_SI_tides_3B2B_PAset
ax2.plot(xp, yp, marker = '+', markersize = 6, linestyle='', alpha = 0.75, color='black', label = r'$\delta t$ Tides')

#labels, limits.. etc:
#ax2.set_xlabel(r'time [sec]')
ax2.set_ylabel(r'$\delta t$ [sec]')
ax2.set_xlim(ptmin, ptmax)
ax2.set_ylim(-3000., 150.)

#fig 3:
xp  = t_SI_obs_3B_R0_PARTset/1000.
yp  = np.log10(abs(frac_dt_obs_3B2B_PARTset))
ax3.plot(xp, yp, marker = 'o', markersize = 8, linestyle='', alpha = 0.75, color='orchid', label = r'$\delta t$ Obs.')

#labels, limits.. etc:
ax3.set_xlabel(r'time [sec]')
ax3.set_ylabel(r'$log\ \delta t/T_{orb}$')
ax3.set_xlim(ptmin, ptmax)
ax3.set_ylim(-4., 1.)

plt.savefig(file_name + '_dt_peak_analysis_1.pdf', bbox_inches='tight')
plt.show()




#PLOT:

fig = plt.figure(figsize=(10, 7))
ax1 = fig.add_subplot(411)
ax2 = fig.add_subplot(412)
ax3 = fig.add_subplot(413)
ax4 = fig.add_subplot(414)



#Fig
xp = t_SI_RF_3B_R0
yp = F_tid_bin
ax1.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'black', rasterized=False)
ax1.set_xlim(-100000., ptmax)
ax1.set_ylim(-1.0, 1.0)


#Fig
xp = data_e_arr
yp = F_tid_bin
ax2.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'black', rasterized=False)
ax2.set_ylim(-1.0, 1.0)


#Fig
xp = data_e_arr[pos_3B_f_peaks_PAset]
yp = dt_SI_tides_3B2B_PAset
ax3.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'black', rasterized=False)


#Fig 4:
xp = data_e_arr
yp = np.log10(data_a_arr)
ax4.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'black', rasterized=False)

xp = data_2B_e_arr
yp = np.log10(data_2B_a_arr)
ax4.plot(xp[:], yp[:], linewidth = 1.0, linestyle='-', marker = '', alpha = 1, color = 'coral', rasterized=False)




#FtidFbin

#a

#e

#dt




plt.show()





exit()
#---------------------------------------------------------------------



#NOTE: assumption about observer. all ok? where do we use obs pos?



















































#compare to Peters:
a0_SI       = data_a_arr[pos_ref]*R_sun_SI 
e0_SI       = data_e_arr[pos_ref] 
t0_SI       = data_t_arr[pos_ref]*tU_tSI
#integrate FORWARDs:
int_SI_t0   = t0_SI
t_eval_SI   = list(output_arr[pos_ref::, 3]*tU_tSI)
sol         = odeint(GW_aet, [a0_SI, e0_SI], t_eval_SI,  args = (mi, mj))
tFW = t_eval_SI
sFW = sol 
#integrate BACKWARDs:
int_SI_t0   = t0_SI
t_eval_SI   = list(reversed(output_arr[0:pos_ref+1, 3]*tU_tSI))
sol         = odeint(GW_aet, [a0_SI, e0_SI], t_eval_SI,  args = (mi, mj))
tBW = t_eval_SI
sBW = sol
#make final arrays: (AND CHANGE UNITS)
sol_t_arr   = np.append(tFW[0::],tBW[1::])/tU_tSI 
sol_a_arr   = np.append(sFW[0::,0],sBW[1::,0])/R_sun_SI
sol_e_arr   = np.append(sFW[0::,1],sBW[1::,1])
#sort
pos_t_sort  = np.argsort(sol_t_arr)
#overwrite:
sol_t_arr   = sol_t_arr[pos_t_sort]
sol_a_arr   = sol_a_arr[pos_t_sort]
sol_e_arr   = sol_e_arr[pos_t_sort]


#calc:
sol_Torb_arr    = 2.*np.pi*np.sqrt((sol_a_arr**3.)/(mi+mj))
data_Torb_arr   = 2.*np.pi*np.sqrt((data_a_arr**3.)/(mi+mj))

#define:
fGW_lim     = 1.0
posf        = np.where(data_f_arr > fGW_lim)[0]
posf10      = np.where(data_f_arr > 10.0)[0]







#PLOT:

fig     = plt.figure(figsize=(10, 8))
ax1     = fig.add_subplot(521)
ax2     = fig.add_subplot(522)
ax3     = fig.add_subplot(523)
ax4     = fig.add_subplot(524)
ax5     = fig.add_subplot(525)
ax6     = fig.add_subplot(526)
ax7     = fig.add_subplot(527)
ax8     = fig.add_subplot(528)
ax9     = fig.add_subplot(529)
ax10    = fig.add_subplot(5,2,10)



#plot ax1:
xp  = data_t_arr - data_t_arr[pos_ref]
yp  = np.log10(data_a_arr)
ax1.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

#xp  = (data_t_arr[posf])
#yp  = np.log10(data_a_arr[posf])
#ax1.plot(xp, yp, markersize=1, linestyle='', marker = 'o', alpha = 0.1, color='red')

xp = sol_t_arr - data_t_arr[pos_ref]
yp = np.log10(sol_a_arr)
ax1.plot(xp, yp, markersize=1, linestyle='-', marker = '', alpha = 0.75, color = 'green')

xp = t_3B_peakvals[pos_3B_p0] - data_t_arr[pos_ref]
yp = np.log10(a_3B_peakvals[pos_3B_p0])
ax1.plot(xp, yp, markersize=1, linestyle='-', marker = 'o', alpha = 0.75, color = 'red')


pos = np.where(data_a_arr[:] > 0.0)[0]
xp  = data_t_arr - data_t_arr[pos_ref]
yp  = np.log10(data_r_arr)
ax1.plot(xp[pos], yp[pos], linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

xp = t_3B_peakvals[pos_3B_p0] - data_t_arr[pos_ref]
yp = np.log10(r_3B_peakvals[pos_3B_p0])
ax1.plot(xp, yp, markersize=1, linestyle='-', marker = 'o', alpha = 0.75, color = 'red')


#plot ax2:
xp  = (data_t_arr)
yp  = 1.-data_e_arr**2
ax2.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

#xp  = (data_t_arr[posf])
#yp  = data_e_arr[posf]
#ax2.plot(xp, yp, markersize=1, linestyle='', marker = 'o', alpha = 0.1, color='red')

xp = (sol_t_arr)
yp = 1.-sol_e_arr**2.
ax2.plot(xp, yp, markersize=1, linestyle='-', marker = '', alpha = 0.75, color = 'green')



#plot ax3:
xp  = (data_t_arr - max_time)*tU_tSI
yp  = (data_f_arr)
ax3.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

#xp  = (data_t_arr[posf]) - max_time)*tU_tSI
#yp  = (data_f_arr[posf])
#ax3.plot(xp, yp, markersize=1, linestyle='', marker = 'o', alpha = 0.1, color='red')



#plot ax4:
xp  = (max_time - data_t_arr)*tU_tSI
yp  = (sol_Torb_arr/data_Torb_arr)
ax4.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

xp  = (max_time - data_t_arr[posf])*tU_tSI
yp  = (sol_Torb_arr[posf]/data_Torb_arr[posf])
ax4.plot(xp, yp, markersize=1, linestyle='', marker = 'o', alpha = 0.1, color='red')

xp  = (max_time - data_t_arr[posf10])*tU_tSI
yp  = (sol_Torb_arr[posf10]/data_Torb_arr[posf10])
ax4.plot(xp, yp, markersize=1, linestyle='', marker = 'x', alpha = 0.1, color='blue')

xp  = (max_time - data_t_arr)*tU_tSI
yp  = (xp*0.0 + 1.)
ax4.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')



#plot ax5:
xp  = np.log10(data_a_arr)
yp  = data_e_arr
ax5.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

xp  = np.log10(data_a_arr[posf])
yp  = data_e_arr[posf]
ax5.plot(xp, yp, markersize=1, linestyle='', marker = 'o', alpha = 0.1, color='red')

xp  = np.log10(sol_a_arr)
yp  = sol_e_arr
ax5.plot(xp, yp, markersize=1, linestyle='-', marker = '', alpha = 0.75, color = 'green')




#plot ax6:
pos = np.where(data_a_arr[:] > 0.0)[0]
xp  = data_t_arr[pos]
yp  = (data_Lz_arr[pos]/(2.*np.pi))*360.
ax6.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='', alpha = 0.75, color='black')



#plot ax7:
pos = np.where(data_a_arr[:] > 0.0)[0]

xp  = data_t_arr
yp  = data_Lx
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

xp  = data_t_arr
yp  = data_Ly
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

xp  = data_t_arr
yp  = data_Lz
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='-', alpha = 0.75, color='black')

xp  = data_t_arr
yp  = data_Lxy
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='-', alpha = 0.75, color='pink')

xp  = data_t_arr
yp  = data_Lt
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=1.5, linestyle='', alpha = 0.75, color='red')

xp  = data_t_arr[pos]
yp  = data_Lij[pos]
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='', alpha = 0.75, color='blue')

data_Lij_xy = data_Lij*np.sin(data_Lz_arr)
xp  = data_t_arr[pos]
yp  = data_Lij_xy[pos]
ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='', alpha = 0.75, color='orange')

data_Lij_z = data_Lij*np.cos(data_Lz_arr)
xp  = data_t_arr[pos]
yp  = data_Lij_z[pos]
#ax7.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='', alpha = 0.75, color='yellow')





#plot ax8:








#plot ax9:
xp  = (data_t_arr - max(data_t_arr))*tU_tSI
yp  = F_tid_bin
ax9.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='', alpha = 0.75, color='black')
ax9.set_ylim(-0.1, 1)


#plot ax10:
xp  = (data_t_arr - max(data_t_arr))*tU_tSI
yp  = np.log10(dist_CMij_k)
ax10.plot(xp, yp, marker = '.', markersize = 0.1, linewidth=0.5, linestyle='', alpha = 0.75, color='black')


plt.show()

#exit()










#---------------------------------------------------------------------
#---------------------------------------------------------------------

#---------------------------------------------------------------------
#FINAL ARRAYS: 2-body and 3-body:
#---------------------------------------------------------------------

#NOTES:
#- we need to match the 2/3 body curves a bit better at R0 (zoom in and look...)
#- how does result depend on e0 for the 2-body solver?
#- check all and make consistent notation, e.g. R0 etc..
#- polish up code...
#-over how many and what peaks should be find: np.mean(Nts_3B_p0 - Nts_2B_p0) ???
#THE PEAK FINDER DOES IT ALWAYS FIND BOTTOM?? what about top??

#REMEMBER:
#- at the moment we can only resolve "burst timing effects". For example, f_GW peak is not correct as it has to be dobbler shiftet. Can we do that? 
#- first paper can be on using "puls timing", but we still need to know if we can "see" the pulses. For that we need a rough estimate for the GW peak f. For a simple estimate, just take the rest frame and dobbler shift it. Should be ok!

#PROVIDE HERE FINAL ARRAYS AND EXPLANATION!!!!
#LAST THING: way to quantify difference and MAKE FIGS FOR ERC!!!
#THEN CLEAN UP AND PREPARE FOR ERC!

#- find a bit more accurate val of a0 and e0
#PROBLEMS:
#1)
#mis match between 2body even when vel is perpendicualr to line of sight. WHY?? Does not depend on: obs_dist_Ua0
#therefore, it has to be something else. Maybe below: t_SI_obs_R0     = t_SI_obs_3B[pos_ref]....
#FIND OUT! I can see that ALL is fine when I put v = 0... but when its not its not mathcing why... I KNOW: v is NOT PERPENDICULAR!!! ONLY TRUE FOR CIRCULAR ORBITS. BUT WE HAVE ECCENTRIC! WHAT I SEE IS JUST NORMAL DOBBLER! CHECK!!!

#2)
#names above has to fit better together...


#---------------------------------------------------------------------




#---------------------------------------------------------------------
#PLOT:
#---------------------------------------------------------------------
ptmin       =  t_SI_obs_3B_R0[pos_ims_0]
ptmax       =  100



fig = plt.figure(figsize=(5, 7))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

#Observed (all effects...)
xp  = t_SI_obs_2B_R0
yp  = data_2B_f_arr
ax1.plot(xp, yp, linewidth=2.0, linestyle='-', alpha = 0.50, color='coral', label = '2-body Obs.')

xp  = t_SI_obs_3B_R0
yp  = data_f_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.85, color='black', label = '3-body Obs.')
ax1.plot(xp[pos_ref], yp[pos_ref], marker = 'o', markersize = 5, alpha = 0.85, color='red')


#labels etc:
ax1.legend(loc='upper left', numpoints = 1, fontsize = 10, frameon = False, ncol=1, framealpha = 1.0)
ax1.set_xlabel(r'time [sec]')
ax1.set_ylabel(r'$f_{GW}$ [Hz]')

ax1.set_xlim(ptmin, ptmax)
ax1.set_ylim(-1., 30.)


#COMPARE PEAKS:
##rest frame
#t_SI_RF_3B_R0_peakvals_PAset    = t_SI_RF_3B_R0_peakvals[pos_3B_peak_PA_3B2B]
#t_SI_RF_2B_R0_peakvals_PAset    = t_SI_RF_2B_R0_peakvals[pos_2B_peak_PA_3B2B]
##obs frame
#t_SI_obs_3B_R0_peakvals_PAset   = t_SI_obs_3B_R0_peakvals[pos_3B_peak_PA_3B2B]
#t_SI_obs_2B_R0_peakvals_PAset   = t_SI_obs_2B_R0_peakvals[pos_2B_peak_PA_3B2B]

#Define linear dopple factor at ref point:
v_LOS_SI    = np.dot(obs_norm_vec, vel_ijCOM_R0)*np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)
Dfac        = (1.-v_LOS_SI/c_SI) 

#total difference between obs 3B and 2B:
xp  = t_SI_obs_3B_R0_peakvals_PAset
yp  = t_SI_obs_3B_R0_peakvals_PAset - t_SI_obs_2B_R0_peakvals_PAset
ax2.plot(xp, yp, marker = 'o', markersize = 8, linestyle='', alpha = 0.75, color='c', label = r'$\delta t$ Obs.')

#contribution from doppler:
#consider now 3B signal moving on linear path as 2B. Comparing this to osb 3D to quitify effect of doppler delay:
xp  = t_SI_obs_3B_R0_peakvals_PAset
yp  = t_SI_obs_3B_R0_peakvals_PAset - Dfac*t_SI_RF_3B_R0_peakvals_PAset
ax2.plot(xp, yp, marker = 'o', markersize = 3, linestyle='', alpha = 0.75, color='black', label = r'$\delta t$ Doppler')

#contribution from 'tides':
xp  = t_SI_obs_3B_R0_peakvals_PAset
yp  = Dfac*(t_SI_RF_3B_R0_peakvals_PAset - t_SI_RF_2B_R0_peakvals_PAset)
ax2.plot(xp, yp, marker = '+', markersize = 6, linestyle='', alpha = 0.75, color='black', label = r'$\delta t$ Tides')

ax2.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = False, ncol=1, framealpha = 1.0)
ax2.set_xlabel(r'time [sec]')
ax2.set_ylabel(r'$\delta t$ [sec]')

ax2.set_xlim(ptmin, ptmax)
ax2.set_ylim(-10., 10000)

plt.savefig(file_name + '_peak_dt_analysis_1.pdf', bbox_inches='tight')

plt.show()










fig = plt.figure(figsize=(5, 7))
ax1 = fig.add_subplot(321)
ax2 = fig.add_subplot(322)
ax3 = fig.add_subplot(323)
ax4 = fig.add_subplot(324)
ax5 = fig.add_subplot(325)
ax6 = fig.add_subplot(326)


#fig 1:
#2-body
xp  = t_SI_RF_2B_R0
yp  = (data_2B_a_arr)
ax1.plot(xp, yp, linewidth=2.0, linestyle='-', alpha = 0.50, color='coral', label = '2-body')
ax1.plot(xp, yp*(1.-data_2B_e_arr), linewidth=2.0, linestyle='-', alpha = 0.50, color='coral', label = '2-body')

#3-body
xp  = t_SI_RF_3B_R0
yp  = (data_a_arr)
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.85, color='black', label = '3-body')
ax1.plot(xp, yp*(1.-data_e_arr), linewidth=1.0, linestyle='-', alpha = 0.85, color='black', label = '3-body')

ax1.plot(xp[pos_ref], yp[pos_ref], marker = 'o', markersize=3, alpha = 0.85, color='black')

#labels etc:
ax1.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = False, ncol=1, framealpha = 1.0)
ax1.set_ylabel(r'log sma $a$ [AU]')
ax1.set_xlim(ptmin, ptmax)
#ax1.set_ylim(-4., -2.5)
ax1.xaxis.set_major_formatter(plt.NullFormatter())


#fig 2:
#2-body
xp  = t_SI_RF_2B_R0
yp  = data_2B_e_arr
ax2.plot(xp, yp, linewidth=2.0, linestyle='-', alpha = 0.50, color='coral', label = '2-body')
#3-body
xp  = t_SI_RF_3B_R0
yp  = data_e_arr
ax2.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.85, color='black', label = '3-body')

ax2.plot(xp[pos_ref], yp[pos_ref], marker = 'o', markersize=3, alpha = 0.85, color='black')

#labels etc:
ax2.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = False, ncol=1, framealpha = 1.0)
ax2.set_ylabel(r'ecc $e$')
ax2.set_xlim(ptmin, ptmax)
#ax2.set_ylim(0.95, 1.0)
ax2.xaxis.set_major_formatter(plt.NullFormatter())


#fig 3:
xp  = t_SI_RF_3B_R0
yp  = np.log10(F_tid_bin)
ax3.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.85, color='black', label = '')
#labels etc:
ax3.set_ylabel(r'$F_{tid}/F_{bin}$')
ax3.set_xlim(ptmin, ptmax)
#ax3.set_ylim(-0.01, 0.25)
ax3.xaxis.set_major_formatter(plt.NullFormatter())


#fig 4:
xp  = t_SI_RF_3B_R0
yp  = angle_vij_vijend1/(2.*np.pi/360.)
ax4.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.85, color='black')
#labels etc:
ax4.set_ylabel(r'$\theta(\Delta v_{ij})$ [deg]')
ax4.set_xlabel(r'time [sec]')
ax4.set_ylim(-10., 150.)
ax4.set_xlim(ptmin, ptmax)


#fig 5:
xp  = t_SI_RF_2B_R0
yp  = np.log10(c0_2B_arr)
ax5.plot(xp[pos_2B_r_peaks], yp[pos_2B_r_peaks], linewidth=1.0, linestyle='-', alpha = 0.85, color='red')

xp  = t_SI_RF_3B_R0
yp  = np.log10(c0_3B_arr)
ax5.plot(xp[pos_3B_r_peaks], yp[pos_3B_r_peaks], linewidth=1.0, linestyle='-', alpha = 0.85, color='black')
ax5.plot(xp[pos_ref], yp[pos_ref], marker = 'o', markersize=3, alpha = 0.85, color='black')



#fig 6:
xp  = t_SI_RF_3B_R0_peakvals_PAset
yp  = t_SI_RF_3B_R0_peakvals_PAset/t_SI_RF_2B_R0_peakvals_PAset
ax6.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.85, color='black')




plt.savefig(file_name + '_analysis_plot_1.pdf', bbox_inches='tight')

plt.show()


exit()












#NO Romer:
xp  = t_SI_RF_3B_R0
yp  = data_f_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='black', label = 'test 1')

xp  = t_SI_RF_2B_R0
yp  = data_2B_f_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'test 2')

xp  = [t_SI_RF_3B_R0[pos_ims_0], t_SI_RF_3B_R0[pos_ims_0]]
yp  = [-0,1e2]
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'test 2')

ax1.set_xlabel(r'time')
ax1.set_ylabel(r'$f_{GW}$')






xp  = t_SI_RF_3B_R0
yp  = data_f_arr
ax3.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='black', label = 'test 1')

xp  = t_SI_obs_3B_R0
yp  = data_f_arr
ax3.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'test 1')
ax3.set_xlabel(r'time')
ax3.set_ylabel(r'$f_{GW}$')



v_LOS_SI    = np.dot(obs_norm_vec, vel_ijCOM_R0)*np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)
Dfac        = (1-v_LOS_SI/c_SI) 
print Dfac

xp  = t_SI_RF_2B_R0
yp  = data_2B_f_arr
ax4.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='black', label = 'test 1')

#YES THIS IS THE SOLUTION: JUST DOPPLER WE CAN MOVE ON!
xp  = Dfac*t_SI_RF_2B_R0
yp  = data_2B_f_arr
ax4.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='blue', label = 'test 1')

xp  = t_SI_obs_2B_R0
yp  = data_2B_f_arr
ax4.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'test 1')

ax4.set_xlabel(r'time')
ax4.set_ylabel(r'$f_{GW}$')



#xp  = (t_3B_peakvals[pos_3B_peak_PA_3B2B] - t_RF_R0)*tU_tSI
#yp  = dt_phs_3B2B*tU_tSI
#ax5.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'test 1')




#compare c0

data_2B_rp_arr  = data_2B_a_arr*(1. - data_2B_e_arr)
data_2B_c0_arr  = data_2B_rp_arr*(1.+data_2B_e_arr)*(data_2B_e_arr**(-12./19.))*((1.+(121./304.)*data_2B_e_arr**2.)**(-870./2299.))

data_rp_arr  = data_a_arr*(1. - data_e_arr)
data_c0_arr  = data_rp_arr*(1.+data_e_arr)*(data_e_arr**(-12./19.))*((1.+(121./304.)*data_e_arr**2.)**(-870./2299.))

xp  = t_SI_RF_2B_R0
yp  = data_2B_c0_arr
ax5.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'test 1')

xp  = t_SI_RF_3B_R0
yp  = data_c0_arr
ax5.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='blue', label = 'test 1')




#LOOK THROUGH THE CODE AGAIN.
#PICK/MAKE FIGS FOR PAPER
#IMPLEMENT STOPPING CRITERIA IN CODE, THEN FIND NEW EXAMPLE FOR a = 1 AU.













exit()










#TEST:
fig = plt.figure(figsize=(5, 3))
ax1 = fig.add_subplot(111)
#ax2 = fig.add_subplot(212)
ptmin       = -2000
ptmax       = 100

xp  = t_SI_obs_3B_R0
yp  = data_r_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='black', label = 'Correct 3-body Path')
#pos = np.where(data_f_arr > 10.)[0]
#ax1.plot(xp[pos], yp[pos], linewidth=1.0, linestyle='-', alpha = 0.5, color='blue')



xp  = t_SI_RF_3B_R0
yp  = data_r_arr
#ax1.plot(xp, yp, linewidth=1.0, linestyle=':', alpha = 0.5, color='black')

xp  = t_SI_obs_2B_R0
yp  = data_2B_r_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red', label = 'Linear 2-body Path')

xp  = t_SI_RF_2B_R0
yp  = data_2B_r_arr
#ax1.plot(xp, yp, linewidth=1.0, linestyle=':', alpha = 0.5, color='red')


ax1.set_xlim(ptmin, ptmax)
ax1.set_ylim(-0.1, 1)
ax1.set_xlabel(r'time')
ax1.set_ylabel(r'distance')
ax1.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = True, ncol=1, framealpha = 1.0)

fig.savefig('rij_time_ex1.pdf', bbox_inches='tight', dpi=500)


plt.show()








#TEST:
fig = plt.figure(figsize=(6, 8))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)
ptmin       = -2000
ptmax       = 100


#ax1.set_xlim(ptmin, ptmax)
#ax1.set_ylim(-0.1, 1)

plt.show()


#exit()







#---------------------------------------------------------------------
#PLOT:
#---------------------------------------------------------------------
plot_t_3B   = t_SI_RF_3B_R0
plot_t_2B   = t_SI_RF_2B_R0

#FIG 1:
fig = plt.figure(figsize=(6, 8))
ax1 = fig.add_subplot(311)
ax2 = fig.add_subplot(312)
ax3 = fig.add_subplot(313)
#plt.setp(ax1, xticklabels=[])
#plt.setp(ax2, xticklabels=[])

ptmin       = -40000#-1750.
ptmax       = 0.0

#plot ax1:
#3-Body
xp  = plot_t_3B
yp  = data_r_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='black', label = r'$3$-body')
#2-Body
xp  = plot_t_2B
yp  = data_2B_r_arr
ax1.plot(xp, yp, linewidth=1.0, linestyle='-', alpha = 0.75, color='red',   label = r's-s GW cap.')
#3B-body fGW limits:
fGW_l           = 1.0
pos_3B_window   = np.where(plot_t_3B > ptmin)[0]
pos_fGW_l       = np.where(data_f_arr > fGW_l)[0]
pos_fGW_l       = min(list(set(pos_fGW_l).intersection(pos_3B_window)))
#ax1.plot([plot_t_3B[pos_fGW_l], plot_t_3B[pos_fGW_l]], [0,10], alpha = 0.5, color='black', linestyle = ':')
fGW_l           = 5.0
pos_3B_window   = np.where(plot_t_3B > ptmin)[0]
pos_fGW_l       = np.where(data_f_arr > fGW_l)[0]
pos_fGW_l       = min(list(set(pos_fGW_l).intersection(pos_3B_window)))
#ax1.plot([plot_t_3B[pos_fGW_l], plot_t_3B[pos_fGW_l]], [0,10], alpha = 0.5, color='black', linestyle = ':')
fGW_l           = 7.0
pos_3B_window   = np.where(plot_t_3B > ptmin)[0]
pos_fGW_l       = np.where(data_f_arr > fGW_l)[0]
pos_fGW_l       = min(list(set(pos_fGW_l).intersection(pos_3B_window)))
#ax1.plot([plot_t_3B[pos_fGW_l], plot_t_3B[pos_fGW_l]], [0,10], alpha = 0.5, color='black', linestyle = ':')
fGW_l           = 9.0
pos_3B_window   = np.where(plot_t_3B > ptmin)[0]
pos_fGW_l       = np.where(data_f_arr > fGW_l)[0]
pos_fGW_l       = min(list(set(pos_fGW_l).intersection(pos_3B_window)))
#ax1.plot([plot_t_3B[pos_fGW_l], plot_t_3B[pos_fGW_l]], [0,10], alpha = 0.5, color='black', linestyle = ':')
fGW_l           = 10.0
pos_3B_window   = np.where(plot_t_3B > ptmin)[0]
pos_fGW_l       = np.where(data_f_arr > fGW_l)[0]
pos_fGW_l       = min(list(set(pos_fGW_l).intersection(pos_3B_window)))
#ax1.plot([plot_t_3B[pos_fGW_l], plot_t_3B[pos_fGW_l]], [0,10], alpha = 0.5, color='black', linestyle = ':')
#axis settings etc:
ax1.legend(loc='upper right', numpoints = 1, fontsize = 10, frameon = True, ncol=1, framealpha = 1.0)
ax1.set_xlim(ptmin, ptmax)
ax1.set_ylim(-0.1, 10)
ax1.set_ylabel(r'dist $r_{ij}$')

#t_SI_RF_3B_R0   = t_SI_RF_3B - t_SI_RF_R0
#t_SI_RF_2B_R0   = t_SI_RF_2B - t_SI_RF_R0
#t_SI_obs_3B_R0  = t_SI_obs_3B - t_SI_obs_R0
#t_SI_obs_2B_R0  = t_SI_obs_2B - t_SI_obs_R0

#plot ax2:
xp  = t_SI_RF_3B_R0
yp  = data_r_arr
ax2.plot(xp, yp, linewidth=1.5, linestyle='-', alpha = 0.5, color='black')

xp  = t_SI_obs_3B_R0
yp  = data_r_arr
ax2.plot(xp, yp, linewidth=1.5, linestyle=':', alpha = 0.5, color='black')

xp  = t_SI_RF_2B_R0
yp  = data_2B_r_arr
ax2.plot(xp, yp, linewidth=1.5, linestyle='-', alpha = 0.5, color='red')

xp  = t_SI_obs_2B_R0
yp  = data_2B_r_arr
ax2.plot(xp, yp, linewidth=1.5, linestyle=':', alpha = 0.5, color='red')

v_LOS_SI    = np.dot(obs_norm_vec, vel_ijCOM_R0)*np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)
Dfac        = (1-v_LOS_SI/c_SI) 
xp  = t_SI_RF_2B_R0*Dfac
yp  = data_2B_r_arr
ax2.plot(xp, yp, linewidth=0.5, linestyle='-', alpha = 0.5, color='blue')

#axis settings etc:
ax2.set_xlim(-10000, 10000)
ax2.set_ylim(-0.1, 10)
ax2.set_ylabel(r'dist $r_{ij}$')


#plot ax3:
xp  = plot_t_3B
yp  = data_e_arr
ax3.plot(xp, yp, linewidth=1.5, linestyle='-', alpha = 0.5, color='purple')
#axis settings etc:
ax3.set_xlim(ptmin, ptmax)
ax3.set_ylim(0.8, 1.0)
ax3.set_xlabel(r'time $t-t_0$')
ax3.set_ylabel(r'eccentricity $e$')

#save and show:
#fig.savefig('shift_2B3B_ecc.pdf', bbox_inches='tight', dpi=500)
plt.show()



#FIG 2:
penr    = 10

fig = plt.figure(figsize=(6, 6))
ax1 = fig.add_subplot(211)
ax2 = fig.add_subplot(212)

#plot ax1:
xp  = plot_t_3B/10000.
yp  = data_f_arr
ax1.plot(xp[::penr], yp[::penr],    linewidth=0.5, linestyle='-', alpha = 0.5, color='black')
ax1.plot([-1e10, 1e10], [1, 1],     linewidth=1.0, linestyle=':', alpha = 0.5, color='black')
ax1.plot([-1e10, 1e10], [5, 5],     linewidth=1.0, linestyle=':', alpha = 0.5, color='black')
ax1.plot([-1e10, 1e10], [10, 10],   linewidth=1.0, linestyle=':', alpha = 0.5, color='black')
#axis settings etc:
ax1.set_yscale('log')
ax1.set_xlim(min(xp), max(xp)+1)
ax1.set_ylim(1e-5, 30.)
ax1.set_xlabel(r'$t-t_0$')
ax1.set_ylabel(r'$f_{GW}$')

#plot ax2:
tfac    = 0.05
xp = b1_posxyz[:,0] + tfac*time[:]
yp = b1_posxyz[:,1]
ax2.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.2, color = 'grey', rasterized=False)
xp = b2_posxyz[:,0] + tfac*time[:]
yp = b2_posxyz[:,1]
ax2.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.2, color = 'grey', rasterized=False)
xp = b3_posxyz[:,0] + tfac*time[:]
yp = b3_posxyz[:,1]
ax2.plot(xp[::penr], yp[::penr], linewidth = 1.0, linestyle='-', marker = '', alpha = 0.2, color = 'red', rasterized=False)

rpij_ae_SI      = (data_a_arr*R_sun_SI)*(1.-data_e_arr)
data_f_ae_arr   = (1./np.pi)*np.sqrt((G_new_SI*(mi+mj)*M_sun_SI)/(rpij_ae_SI**3.))
pos = np.where(data_f_ae_arr > 1.0)[0]
xp = b1_posxyz[pos,0] + tfac*time[pos]
yp = b1_posxyz[pos,1]
#ax2.plot(xp, yp, markersize=0.05, linestyle='', marker = '.', alpha = 0.75, color = 'black')
xp = b2_posxyz[pos,0] + tfac*time[pos]
yp = b2_posxyz[pos,1]
ax2.plot(xp[::penr], yp[::penr], markersize=1.0, linestyle='', marker = '.', alpha = 0.75, color = 'black', rasterized=False)
xp = b3_posxyz[pos,0] + tfac*time[pos]
yp = b3_posxyz[pos,1]
ax2.plot(xp[::penr], yp[::penr], markersize=1.0, linestyle='', marker = '.', alpha = 0.75, color = 'red', rasterized=False)

#axis settings etc:
ax2.set_xlim(-5, 12)
ax2.set_ylim(-4, 4)
ax2.set_xlabel(r'$x$ pos')
ax2.set_ylabel(r'$y$ pos')

#save and show:
#fig.savefig('fGW_3borbits.pdf', bbox_inches='tight', dpi=500)
plt.show()







#FIG 3:
#GW190814: 2.5, 25.
#GW190412: 8., 30.
mb1_cp  = 8.0
mb2_cp  = 30.
nr_Ms   = 1000
nr_as   = 1000
Ms_pl   = np.linspace(10,100,nr_Ms)
as_pl   = 10.**(np.linspace(-2,1,nr_as))
fp_arr  = np.zeros((nr_Ms, nr_as), dtype=np.float64)
kfac    = (6.*np.sqrt(2.)/(85.*np.pi))*((c_SI**5.)/(G_new_SI**(5./2.)))*((G_new_SI**(7./6.))/(np.pi**(7./3.)))
tgfac   = 1./32.
for mc in range(0,nr_Ms):
    for ac in range(0,nr_as):
        fp_arr[mc,ac]   = (((kfac**3.)*((mb1_cp+mb2_cp)*M_sun_SI)/((mb1_cp*M_sun_SI*mb2_cp*M_sun_SI)**3.))*(1./tgfac)*((Ms_pl[mc]*M_sun_SI)/((as_pl[ac]*AU_SI)**3.)))**(1./7.)
X, Y    = np.meshgrid(Ms_pl, as_pl)
Z       = np.transpose(fp_arr[:,:]) 
#plot colors:
fig = plt.figure(figsize=(6, 5))
ax  = plt.subplot(111)
quad = plt.pcolormesh(X, Y, Z, cmap='gist_earth', linewidth=0, rasterized=True, vmin=1.0, vmax=10.)
plt.colorbar()
#oplot contours:
mpl.rcParams['contour.negative_linestyle'] = 'solid'
CS = plt.contour(X, Y, Z, list([1,2,3,4,5,6,7,8,9,10]), colors='red')
plt.clabel(CS, fontsize=9, inline=1, fmt='%1.1f')

#oplot line:
fEbeta          = 2./21.
mb2Ms_mu        = (mb2_cp*Ms_pl/(mb2_cp + Ms_pl))*M_sun_SI  #SI
v_esc           = 100.0*(1000.)  #m/s                        #SI
a_esc_AU        = ((fEbeta*G_new_SI*mb2Ms_mu)/(v_esc**2.))/AU_SI
ax.plot(Ms_pl, a_esc_AU, linestyle = '--', linewidth = 2, color='yellow')

fEbeta          = 2./21.
mb2Ms_mu        = (mb2_cp*Ms_pl/(mb2_cp + Ms_pl))*M_sun_SI  #SI
v_esc           = 250.0*(1000.)  #m/s                        #SI
a_esc_AU        = ((fEbeta*G_new_SI*mb2Ms_mu)/(v_esc**2.))/AU_SI
ax.plot(Ms_pl, a_esc_AU, linestyle = '--', linewidth = 2, color='yellow')

ax.set_yscale('log')
ax.set_xlabel(r'$m_{p}$ [$M_{\odot}$]')
ax.set_ylabel(r'SMA $a_{p}$ $[AU]$')
ax.set_title(r'$m_1 = 8M_{\odot}, m_2 = 30M_{\odot}$ ($\sim GW190412$)')

#save and show:
#plt.savefig('fGWp_mp_ap.pdf', bbox_inches='tight')     
plt.show()




#-----------------------------------------------------------------
exit()







#---------------------------------------
#2-D Plot:
#---------------------------------------
#make plot window:
fig = plt.figure(figsize=(7, 7))
#object 1:
xp1 = b1_posxyz[0::10,0]/sma_a0 + 0.*time[0::10]
yp1 = b1_posxyz[0::10,1]/sma_a0
zp1 = b1_posxyz[0::10,2]/sma_a0
fig.add_subplot(111).plot(xp1, yp1, linewidth=2, linestyle='-', alpha = 0.75, color='dodgerblue')
#object 2:
xp2 = b2_posxyz[0::10,0]/sma_a0 + 0.*time[0::10]
yp2 = b2_posxyz[0::10,1]/sma_a0
zp2 = b2_posxyz[0::10,2]/sma_a0
fig.add_subplot(111).plot(xp2, yp2, linewidth=2, linestyle='-', alpha = 0.75, color='green')
#object 3:
xp3 = b3_posxyz[0::10,0]/sma_a0 + 0.*time[0::10]
yp3 = b3_posxyz[0::10,1]/sma_a0
zp3 = b3_posxyz[0::10,2]/sma_a0
fig.add_subplot(111).plot(xp3, yp3, linewidth=2, linestyle='-', alpha = 0.75, color='purple')
#---------------------------------------
#axis settings/labels/ranges, etc:
#---------------------------------------
ax = fig.add_subplot(111)
ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
#ax.set_xlim(-8.5, 8.5)
#ax.set_ylim(-8.5, 8.5)
ax.set_xlim(1.9, 1.98)
ax.set_ylim(-3.2, -2.7)

#ax.text(-10.0, 1.9, r'$\rm{Binary}$ $\rightarrow$', size = 20,
#        horizontalalignment='left',
#       verticalalignment='center',
#        rotation = -22.5)
#ax.text(10.5, -0.2, r'$\leftarrow$ $\rm{Single}$', size = 20,
#        horizontalalignment='right',
#        verticalalignment='center',
#        rotation = - 21.0)
#ax.text(-3.5, -1.5, r'$\rm{GW\ capture}$', size = 15,
#        horizontalalignment='center', fontname='serif',
#        verticalalignment='center',
#        rotation = 0.0)
#show:    
plt.show()
#fig.savefig('orbit3body_ex_1.eps', bbox_inches='tight')
#fig.savefig('AGNorb1_all.eps', bbox_inches='tight')
#fig.savefig('AGNorb1_zoom.eps', bbox_inches='tight')
#-----------------------------------------------------------------    





#-----------------------------------------------------------------
#fGW plot:
#-----------------------------------------------------------------
m1_SI   = 50.*M_sun_SI
m2_SI   = 50.*M_sun_SI
m3_SI   = 50.*M_sun_SI
T0      = (2.*np.pi*np.sqrt(((0.1*AU_U)**3.)/((m1_SI/M_sun_SI+m2_SI/M_sun_SI))))
Tminp   = 0
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
#r12:
r_12    = np.sqrt((xp1 - xp2)**(2.) + (yp1 - yp2)**(2.) + (zp1 - zp2)**(2.))*R_sun_SI
f_12    = (1./np.pi)*np.sqrt(G_new_SI*(m1_SI+m2_SI)/(r_12**3.))
#r13
r_13    = np.sqrt((xp1 - xp3)**(2.) + (yp1 - yp3)**(2.) + (zp1 - zp3)**(2.))*R_sun_SI
f_13    = (1./np.pi)*np.sqrt(G_new_SI*(m1_SI+m3_SI)/(r_13**3.))
#r23
r_23    = np.sqrt((xp2 - xp3)**(2.) + (yp2 - yp3)**(2.) + (zp2 - zp3)**(2.))*R_sun_SI
f_23    = (1./np.pi)*np.sqrt(G_new_SI*(m2_SI+m3_SI)/(r_23**3.))
#make plot window:
fig = plt.figure(figsize=(7, 7))
ax = fig.add_subplot(111)
ax.plot(time/T0 - Tminp, f_12, linewidth=2, linestyle='-', alpha = 0.75, color='grey')
ax.plot(time/T0 - Tminp, f_13, linewidth=2, linestyle='-', alpha = 0.75, color='grey')
ax.plot(time/T0 - Tminp, f_23, linewidth=2, linestyle='-', alpha = 0.75, color='red')

plt.yscale('log')

ax.set_xlim(Tminp, 60)
ax.set_ylim(1e-7, 20)

plt.show()
#fig.savefig('AGNorb1_fGW.eps', bbox_inches='tight')
#-----------------------------------------------------------------
exit()




   
#---------------------------------------
#3-D Plot:
#---------------------------------------
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#make plot window:
fig = plt.figure(figsize=(6, 6))
ax  = fig.add_subplot(111, projection='3d')
nrp = len(b1_posxyz[:,0])
ep = 100
for i in range(0,nrp-1,ep):
    veldist_x = 0.01*time[i]
    #object 1:
    xp1 = b1_posxyz[i:i+ep,0]/sma_a0 + veldist_x
    yp1 = b1_posxyz[i:i+ep,1]/sma_a0
    zp1 = b1_posxyz[i:i+ep,2]/sma_a0
    ax.plot(xp1, yp1, zp1,      linewidth = 1.5, color='black',  alpha = 1)
    #object 2:
    xp2 = b2_posxyz[i:i+ep,0]/sma_a0 + veldist_x
    yp2 = b2_posxyz[i:i+ep,1]/sma_a0
    zp2 = b2_posxyz[i:i+ep,2]/sma_a0
    ax.plot(xp2, yp2, zp2,      linewidth = 1.5, color='black',  alpha = 1)
    #object 3:
    xp3 = b3_posxyz[i:i+ep,0]/sma_a0 + veldist_x
    yp3 = b3_posxyz[i:i+ep,1]/sma_a0
    zp3 = b3_posxyz[i:i+ep,2]/sma_a0
    ax.plot(xp3, yp3, zp3,      linewidth = 1.5, color='orange',  alpha = ((1.*i)/(1.*nrp))**0.25)

p3 = b3_posxyz[ep,:]/sma_a0
#ax.text(1.0*p3[0] + veldist_x, 1.2*p3[1], 1.2*p3[2], "Incoming planet", color='dodgerblue', fontsize = 12)
ax.plot([p3[0]], [p3[1]], [p3[2]],  marker = 'o', markersize = 5, color='orange', markeredgecolor = 'black')
p3 = b3_posxyz[nrp-1,:]/sma_a0
#ax.text(1.0*p3[0] + veldist_x, 2.0*p3[1], 1.0*p3[2], r"Disruption",   color='dodgerblue', fontsize = 12)
ax.plot([p3[0]] + veldist_x, [p3[1]], [p3[2]],  marker = '*', markersize = 10, color='red', markeredgecolor = 'black')

WR = 1.0
#ax.set_xlim3d(-WR, WR)
ax.set_ylim3d(-3.*WR,  3.*WR)
ax.set_zlim3d(-1.25*WR, 1.25*WR)
ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
ax.set_zlabel(r'$z/a$')
ax.view_init(azim=40,elev=19)

#-----------------------------------------------------------------    
#Finalize plot and save:
#-----------------------------------------------------------------
plt.show()
fig.savefig('KK4.pdf', bbox_inches='tight')
exit()
#-----------------------------------------------------------------    

 
    
    
    
    
    
    
ax  = fig.add_subplot(111, projection='3d')
nrp = len(b1_posxyz[:,0])
for i in range(0,nrp,100):
    #object 1:
    xp1 = b1_posxyz[i,0]/sma_a0
    yp1 = b1_posxyz[i,1]/sma_a0
    zp1 = b1_posxyz[i,2]/sma_a0
    ax.scatter(xp1, yp1, zp1,     marker='o', c='black', color = 'none', alpha = 0.2)
    #object 2:
    xp2 = b2_posxyz[i,0]/sma_a0
    yp2 = b2_posxyz[i,1]/sma_a0
    zp2 = b2_posxyz[i,2]/sma_a0
    ax.scatter(xp2, yp2, zp2,     marker='o', c='blue', color = 'none', alpha = 0.2)
    #object 3:
    xp3 = b3_posxyz[i,0]/sma_a0
    yp3 = b3_posxyz[i,1]/sma_a0
    zp3 = b3_posxyz[i,2]/sma_a0
    ax.scatter(xp3, yp3, zp3,     marker='o', c='orange', color = 'none')

    
    
    
    
    
    
    
    

