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
import matplotlib.gridspec as gridspec

#Input IC filename:
IC_filename = 'MCinput_Nbody.txt'#'paper1_ex_tidalinsp.txt'#'tidalinsp_ex1_WD12NS12.txt'#'WD06Tidalinsp_1_3body.txt'#'insp_ill_1_inputparams.txt'#'test10.txt'#'MCinput_Nbody.txt'#'MCinput_Nbody.txt'#'3body_Tidal_insp.txt'



#------------------------------------------------------------------------------------------
#Read IC data and info from input file:
#------------------------------------------------------------------------------------------
tf = open(IC_filename, "r")
fline_split = tf.readline().split()
n_particles = int(fline_split[0])
#Nbody params:
nbody_params_arr_1_INT  = np.array([int(i)        for i in tf.readline().split()[0:10]])
nbody_params_arr_2_REAL = np.array([float(i)      for i in tf.readline().split()[0:10]])
#Define:
IC_allobj_const_arr    = np.zeros((n_particles, 10),    dtype=np.float64)
IC_allobj_posxyz_CM    = np.zeros((n_particles, 3),     dtype=np.float64)  
IC_allobj_velxyz_CM    = np.zeros((n_particles, 3),     dtype=np.float64)
IC_allobj_q            = np.zeros((n_particles, 3,3),   dtype=np.float64)
IC_allobj_qdot         = np.zeros((n_particles, 3,3),   dtype=np.float64)
#Read IC infor per obj:
for pc in range(0,n_particles):
        #const arr:
        IC_allobj_const_arr[pc,:]   = np.array([float(i)      for i in tf.readline().split()[0:10]])
        #pos,vel:
        IC_allobj_posxyz_CM[pc,:]   = np.array([float(i)      for i in tf.readline().split()[0:3]])
        IC_allobj_velxyz_CM[pc,:]   = np.array([float(i)      for i in tf.readline().split()[0:3]])
        #q 
        IC_allobj_q[pc,0,:]         = np.array([float(i)      for i in tf.readline().split()[0:3]])
        IC_allobj_q[pc,1,:]         = np.array([float(i)      for i in tf.readline().split()[0:3]])
        IC_allobj_q[pc,2,:]         = np.array([float(i)      for i in tf.readline().split()[0:3]])
        #qdot: 
        IC_allobj_qdot[pc,0,:]      = np.array([float(i)      for i in tf.readline().split()[0:3]])
        IC_allobj_qdot[pc,1,:]      = np.array([float(i)      for i in tf.readline().split()[0:3]])
        IC_allobj_qdot[pc,2,:]      = np.array([float(i)      for i in tf.readline().split()[0:3]])
#Define:
allobj_mass     = IC_allobj_const_arr[:,0]          
allobj_radius   = IC_allobj_const_arr[:,1]
#calc ini E of single wrt. bin:     
vel_CM_0_bin12  = (allobj_mass[0]*IC_allobj_velxyz_CM[0,:] + allobj_mass[1]*IC_allobj_velxyz_CM[1,:])/(allobj_mass[0]+allobj_mass[1])
pos_CM_0_bin12  = (allobj_mass[0]*IC_allobj_posxyz_CM[0,:] + allobj_mass[1]*IC_allobj_posxyz_CM[1,:])/(allobj_mass[0]+allobj_mass[1])
mbin = allobj_mass[0]+allobj_mass[1]
msin = allobj_mass[2]
mtot    = mbin+msin
mred    = mbin*msin/mtot
v_binsin    = np.linalg.norm(IC_allobj_velxyz_CM[2,:] - vel_CM_0_bin12[:])
r_binsin    = np.linalg.norm(IC_allobj_posxyz_CM[2,:] - pos_CM_0_bin12[:])
E_CM_0_sin  = (1./2.)*mred*v_binsin**2. - mred*mtot/r_binsin   
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#Run Nbody code:
#------------------------------------------------------------------------------------------
#subprocess.call('./Nbody_AffineTides_solver_6.exe' + '<' + IC_filename, shell=True)
simyesno = raw_input('simulate (again)? (yes=1, no=0): ')
if (int(simyesno) == 1):
    print 'TEST'
    subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + IC_filename, shell=True)
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#Read in data and define:
#------------------------------------------------------------------------------------------
#---------------------------------------
#Read in data from Nbody solver:
#---------------------------------------
#File names:
fn_NbodyTides_dataout_pos           = 'NbodyTides_dataout_pos.dat'
fn_NbodyTides_dataout_a1a2a3        = 'NbodyTides_dataout_a1a2a3.dat' 
fn_NbodyTides_dataout_Eself         = 'NbodyTides_dataout_Eself.dat'
fn_NbodyTides_dataout_Etot          = 'NbodyTides_dataout_Etot.dat'
fn_NbodyTides_dataout_binij_info    = 'NbodyTides_dataout_binij_info.dat' 
fn_NbodyTides_dataout_full_q        = 'NbodyTides_dataout_full_q.dat' 
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
tf = open(fn_NbodyTides_dataout_full_q,  "r")
NbodyTides_dataout_full_q       = np.loadtxt(tf, dtype=float)
tf.close()
#---------------------------------------
#Define (3 objects:)
#---------------------------------------
#---------------------
#Format (for 3 objs):
#---------------------
#NbodyTides_dataout_pos:        pos(1,:),  pos(2,:), pos(3,:)
#NbodyTides_dataout_a1a2a3:     time_t, body_all_a1a2a3(1,:),  body_all_a1a2a3(2,:), body_all_a1a2a3(3,:)      
#NbodyTides_dataout_Eself:      time_t, body_all_Eselftot(1),  body_all_Eselftot(2), body_all_Eselftot(3)
#NbodyTides_dataout_Etot:       time_t, E_tot_kin, E_tot_pot_pmass, E_tot_pot_tides, E_tot_internal_terms, E_tot_external_terms, E_tot_system
#NbodyTides_dataout_binij_info: time_t, bin_i, bin_j, a_bin_ij, e_bin_ij, rperi_ij, T_ij, T_ji
#---------------------
#pos:
#---------------------
b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
b3_posxyz   = NbodyTides_dataout_pos[:,6:9]
#---------------------
#time:
#---------------------
time = NbodyTides_dataout_a1a2a3[:,0]
#---------------------
#stellar axes:
#---------------------
b1_a1a2a3   =  NbodyTides_dataout_a1a2a3[:,1:4]
b2_a1a2a3   =  NbodyTides_dataout_a1a2a3[:,4:7]
b3_a1a2a3   =  NbodyTides_dataout_a1a2a3[:,7:10]
#full q (for 3D plot):
nr_time = len(time)
RS_NbodyTides_dataout_full_q = NbodyTides_dataout_full_q.reshape(nr_time,3,3,3) #time, nr obj, q(3,3)
#---------------------
#tot self energy each star:
#---------------------
b1_Eselftot = NbodyTides_dataout_Eself[:,1]
b2_Eselftot = NbodyTides_dataout_Eself[:,2]
b3_Eselftot = NbodyTides_dataout_Eself[:,3]
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
#define:
#---------------------
r_coll  = allobj_radius[0]+allobj_radius[1]#2*allobj_radius[0]    #assume all objs have the same radius.
a_0     = a_bin_ij[0]           #assume system is in well defined bin-sin state from start (t0).  
Torb_0  = 2.*np.pi*np.sqrt((a_0**3.)/(allobj_mass[0]+allobj_mass[1])) 
#---------------------
#---------------------------------------
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#PLOT:
#------------------------------------------------------------------------------------------
#--------------------------------------------------------------
#Set general plot settings incl. font size etc.:
#--------------------------------------------------------------
font = {'family' : 'serif'}
mpl.rc('font', **font)
#--------------------------------------------------------------



#--------------------------------------------------------------
#set time range for last insp:
#--------------------------------------------------------------
timeT0  = time/Torb_0                       #time in units of initial orbital time Torb_0
ti_l    = min(np.where(timeT0 > 0.0)[0])    #lower time limit index
ti_u    = len(timeT0)-1                     #upper time limit index
ti_p    = ti_l                              #time for 3d obj plot
ptimeT0 = timeT0[ti_l:ti_u]                 #time range
#--------------------------------------------------------------

print timeT0

#--------------------------------------------------------------
#Make 3D plot:
#--------------------------------------------------------------
#Define:
posxyz_all_arr  = np.zeros((nr_time,3,3),   dtype=np.float64)
posxyz_all_arr[:,0,:] = b1_posxyz
posxyz_all_arr[:,1,:] = b2_posxyz
posxyz_all_arr[:,2,:] = b3_posxyz

print 'TIME: ', timeT0[ti_p], nr_time

#Make figure:
fig = plt.figure(2,figsize=(5,8))

#------------------------------------------
#loop over plot windows:
#------------------------------------------
for wc in range(0,2):

    #Make plot:
    ax = fig.add_subplot(2, 1, wc+1, projection='3d')
    #ax = fig.add_subplot(111, projection='3d')

    #--------------------------------------
    #window range/view angles:
    #--------------------------------------
    #calc 'COM' for plot window:
    bi  = bin_i[ti_u]-1   #-1 to get python index
    bj  = bin_j[ti_u]-1   #-1 to get python index
    pos_CM_binij_t      = (allobj_mass[bi]*posxyz_all_arr[ti_u,bi,:]+allobj_mass[bj]*posxyz_all_arr[ti_u,bj,:])/(allobj_mass[bi]+allobj_mass[bj])
    pos_CM_binij_all_t  = (allobj_mass[bi]*posxyz_all_arr[:,bi,:]+allobj_mass[bj]*posxyz_all_arr[:,bj,:])/(allobj_mass[bi]+allobj_mass[bj])
    #window 1:
    if (wc == 0):
        WR = 10.5*a_0
        ax.set_xlim3d(-5.0*WR, 40.*WR)
        ax.set_ylim3d(-WR, WR)
        ax.set_zlim3d(-WR, WR)
    #window 2:
    if (wc == 1):
        WR = 0.075*a_0
        #ax.set_xlim3d(pos_CM_binij_t[0]-WR, pos_CM_binij_t[0]+WR)
        #ax.set_ylim3d(pos_CM_binij_t[1]-WR, pos_CM_binij_t[1]+WR)
        #ax.set_zlim3d(pos_CM_binij_t[2]-WR, pos_CM_binij_t[2]+WR)
        ax.set_xlim3d(-1.0*WR, 1.0*WR)
        ax.set_ylim3d(-1.0*WR, 1.0*WR)
        ax.set_zlim3d(-0.7*WR, 0.7*WR)

    if (wc == 0): ax.view_init(azim=-65,elev=38)
    if (wc == 1): ax.view_init(azim=35,elev=30)      
    #--------------------------------------
    
    #--------------------------------------
    #window labels/titles:
    #--------------------------------------
    if (wc == 0):
        ax.text2D(0.12, 0.95, "Binary-Single Interaction (Full)", transform=ax.transAxes, fontsize = 16.5)       
    if (wc == 1):
        ax.text2D(0.12, 0.95, "Final Tidal Interaction (Zoom in)", transform=ax.transAxes, fontsize = 16.5)
    #--------------------------------------
 
    #--------------------------------------
    #plot object trajectories:
    #--------------------------------------
    colornr=[25,100,210]
    vxcom = 6.5*timeT0[:] #create a COM vel for better visualization 
    
    if (wc == 0):
        #2d projections:
        for i in range(0,3):
            xp = posxyz_all_arr[:,i,0] + vxcom[:]
            yp = posxyz_all_arr[:,i,1]
            zp = posxyz_all_arr[:,i,2]
            pp = WR+0.0*xp
            ax.plot(xp, yp, -pp,    linewidth = 0.5, alpha=0.75, color='black')
            ax.plot(xp, pp, zp,     linewidth = 0.5, alpha=0.75, color='black')
            ax.plot(-5.0*pp, yp, zp,    linewidth = 0.5, alpha=0.75, color='black')    
        #3d plot:        
        for i in range(0,3):
            xp = posxyz_all_arr[:,i,0] + vxcom[:]
            yp = posxyz_all_arr[:,i,1]
            zp = posxyz_all_arr[:,i,2]
            ax.plot(xp, yp, zp,     linewidth = 1.0, color=plt.cm.gist_earth(colornr[i]))
            
    if (wc == 1):
        for i in range(0,3):
            #3d plot:
            #---------
            xp = posxyz_all_arr[ti_l:ti_u,i,0] - pos_CM_binij_all_t[ti_l:ti_u,0]
            yp = posxyz_all_arr[ti_l:ti_u,i,1] - pos_CM_binij_all_t[ti_l:ti_u,1]
            zp = posxyz_all_arr[ti_l:ti_u,i,2] - pos_CM_binij_all_t[ti_l:ti_u,2]
            #2d projections:
            #---------
            #pp = WR+0.0*xp
            #ax.plot(xp, yp, -pp, linewidth = 0.75, alpha=0.25, color='grey')
            #ax.plot(xp, pp, zp,  linewidth = 0.75, alpha=0.5, color='grey')
            #ax.plot(pp, yp, zp,  linewidth = 0.75, alpha=0.25, color='grey')
            #---------
            ax.plot(xp, yp, zp,  linewidth = 1.0, alpha=0.75, color=plt.cm.gist_earth(colornr[i]))            
    #--------------------------------------       


    #--------------------------------------
    #plot object shapes:
    #--------------------------------------    
    if (wc == 1):
        for i in range(0,3):
            #----------------------------------
            #Object i info:
            #----------------------------------
            #const:
            q_t     = np.matrix(RS_NbodyTides_dataout_full_q[ti_p,i,:,:])
            S_t     = q_t*q_t.T        
            Eigen_vals, Eigen_vecs = np.linalg.eigh(S_t)
            Robj = allobj_radius[i]
            if (Robj < 1e-4): Robj = 0.0025*a_0
            [R1,R2,R3]  = Robj*Eigen_vals**(1./2.)#[a1,a2,a3]
            #[R1,R2,R3]  = 0.1*Eigen_vals**(1./2.)#[a1,a2,a3] 
            #----------------------------------
            #----------------------------------
            #Draw Ellipse of Obj i:
            #----------------------------------
            # Set of all spherical angles:
            plot_ell_res = 20
            u = np.linspace(0, 2*np.pi, plot_ell_res)
            v = np.linspace(0, np.pi, plot_ell_res)
            # Cartesian coordinates that correspond to the spherical angles:
            # (this is the equation of an ellipsoid):
            x_ell = (np.array(R1*np.outer(np.cos(u), np.sin(v))))
            y_ell = (np.array(R2*np.outer(np.sin(u), np.sin(v))))
            z_ell = (np.array(R3*np.outer(np.ones_like(u), np.cos(v))))
            R_ell = np.sqrt(x_ell**2 + y_ell**2 + z_ell**2)
            #rotate in CM:
            xyz_array       = np.reshape([x_ell,y_ell,z_ell], (3,plot_ell_res*plot_ell_res))
            xyz_array_ROT   = np.matrix(Eigen_vecs)*np.matrix(xyz_array)    
            x_ell_ROT       = np.reshape(xyz_array_ROT[0,:], (plot_ell_res,plot_ell_res))
            y_ell_ROT       = np.reshape(xyz_array_ROT[1,:], (plot_ell_res,plot_ell_res))
            z_ell_ROT       = np.reshape(xyz_array_ROT[2,:], (plot_ell_res,plot_ell_res))  
            #----------------------------------
            #----------------------------------
            #pos of CM for obj i:
            #----------------------------------
            posCM_x     = posxyz_all_arr[ti_p,i,0] - pos_CM_binij_all_t[ti_p,0]
            posCM_y     = posxyz_all_arr[ti_p,i,1] - pos_CM_binij_all_t[ti_p,1]
            posCM_z     = posxyz_all_arr[ti_p,i,2] - pos_CM_binij_all_t[ti_p,2]
            #----------------------------------
            #shift ellip coord to CM pos:
            #----------------------------------
            xp = x_ell_ROT + posCM_x
            yp = y_ell_ROT + posCM_y
            zp = z_ell_ROT + posCM_z
            #----------------------------------
            #plot:
            #----------------------------------            
            #draw object shapes:
            #color_R_fac = 0.5
            #color_xyz = ((R_ell/allobj_radius[i]-1.0)/color_R_fac)*0.5 + 0.5
            ##ax.plot_surface(xp, yp, zp,  rstride=1, cstride=1, facecolors=cm.bone(color_xyz), linewidth=0.1, antialiased=False, shade=False)
            #ax.plot_wireframe(xp, yp, zp, linewidth=1.0, color='black', rstride=1, cstride=1, alpha=0.5) 
            ##ax.view_init(elev=34, azim=-165)
            if (i == 0):
                color_R_fac = 0.5
                color_xyz = ((R_ell/allobj_radius[i]-1.0)/color_R_fac)*0.5 + 0.5
                #ax.plot_surface(xp, yp, zp,  rstride=1, cstride=1, facecolors=cm.bone(color_xyz), linewidth=0.1, antialiased=False, shade=False)
                ax.plot_wireframe(xp, yp, zp, linewidth=1.0, color='black', rstride=1, cstride=1, alpha=0.5) 
                #ax.view_init(elev=34, azim=-165)
            if (i == 1):
                color_R_fac = 0.5
                color_xyz = 0.0*((R_ell/allobj_radius[i]-1.0)/color_R_fac)
                ax.plot_surface(xp, yp, zp,  rstride=1, cstride=1, facecolors=cm.bone(color_xyz), linewidth=0.1, antialiased=False, shade=False)
                #ax.plot_wireframe(xp, yp, zp, linewidth=1.0, color='black', rstride=1, cstride=1, alpha=0.5) 
                #ax.view_init(elev=34, azim=-165)
            #----------------------------------        
    #--------------------------------------
    
    #--------------------------------------
    #Win figure settings:
    #--------------------------------------
    ax.w_xaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_yaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    ax.w_zaxis.set_pane_color((1.0, 1.0, 1.0, 1.0))
    
    ax.w_xaxis.set_ticklabels([])
    ax.w_yaxis.set_ticklabels([])
    ax.w_zaxis.set_ticklabels([])

    ax.w_xaxis.gridlines.set_lw(1.)
    ax.w_yaxis.gridlines.set_lw(1.)
    ax.w_zaxis.gridlines.set_lw(1.)
    #--------------------------------------
#------------------------------------------

plt.subplots_adjust(bottom=0.0, left=-0.12, right=1.12, top=1.0)
plt.subplots_adjust(hspace=0.001)
plt.subplots_adjust(wspace=0.001)

#margins(0,0)

print 'SAVING:'
#plt.savefig('3bodyex_test_3D_ill.eps')#, bbox_inches='tight', pad_inches = 0)

plt.show()
#--------------------------------------------------------------



#--------------------------------------------------------------
#Orbits and Energy:
#--------------------------------------------------------------
f = plt.figure(figsize=(5,3))
gs = gridspec.GridSpec(2, 1,height_ratios=[1,1], hspace=0.0)
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])

#---------------------------------------
#plt.subplot(211)
#---------------------------------------
#---------------------
ax1.plot(ptimeT0, b1_a1a2a3[ti_l:ti_u,0], linewidth=1.5, alpha=0.5, color='red')#, label=r'Star 1')
ax1.plot(ptimeT0, b1_a1a2a3[ti_l:ti_u,1], linewidth=1.5, alpha=0.5, color='orange')
ax1.plot(ptimeT0, b1_a1a2a3[ti_l:ti_u,2], linewidth=1.5, alpha=0.5, color='blue')
#---------------------
ax1.set_xlim([min(ptimeT0),  max(ptimeT0)])
ax1.set_ylim([0.6,2.0])
ax1.set_yticks([1.0, 1.5, 2.0])
ax1.set_ylabel(r'$a_{1}$, $a_{2}$, $a_{3}$ $[R_{WD}]$')
ax1.set_title('Dynamical evolution of final inspiral')
#---------------------
#---------------------------------------
#plt.subplot(212)
#---------------------------------------
#---------------------
E_scale = (allobj_mass[0]**2.)/allobj_radius[0]
ax2.plot(ptimeT0, a_bin_ij[ti_l:ti_u]/a_0,        linewidth=1.5, color='black', linestyle='-',  alpha=1, label=r'$a/a_0$')
ax2.plot(ptimeT0, rperi_ij[ti_l:ti_u]/a_0,        linewidth=1.5, color='black', linestyle='--',  alpha=1, label=r'$r_p/a_0$')
ax2.plot(ptimeT0, e_bin_ij[ti_l:ti_u],            linewidth=1.5, color='black', linestyle=':', alpha=1, label=r'$e$')
ax2.plot(ptimeT0, (1./E_scale)*(E_tot_external_terms[0]-E_tot_external_terms[ti_l:ti_u]),    linewidth=0.75, alpha=1, color='grey', label=r'$\Delta{E}/E_{WD}$')
#---------------------
ax2.set_xlim([min(ptimeT0),  max(ptimeT0)])
ax2.set_ylim([-0.1,  1.5])
ax2.set_yticks([0.0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.4])
ax2.set_xlabel('time $[T_{orb,0}]$')
ax2.set_ylabel(r'$a$, $r_{p}$, $e$, $\Delta{E}$')
ax2.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False, ncol=4)
#---------------------
#---------------------------------------
ax1.tick_params(
    axis='x',           # changes apply to the x,y-axis
    which='both',       # both major and minor ticks are affected
    bottom='on',        # ticks along the bottom edge are off
    top='off',           # ticks along the top edge are off
    labelbottom='off',  # labels along the bottom edge are off
    right='off',
    left='off',
    labelleft='off')
#plt.subplots_adjust(hspace=0.3)
#plt.subplots_adjust(wspace=0.3)
plt.savefig('3bodyex_test.eps', bbox_inches='tight')
#--------------------------------------------------------------

plt.show()
#exit()
#--------------------------------------------------------------
#--------------------------------------------------------------





exit()





#--------------------------------------------------------------
#Orbital, Energy Analysis:
#--------------------------------------------------------------
plt.figure(1, figsize=(10,5.2))                # the first figure

#---------------------------------------
plt.subplot(232)
#---------------------------------------
#---------------------
plt.plot(time, b1_a1a2a3[:,0], linewidth=0.5, alpha=0.75, color=plt.cm.brg(25), label=r'Star 1')
plt.plot(time, b1_a1a2a3[:,1], linewidth=0.5, alpha=0.75, color=plt.cm.brg(25))
plt.plot(time, b1_a1a2a3[:,2], linewidth=0.5, alpha=0.75, color=plt.cm.brg(25))
#---------------------
plt.plot(time, b2_a1a2a3[:,0], linewidth=0.5, alpha=0.75, color=plt.cm.brg(100), label=r'Star 2')
plt.plot(time, b2_a1a2a3[:,1], linewidth=0.5, alpha=0.75, color=plt.cm.brg(100))
plt.plot(time, b2_a1a2a3[:,2], linewidth=0.5, alpha=0.75, color=plt.cm.brg(100))
#---------------------
plt.plot(time, b3_a1a2a3[:,0], linewidth=0.5, alpha=0.5, color=plt.cm.brg(210), label='Star 3')
plt.plot(time, b3_a1a2a3[:,1], linewidth=0.5, alpha=0.5, color=plt.cm.brg(210))
plt.plot(time, b3_a1a2a3[:,2], linewidth=0.5, alpha=0.5, color=plt.cm.brg(210))
#---------------------
plt.xlim(min(time),  max(time))
plt.ylim(0.75,  10.5)
plt.locator_params(nbins=8, axis='y')
plt.locator_params(nbins=8, axis='x')
plt.xlabel(r'Time')#'/$\sqrt{R^3/GM}$')
plt.ylabel(r'$a_1,a_2,a_3$')
plt.legend(loc='upper left', numpoints = 1, fontsize = 6.0, frameon = False)
plt.title('Principal Axes')
#---------------------
#---------------------------------------
plt.subplot(233)
#---------------------------------------
#---------------------
nr_colors   = 6
color_vals  = np.arange(0,255,255/nr_colors)
plt.plot(time, E_tot_kin,               linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(color_vals[0]), label=r'$E_{kin}$')
plt.plot(time, E_tot_pot_pmass,         linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(color_vals[1]), label=r'$E_{pot} (p)$')
plt.plot(time, E_tot_pot_tides,         linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(color_vals[2]), label=r'$E_{pot} (t)$')
plt.plot(time, E_tot_internal_terms,    linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(color_vals[3]), label=r'$E_{int}$')
plt.plot(time, E_tot_external_terms,    linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(color_vals[4]), label=r'$E_{ext}$')
plt.plot(time, E_tot_system,            linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(color_vals[5]), label=r'$E_{tot}$')
#---------------------
print E_tot_system
print E_tot_pot_pmass
plt.xlim(min(time),  max(time))
plt.ylim(-2.5,  2.0)
plt.locator_params(nbins=8, axis='y')
plt.locator_params(nbins=8, axis='x')
plt.xlabel('Time')
plt.ylabel('Energy')
plt.legend(loc='upper left', numpoints = 1, fontsize = 6.0, ncol=2, frameon = False)
plt.title('Total Energy')
#---------------------
#---------------------------------------
plt.subplot(234)
#---------------------------------------
#---------------------
nr_colors   = 3
color_vals  = np.arange(0,255,255/nr_colors)
plt.plot(time, a_bin_ij/a_0,        linewidth=0.75, color=plt.cm.cubehelix(color_vals[0]), alpha=0.85, label=r'$a/a_0$')
plt.plot(time, e_bin_ij,            linewidth=0.75, color=plt.cm.cubehelix(color_vals[1]), alpha=0.85, label=r'$e$')
plt.plot(time, rperi_ij/a_0,        linewidth=0.75, color=plt.cm.cubehelix(color_vals[2]), alpha=0.85, label=r'$r_p/a_0$')
#plt.plot(time, T_ij,                linewidth=1.0, alpha=0.5)
#plt.plot(time, T_ji,                linewidth=1.0, alpha=0.5)
#---------------------
plt.xlim(min(time),  max(time))
plt.ylim(0.0,  1.75)
plt.locator_params(nbins=8, axis='y')
plt.locator_params(nbins=8, axis='x')
plt.xlabel('Time')
plt.ylabel(r'$a,e,r_p$')
plt.legend(loc='upper right', numpoints = 1, fontsize = 6.0, frameon = False)
plt.title('Binary parameters')
#---------------------
#---------------------------------------
plt.subplot(235)
#---------------------------------------
#---------------------
a_arr       = np.arange(0.001, 3*a_0, 0.01)
e_arr_coll  = (1. - r_coll/a_arr)
plt.plot(a_arr/a_0, e_arr_coll, alpha=0.25, color=plt.cm.Greys(200))
#plt.fill_between(a_arr/a_0, e_arr_coll, 1, alpha=0.25, color=plt.cm.Reds(150))
#plt.fill_between(a_arr/a_0, 0, e_arr_coll, alpha=0.25, color=plt.cm.Greys(150))
plt.plot([1.0,1.0], [0,1], alpha=0.25, linewidth=0.5, linestyle='-', color=plt.cm.Greys(250))
plt.scatter(a_bin_ij/a_0, e_bin_ij, s=0.75, c=plt.cm.gist_earth(0), lw=0)
#---------------------
plt.xlim(0.0,  5.0)
plt.ylim(0.0,  1.0)
plt.locator_params(nbins=8, axis='y')
plt.locator_params(nbins=8, axis='x')
plt.xlabel(r'semi-major axis $a/a_0$')
plt.ylabel(r'eccentricity $e$')
plt.title(r'Evolution in $(a,e)$')
#---------------------
#---------------------------------------
plt.subplot(231)
#---------------------------------------
#---------------------
plt.plot(time, bin_i,        linewidth=0.75, color=plt.cm.gray(0))
plt.plot(time, bin_j,        linewidth=0.75, color=plt.cm.gray(0))
#---------------------
plt.xlim(min(time), max(time))
plt.ylim(0, 4)
plt.yticks([1,2,3])
plt.locator_params(nbins=8, axis='x')
plt.xlabel('Time')
plt.ylabel('Object index: 1,2,3')
plt.title('Binary pair indices')
#---------------------
#---------------------------------------
plt.subplot(236)
#---------------------------------------
#---------------------
E_bin_eff   = E_tot_external_terms[:]-E_CM_0_sin
a_eff   = -allobj_mass[0]*allobj_mass[1]/(2.*E_bin_eff)
plt.plot(time, a_eff/a_0,   linewidth=0.75, alpha=0.85, color=plt.cm.cubehelix(0), label='test')
#---------------------
plt.xlim(min(time),  max(time))
plt.ylim(0.0,  1.1)
plt.locator_params(nbins=8, axis='y')
plt.locator_params(nbins=8, axis='x')
plt.xlabel('Time')
plt.ylabel(r'$a_{eff}$')
plt.title('Effective semi-major axis')
#---------------------
#---------------------------------------
plt.subplots_adjust(hspace=0.3)
plt.subplots_adjust(wspace=0.3)

#plt.savefig('3bodyex_test.pdf', bbox_inches='tight')
#--------------------------------------------------------------



plt.show()
#exit()









#------------------------------------------------------------------------------------------












































