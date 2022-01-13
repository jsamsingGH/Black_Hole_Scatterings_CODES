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
from matplotlib.patches import Ellipse


#input names of Nbody IC files:
#list_Nbody_IC_file_names    = ['quick3body.txt', 'quick3body.txt']#, 'quick3body.txt', 'quick3body.txt', 'quick3body.txt', 'quick3body.txt']  #we try optimize the plot for 6 plots.
#list_Nbody_IC_file_names    = ['TEST_t0.txt', 'TEST_t1.txt', 'TEST_t2.txt', 'TEST_t3.txt']
list_Nbody_IC_file_names    = ['TEST_t0.txt', 'TEST_t1.txt', 'TEST_t3.txt']
#list_Nbody_IC_file_names    = ['t1_ex_AT.txt', 't1_ex_DTn4.txt', 't1_ex_DTn10.txt']



#['paper1_ex_tidalinsp.txt', 'paper1_ex_tidalinsp_FTM.txt', 'paper1_ex_tidalinspNT.txt'] #['MCinput_Nbody.txt'] #['MCinput_Nbody_102020BH_paper2.txt']
#['GWinsp_14WD_E1.txt','GWinsp_14WD_E1.txt']         #PAPER1 :['paper1_ex_tidalinsp.txt', 'paper1_ex_tidalinspNT.txt']
#['MCinput_Nbody.txt', 'MCinput_Nbody.txt']#['GWinsp_14WD_E1.txt', 'GWinsp_04WD_E10.txt']#['insp_ill_4_WT_inputparams.txt', 'insp_ill_4_NT_inputparams.txt']
list_plot_annotations       = [r'No Tides', r'AFT($-$damping)', r'DFT($+$damping, $n=10$)']
#list_plot_annotations       = [r'AFT($-$damping)', r'DFT($+$damping, $n=4$)', r'DFT($+$damping, $n=10$)']

nr_file_names               = len(list_Nbody_IC_file_names)

#CHECK: GWinsp_Tinsp_longer_Tiso.txt
#paper 2: MCinput_Nbody_102020BH_paper2.txt











#-------------------------------------------------------------------------------------------------------
#3-body tidal capture example:
#-------------------------------------------------------------------------------------------------------
#make plot window:
fig = plt.figure(figsize=(6, 12))
nrfigx = 3
nrfigy = 1
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)
#-----------------------------------------------------------------
#loop over file-names (figures we sim and plot):
#-----------------------------------------------------------------
for i in range(0,nr_file_names):
    #---------------------------------------
    #input file name:
    #---------------------------------------
    file_name = list_Nbody_IC_file_names[i]
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
    #get time:
    tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
    NbodyTides_dataout_a1a2a3       = np.loadtxt(tf, dtype=float)
    tf.close()
    time = NbodyTides_dataout_a1a2a3[:,0]
    nr_sim_steps    = len(time)
    #calc characteristic vel:
    r0      = np.sqrt(sum((b1_posxyz[0,:]-b2_posxyz[0,:])**2))
    velfac  = 0.0/r0   
    print velfac
    #---------------------------------------
    #Plot:
    #---------------------------------------
    colornr=[25,100,210]    
    #object 1:
    xp1 = b1_posxyz[:,0] + velfac*time
    yp1 = b1_posxyz[:,1]
    zp1 = b1_posxyz[:,2]
    #if (i==0):  fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp1, yp1, linewidth=0.75, color=plt.cm.gist_earth(colornr[0]), label=r'WD ($1M_{\odot}, \gamma = 5/3, n=1.5$)')
    #if (i==1):  fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp1, yp1, linewidth=0.75, color=plt.cm.gist_earth(colornr[0]), label=r'WD ($1M_{\odot}$)')
    fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp1, yp1, linewidth=0.5, linestyle='-', color='black', label=r'WD (0.6M$_{\odot}$)')
    #object 2:
    xp2 = b2_posxyz[:,0] + velfac*time
    yp2 = b2_posxyz[:,1]
    zp2 = b2_posxyz[:,2]
    #fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp2, yp2, linewidth=0.75, color=plt.cm.gist_earth(colornr[1]), label=r'CO ($1M_{\odot}$)')
    fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp2, yp2, linewidth=0.5, linestyle='-', color='purple', label=r'CO (1M$_{\odot}$)')
    #object 3:
    xp3 = b3_posxyz[:,0] + velfac*time
    yp3 = b3_posxyz[:,1]
    zp3 = b3_posxyz[:,2]
    #fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp3, yp3, linewidth=0.75, color=plt.cm.gist_earth(colornr[2]), label=r'CO ($1M_{\odot}$)')
    fig.add_subplot(nrfigx,nrfigy,i+1).plot(xp3, yp3, linewidth=0.5, linestyle='-', color='orange', label=r'CO (1M$_{\odot}$)')
    #---------------------------------------
    #axis settings/labels/ranges, etc:
    #---------------------------------------
    #plt.tick_params(
    #    axis='both',        # changes apply to the x,y-axis
    #    which='both',       # both major and minor ticks are affected
    #    bottom='off',       # ticks along the bottom edge are off
    #    top='off',          # ticks along the top edge are off
    #    labelbottom='off',  # labels along the bottom edge are off
    #    right='off',
    #    left='off',
    #    labelleft='off')
    ax = fig.add_subplot(nrfigx,nrfigy,i+1)   
    #ax.axis('equal')
    ax.set_xlabel(r'pos $x/R_{\odot}$')
    ax.set_ylabel(r'pos $y/R_{\odot}$')
    ax.set_xlim(-35,20)
    ax.set_ylim(-25,30)
    #plt.title(list_plot_annotations[i])
    plt.grid()
    
    plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False, ncol=1)
        
    #currentAxis = plt.gca()
    #deltabox_x    = 100.*r0
    #deltabox_y    = 1.*r0
    #currentAxis.add_patch(Rectangle((xp2[nr_sim_steps-1]-deltabox_x/2., yp2[nr_sim_steps-1]-deltabox_y/2.), deltabox_x, deltabox_y, fill=None))
    
    #inset zoomin box:
    if (i > -1):
        box_dxdydz = [3.*0.085, 0.085]
        dc = 0.1
        if (i == 0): subfig = plt.axes([0.2, 2./3.+1.*0.17 - 3.*0.051, box_dxdydz[0], box_dxdydz[1]], axisbg='white')
        if (i == 1): subfig = plt.axes([0.2, 1./3.+1.*0.17 - 2.*0.051, box_dxdydz[0], box_dxdydz[1]], axisbg='white')
        if (i == 2): subfig = plt.axes([0.2, 0.0  +1.*0.17 - 1.*0.051, box_dxdydz[0], box_dxdydz[1]], axisbg='white')
        ##subfig.plot(xp1, yp1, linewidth=0.5, color=plt.cm.gist_earth(colornr[0]))
        ##subfig.plot(xp2, yp2, linewidth=0.5, color=plt.cm.gist_earth(colornr[1]))
        subfig.plot(xp1, yp1, linewidth=0.5, color='black')
        subfig.plot(xp3, yp3, linewidth=0.5, color='orange')
        #test COM:
        m1 = 0.6
        m3 = 1.0
        xp13 = (m1*xp1+m3*xp3)/(m1 + m3) 
        yp13 = (m1*yp1+m3*yp3)/(m1 + m3)
        subfig.plot(xp13, yp13, linewidth=0.75, linestyle = '--', color='green')
        #subfig.set_xlim(-1.5*r0+xp1[nr_sim_steps-1], 0.5*r0+xp1[nr_sim_steps-1])
        #subfig.set_ylim(-0.35*r0+yp1[nr_sim_steps-1], 0.25*r0+yp1[nr_sim_steps-1])
        #print -1.5*r0+xp1[nr_sim_steps-1], 0.5*r0+xp1[nr_sim_steps-1]
        #print -0.35*r0+yp1[nr_sim_steps-1], 0.25*r0+yp1[nr_sim_steps-1]
        subfig.set_xlim(-6.65, -2.30)
        subfig.set_ylim(-7.55, -6.25)        
        #plt.setp(subfig, xticks=[], yticks=[])            
        plt.grid()
        plt.title('zoom box')
        ##text:
        #subfig.plot([0], [0], linewidth=0.0, color='black', label='Zoom in GW inspiral')
        #subfig.legend(loc='upper left', numpoints = 1, fontsize = 7.0, frameon = False)
        im = nr_sim_steps-1    
        if (i == 0): ax.annotate(r'Exchange',       xy=(xp1[im]+1.0*r0, yp1[im]-2.0*r0), fontsize = 10.0, horizontalalignment='center', verticalalignment='center', color='red')
        if (i >= 1): ax.annotate(r'Tidal Capture',  xy=(xp1[im]+1.0*r0, yp1[im]-2.0*r0), fontsize = 10.0, horizontalalignment='center', verticalalignment='center', color='red')
        ax.annotate(list_plot_annotations[i],       xy=(-30.0, 25.0), fontsize = 15.0, horizontalalignment='left', verticalalignment='center', color='blue')
        
        #fig.add_subplot(nrfigx,nrfigy,i+1).plot([xp2[im], xp2[im]+50.*r0], [yp2[im]-0.2*r0, yp2[im]-1.0*r0], linewidth=1.0, color='black')        
        
    
#-----------------------------------------------------------------    
#Finalize plot and save:
#-----------------------------------------------------------------
#plt.subplots_adjust(bottom=0.04, left=0.04, right=1.0-0.04, top=1.0-0.04)
#plt.subplots_adjust(hspace=0.05)
#plt.subplots_adjust(wspace=0.05)
#show:    
plt.show()
#save pdf:
#fig.savefig('compare3body_interactions.eps')#, bbox_inches='tight')
fig.savefig('threebody_ex_DFT.pdf', bbox_inches='tight')
#-----------------------------------------------------------------    
exit()
#-------------------------------------------------------------------------------------------------------    
#-------------------------------------------------------------------------------------------------------    
    
  








#-------------------------------------------------------------------------------------------------------
#2-body tidal capture example:
#-------------------------------------------------------------------------------------------------------
#make plot window:
fig = plt.figure(figsize=(6, 8))
nrfigx = 2
nrfigy = 1
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)
#-----------------------------------------------------------------
#loop over file-names (figures we sim and plot):
#-----------------------------------------------------------------
for i in range(0,nr_file_names):
    #---------------------------------------
    #input file name:
    #---------------------------------------
    file_name = list_Nbody_IC_file_names[i]
    #---------------------------------------
    #run Nbody code with that input file:
    #---------------------------------------
    subprocess.call('./TEST_main_Nbody_AffineTides_solver.exe' + '<' + file_name, shell=True)
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
    
    #define:
    b1_posxyz   = NbodyTides_dataout_pos[:,0:3]
    b2_posxyz   = NbodyTides_dataout_pos[:,3:6]
    time        = NbodyTides_dataout_a1a2a3[:,0]
    r0          = np.sqrt(sum((b1_posxyz[0,:]-b2_posxyz[0,:])**2))
    bin_i       = NbodyTides_dataout_binij_info[:,1]
    bin_j       = NbodyTides_dataout_binij_info[:,2]
    a_bin_ij    = NbodyTides_dataout_binij_info[:,3]
    e_bin_ij    = NbodyTides_dataout_binij_info[:,4]
    rperi_ij    = NbodyTides_dataout_binij_info[:,5]
    T_ij        = NbodyTides_dataout_binij_info[:,6]
    T_ji        = NbodyTides_dataout_binij_info[:,7]
    #---------------------------------------
    #Plot:
    #---------------------------------------
    
    #plot 1:
    
    ax = fig.add_subplot(nrfigx,nrfigy,1)   
    colorname=['grey','dodgerblue','peru']    
    xp = b1_posxyz[:,0] - b2_posxyz[:,0]
    yp = b1_posxyz[:,1] - b2_posxyz[:,1]   
    ax.plot(xp, yp, linewidth=1.0, linestyle='-', color=colorname[i], alpha = 1., label = list_plot_annotations[i])
    ax.set_title('Two-body Tidal Evolution')
    ax.set_xlabel(r'pos $x/{R_{\odot}}$')
    ax.set_ylabel(r'pos $y/{R_{\odot}}$')
    ax.set_xlim(-5,45)
    ax.set_ylim(-15,15)
    plt.grid()
    plt.legend(loc='lower right', numpoints = 1, fontsize = 10.0, frameon = False, ncol=1)
    ax.add_patch(Ellipse((0,0), 1.0, 1.0, color='black'))

    #plot 2:
    
    #first analytical estimate for SMA a(t):
    b1_mass = 1.0
    b1_rad  = 1.0
    b2_mass = 1.4
    #calc PT dE (orbit) from initial rp: (ASSUMES obj1 IS THE ONLY TIDAL OBJECT - but works for any mass ratio)
    rp_ini  = rperi_ij[0]
    eta_imp = ((b1_mass/(b1_mass+b2_mass))**(1./2.))*((rp_ini/b1_rad)**(3./2.))
    log10eta    = np.log10(eta_imp)
    #For n=3.0:
    fA = -1.124
    fB = 0.877
    fC = -13.37
    fD = 21.55
    fE = -16.48
    fF = 4.124
    log10T2         = fA + fB*(log10eta**(1.)) + fC*(log10eta**(2.)) + fD*(log10eta**(3.)) + fE*(log10eta**(4.)) + fF*(log10eta**(5))
    T2_PT           = 10.**(log10T2)
    dEorb_PT_tides  = ((b2_mass**2.)/b1_rad)*((b1_rad/rp_ini)**(6.))*T2_PT
    Mtot12  = b1_mass + b2_mass
    mu12    = b1_mass*b2_mass/Mtot12
    gamma   = (1./2.)*np.sqrt(Mtot12*mu12)*(dEorb_PT_tides/Mtot12)*(1./np.pi)*(mu12**(-3./2.))   #(dEorb_PT_tides/(io_mass**(3./2.)))*(1./(np.pi*np.sqrt(2.)))
    SMA_a_t_analesti = (np.sqrt(a_bin_ij[0])-gamma*time)**2.        
    
    ax = fig.add_subplot(nrfigx,nrfigy,2)   
    a0      = a_bin_ij[0]
    rp0     = rperi_ij[0]
    Torb0   = 2.*np.pi*np.sqrt(a0**(3.)/(b1_mass + b2_mass)) 
    if (i == 2): ax.plot(time/Torb0, SMA_a_t_analesti/a0,  linewidth=2.0, linestyle='-.',  color='black', alpha=0.75)
    ax.plot(time/Torb0, a_bin_ij/a0,          linewidth=2.0, linestyle='-',   color=colorname[i], alpha = 0.75)
    ax.plot(time/Torb0, rperi_ij/rp0,         linewidth=2.0, linestyle='--',  color=colorname[i], alpha = 0.75)
    ax.plot(time/Torb0, e_bin_ij,             linewidth=2.0, linestyle=':',   color=colorname[i], alpha = 0.75)
    #dummy plot for legend:
    if (i == 0):
        ax.plot(0.0, 0.0,   linewidth=2.0, linestyle='-',   color='black', alpha = 0.75, label=r'semi-major axis $a(t)/a_{0}$')
        ax.plot(0.0, 0.0,   linewidth=2.0, linestyle='--',  color='black', alpha = 0.75, label=r'pericenter dist. $r_{\rm p}(t)/r_{\rm p,0}$')
        ax.plot(0.0, 0.0,   linewidth=2.0, linestyle=':',   color='black', alpha = 0.75, label=r'eccentricity $e(t)$')
        ax.plot(0.0, 0.0,   linewidth=2.0, linestyle='-.',  color='black',  alpha = 0.75, label=r'analytical sol. $a(t)/a_{0}$')
    ax.set_xlabel(r'time $t/T_{0}$')
    ax.set_ylabel(r'orbital parameters $a(t)$, $e(t)$, $r_{\rm p}(t)$')
    ax.set_ylim(0.0,1.5)
    plt.grid()
    plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False, ncol=2)
    



#-----------------------------------------------------------------    
#Finalize plot and save:
#-----------------------------------------------------------------
#plt.subplots_adjust(bottom=0.04, left=0.04, right=1.0-0.04, top=1.0-0.04)
#plt.subplots_adjust(hspace=0.05)
#plt.subplots_adjust(wspace=0.05)
#show:    
plt.show()
#save pdf:
#fig.savefig('compare3body_interactions.eps')#, bbox_inches='tight')
fig.savefig('2body_ex_1.pdf')#, bbox_inches='tight')
#-----------------------------------------------------------------
exit()
#-------------------------------------------------------------------------------------------------------
#-------------------------------------------------------------------------------------------------------






    
    
    
    
    

