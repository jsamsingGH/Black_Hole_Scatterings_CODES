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



#---------------------------------------
#input names of Nbody IC files:
#---------------------------------------
file_name    = 'fig_paper_eccBH.txt'
#file_name    = 'fig_paper_eccBH_NOGR.txt'
#---------------------------------------


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

#open time info:
tf = open('NbodyTides_dataout_a1a2a3.dat',  "r")
NbodyTides_dataout_a1a2a3       = np.loadtxt(tf, dtype=float)
tf.close()
time = NbodyTides_dataout_a1a2a3[:,0]
nr_sim_steps    = len(time)
max_time        = 1.0*max(time)

nr_eqt_steps    = 1000
dt_eqt          = max_time/(1.0*nr_eqt_steps)
eqt_arr_time    = np.arange(0.0, max_time, dt_eqt)

#calc info:
sma_a0      = np.sqrt(sum((b1_posxyz[0,:]-b2_posxyz[0,:])**2.))
print sma_a0/AU_U
#---------------------------------------

#---------------------------------------
#Plot:
#---------------------------------------
fc = 0

for tc in range(0,nr_sim_steps):
    sim_t   = time[tc]
    eq_t    = eqt_arr_time[fc+1]
    
    if (sim_t > eq_t):
        #update frame counter:
        fc = fc+1
        print fc,nr_eqt_steps
        
        #make plot window:        
        fig = plt.figure(figsize=(6, 5))

        plt.tick_params(
            axis='both',        # changes apply to the x,y-axis
            which='both',       # both major and minor ticks are affected
            bottom='off',       # ticks along the bottom edge are off
            top='off',          # ticks along the top edge are off
            labelbottom='off',  # labels along the bottom edge are off
            right='off',
            left='off',
            labelleft='off')

        #object 1:
        xp1 = b1_posxyz[0:tc,0]/sma_a0
        yp1 = b1_posxyz[0:tc,1]/sma_a0
        zp1 = b1_posxyz[0:tc,2]/sma_a0
        fig.add_subplot(111).plot(xp1, yp1, linewidth=0.25, linestyle='-', color='black')
        fig.add_subplot(111).scatter(xp1[tc-1], yp1[tc-1], linewidth=0.01, c='black')

        #object 2:
        xp2 = b2_posxyz[0:tc,0]/sma_a0
        yp2 = b2_posxyz[0:tc,1]/sma_a0
        zp2 = b2_posxyz[0:tc,2]/sma_a0
        fig.add_subplot(111).plot(xp2, yp2, linewidth=0.25, linestyle='-', color='black')
        fig.add_subplot(111).scatter(xp2[tc-1], yp2[tc-1], linewidth=0.01, c='black')
        
        #object 3:
        xp3 = b3_posxyz[0:tc,0]/sma_a0
        yp3 = b3_posxyz[0:tc,1]/sma_a0
        zp3 = b3_posxyz[0:tc,2]/sma_a0
        fig.add_subplot(111).plot(xp3, yp3, linewidth=0.25, linestyle='-', color='black')
        fig.add_subplot(111).scatter(xp3[tc-1], yp3[tc-1], linewidth=0.01, c='black')
        #---------------------------------------
        #axis settings/labels/ranges, etc:
        #---------------------------------------
        ax = fig.add_subplot(111)
        ax.set_xlim(-11, 11)
        ax.set_ylim(-4, 5)
        #ax.set_xlim(-20, 20)
        #ax.set_ylim(-20, 20)


        #-----------------------------------------------------------------    
        #Finalize plot and save:
        #-----------------------------------------------------------------
        fig_name = 'frame_' + str(fc) + '.jpg'
        plt.savefig('MOVIE_FRAMES_2/' + fig_name)    
        plt.close()        
        #-----------------------------------------------------------------    
        
        
        
    
    
    
    
exit()












#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#make plot window:
fig = plt.figure(figsize=(6, 5))

#object 1:
xp1 = b1_posxyz[:,0]/sma_a0
yp1 = b1_posxyz[:,1]/sma_a0
zp1 = b1_posxyz[:,2]/sma_a0
fig.add_subplot(111).plot(xp1, yp1, linewidth=0.5, linestyle='-', color='black')
#object 2:
xp2 = b2_posxyz[:,0]/sma_a0
yp2 = b2_posxyz[:,1]/sma_a0
zp2 = b2_posxyz[:,2]/sma_a0
fig.add_subplot(111).plot(xp2, yp2, linewidth=0.5, linestyle='-', color='black')
#object 3:
xp3 = b3_posxyz[:,0]/sma_a0
yp3 = b3_posxyz[:,1]/sma_a0
zp3 = b3_posxyz[:,2]/sma_a0
fig.add_subplot(111).plot(xp3, yp3, linewidth=0.5, linestyle='-', color='black')
#---------------------------------------
#axis settings/labels/ranges, etc:
#---------------------------------------
ax = fig.add_subplot(111)
ax.set_xlabel(r'$x/a$')
ax.set_ylabel(r'$y/a$')
ax.set_xlim(-11, 11)
ax.set_ylim(-4, 5)
#ax.set_xlim(-20, 20)
#ax.set_ylim(-20, 20)

ax.text(-10.0, 1.9, r'$\rm{Binary}$ $\rightarrow$', size = 20,
        horizontalalignment='left',
       verticalalignment='center',
        rotation = -22.5)

ax.text(10.5, -0.2, r'$\leftarrow$ $\rm{Single}$', size = 20,
        horizontalalignment='right',
        verticalalignment='center',
        rotation = - 21.0)

#ax.text(-3.5, -1.5, r'$\rm{GW\ capture}$', size = 15,
#        horizontalalignment='center', fontname='serif',
#        verticalalignment='center',
#        rotation = 0.0)

#-----------------------------------------------------------------    
#Finalize plot and save:
#-----------------------------------------------------------------
#show:    
plt.show()
#save pdf:
#fig.savefig('orbit3body_ex_1.eps', bbox_inches='tight')
fig.savefig('orbit3body_ex_1_NOGR.eps', bbox_inches='tight')
#-----------------------------------------------------------------    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

