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

import astropy.units as u
from astropy.cosmology import Planck13, z_at_value

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


#open data:
#filename = 'BBHmergers3_with_age_dist.dat'
#print '1'
#tf = open(filename, "r")
#output_arr      = np.loadtxt(tf, dtype=float)
#tf.close()
#print '2'
#np.savez('BBHmergers3_with_age_dist_NPZ', output_arr)
info_data_arr   = np.load('BBHmergers3_with_age_dist_NPZ' + '.npz')['arr_0']

#0)model 1)rv 2)Rgc 3)Z 4)N(x10^5) 5)t_merge(Myr) 6)id1 7)id2 8)m1 9)m2 10)spin1 11)spin2 12)v_kick 13)v_esc 14)channel 15)id_merger 16)m_merger 17)spin_merger 18)a_final (AU) 19)e_final 20)cluster age (Myr) 21)t_merger_actual(Myr)
#Merger channels: 1)Ejected 2)In-cluster 2-body 3)In-cluster 3-body 4)In-cluster 4-body 5)In-cluster single-single capture


analyze_YN  = 0


#-----------------------------------------------------------------------------
#ANALYZE:
#-----------------------------------------------------------------------------
if (analyze_YN == 1):
    
    #input dataname:
    save_data_filename = raw_input('Save data filename: ')
    
    #define:
    nrm         = len(info_data_arr[:,0])
    save_arr    = np.zeros((nrm, 5), dtype=np.float64)

    for nrc in range(0,nrm):
        print nrc, nrm

        #all info about merger:
        info_m      = info_data_arr[nrc, :]
        #define:
        m_id        = info_m[14]
    
        #select based on ID (to speed up loop for now ...)
        #if (m_id == 3 or m_id == 4 or m_id == 5):   #only GW captures
        if (m_id > -1):                             #all mergers.
    
            #calc z:
            tmerg_Gyr   = info_m[21]*(1e-3)
            z_merg      = z_at_value(Planck13.age,  tmerg_Gyr*u.Gyr) 
        
            #calc ecc at XHz:
            fHz_lim_OF  = 10.0                      #OF = observer frame
            fHz_lim_RF  = fHz_lim_OF*(1.+z_merg)    #RF = rest frame  
            
            #BBH parameters:
            a0          = info_m[18]*AU_SI
            e0          = info_m[19]
            rp0         = a0*(1.-e0)
            mi          = info_m[8]*M_sun_SI
            mj          = info_m[9]*M_sun_SI
            #calc/define:
            c0          = rp0*(1.+e0)*(e0**(-12./19.))*((1.+(121./304.)*e0**2.)**(-870./2299.))
            rpmax       = (c0/2.)*(425./304.)**(870./2299.)
            fmin        = (1./np.pi)*np.sqrt(G_new_SI*(mi+mj)/rpmax**3.)     

            #OBS-frame (OF):
            fGW_limit   = fHz_lim_RF
            eccM_OF     = -1
            #check if ecc or not at fGW_limit:
            if (fmin > fGW_limit):
                eccM_OF     = 1.0
            if (fmin < fGW_limit):        
                func        = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mi+mj)/(((c0*ecc**(12./19.))/(1.+ecc))*((1.+(121./304.)*ecc**2.)**(870./2299.)))**3.)
                ecc_initial_guess   = 1e-8
                ecc_fGW_limit       = fsolve(func, ecc_initial_guess)
                eccM_OF     = ecc_fGW_limit 
    
            #REST-frame (RF):
            fGW_limit   = fHz_lim_OF
            eccM_RF     = -1
            #check if ecc or not at fGW_limit:
            if (fmin > fGW_limit):
                eccM_RF     = 1.0
            if (fmin < fGW_limit):        
                func        = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mi+mj)/(((c0*ecc**(12./19.))/(1.+ecc))*((1.+(121./304.)*ecc**2.)**(870./2299.)))**3.)
                ecc_initial_guess   = 1e-8
                ecc_fGW_limit       = fsolve(func, ecc_initial_guess)
                eccM_RF     = ecc_fGW_limit 
        
            #save info:
            save_arr[nrc,0] = m_id 
            save_arr[nrc,1] = z_merg 
            save_arr[nrc,2] = eccM_OF 
            save_arr[nrc,3] = eccM_RF 
                        
    #save data:
    np.savez(save_data_filename, save_arr)
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------    
    
    

#-----------------------------------------------------------------------------    
#PLOT etc:
#-----------------------------------------------------------------------------
#load data:
save_data_filename  = 'test1'    
save_arr            = np.load(save_data_filename + '.npz')['arr_0']

#merge with other datafile:
#K_data_arr  = np.load('BBHmergers3_with_age_dist_NPZ' + '.npz')['arr_0']
#J_data_arr  = save_arr[:,:]
#KJ_data_arr = np.append(K_data_arr, J_data_arr, axis=1)
#fn = 'KJ_BBHmergers3_with_age_dist.txt'
#text_file = open(fn, "w")
#np.savetxt(text_file, KJ_data_arr,   fmt='%1.10e')
#text_file.close()    
#exit()

#save_arr[nrc,0] = m_id 
#save_arr[nrc,1] = z_merg 
#save_arr[nrc,2] = eccM_OF 
#save_arr[nrc,3] = eccM_RF 

#PLOT 1:
fig = plt.figure(figsize=(5.0, 4.0))    
ax1 = fig.add_subplot(111)    

#define:
pos_id3             = np.where(save_arr[:,0] == 3)[0]
pos_id4             = np.where(save_arr[:,0] == 4)[0]
pos_EMXhz_OF        = np.where(save_arr[:,2] > 0.1)[0]
pos_EMXhz_RF        = np.where(save_arr[:,3] > 0.1)[0]
pos_id3_EMXhz_OF    = list(set(pos_id3).intersection(pos_EMXhz_OF))
pos_id3_EMXhz_RF    = list(set(pos_id3).intersection(pos_EMXhz_RF))
pos_id4_EMXhz_OF    = list(set(pos_id4).intersection(pos_EMXhz_OF))
pos_id4_EMXhz_RF    = list(set(pos_id4).intersection(pos_EMXhz_RF))
z_merg_arr          = save_arr[:,1]

#plot settings:
nr_Hbins = 50
#compare mergers (id3):
#OF:
D_arr               = z_merg_arr[pos_id3_EMXhz_OF]
Hy_A, bin_edges     = np.histogram(D_arr, bins = nr_Hbins, range=[0, 10], density=False)
Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
#RF:
D_arr               = z_merg_arr[pos_id3_EMXhz_RF]
Hy_B, bin_edges     = np.histogram(D_arr, bins = nr_Hbins, range=[0, 10], density=False)
Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
#plot:
ax1.step(Hx, (1.*Hy_A)/(1.*Hy_B))
#compare mergers (id4):
#OF:
D_arr               = z_merg_arr[pos_id4_EMXhz_OF]
Hy_A, bin_edges     = np.histogram(D_arr, bins = nr_Hbins, range=[0, 10], density=False)
Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
#RF:
D_arr               = z_merg_arr[pos_id4_EMXhz_RF]
Hy_B, bin_edges     = np.histogram(D_arr, bins = nr_Hbins, range=[0, 10], density=False)
Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
#plot:
ax1.step(Hx, (1.*Hy_A)/(1.*Hy_B))

#analytical solutions:
ax1.plot(Hx, (1.+Hx)**(-2./3.), linestyle = '--', color = 'black')
ax1.plot(Hx, (1.+Hx)**(-1./3.), linestyle = ':', color = 'black')

plt.show()
    
  
  
    
#PLOT 2:
fig = plt.figure(figsize=(5.0, 4.0))    
ax1 = fig.add_subplot(111)    

#save_arr[nrc,0] = m_id 
#save_arr[nrc,1] = z_merg 
#save_arr[nrc,2] = eccM_OF 
#save_arr[nrc,3] = eccM_RF 

#define:
eccM_OF_arr         = save_arr[:,2]
eccM_RF_arr         = save_arr[:,3]
pos_zGT1            = np.where(save_arr[:,1] > 1.0)[0]
pos_zLT1            = np.where(save_arr[:,1] < 1.0)[0]

#plot settings:
nr_Hbins = 100

#OF:
D_arr               = np.log10(eccM_OF_arr[:])
Hy, bin_edges       = np.histogram(D_arr, bins = nr_Hbins, range=[-8, 0], density=False)
Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
#plot:
ax1.step(Hx, Hy)

#RF:
D_arr               = np.log10(eccM_RF_arr[:])
Hy, bin_edges       = np.histogram(D_arr, bins = nr_Hbins, range=[-8, 0], density=False)
Hx                  = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
#plot:
ax1.step(Hx, Hy)


plt.show()


    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    