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
from pylab import *
from scipy.interpolate import interp1d
import h5py
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.ndimage.filters import gaussian_filter

def UNIT_system(mass):
    
    m_unit  = sum(mass)        # [M_sol]  # The total mass of 1 inital cluster (the mass unit)
    r_unit  = 1                # [pc]     # The core radius of 1 initial cluster (the length unit)
    G       = 0.0043009       # [pc/M_sol*km^2/s^2] # Gravtitational constant
    pc      = 3.0860e+13      # [km]     # 1 pc in km-s
    yr      = 31556926        # [s]      # 1 year in seconds

    t_unit  = sqrt(r_unit**2/(G*m_unit))*pc/yr*10**(-6) # [Myr]
    v_unit  = sqrt(G*m_unit/r_unit)                     # [km/s]

    print('The length unit is   ', r_unit,'       pc')
    print('The mass unit is     ', around(m_unit,decimals=1),'  M_sol')
    print('The time unit is     ',around(t_unit,decimals=4),'  Myr')
    print('The velocity unit is ',around(v_unit,decimals=2),'    km/s')
    
    return m_unit, r_unit, t_unit, v_unit


mpl.rcParams['font.family'] = 'serif'
mpl.rcParams['mathtext.fontset'] = 'cm'
mpl.rcParams['mathtext.rm'] = 'serif'



#OPEN DATA:
#[u'Mass', u'NAM', u'V1', u'V2', u'V3', u'X1', u'X2', u'X3', u'tg']
filefolder  = '/Users/jsamsing/Desktop/TIDES_PROJ/ROT_cluster_AB/'
#fname       = 'M64.6.2.h5part'
fname       = 'M64.0.2.h5part'

read_analyze_data_YN    = 0



nres_Lr     = 50
nres_cost   = 50
nres_massb  = 3*3
#--------------------------------------------------------------
#READ DATA ETC:
#--------------------------------------------------------------
if (read_analyze_data_YN == 1):

    filename    = filefolder + fname
    f = h5py.File(filename, 'r')
    list_fnames = f.keys()
    nr_files    = len(list_fnames)
    #define empty lists:
    m_arr   = []
    x_arr   = []
    y_arr   = []
    z_arr   = []
    vx_arr  = []
    vy_arr  = []
    vz_arr  = []
    #open files and save in lists:
    for fc in range(0,nr_files):
        #read data:
        m_arr_fc    = f[list(f.keys())[fc]]['Mass'][:]
        vx_arr_fc   = f[list(f.keys())[fc]]['V1'][:]
        vy_arr_fc   = f[list(f.keys())[fc]]['V2'][:]
        vz_arr_fc   = f[list(f.keys())[fc]]['V3'][:]
        x_arr_fc    = f[list(f.keys())[fc]]['X1'][:]
        y_arr_fc    = f[list(f.keys())[fc]]['X2'][:]
        z_arr_fc    = f[list(f.keys())[fc]]['X3'][:]
        #append data to lists:    
        m_arr.append(m_arr_fc)
        vx_arr.append(vx_arr_fc)
        vy_arr.append(vy_arr_fc)
        vz_arr.append(vz_arr_fc)
        x_arr.append(x_arr_fc)
        y_arr.append(y_arr_fc)
        z_arr.append(z_arr_fc)
    #flatten lists:
    m_arr   = [item for sublist in m_arr for item in sublist]
    vx_arr  = [item for sublist in vx_arr for item in sublist]
    vy_arr  = [item for sublist in vy_arr for item in sublist]
    vz_arr  = [item for sublist in vz_arr for item in sublist]
    x_arr   = [item for sublist in x_arr for item in sublist]
    y_arr   = [item for sublist in y_arr for item in sublist]
    z_arr   = [item for sublist in z_arr for item in sublist]
    #turn into arrays:
    m_arr   = np.asarray(m_arr)
    vx_arr  = np.asarray(vx_arr)
    vy_arr  = np.asarray(vy_arr)
    vz_arr  = np.asarray(vz_arr)
    x_arr   = np.asarray(x_arr)
    y_arr   = np.asarray(y_arr)
    z_arr   = np.asarray(z_arr)

    #sort all arrays according to 'r':
    r_arr   = (x_arr**2. + y_arr**2. + z_arr**2.)**(1./2.)
    pos     = np.argsort(r_arr)
    m_arr   = m_arr[pos]
    x_arr   = x_arr[pos]
    y_arr   = y_arr[pos]
    z_arr   = z_arr[pos]
    vx_arr  = vx_arr[pos]
    vy_arr  = vy_arr[pos]
    vz_arr  = vz_arr[pos]

    nr_obj  = len(m_arr)
    Mtot    = sum(m_arr)

    print nr_obj, 'nr_obj'

    #coordinate: r
    r_arr       = (x_arr**2. + y_arr**2. + z_arr**2.)**(1./2.)
    #coordinate: phi: WE DONT USE THIS.
    phi_arr     = np.zeros(nr_obj, dtype=np.float64)
    pos_xG0     = np.where(x_arr[:] > 0)[0] 
    pos_xL0     = np.where(x_arr[:] < 0)[0] 
    pos_yG0     = np.where(y_arr[:] > 0)[0] 
    pos_yL0     = np.where(y_arr[:] < 0)[0] 
    pos_xG0yG0  = list(set(pos_xG0).intersection(pos_yG0))
    pos_xL0yG0  = list(set(pos_xL0).intersection(pos_yG0))
    pos_xL0yL0  = list(set(pos_xL0).intersection(pos_yL0))
    pos_xG0yL0  = list(set(pos_xG0).intersection(pos_yL0))
    phi_arr[pos_xG0yG0] = 0.*np.pi  + np.arctan(abs(y_arr[pos_xG0yG0]/x_arr[pos_xG0yG0]))
    phi_arr[pos_xL0yG0] = 1.*np.pi  - np.arctan(abs(y_arr[pos_xL0yG0]/x_arr[pos_xL0yG0]))
    phi_arr[pos_xL0yL0] = 1.*np.pi  + np.arctan(abs(y_arr[pos_xL0yL0]/x_arr[pos_xL0yL0]))
    phi_arr[pos_xG0yL0] = 2.*np.pi  - np.arctan(abs(y_arr[pos_xG0yL0]/x_arr[pos_xG0yL0]))
    #coordinate: theta
    theta_arr       = np.arccos(z_arr/r_arr)

    #define:
    costheta_arr    = np.cos(theta_arr)
    logm_arr        = np.log10(m_arr)

    #nr_Hbins        = nres_massb
    #Hdata           = logm_arr
    #hist, bin_edges = np.histogram(Hdata, bins = nr_Hbins, range=[min(logm_arr),max(logm_arr)], density=False)
    #bin_centers     = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    #saveoutput_arr  = [bin_centers, np.log10(hist)]
    #np.savez(filefolder + fname + '-' + '_mbdist', saveoutput_arr)
    #exit()

    #find v_r_arr, v_theta_arr, v_phi_arr:
    v_r_arr     = (x_arr*vx_arr + y_arr*vy_arr + z_arr*vz_arr)/r_arr
    v_theta_arr = (-r_arr/np.sqrt(1.-(z_arr/r_arr)**2.))*((vz_arr/r_arr) - (z_arr*v_r_arr/(r_arr**2.))) 
    v_phi_arr   = (np.sqrt(x_arr**2. + y_arr**2)/((y_arr/x_arr)**2. + 1.))*((vy_arr/x_arr) - (y_arr*vx_arr/(x_arr**2.)))

    #find lagrangian-radii:
    m_arr_CS        = np.cumsum(m_arr)
    m_arr_CS_NM     = m_arr_CS/Mtot   
    Lr_arr          = m_arr_CS_NM       #from 0 - 1 (in mass fraction)
    #find nr-particle-radii:
    PN_arr_CS       = np.cumsum(np.ones(nr_obj))
    PN_arr_CS_NM    = PN_arr_CS/(1.*nr_obj)   
    PNr_arr         = PN_arr_CS_NM       #from 0 - 1 (in particle nr (PN) fraction)
    #MAKE THIS WAY OF CHOOSING BETTER!!! easy
    Lr_arr = PNr_arr
    #make sure last entry is = 1 (otherwise interpol complains)
    Lr_arr[nr_obj-1] = 1

    #make interpolation function: go from Sr (scaled r) to r (real r):
    f_intp_Sr_r = interp1d(np.insert(Lr_arr,0,0), np.insert(r_arr,0,0))

    #for hver af disse beregn: v_i, sigma_i (for i - r,t,p)
    Lr_ct_info          = np.zeros((nres_Lr, nres_cost, 10), dtype=np.float64)
    Lr_ct_m_bin_info    = np.zeros((nres_Lr, nres_cost, nres_massb, 10), dtype=np.float64)
    Lr_arr_edges        = np.linspace(0,  1, num=nres_Lr+1)
    cost_arr_edges      = np.linspace(-1, 1, num=nres_cost+1)
    logmass_arr_edges   = np.linspace(min(logm_arr), max(logm_arr), num=nres_massb+1) 
    r_arr_edges         = f_intp_Sr_r(Lr_arr_edges)

    for rc in range(0,nres_Lr):
        print rc, nres_Lr
        pos_r   = np.where((Lr_arr > Lr_arr_edges[rc+0]) & (Lr_arr < Lr_arr_edges[rc+1]))[0]
        for tc in range(0,nres_cost):
            pos_ct      = np.where((costheta_arr > cost_arr_edges[tc + 0]) & (costheta_arr < cost_arr_edges[tc + 1]))[0]
            pos_rct     = list(set(pos_r).intersection(pos_ct))
            dV_bin      = abs((2.*np.pi/3.)*(r_arr_edges[rc+1]**3. - r_arr_edges[rc+0]**3.)*(cost_arr_edges[tc+1] - cost_arr_edges[tc+0]))
            meanM_dV    = np.mean(m_arr[pos_rct]) 
            Lr_ct_info[rc, tc, 0] = meanM_dV
            for mc in range(0,nres_massb):
                pos_m       = np.where((logm_arr > logmass_arr_edges[mc+0]) & (logm_arr < logmass_arr_edges[mc+1]))[0]
                pos_rct_m   = list(set(pos_rct).intersection(pos_m))
                #calc:
                v_r_arr_rct_m       = v_r_arr[pos_rct_m]
                v_theta_arr_rct_m   = v_theta_arr[pos_rct_m]
                v_phi_arr_rct_m     = v_phi_arr[pos_rct_m]
                #vr:
                vr_mean     = np.mean(v_r_arr_rct_m)                #mean(x)
                vr_sigma    = np.sqrt(np.var(v_r_arr_rct_m))        #RMS(x - mean(x))
                #vt:
                vt_mean     = np.mean(v_theta_arr_rct_m)            #mean(x)
                vt_sigma    = np.sqrt(np.var(v_theta_arr_rct_m))    #RMS(x - mean(x))
                #vp:
                vp_mean     = np.mean(v_phi_arr_rct_m)              #mean(x)
                vp_sigma    = np.sqrt(np.var(v_phi_arr_rct_m))      #RMS(x - mean(x))
                #define:
                nr_obj_dV   = 1.*len(pos_rct_m)
                nrn_dV      = nr_obj_dV/dV_bin
                #save:
                Lr_ct_m_bin_info[rc,tc,mc, 0]   = nr_obj_dV   
                Lr_ct_m_bin_info[rc,tc,mc, 1]   = nrn_dV
                #Lr_ct_m_bin_info[rc,tc,mc, 2]   = ...  
                Lr_ct_m_bin_info[rc,tc,mc, 3:5] = [vr_mean,vr_sigma]    
                Lr_ct_m_bin_info[rc,tc,mc, 5:7] = [vt_mean,vt_sigma]    
                Lr_ct_m_bin_info[rc,tc,mc, 7:9] = [vp_mean,vp_sigma]
            
                                
                                
    #make encounter matrix:
    Lr_ct_enc_arr       = np.zeros((nres_Lr, nres_cost, nres_massb, nres_massb, 5), dtype=np.float64)
    mass_arr_logcenters = 10.**((logmass_arr_edges[1:nres_massb+1] + logmass_arr_edges[0:nres_massb])/2.)

    for rc in range(0,nres_Lr):
        for tc in range(0,nres_cost):
            for mc_t in range(0,nres_massb):
                for mc_i in range(0,nres_massb):
                
                    m_t = mass_arr_logcenters[mc_t]
                    m_i = mass_arr_logcenters[mc_i]
                
                    nr_obj_dV_t             = Lr_ct_m_bin_info[rc,tc,mc_t, 0]
                    nrn_dV_t                = Lr_ct_m_bin_info[rc,tc,mc_t, 1]
                    [vr_mean_t,vr_sigma_t]  = Lr_ct_m_bin_info[rc,tc,mc_t, 3:5]
                    [vt_mean_t,vt_sigma_t]  = Lr_ct_m_bin_info[rc,tc,mc_t, 5:7]   
                    [vp_mean_t,vp_sigma_t]  = Lr_ct_m_bin_info[rc,tc,mc_t, 7:9]
                
                    nr_obj_dV_i             = Lr_ct_m_bin_info[rc,tc,mc_i, 0]
                    nrn_dV_i                = Lr_ct_m_bin_info[rc,tc,mc_i, 1]
                    [vr_mean_i,vr_sigma_i]  = Lr_ct_m_bin_info[rc,tc,mc_i, 3:5]
                    [vt_mean_i,vt_sigma_i]  = Lr_ct_m_bin_info[rc,tc,mc_i, 5:7]   
                    [vp_mean_i,vp_sigma_i]  = Lr_ct_m_bin_info[rc,tc,mc_i, 7:9]                                    
                
                    if ((nr_obj_dV_t >= 2) & (nr_obj_dV_i >= 2)):
                        v_rel       = ((vr_mean_t - vr_mean_i)**2. +  (vt_mean_t - vt_mean_i)**2. + (vp_mean_t - vp_mean_i)**2. + (vr_sigma_t**2. + vt_sigma_t**2. + vp_sigma_t**2.) + (vr_sigma_i**2. + vt_sigma_i**2. + vp_sigma_i**2.))**(1./2.)
                        sigma_enc   = (m_t+m_i)/(v_rel**2.)
                        sigma_GWc   = (m_t**(2./7.))*(m_i**(2./7.))*((m_t+m_i)**(10./7.))/(v_rel**(18./7.))
                        Gamma_enc   = nr_obj_dV_t*(nrn_dV_i*sigma_enc*v_rel) #correct up to a factor of order unity
                        Gamma_GWcap = nr_obj_dV_t*(nrn_dV_i*sigma_GWc*v_rel) #correct up to a factor of order unity
                        flag_param1 = 1
                                
                    if ((nr_obj_dV_t <= 1) or (nr_obj_dV_i <= 1)):
                        Gamma_enc   = 0.0
                        Gamma_GWcap = 0.0
                        flag_param1 = -1
                     
                    Lr_ct_enc_arr[rc, tc, mc_t, mc_i, 0]    = Gamma_enc  
                    Lr_ct_enc_arr[rc, tc, mc_t, mc_i, 1]    = Gamma_GWcap  


    Gammatot_arr    = np.zeros((nres_massb, nres_massb, 5), dtype=np.float64)
    for mc_t in range(0,nres_massb):
        for mc_i in range(0,nres_massb):
            Gammatot_arr[mc_t, mc_i,0] = np.sum(Lr_ct_enc_arr[:,:,mc_t, mc_i, 0])
            Gammatot_arr[mc_t, mc_i,1] = np.sum(Lr_ct_enc_arr[:,:,mc_t, mc_i, 1])
        
    #save results:
    saveoutput_arr  = Lr_ct_info
    np.savez(filefolder + fname + '-' + 'Lr_ct_info', saveoutput_arr)
    
    saveoutput_arr  = Lr_ct_m_bin_info  
    np.savez(filefolder + fname + '-' + 'Lr_ct_m_bin_info', saveoutput_arr)
    
    saveoutput_arr  = Lr_ct_enc_arr  
    np.savez(filefolder + fname + '-' + 'Lr_ct_enc_arr', saveoutput_arr)

    saveoutput_arr  = Gammatot_arr  
    np.savez(filefolder + fname + '-' + 'Gammatot_arr', saveoutput_arr)
    
    #Lr_arr_edges, cost_arr_edges, logmass_arr_edges, r_arr_edges...
    
    
    #PLOT TEST/CHECK:
    fig = plt.figure(figsize=(4, 4))
    ax      = fig.add_subplot(111)
    data    = np.log10(Gammatot_arr[:,:,0])
    ax.imshow(data.T, interpolation='none', origin='low', extent=[logmass_arr_edges[0], logmass_arr_edges[-1], logmass_arr_edges[0], logmass_arr_edges[-1]], cmap = 'rainbow')
    plt.show()    

    fig = plt.figure(figsize=(4, 4))
    ax      = fig.add_subplot(111)
    data    = np.log10(Gammatot_arr[:,:,1])
    ax.imshow(data.T, interpolation='none', origin='low', extent=[logmass_arr_edges[0], logmass_arr_edges[-1], logmass_arr_edges[0], logmass_arr_edges[-1]], cmap = 'rainbow')
    plt.show()    

    nr_Hbins        = 100
    Hdata           = logm_arr
    hist, bin_edges = np.histogram(Hdata, bins = nr_Hbins, range=[-6,-3], density=False)
    bin_centers     = (bin_edges[0:nr_Hbins]+bin_edges[1:nr_Hbins+1])/2.
    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    ax1.step(bin_centers, np.log10(hist))
    plt.show()    

    fig = plt.figure(figsize=(7, 7))
    ngrid = np.sqrt(nres_massb)
    for i in range(0,nres_massb):
        ax      = fig.add_subplot(ngrid,ngrid,i+1)
        n_den   = Lr_ct_m_bin_info[0:nres_Lr-1,:,i,1]
        data    = np.log10(n_den)
        ax.imshow(data.T, interpolation='none', origin='low',
            extent=[Lr_arr_edges[0], Lr_arr_edges[-1], cost_arr_edges[0], cost_arr_edges[-1]], cmap = 'rainbow', aspect='auto')
    plt.show()    
    
    fig = plt.figure(figsize=(7, 7))
    for i in range(0,nres_massb):
        ax      = fig.add_subplot(ngrid,ngrid,i+1)
        vp_mean = Lr_ct_m_bin_info[:,:,i,7]
        data    = vp_mean
        ax.imshow(data.T, interpolation='none', origin='low',
            extent=[Lr_arr_edges[0], Lr_arr_edges[-1], cost_arr_edges[0], cost_arr_edges[-1]], cmap = 'rainbow', aspect='auto')
    plt.show()    

    fig = plt.figure(figsize=(5.0, 4.0))
    ax1 = fig.add_subplot(111)
    ax1.plot(Lr_arr, r_arr)
    ax1.plot(PNr_arr, r_arr)
    ax1.set_yscale('log')
    plt.show()    





    exit()    
#--------------------------------------------------------------



#--------------------------------------------------------------
#ANALYZE DATA sets:
#--------------------------------------------------------------
dataset_arr = ['M64.6.2.h5part', 'M64.0.2.h5part']
nr_dataset  = len(dataset_arr)

Dset_Lr_ct_info         = []
Dset_Lr_ct_m_bin_info   = []
Dset_Lr_ct_enc_arr      = []
Dset_Gammatot_arr       = []
Dset_mbdist_arr         = []

for i in range(0,nr_dataset):
    #open data:
    Dset_fname              = dataset_arr[i]
    dsi_Lr_ct_info          = np.load(filefolder + Dset_fname + '-' + 'Lr_ct_info' + '.npz')['arr_0']
    dsi_Lr_ct_m_bin_info    = np.load(filefolder + Dset_fname + '-' + 'Lr_ct_m_bin_info' + '.npz')['arr_0']
    dsi_Lr_ct_enc_arr       = np.load(filefolder + Dset_fname + '-' + 'Lr_ct_enc_arr' + '.npz')['arr_0']
    dsi_Gammatot_arr        = np.load(filefolder + Dset_fname + '-' + 'Gammatot_arr' + '.npz')['arr_0']
    dsi_mbdist_arr          = np.load(filefolder + Dset_fname + '-' + '_mbdist' + '.npz')['arr_0']
    
    #save data:
    Dset_Lr_ct_info.append(dsi_Lr_ct_info)
    Dset_Lr_ct_m_bin_info.append(dsi_Lr_ct_m_bin_info)
    Dset_Lr_ct_enc_arr.append(dsi_Lr_ct_enc_arr)
    Dset_Gammatot_arr.append(dsi_Gammatot_arr)
    Dset_mbdist_arr.append(dsi_mbdist_arr)
#define:
nr_rc       = len(Dset_Lr_ct_m_bin_info[0][:,0,0,0])  
nr_tc       = len(Dset_Lr_ct_m_bin_info[0][0,:,0,0])  
nres_massb  = len(Dset_Lr_ct_m_bin_info[0][0,0,:,0])  

rp_arr      = np.linspace(0,  1,    num=nr_rc)
ct_arr      = np.linspace(-1,  1,   num=nr_tc)

#Lr_ct_info[rc, tc, 0]  = meanM_dV

#Lr_ct_m_bin_info[rc,tc,mc, 0]   = nr_obj_dV   
#Lr_ct_m_bin_info[rc,tc,mc, 1]   = nrn_dV
#Lr_ct_m_bin_info[rc,tc,mc, 3:5] = [vr_mean,vr_sigma]    
#Lr_ct_m_bin_info[rc,tc,mc, 5:7] = [vt_mean,vt_sigma]    
#Lr_ct_m_bin_info[rc,tc,mc, 7:9] = [vp_mean,vp_sigma]

#Lr_ct_enc_arr[rc, tc, mc_t, mc_i, 0]    = Gamma_enc  
#Lr_ct_enc_arr[rc, tc, mc_t, mc_i, 1]    = Gamma_GWcap  

#Lr_ct_enc_arr[rc, tc, mc_t, mc_i, 0]    = Gamma_enc  
#Lr_ct_enc_arr[rc, tc, mc_t, mc_i, 1]    = Gamma_GWcap  
#Gammatot_arr[mc_t, mc_i,0] = np.sum(Lr_ct_enc_arr[:,:,mc_t, mc_i, 0])
#Gammatot_arr[mc_t, mc_i,1] = np.sum(Lr_ct_enc_arr[:,:,mc_t, mc_i, 1])






#NOTES:
#CHECK AND FINISH fig. 3 BELOW:
#fig 1,2 SHOULD BE OK! WAIT FOR AKOS TO REPORT BACK!!
#remove notation _dV below...
#-NEED TO FIND OUT WHAT THE RIGHT DISPERSION IS!!! CROSS TERMS?

fig = plt.figure(figsize=(4, 7))
ax  = fig.add_subplot(2,1,1)
pbin    = Dset_mbdist_arr[0][0][:]
pdata   = Dset_mbdist_arr[0][1][:]/Dset_mbdist_arr[1][1][:]
ax.plot(pdata)
plt.show()    
#FIG 3:
fig = plt.figure(figsize=(4, 7))
X       = np.arange(1,nres_massb+2, dtype=float)
Y       = np.arange(1,nres_massb+2, dtype=float)
#fig ENC:
ax      = fig.add_subplot(2,1,1)
s1      = Dset_Gammatot_arr[0][:,:,0]
s2      = Dset_Gammatot_arr[1][:,:,0]
Zdata   = 100*((s1-s2)/s2)
Z       = Zdata.T
quad    = plt.pcolormesh(X, Y, Z, cmap='RdBu', linewidth=0, rasterized=False, vmin=-10, vmax=10)
fig.colorbar(quad, ax=ax)
#fig GWC:
ax      = fig.add_subplot(2,1,2)
s1      = Dset_Gammatot_arr[0][:,:,1]
s2      = Dset_Gammatot_arr[1][:,:,1]
Zdata   = 100*((s1-s2)/s2)
Z       = Zdata.T
quad    = plt.pcolormesh(X, Y, Z, cmap='RdBu', linewidth=0, rasterized=False, vmin=-10, vmax=10)
fig.colorbar(quad, ax=ax)
plt.savefig('test.pdf', rasterized=False, bbox_inches='tight')
plt.show()    
exit()


#FIG 1:
fig = plt.figure(figsize=(4, 7))
X       = rp_arr 
Y       = ct_arr
dataset = Dset_Lr_ct_m_bin_info
mbc1    = 0
mbc2    = 6
#loop over dataset:
for dc in range(0,2):
    #plot settings:
    ax      = fig.add_subplot(2,1,dc+1)
    
    #plot 1:
    s1      = dataset[dc][:,:,mbc1,1]
    Zdata   = np.log10(s1/np.nanmedian(s1))
    Zdata   = gaussian_filter(Zdata, 2)
    Z       = Zdata.T
    quad    = plt.pcolormesh(X, Y, Z, cmap='gnuplot', linewidth=0, rasterized=False, vmin=-1, vmax=1)
    CS      = plt.contour(X, Y, Z, list([np.linspace(-3, 1, num=20)]), colors='black', linestyles='solid')
    ax.plot(0,0, color = 'black', linestyle='solid', label = r'light')
    
    #plot 2:
    s2      = dataset[dc][:,:,mbc2,1]
    Zdata   = np.log10(s2/np.nanmedian(s2))
    Zdata   = gaussian_filter(Zdata, 2)
    Z       = Zdata.T
    CS      = plt.contour(X, Y, Z, list([np.linspace(-3, 2, num=20)]), colors='limegreen', linestyles='solid')
    ax.plot(0,0, color = 'limegreen', linestyle='solid', label = r'heavy')
    
    #plot axis etc:
    ax.set_xlim(0.0, 0.75)
    ax.set_ylim(-1, 1)
    ax.set_xlabel(r'Radius $R_s$')
    ax.set_ylabel(r'cos($\theta$)')
    if (dc == 0):
        ax.text(0.7, 0.8, 'Rotating',       color = 'white', fontsize=12, horizontalalignment = 'right')
        ax.set_title(r'Density $n_{i}$')
    if (dc == 1):
        ax.text(0.7, 0.8, 'Non-Rotating',   color = 'white', fontsize=12, horizontalalignment = 'right')

#legend, etc.,...:
ax.legend(loc='lower right', numpoints = 1, fontsize = 10, frameon = True, facecolor = 'white', ncol=1, framealpha = 0.85)
plt.savefig('ill_ndist.pdf', rasterized=False, bbox_inches='tight')
plt.show()    
exit()



#FIG 2:
fig = plt.figure(figsize=(4, 7))
X       = rp_arr 
Y       = ct_arr
dataset = Dset_Lr_ct_m_bin_info
mbc1    = 0
mbc2    = 6
#loop over dataset:
for dc in range(0,2):
    #plot settings:
    ax      = fig.add_subplot(2,1,dc+1)
    
    #plot 1:
    Nnum_dV = dataset[dc][:,:,mbc1,0] 
    nden_dV = dataset[dc][:,:,mbc1,1]
    vdis_dV = np.sqrt(dataset[dc][:,:,mbc1,4]**2. + dataset[dc][:,:,mbc1,6]**2. + dataset[dc][:,:,mbc1,8]**2.)
    Rate_dV = Nnum_dV*nden_dV*(vdis_dV**(-11./7.))    
    Zdata   = np.log10(Rate_dV/np.nanmedian(Rate_dV))
    Zdata   = gaussian_filter(Zdata, 2)
    Z       = Zdata.T
    quad    = plt.pcolormesh(X, Y, Z, cmap='gnuplot', linewidth=0, rasterized=False, vmin=-1, vmax=1)
    CS      = plt.contour(X, Y, Z, list([np.linspace(-1, 1, num=10)]), colors='black', linestyles='solid')
    ax.plot(0,0, color = 'black', linestyle='solid', label = r'light')
    print sum(Rate_dV), dc, mbc1
   
    #plot 2:
    Nnum_dV = dataset[dc][:,:,mbc2,0] 
    nden_dV = dataset[dc][:,:,mbc2,1]
    vdis_dV = np.sqrt(dataset[dc][:,:,mbc2,4]**2. + dataset[dc][:,:,mbc2,6]**2. + dataset[dc][:,:,mbc2,8]**2.)
    Rate_dV = Nnum_dV*nden_dV*(vdis_dV**(-11./7.))    
    Zdata   = np.log10(Rate_dV/np.nanmedian(Rate_dV))
    Zdata   = gaussian_filter(Zdata, 2)
    Z       = Zdata.T
    CS      = plt.contour(X, Y, Z, list([np.linspace(-1, 2, num=10)]), colors='limegreen', linestyles='solid')
    ax.plot(0,0, color = 'limegreen', linestyle='solid', label = r'heavy')
    print sum(Rate_dV), dc, mbc2
    
    #plot axis etc:
    if (dc == 0):
        ax.text(0.7, 0.8, 'Rotating',       color = 'white', fontsize=12, horizontalalignment = 'right')
        ax.set_title(r'Rate $\Gamma_{ii}^{GW}$')
    if (dc == 1):
        ax.text(0.7, 0.8, 'Non-Rotating',   color = 'white', fontsize=12, horizontalalignment = 'right')
    ax.set_xlim(0.0, 0.75)
    ax.set_ylim(-1, 1)
    ax.set_xlabel(r'Radius $R_s$')
    ax.set_ylabel(r'cos($\theta$)')
        
#legend, etc.,...:
ax.legend(loc='lower right', numpoints = 1, fontsize = 10, frameon = True, facecolor = 'white', ncol=1, framealpha = 0.85)
plt.savefig('illmRate.pdf', rasterized=False, bbox_inches='tight')
plt.show()    
exit()


#NOTES:
# first row in r,ct: wrong?
# hvorfor er raterne altid en faktor 2 forskellige?

exit()
#--------------------------------------------------------------








#########################################################
# In order to convert data to physical units one should 
# fix the mass and distance units.
# For example let the lightest object have the mass of 0.1 M_sol.
# Then one can rescale all the masses by:

M = 0.1/min(m)*m

# One should also fix the scale of distances. 
# Let's say what is 1 in the data set is 1 pc in physical units.
# So it is just multipling x,y,z by one.

# The function UNIT_system(mass) calculates what is velocity 
# unit in km/s and the time unit in Myr. 

m_unit, r_unit, t_unit, v_unit = UNIT_system(M)

# Then: 

# m = m_unit*m   # M_sol
# x = r_unit*x   # pc
# y = r_unit*y   # pc
# z = r_unit*z   # pc
# vx = v_unit*vx # km/s
# vy = v_unit*vy # km/s
# vz = v_unit*vz # km/s
