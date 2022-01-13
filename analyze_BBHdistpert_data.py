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
from os import listdir
from os.path import isfile, join



#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


#----------------------------------------------------------
#Units and conversions:
#----------------------------------------------------------
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
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
tH_SI       = (10.**10.)*sec_year
#----------------------------------------------------------

file_folder     = '/Users/jsamsing/Desktop/TIDES_PROJ/BBHpert_files/'

#file_folder     = '/Users/jsamsing/Desktop/TIDES_PROJ/BBHpert_files/testfiles/'
#file_folder     = '/Users/jsamsing/Desktop/TIDES_PROJ/BBHpert_files/TC_files/'

#tan4: R=2. tan6: R=50, de=0. tan3: R=50, de=HS. tan5: R=50, de=HR.
#TC_N50000_m20_PNY_v1050n105_rp252_a05e2e_
#datasets_list   = ['T1000_2']#['tan4', 'tan6', 'tan1', 'tan3', 'tan5']#['TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp255_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp2510_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp2525_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp2550_a05e2e_']#['T1000_2', 'T1000_50', 'T1000_HS50', 'TNPN1000_HS50', 'TNPN1000_2']

#plot settings:


#convergence nr incluster vs. ejected plot:
#fig3type = -1
#fig4type = -1
##datasets_list   = ['TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp255_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp2510_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp2525_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp2550_a05e2e_HB']
#datasets_list   = ['TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp22_a05e2e_HS19', 'TC_N5000_m20_PNY_v1050n105_rp255_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp25_a05e2e_HS19', 'TC_N5000_m20_PNY_v1050n105_rp2510_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp210_a05e2e_HS19', 'TC_N5000_m20_PNY_v1050n105_rp2525_a05e2e_HB' ,'TC_N5000_m20_PNY_v1050n105_rp225_a05e2e_HS19', 'TC_N5000_m20_PNY_v1050n105_rp2550_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp250_a05e2e_HS19']
#colors_datasets = ['black', 'dodgerblue', 'black', 'dodgerblue', 'black', 'dodgerblue', 'black', 'dodgerblue', 'black', 'dodgerblue']
#xvalplot        = [1,1,2,2,3,3,4,4,5,5]
#msize           = [10,8,10,8,10,8,10,8,10,8]
#alpha_datasets  = [0.85, 0.75, 0.85, 0.75, 0.85, 0.75, 0.85, 0.75, 0.85, 0.75]

#convergence e>ec plot:
#datasets_list   = ['TC_N5000_m20_PNY_v1050n105_rp2550_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp2525_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp2510_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp255_a05e2e_', 'TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_']
#colors_datasets = ['black', 'blue', 'purple', 'red', 'orange']
#linewi_datasets = [2.0, 2.0, 2.0, 2.0, 2.0]
#alpha_datasets  = [0.5, 0.5, 0.5, 0.5, 0.5]
#linest_datasets = ['-','-','-','-','-']
#legend_datasets = [r'$R/a = 50$', r'$R/a = 25$', r'$R/a = 10$', r'$R/a = 5$', r'$R/a = 2$']


#starting sim:
#hybrid 5000. PNY.
# 252  : 4 : OK
# 255  : 4 : OK
# 2510 : 4 : OK
# 2525 : 4 : OK
# 2550 : 4 : OK 
#hybrid 5000. PNN.
# 252  : 0
#HS19 5000. PNY.
# 250  : 0

#last: dowload them all to make sure we have the complete dataset.

#MULTIPLE DATASETS:
#fig3type = 2
#fig4type = 2
#datasets_list   = ['TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp2550_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp250_a05e2e_HS19']#,'TC_N5000_m20_PNN_v1050n105_rp252_a05e2e_HB']    #['TC_N5000_m20_PNY_v1050n105_rp250_a05e2e_HR96']
#colors_datasets = ['red', 'black', 'dodgerblue', 'forestgreen', 'yellow']
#linewi_datasets = [2.5, 2.5, 2.0, 2.0, 2.5]
#alpha_datasets  = [0.1, 0.1, 0.0, 0.0, 1]
#linest_datasets = ['-','-',':',':','-']
#legend_datasets = [r'$R/a = 2$',r'$R/a = 50$ (hybrid)',r'$R/a = 50$ (HS19)', r'$R/a = 50$ (HR96)', r'$R/a = 2$ (no 2.5PN)', r'$R/a = 50$']

#MULTIPLE DATASETS:
#fig3type = 2
#fig4type = 2
#datasets_list   = ['TC_N5000_m20_PNY_v1050n105_rp252_a05e2e_HB', 'TC_N5000_m20_PNY_v1050n105_rp2550_a05e2e_HB', 'TC_N5000_m20_PNN_v1050n105_rp252_a05e2e_HB']#,'TC_N5000_m20_PNN_v1050n105_rp252_a05e2e_HB']    #['TC_N5000_m20_PNY_v1050n105_rp250_a05e2e_HR96']
#colors_datasets = ['red', 'black', 'orange', 'forestgreen', 'yellow']
#linewi_datasets = [2.5, 2.5, 2.5, 2.0, 2.5]
#alpha_datasets  = [0.1, 0.1, 0.0, 0.0, 1]
#linest_datasets = ['-','-',':',':','-']
#legend_datasets = [r'$R/a = 2$',r'$R/a = 50$ (hybrid)', r'$R/a = 2$ (no 2.5PN)', r'$R/a = 50$']

#ONE DATASET:
fig3type = 1
fig4type = 1
datasets_list   = ['TC_N15000_m20_PNY_v1050n105_rp2100_a05e2e_HS19']
colors_datasets = ['red', 'black', 'black', 'orange', 'yellow']
linewi_datasets = [2.5, 2.5, 2.0, 1.0, 2.5]
alpha_datasets  = [0.1, 0.1, 0.0, 0.1, 1]
linest_datasets = ['-','-',':','-','-']
legend_datasets = [r'$R/a = 2$',r'$R/a = 50$ (hybrid)',r'$R/a = 50$ (HS19)', r'$R/a = 2$ (no 2.5PN)', r'$R/a = 50$']

#input:
mBH     = 20.*M_sun_SI            

#--------------------------------------------------------------
#loop over input data file sets:
#--------------------------------------------------------------
#plot settings:
fig1, ax1 = plt.subplots(figsize=(5, 4))
fig2, ax2 = plt.subplots(figsize=(5, 4))
fig3, ax3 = plt.subplots(figsize=(5, 4))
fig4, ax4 = plt.subplots(figsize=(5, 4))
fig5, ax5 = plt.subplots(figsize=(5, 4))
fig6, ax6 = plt.subplots(figsize=(5, 4))
fig7, ax7 = plt.subplots(figsize=(5, 4))
fig8, ax8 = plt.subplots(figsize=(5, 4))
fig9, ax9 = plt.subplots(figsize=(5, 4))
fig10, ax10 = plt.subplots(figsize=(5, 4))
fig11, ax11 = plt.subplots(figsize=(5, 4))

nr_datasets     = len(datasets_list)

for dsc in range(0,nr_datasets):


    #----------------------------------------------------------
    #create file set etc.:
    #----------------------------------------------------------
    fileSW          = datasets_list[dsc]
    filenames_list  = [f for f in listdir(file_folder) if f.startswith(fileSW)]
    nr_files        = len(filenames_list)
    endstates_arr   = np.zeros((nr_files,20), dtype=np.float64)
    distdata_wfcs   = np.zeros((nr_files*10000,20), dtype=np.float64)
    #----------------------------------------------------------
    
    
    #----------------------------------------------------------         
    #loop through data files:        
    #----------------------------------------------------------        
    #initialize:
    tc  = 0
    for i in range(0,10000):#range(0,nr_files):
    
        #print info:
        print i+1,nr_files


        #read data file 'i':
        filename_i  = filenames_list[i]    
        tf = open(file_folder + filename_i,  "r")
        data = np.loadtxt(tf, dtype=float)
        tf.close()
        #define:
        nl  = len(data)
    
    
        #endstate data from file 'i':
        endstates_arr[i,0:13] = data[nl-1,:]    #[ [out_breakID, out_encID, out_a, out_e, out_rpENC] ]
    
    
        #evolve eccentricity to fGW:
        am      = endstates_arr[i,2]            #[ [out_breakID, out_encID, out_a, out_e, out_rpENC] ]
        em      = endstates_arr[i,3]            #[ [out_breakID, out_encID, out_a, out_e, out_rpENC] ]
        rp      = am*(1.-em)                                    #rp at formation of merger
        fGW0    = (1./np.pi)*np.sqrt(2.*G_new_SI*mBH/(rp**3.))  #fGW at formation of merger
        #initialize:
        ecc_A       = -1
        ecc_B       = -1
        #--------------
        #fGW cut set A:
        #--------------
        fGW_limit   = 10**(-2.0) #fGW CUT: A
        if (fGW0 < fGW_limit):
            #propagating ini a,e to a,e(fGW_limit)
            c0      = am/((em**(12./19.)/(1.-em**2.))*((1.+(121./304.)*em**2.)**(870./2299.)))
            func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mBH+mBH)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 1e-8
            ecc_fGWlimit        = fsolve(func, ecc_initial_guess)   #eccentricity at fGW_limit
            ecc_A               = ecc_fGWlimit      
        endstates_arr[i,13]     = ecc_A  
        #--------------    
        #fGW cut set B: (REMEMBER: too low a cut might not be ok when weak scatterings are included: scatterings can occur at low fGW...)
        #--------------
        fGW_limit   = 10**(-2.0) #fGW CUT: B
        if (fGW0 < fGW_limit):
            #propagating ini a,e to a,e(fGW_limit)
            c0      = am/((em**(12./19.)/(1.-em**2.))*((1.+(121./304.)*em**2.)**(870./2299.)))
            func    = lambda ecc : fGW_limit - (1./np.pi)*np.sqrt(G_new_SI*(mBH+mBH)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 1e-8
            ecc_fGWlimit        = fsolve(func, ecc_initial_guess)   #eccentricity at fGW_limit
            ecc_B               = ecc_fGWlimit      
        endstates_arr[i,14]     = ecc_B                  
        #--------------
        #calc and save t_insp:
        #--------------
        dErp        = (85.*np.pi/12.)*(G_new_SI**(7./2.))*(c_SI**(-5.))*(mBH**(9./2.))*(rp**(-7./2.))
        t_insp_UHt  = ((np.pi*np.sqrt(2.*G_new_SI)*(mBH**(3./2.))*np.sqrt(am)/dErp)/tH_SI)  #I have checked: this is equivalent to Peters high ecc limit with (1+e) = 2. so OK!
        endstates_arr[i,15] = t_insp_UHt     
        #save fGW:
        endstates_arr[i,16] = fGW0     
        #--------------
    
        #format:
        #in code:       [ [a_bin_0, e_bin_0, dtime, t_insp], [a_bin_dt, e_bin_dt], [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN] ]
        #after break:   [ [out_breakID, out_encID, out_a, out_e, out_rpENC] ]
    
        #additional info:
        data_noendline  = data[0:nl-1,:]
        pos_enctype1    = np.where(data_noendline[:,7] == 1)[0] #enctype = 1 = strong encounter.
        if (len(pos_enctype1) >  0): pos_lastenc1    = max(pos_enctype1)
        if (len(pos_enctype1) == 0): pos_lastenc1    = -1
        
        if (data_noendline[pos_lastenc1,12] == 1): endstates_arr[i,17] =  1     
        if (data_noendline[pos_lastenc1,12] != 1): endstates_arr[i,17] = -1     
        
        #read each line 'j' from the file exept the last (endstate) line:
        flag_lastSEmerge_YN = 0
        for j in range(0,nl-1):
            #define/calc:
            a_bin_0         = data[j,0]
            e_bin_0         = data[j,1]
            a_bin_enc       = data[j,9]
            e_bin_enc       = data[j,10]
            BBHmerge_SE_YN  = data[j,12]
            rpENC_Ua        = data[j,6]
            rp0     = a_bin_0*(1. - e_bin_0)
            fGWp0   = (1./np.pi)*np.sqrt(2.*G_new_SI*mBH/(rp0**3.))
            t_int   = data[j,2]/tH_SI
            t_insp  = data[j,3]/tH_SI
            wfc     = min([t_int,t_insp])
            #save:
            distdata_wfcs[tc,0]  = fGWp0    
            distdata_wfcs[tc,1]  = wfc
            distdata_wfcs[tc,2]  = endstates_arr[i,0]   
            distdata_wfcs[tc,3]  = e_bin_0
            
            #TEST:
            de = abs(data[j,5] - data[j,10])
            if (de == 0.0):
                distdata_wfcs[tc,4]  = data[j,6]
                distdata_wfcs[tc,5]  = data[j,5]

            #pos flags:
            if (j == pos_lastenc1): distdata_wfcs[tc,6] = 1
            if (j >  pos_lastenc1): distdata_wfcs[tc,6] = 2            
            if (j == pos_lastenc1 and BBHmerge_SE_YN == 0): flag_lastSEmerge_YN = -1
            if (j == pos_lastenc1 and BBHmerge_SE_YN == 1): flag_lastSEmerge_YN = 1

            t_inspc = t_insp/((768./425.)*(1.-e_bin_0**2)**(7./2.))
            ecr_me  = (1.-(t_int/((768./425.)*t_inspc))**(2./7.))**(1./2.)
            

            distdata_wfcs[tc,7]     = a_bin_0
            distdata_wfcs[tc,8]     = a_bin_enc
            distdata_wfcs[tc,9]     = e_bin_enc
            distdata_wfcs[tc,10]    = j-pos_lastenc1
            distdata_wfcs[tc,11]    = BBHmerge_SE_YN
            distdata_wfcs[tc,12]    = flag_lastSEmerge_YN
            distdata_wfcs[tc,13]    = e_bin_0/ecr_me
            distdata_wfcs[tc,14]    = i
            distdata_wfcs[tc,15]    = rpENC_Ua
                            
            
                
            #update counter:
            tc = tc + 1
        
    #trim output arrays etc.:
    distdata_wfcs   = distdata_wfcs[0:tc,:]
    #----------------------------------------------------------


    #NOTES:
    #make a better track of when the analyt model breaks down (put in main code)!!!!

    #----------------------------------------------------------
    #FIGURES:
    #----------------------------------------------------------
    

    #------------------------------------------------
    #sorting endstates data:
    #------------------------------------------------
    #in code:       [ [a_bin_0, e_bin_0, dtime, t_insp], [a_bin_dt, e_bin_dt], [rpENC_Ua, enc_type, endstate_ID, a_bin_enc, e_bin_enc, Lz_sign, BBHmerge_SE_YN] ]
    #after break:   [ [out_breakID, out_encID, out_a, out_e, out_rpENC] ]

    #define:
    out_breakID_arr = endstates_arr[:,0]
    out_a_arr       = endstates_arr[:,2]
    out_e_arr       = endstates_arr[:,3]
    ecc_fGW_A_arr   = endstates_arr[:,13]
    ecc_fGW_B_arr   = endstates_arr[:,14]
    out_tinsp_arr   = endstates_arr[:,15]   #t_insp in units of tH    
    out_fGW_arr     = endstates_arr[:,16]
    lastSEmerge_YN  = endstates_arr[:,17]
    pos1            = np.where(out_breakID_arr[:] == 1)[0]
    pos2            = np.where(out_breakID_arr[:] == 2)[0]
    pos3            = np.where(out_breakID_arr[:] == 3)[0]
    pos5            = np.where(out_breakID_arr[:] == 5)[0]
    pos6            = np.where(out_breakID_arr[:] == 6)[0]
    pos_tLTtH       = np.where(out_tinsp_arr[:] < 1.0)[0]
    pos_lastSEmerge_N       = np.where(lastSEmerge_YN[:] == -1)[0]  
    pos_lastSEmerge_Y       = np.where(lastSEmerge_YN[:] ==  1)[0]  
    pos_ej          = list(pos2)  
    pos_ej_tLTtH    = list(set(pos_ej).intersection(pos_tLTtH))
    pos_2body       = list(pos1) + list(pos5) + list(pos6)
    pos_2b_tLTtH    = list(set(pos_2body).intersection(pos_tLTtH))
    pos_3body       = list(pos3)
    pos_2body_lastSEmerge_N = list(set(pos_2body).intersection(pos_lastSEmerge_N))
    pos_2body_lastSEmerge_Y = list(set(pos_2body).intersection(pos_lastSEmerge_Y))
    #print info:
    print 1.*(len(pos_2body))/(1.*len(pos_ej))
    print 1.*(len(pos_2body))/(1.*len(pos_ej_tLTtH))
    print len(pos_2body), len(pos_3body), len(pos_ej), len(pos_ej_tLTtH) 
    print len(pos5)
    print len(pos3)
    print len(np.where(out_breakID_arr[:] == 10)[0])
    print len(np.where(out_breakID_arr[:] == 11)[0])
    print len(np.where(out_breakID_arr[:] == 12)[0])
    #------------------------------------------------
    #------------------------------------------------
    #fGW distribution:
    #------------------------------------------------
    #fig, ax1 = plt.subplots(figsize=(5, 5))
    Hmin    = -5
    Hmax    = 3   
    Hnrbins = 40
    
    #JUST ONE dataset:
    if (fig3type == 1):
        ax3addfilename = 'figt1_'
        #all:
        pos = pos_ej_tLTtH + pos_2body + pos_3body
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=True,                                 color = 'grey',   linewidth = 1, alpha = 0.025)#hatch = '\\\\\\\\\\')
        #ejected:
        pos = pos_ej_tLTtH
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = 'ejected',             color = 'blue',    linewidth = 1, hatch = '/////////')
        #2-body:        
        pos = pos_2body_lastSEmerge_N
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (subset A)',   color = 'brown',   linewidth = 0.5)#, hatch = '\\\\\\\\')
        pos = pos_2body_lastSEmerge_Y
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (subset B)',   color = 'orange',  linewidth = 0.5)#, hatch = '\\\\\\')
        pos = pos5
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (subset C)',   color = 'yellow',  linewidth = 0.5)#, hatch = '//////')
        pos = pos_2body
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (all)',        color = 'green',   linewidth = 1, hatch = '\\\\\\\\\\\\')
        #3-body:        
        pos = pos_3body
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '3-body (all)',        color = 'red',     linewidth = 1, hatch = '///////////')
        #all:
        pos = pos_ej_tLTtH + pos_2body + pos_3body
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = 'all',                 color = 'black',   linewidth = 1, alpha = 1.0)#hatch = '\\\\\\\\\\')
        #plot settings: 
        yaxfac = 500.0
        ax3.text(-2.15, 0.85*yaxfac, 'LISA', fontsize=10, color = 'black', alpha = 0.5)
        ax3.text(-3.0,  0.72*yaxfac, '/////////////////////', fontsize=10, color = 'black', alpha = 0.5)
        ax3.text(-0.45, 0.85*yaxfac, 'DECIGO', fontsize=10, color = 'black', alpha = 0.5)
        ax3.text(-0.5,  0.72*yaxfac, '/////////////', fontsize=10, color = 'black', alpha = 0.5)
        ax3.text(1.65,  0.85*yaxfac, 'LIGO', fontsize=10, color = 'black', alpha = 0.5)
        ax3.text(1.0,   0.72*yaxfac, '////////////////////', fontsize=10, color = 'black', alpha = 0.5)
        ax3.set_title(r'peak GW freq. at formation')
        ax3.set_xlabel(r'log $f_{\rm GW}$ [Hz]')
        ax3.set_ylabel(r'nr events [rnd norm]')
        ax3.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)
    
    #MULTIPLE dataset:
    if (fig3type == 2):
        ax3addfilename = 'figt2_'
        #all:
        pos = pos_ej_tLTtH + pos_2body + pos_3body
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=True,  normed = True, facecolor = colors_datasets[dsc],   color = colors_datasets[dsc], linewidth = 0, alpha = alpha_datasets[dsc])    
        ax3.hist(np.log10(out_fGW_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, normed = True, label = legend_datasets[dsc],       color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = 1.0, linestyle = linest_datasets[dsc])    
        #plot settings: 
        yaxfac = 0.4
        #ax3.text(-2.25, 0.85*yaxfac, 'LISA', fontsize=10, color = 'black', alpha = 0.5)
        #ax3.text(-3.0,  0.75*yaxfac, '/////////////////////', fontsize=10, color = 'black', alpha = 0.5)
        #ax3.text(-0.45,  0.85*yaxfac, 'DECIGO', fontsize=10, color = 'black', alpha = 0.5)
        #ax3.text(-0.5,  0.75*yaxfac, '/////////////', fontsize=10, color = 'black', alpha = 0.5)
        #ax3.text(1.65,  0.85*yaxfac, 'LIGO', fontsize=10, color = 'black', alpha = 0.5)
        #ax3.text(1.0,   0.75*yaxfac, '/////////////////////', fontsize=10, color = 'black', alpha = 0.5)
        ax3.set_title(r'peak GW freq. at formation')
        ax3.set_xlabel(r'log $f_{\rm GW}$ [Hz]')
        ax3.set_ylabel(r'nr events [normed]')
        ax3.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)
    
    #plt.show()
    #exit()
    #------------------------------------------------
    #------------------------------------------------
    #ecc dist at X Hz:
    #------------------------------------------------
    #fig, ax1 = plt.subplots(figsize=(5, 5))
    Hmin    = -3.5
    Hmax    = 0.0   
    Hnrbins = 30
    
    #JUST ONE dataset:
    if (fig4type == 1):
        ax4addfilename = 'figt1_'
        #all:
        pos = pos_ej_tLTtH + pos_2body + pos_3body
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=True,                                 color = 'grey',   linewidth = 1, alpha = 0.025)#hatch = '\\\\\\\\\\')
        #ejected:
        pos = pos_ej_tLTtH
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = 'ejected',             color = 'blue',    linewidth = 1, hatch = '/////////')
        #2-body:        
        pos = pos_2body_lastSEmerge_N
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (subset A)',   color = 'brown',   linewidth = 0.5)#, hatch = '\\\\\\\\')
        pos = pos_2body_lastSEmerge_Y
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (subset B)',   color = 'orange',  linewidth = 0.5)#, hatch = '\\\\\\')
        #pos = pos5
        #ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (subset C)',   color = 'yellow',  linewidth = 0.5)#, hatch = '//////')
        pos = pos_2body
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '2-body (all)',        color = 'green',   linewidth = 1, hatch = '\\\\\\\\\\\\')
        #3-body:        
        #pos = pos_3body
        #ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = '3-body (all)',        color = 'red',     linewidth = 1, hatch = '///////////')
        #all:
        pos = pos_ej_tLTtH + pos_2body + pos_3body
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, label = 'all',                 color = 'black',   linewidth = 1, alpha = 1.0)#hatch = '\\\\\\\\\\')
        #plot settings:
        ax4.set_xlim([Hmin, 1.0]) 
        ax4.set_title(r'Eccentricity at $10^{-2}$ Hz')
        ax4.set_xlabel(r'log $e$ [$10^{-2}$ Hz]')
        ax4.set_ylabel(r'nr events [rnd norm]')
        ax4.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)
    
        
    #MULTIPLE dataset:
    if (fig4type == 2):
        ax4addfilename = 'figt2_'
        pos = pos_ej_tLTtH + pos_2body
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=True,  normed = True, facecolor = colors_datasets[dsc],   color = colors_datasets[dsc], linewidth = 0, alpha = alpha_datasets[dsc])    
        ax4.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, normed = True, label = legend_datasets[dsc],       color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = 1.0, linestyle = linest_datasets[dsc])    
        ax4.set_xlim([Hmin, 0.0]) 
        ax4.set_title(r'Eccentricity $e$ at $10^{-2}$ Hz')
        ax4.set_xlabel(r'log $e$ [$10^{-2}$ Hz]')
        ax4.set_ylabel(r'nr events [normed]')
        ax4.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)

        ax9addfilename = 'figt2_'
        pos = pos_ej_tLTtH + pos_2body
        ax9.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', cumulative=-1, stacked=True, fill=True,  normed = True, facecolor = colors_datasets[dsc],   color = colors_datasets[dsc], linewidth = 0, alpha = alpha_datasets[dsc])    
        ax9.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', cumulative=-1, stacked=True, fill=False, normed = True, label = legend_datasets[dsc],       color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = 1.0, linestyle = linest_datasets[dsc])    
        ax9.set_xlim([-2.0, 0.0]) 
        ax9.set_ylim([0.0,  0.6])  
        ax9.set_xlabel(r'log $e$ [$10^{-2}$ Hz]')
        ax9.set_ylabel(r'cumulative dist.')
        ax9.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)

        ax10addfilename = 'figt3_'
        pos = pos_ej_tLTtH + pos_2body
        yvals, binvals, patches = ax10.hist(np.log10(ecc_fGW_A_arr[pos]), range=(Hmin,Hmax), bins=Hnrbins, histtype='step', cumulative=-1, stacked=True, fill=False, normed = True, linewidth = 0.0)    
        xvals = [(binvals[bc]+binvals[bc+1])/2. for bc in range(0,len(binvals)-1)]
        #use the first dataset as a norm:
        if (dsc == 0): norm_yvals = yvals
        ax10.plot(xvals, yvals-norm_yvals, color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = alpha_datasets[dsc], linestyle = linest_datasets[dsc], label = legend_datasets[dsc])
        ax10.set_xlim([Hmin, 0.0])
        ax10.set_ylim([-0.2, 0.2])  
        ax10.set_title(r'Convergence R/a')
        ax10.set_xlabel(r'log $e$ [$10^{-2}$ Hz]')
        ax10.set_ylabel(r'Frac$(e>e_{-2})$-Frac$(e>e_{-2})_{50}$')
        ax10.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 2, frameon = False)


    #plt.show()
    #exit()
    #------------------------------------------------
    #nr in-cluster vs. ejected mergers:
    #------------------------------------------------
    totnr23     = 1.*len(pos_2body) + 1.*len(pos_3body) #nr in-cluster
    totnrejtH   = 1.*len(pos_ej_tLTtH)                  #nr ejected t<tH
    totnr23ejtH = totnr23 + totnrejtH                   #nr tot mergs.
    xp          = xvalplot[dsc]
    msizeplot   = msize[dsc]
    #set in-cluster:
    yp          = totnr23/totnr23ejtH
    yp_err      = np.sqrt(totnr23)/totnr23ejtH
    ax11.errorbar(xp, yp, yp_err, markersize=msizeplot, marker='o', alpha = alpha_datasets[dsc], markeredgewidth=1, color = colors_datasets[dsc], linewidth = 1.0)
    #dummy for legend:
    if (dsc == 0):  ax11.plot(xp-1000, yp-1000, markersize=10, marker='o', alpha = 0.85, markeredgewidth=1, color = 'black', linewidth = 0.0, label = r'$N_{\rm IC}/N_{\rm tot}^{<t_{\rm H}}$ (hybrid)')
    if (dsc == 0):  ax11.plot(xp-1000, yp-1000, markersize=8, marker='o', alpha = 0.75, markeredgewidth=1, color = 'dodgerblue', linewidth = 0.0, label = r'$N_{\rm IC}/N_{\rm tot}^{<t_{\rm H}}$ (HS19)')
    #set ejected:
    yp          = totnrejtH/totnr23ejtH
    yp_err      = np.sqrt(totnrejtH)/totnr23ejtH
    ax11.errorbar(xp, yp, yp_err, markersize=msizeplot, marker='x', alpha = alpha_datasets[dsc], markeredgewidth=1, color = colors_datasets[dsc], linewidth = 1.0)
    #dummy for legend:
    if (dsc == 0):  ax11.plot(xp-1000, yp-1000, markersize=10, marker='x', alpha = 0.85, markeredgewidth=1, color = 'black', linewidth = 0.0, label = r'$N_{\rm EJ}^{<t_{\rm H}}/N_{\rm tot}^{<t_{\rm H}}$ (hybrid)')
    if (dsc == 0):  ax11.plot(xp-1000, yp-1000, markersize=8, marker='x', alpha = 0.75, markeredgewidth=1, color = 'dodgerblue', linewidth = 0.0, label = r'$N_{\rm EJ}^{<t_{\rm H}}/N_{\rm tot}^{<t_{\rm H}}$ (HS19)')
    #set axes label names:
    legend_names  = [r'$R/a = 2$', r'$R/a = 5$', r'$R/a = 10$', r'$R/a = 25$', r'$R/a = 50$']
    plt.xticks(range(1,len(legend_names)+1,1), legend_names[:], rotation = -25)
    ax11.set_ylabel(r'$N/N_{\rm tot}^{<t_{\rm H}}$')
    ax11.set_title(r'in-cluster vs. ejected mergers')
    #guide lines:
    ax11.plot([0,10], [0.5,0.5], color = 'black', linestyle = ':', linewidth = 1.0)
    ax11.set_xlim([1-0.5, 5+0.5])
    ax11.set_ylim([0.35,0.65])  
    ax11.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 2, frameon = False)
    #----------------------------------------------------------
    #----------------------------------------------------------


#--------------------------------------------------------------
fig1.savefig('BBHpert_eccdist.pdf',     bbox_inches='tight')        
fig2.savefig('BBHpert_obsFGWdist.pdf',  bbox_inches='tight')        
fig11.savefig('BBHpert_nrBBHmergsinout.pdf',    bbox_inches='tight')        

if (fig3type > 0 and fig4type > 0):
    fig3.savefig(ax3addfilename + 'BBHpert_fGWdist.pdf',     bbox_inches='tight')        
    fig4.savefig(ax4addfilename + 'BBHpert_edistXHz.pdf',    bbox_inches='tight')        

if (fig4type == 2):
    fig9.savefig(ax9addfilename + 'BBHpert_edistXHzcumu.pdf',    bbox_inches='tight')        
    fig10.savefig(ax10addfilename + 'BBHpert_edistXHzcumu.pdf',    bbox_inches='tight')        

#TESTS:
fig5.savefig('BBHpert_Radistbreakdown.pdf',    bbox_inches='tight')        
fig6.savefig('BBHpert_Ra1ebreakdown.pdf',    bbox_inches='tight')        




plt.show()

exit()

#dist of ecc at 10-2 Hz
#nr BBH mergers
#fGW can it be used?
#all seems ok compared to my previous codes: results very sensible to vesc (amin), so I think its fine.
#chekc this script again.
#can we prove that there is a flow towards high ecc in the ecc limit?
#why do we see more mergers? is it a mistake?
#maybe make one mode check against previous results.











#TESTS:




#distdata_wfcs[tc,0]  = fGWp0    
#distdata_wfcs[tc,1]  = wfc
#distdata_wfcs[tc,2]  = endstates_arr[i,0]   
#distdata_wfcs[tc,3]  = e_bin_0
#pos2 = np.where(distdata_wfcs[:,2] == 2)[0]
#distdata_id2    = distdata_wfcs[pos2,:]
     
#sort data:
endID_arr               = distdata_wfcs[:,2]
flagenc_arr             = distdata_wfcs[:,6]
flag_lastSEmerge_YN_arr = distdata_wfcs[:,12]
rpENC_Ua_arr            = distdata_wfcs[:,15]
pos_2body       = np.where(endID_arr[:] == 1)[0]
pos_encflag1    = np.where(flagenc_arr[:] == 1)[0]
pos_encflag2    = np.where(flagenc_arr[:] == 2)[0]
pos_2body_encflag1  = list(set(pos_2body).intersection(pos_encflag1)) 
pos_2body_encflag2  = list(set(pos_2body).intersection(pos_encflag2))
pos_lastSEmerge_N  = np.where(flag_lastSEmerge_YN_arr[:] == -1)[0]
pos_lastSEmerge_Y  = np.where(flag_lastSEmerge_YN_arr[:] ==  1)[0]
pos_2body_encflag2_lastSEmergeN = list(set(pos_2body_encflag2).intersection(pos_lastSEmerge_N))
pos_2body_encflag2_lastSEmergeY = list(set(pos_2body_encflag2).intersection(pos_lastSEmerge_Y))

datasetsID_arr  = distdata_wfcs[:,14]
maxdatasID      = int(max(datasetsID_arr))

ax7.plot([0,100], [1,1], color = 'black', alpha = 1, linestyle = ':')   
for dsIDc in range(0,maxdatasID+1):
    pos_dataset = np.where(datasetsID_arr[:] == dsIDc)[0]
    posplot     = list(set(pos_dataset).intersection(pos_2body_encflag2_lastSEmergeN))
    posplot.sort()
    if (len(posplot) > 1):
        xp    = distdata_wfcs[posplot,10]
        yp    = distdata_wfcs[posplot,13]
        #ax7.plot(xp, yp,                    markersize=2, marker='o', color = colors_datasets[dsc], alpha = 0.5, markeredgewidth=0, linestyle = '-')#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    
        ax7.plot(xp, yp,                    markersize=2, marker='o', alpha = 0.5, markeredgewidth=0, linestyle = '-')#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    
        ax7.plot(xp[len(posplot)-1], yp[len(posplot)-1],    markersize=3, marker='o', color = 'black',              alpha = 0.5, markeredgewidth=0)#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    

        #CHECK BELOW AGAIN!!!
        xp    = distdata_wfcs[posplot,10]
        posploty = posplot
        posploty = [x - 1 for x in posploty]            
        yp    = distdata_wfcs[posploty,15]
        ax8.plot(xp, yp,                    markersize=2, marker='o', color = colors_datasets[dsc], alpha = 0.5, markeredgewidth=0, linestyle = '')#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    
        ax8.plot(xp[len(posplot)-1], yp[len(posplot)-1],    markersize=6, marker='o', color = 'black',              alpha = 0.5, markeredgewidth=0)#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    

ax7.set_xlim(0,80)
#xp    = distdata_wfcs[pos_2body_encflag2_lastSEmergeY,10]
#yp    = distdata_wfcs[pos_2body_encflag2_lastSEmergeY,13]
#ax7.plot(xp, yp, markersize=2, marker='o', color = 'black', alpha = 0.5, markeredgewidth=0, linestyle = '')#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    

#plt.show()
#exit()

#FINISH THIS FIG: put it below the figure below
Hmin    = 2
Hmax    = 10   
Hnrbins = 50
Ra  = distdata_wfcs[:,4]
data_arr    = Ra
#ax5.hist(data_arr[:], range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, normed = False, label = legend_datasets[dsc], color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = alpha_datasets[dsc], linestyle = linest_datasets[dsc])    
ax5.set_xlim([Hmin, Hmax])
ax5.set_xlabel(r'$E$')
ax5.set_ylabel(r'nr events [normed]')
ax5.legend(loc='upper left', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)
#FINISH THIS FIG: e.g. below it add hist of R/a.
Ra  = distdata_wfcs[:,4]
ec  = distdata_wfcs[:,5]
xp    = Ra
yp    = 1.-ec
#ax6.plot(xp, yp, markersize=4, marker='o', color = colors_datasets[dsc], alpha = 0.5, markeredgewidth=0, linestyle = '')#marker='s', markeredgecolor = 'green',  markeredgewidth=1, markersize=1.5*markersize, fillstyle = 'none')    
deFAC           = -(15.*np.pi/16.)*(((2.*(mBH**2.))/((3.*mBH)*(2.*mBH)))**(1./2.))*(ec*np.sqrt(1.0-ec**2.))
rpDEeq1E_Ua     = (abs(deFAC)/(1.0-ec))**(2./3.)         #rp where abs(de) = 1-e (important for high ecc)
#ax6.scatter(rpDEeq1E_Ua, 1-ec, s=4, marker='o', color = 'black')
ax6.set_xlim(1,20)
ax6.set_ylim(1e-6,1e0)
ax6.set_xscale('log')
ax6.set_yscale('log')




#------------------------------------------------
#eccentricity:
#------------------------------------------------
Hmin    = 0
Hmax    = 1   
Hnrbins = 50

#endID_arr       = distdata_wfcs[:,2]
#pos_2body       = list(np.where(endID_arr[:] == 1)[0]) + list(np.where(endID_arr[:] == 6)[0])
#pos_ej          = list(np.where(endID_arr[:] == 2)[0])
#pos_2body_ej    = pos_2body + pos_ej

data_arr    = distdata_wfcs[:,3]
wf          = distdata_wfcs[:,1]
ax1.hist(data_arr[:], weights = wf, range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, normed = True, label = legend_datasets[dsc], color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = alpha_datasets[dsc], linestyle = linest_datasets[dsc])    
ax1.set_xlim([Hmin, Hmax])
#ax1.set_ylim([0.0, 2.0])
ax1.set_xlabel(r'$E$')
ax1.set_ylabel(r'nr events [normed]')
ax1.legend(loc='upper left', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)

#plt.show()
#exit()
#------------------------------------------------
#obs in-cluster fGW dist:
#------------------------------------------------
Hmin    = -6
Hmax    = -1   
Hnrbins = 50

endID_arr       = distdata_wfcs[:,2]
pos_2body       = list(np.where(endID_arr[:] == 1)[0]) + list(np.where(endID_arr[:] == 6)[0])
pos_ej          = list(np.where(endID_arr[:] == 2)[0])
pos_2body_ej    = pos_2body + pos_ej

#xp  = distdata_wfcs[pos_2body_ej,0]
#wf  = distdata_wfcs[pos_2body_ej,1]
#ax2.hist(np.log10(xp),    weights = wf,    range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, color = 'black', linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')

#xp  = distdata_wfcs[pos_2body,0]
#wf  = distdata_wfcs[pos_2body,1]
#ax2.hist(np.log10(xp),    weights = wf,    range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, color = 'red', linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')

#xp  = distdata_wfcs[pos_ej,0]
#wf  = distdata_wfcs[pos_ej,1]
#ax2.hist(np.log10(xp),    weights = wf,    range=(Hmin,Hmax), bins=Hnrbins, histtype='step', fill=False, color = 'blue', linestyle = '-',    linewidth = 1, alpha = 1.0, label = r'Newt., all.')

#MULTIPLE dataset:
data_arr    = np.log10(distdata_wfcs[:,0])
wf          = distdata_wfcs[:,1]
ax2.hist(data_arr[:], weights = wf, range=(Hmin,Hmax), bins=Hnrbins, histtype='step', stacked=True, fill=False, normed = True, label = legend_datasets[dsc], color = colors_datasets[dsc], linewidth = linewi_datasets[dsc], alpha = alpha_datasets[dsc], linestyle = linest_datasets[dsc])    
ax2.set_xlim([Hmin, Hmax])
ax2.set_ylim([1e-6, 1e0])
ax2.set_yscale('log')
ax2.set_xlabel(r'log $F_{\rm GW}$ [Hz]')
ax2.set_ylabel(r'nr events [normed]')
ax2.legend(loc='upper right', numpoints = 1, fontsize = 10.0, ncol = 1, frameon = False)

#save figure:
#plt.savefig('BBHpert_obsFGWdist.pdf', bbox_inches='tight')        
#plt.show()
#plt.show()
#exit()
#------------------------------------------------


        








