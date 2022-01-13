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
import csv



#----------------------------------------------------------
#Units and conversions:
#----------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
c_SI        = 299792458.0        #m/s
M_sun_SI    = 1.989*(10.**30.)   #kg
R_sun_SI    = 695800000.         #m
AU_SI       = 149597871000.      #m 
G_new_SI    = 6.67*(10.**(-11.))
AU_U        = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U     = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
sec_year    = 31536000.
#----------------------------------------------------------

#GW inspiral life time:
tf = open('a_e_tlife_yrs_interpoltable.txt', "r")
a_e_tlife_yrs_arr   = np.loadtxt(tf, dtype=float)
tf.close()
#read data and make interpolation function:
tdata_e_arr                 = a_e_tlife_yrs_arr[:,1]
tdata_t_arr                 = a_e_tlife_yrs_arr[:,2]
ecc_to_tdata_interpol_func  = interpolate.interp1d(tdata_e_arr, tdata_t_arr)  

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#----------------------------------------------------------
#Read in data from files:
#----------------------------------------------------------
#data name:

data_name = 'SIM1_BH202020_NT_WGR_01_02_5_10kmsec_T2500'#'SIM1_BH202020_NT_WGR_001_1_10kmsec_T2500'#'SIM1_BH101010_NT_WGR_0001_10_10kmsec_T2500'#'TEST2_IMSae'

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
#----------------------------------------------------------
#Define:
#----------------------------------------------------------
nr_SMA                          = MC_settings_list_INT[0]
nr_vinf                         = MC_settings_list_INT[1]
nr_scatterings_per_paramcomb    = MC_settings_list_INT[2]
#----------------------------------------------------------
#reshape input arrays for easy analysis:
#----------------------------------------------------------
RS_avs_outlist_MC_info_INT          = outlist_MC_info_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)                #[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0]
RS_avs_outlist_MC_info_REAL         = outlist_MC_info_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)               #[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]
RS_avs_output_Nbody_endstate_INT    = output_Nbody_endstate_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)          #[out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]   
RS_avs_output_Nbody_endstate_REAL   = output_Nbody_endstate_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)         #[1 E_kin(ij), 2 E_pot(ij), 3 E_tot(ij), 4 a_bin(ij), 5 e_bin(ij), 6 E_kin(ijk), 7 E_pot(ijk), 8 E_tot(ijk), 9 a_bin(ijk), 10 e_bin(ijk)]
RS_avs_output_Nbody_xtra_info_INT   = output_Nbody_xtra_info_INT.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)
RS_avs_output_Nbody_xtra_info_REAL  = output_Nbody_xtra_info_REAL.reshape(nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, 10)        #[rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
#From data define:
MCvalues_SMA_arr                    = RS_avs_outlist_MC_info_REAL[:,0,0,0]
MCvalues_vinf_arr                   = RS_avs_outlist_MC_info_REAL[0,:,0,1]
#----------------------------------------------------------

#----------------------------------------------------------
#Allocate arrays:
#----------------------------------------------------------
#below: 'sim' is the direct count from the sim. This can be slightly biased since some interactions may not have been completed wihtin the time limit. 
nr_eachstate_sim_DI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
nr_eachstate_sim_RI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
#below: 'fin' is our (best try) corrected nr of endstates, where we try to correct for the set of non-completed interactions.
nr_eachstate_fin_DI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
nr_eachstate_fin_RI                 = np.zeros((3,10), dtype=np.float64)  #temp arr for a given a,v comb.
#data save:
tot_nr_scat = nr_SMA*nr_scatterings_per_paramcomb
IMS_GWinsp_123_aeINFO               = np.zeros((tot_nr_scat, 10), dtype=np.float64)
IMS_BSinsp_123_aeINFO               = np.zeros((tot_nr_scat, 10), dtype=np.float64)
IMS_GWinsp_123_nrc                  = np.zeros((nr_SMA, 10), dtype=np.float64)
IMS_GWinsp_123_crosssecAU2          = np.zeros((nr_SMA, 10), dtype=np.float64)
#----------------------------------------------------------

#----------------------------------------------------------
#loop over a,v,s:
#----------------------------------------------------------
vc      = 0     #in this analysis we keep vc fixed.
totc    = -1    #initialize tot counter     

for ac in range(0,nr_SMA):                              # ac-loop over SMA
   
    #---------------------------------------
    #count nr of outcomes:
    #---------------------------------------
    for sc in range(0,nr_scatterings_per_paramcomb):    # sc-loop over nr scatterings per param combination (a(ac),v(vc),...)          
        
        #1 = TDE, 2 = COLLISION, 3 = Bin-Sin, 4 = Ionization, 5 = Inspiral, 10 = time limit (Nbody time), 11 = steps, 12 = time limit (wall clock in sec), 13 = Solver not succesfull, ...
		#Return_3body_Info_REAL(1:10)		= [E_kin(ij), E_pot(ij), E_tot(ij), a_bin(ij), e_bin(ij), E_kin(ijk), E_pot(ijk), E_tot(ijk), a_bin(ijk), e_bin(ijk)]
        #RS_avs_output_Nbody_endstate_INT: [out_end_state_flag, out_bin_i, out_bin_j, out_sin_k, out_IMS_bin_yesno, IMS_rp_counter, IMS_binsin_counter, ...]
        #RS_avs_output_Nbody_xtra_info_REAL: [rmin_12, rmin_13, rmin_23, MRini_IMS_a, MRini_IMS_e]
        
        #update total scattering counter:
        totc = totc+1
        print totc
        
        #calc weight fac:
        bmax_Rsun   = RS_avs_outlist_MC_info_REAL[ac,vc,0,2]    #in Rsun
        Rsun2_AU2   = (R_sun_SI/AU_SI)**2.
        cs_wfac     = (1./(1.*nr_scatterings_per_paramcomb))*(np.pi*(bmax_Rsun**2.))*Rsun2_AU2
        
        #endstate:
        endstate_id = RS_avs_output_Nbody_endstate_INT[ac,vc,sc,0]
        
        #-----------------------------------
        #EQUAL MASS GW MERGER SCENARIO:
        #-----------------------------------
        #input:
        f_GW    = 10.0 #in Hz
        mMsun   = 20.
        mSI     = mMsun*M_sun_SI        
        
        #-----------------------------------
        if (endstate_id == 3):
        #-----------------------------------    
            #initial a,e after BS endstate:
            iniRsun_BS_a    = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,3]
            iniSI_BS_a      = iniRsun_BS_a*R_sun_SI
            iniSI_BS_e      = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,4]
            #calc insp time:
            #tfac    = (iniRsun_BS_a**(4.))/(mMsun*mMsun*(mMsun+mMsun))
            #read ecc:
            #e_bin   = iniSI_BS_e
            #calc tlife: from interpolation:
            #if (e_bin <  max(tdata_e_arr)):
            #    tbin_yr = tfac*ecc_to_tdata_interpol_func(e_bin)
            #if (e_bin >= max(tdata_e_arr)):
            #    tbin_yr_circular = tfac*ecc_to_tdata_interpol_func(0.0)
            #    tbin_yr = (768./425.)*tbin_yr_circular*((1.-e_bin**2.)**(7./2.))
            
            #calc tlife: fast procedure (assumes e>>0 limit):
            tbin_yr_circular    = ((iniSI_BS_a**4.)/(4.*(64./5.)*(G_new_SI**3.)*mSI*mSI*(mSI+mSI)/(c_SI**5.)))/sec_year
            tbin_yr             = ((768./425.)*tbin_yr_circular*((1.-iniSI_BS_e**2.)**(7./2.)))            
            #calc v_kick:
            Etot_binsin         = RS_avs_output_Nbody_endstate_REAL[ac,vc,sc,7]
            ESI_fac             = G_new_SI*(M_sun_SI**2.)/R_sun_SI
            Etot_SI             = Etot_binsin*ESI_fac
            mbin                = mSI+mSI
            msin                = mSI
            v_inf_bin           = np.sqrt(2.*Etot_SI*((mbin*(1.+mbin/msin))**(-1.)))   #SI units
            v_kick_kms          = v_inf_bin/1000.
            
            
            #print (((0.1*AU_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*mSI*mSI*(mSI+mSI)/(c_SI**5.)))/sec_year
            #exit()
            #print v_kick_kms
            
            
            if (tbin_yr < 10.**(10.) and v_kick_kms > 50.0):
                #propagating ini a,e to a,e(f_GW)
                c0      = iniSI_BS_a/((iniSI_BS_e**(12./19.)/(1.-iniSI_BS_e**2.))*((1.+(121./304.)*iniSI_BS_e**2.)**(870./2299.)))
                #a(ecc)  = (c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))
                func    = lambda ecc : f_GW - (1./np.pi)*np.sqrt(G_new_SI*(mSI+mSI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
                ecc_initial_guess   = 1e-8
                fGW_BS_e            = fsolve(func, ecc_initial_guess)   #value of bin ecc when the bin is at f_GW            
                #save ini e, fin e, weightfac, ...
                IMS_BSinsp_123_aeINFO[totc,0] = iniSI_BS_e
                IMS_BSinsp_123_aeINFO[totc,1] = fGW_BS_e
                IMS_BSinsp_123_aeINFO[totc,2] = cs_wfac
                print iniSI_BS_e, fGW_BS_e
        #-----------------------------------           
        
        #-----------------------------------
        if (endstate_id == 5):
        #-----------------------------------    
            #initial a,e at formation of IMS
            iniSI_IMS_a     = RS_avs_output_Nbody_xtra_info_REAL[ac,vc,sc,3]*R_sun_SI
            iniSI_IMS_e     = RS_avs_output_Nbody_xtra_info_REAL[ac,vc,sc,4]
            #propagating ini a,e to a,e(f_GW)
            c0      = iniSI_IMS_a/((iniSI_IMS_e**(12./19.)/(1.-iniSI_IMS_e**2.))*((1.+(121./304.)*iniSI_IMS_e**2.)**(870./2299.)))
            #a(ecc)  = (c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))
            func    = lambda ecc : f_GW - (1./np.pi)*np.sqrt(G_new_SI*(mSI+mSI)/((c0*((ecc**(12./19.)/(1.-ecc**2.))*((1.+(121./304.)*ecc**2.)**(870./2299.))))**3.0))*((1.+ecc)**1.1954)/((1.-ecc**2.)**1.5)
            ecc_initial_guess   = 0.001
            fGW_IMS_e           = fsolve(func, ecc_initial_guess)   #value of bin ecc when the bin is at f_GW
            #save counters:
            #cut on ecc (1):
            if (fGW_IMS_e > 0.00 and fGW_IMS_e < 0.25):
                IMS_GWinsp_123_nrc[ac,0] = IMS_GWinsp_123_nrc[ac,0]+1
            #cut on ecc (2):
            if (fGW_IMS_e > 0.25 and fGW_IMS_e < 0.50):
                IMS_GWinsp_123_nrc[ac,1] = IMS_GWinsp_123_nrc[ac,1]+1
            #cut on ecc (3):
            if (fGW_IMS_e > 0.50 and fGW_IMS_e < 0.75):
                IMS_GWinsp_123_nrc[ac,2] = IMS_GWinsp_123_nrc[ac,2]+1
            #cut on ecc (4):
            if (fGW_IMS_e > 0.75 and fGW_IMS_e < 1.00):
                IMS_GWinsp_123_nrc[ac,3] = IMS_GWinsp_123_nrc[ac,3]+1
            #no cut:
            IMS_GWinsp_123_nrc[ac,4] = IMS_GWinsp_123_nrc[ac,4]+1
            #save ini e, fin e, weightfac, ...
            IMS_GWinsp_123_aeINFO[totc,0] = iniSI_IMS_e
            IMS_GWinsp_123_aeINFO[totc,1] = fGW_IMS_e
            IMS_GWinsp_123_aeINFO[totc,2] = cs_wfac    
            #print info:
            #print np.log(MCvalues_SMA_arr[ac]), iniSI_IMS_e, fGW_IMS_e, func(fGW_IMS_e), (1./np.pi)*np.sqrt(G_new_SI*(mSI+mSI)/((iniSI_IMS_a)**3.0))*((1.+iniSI_IMS_e)**1.1954)/((1.-iniSI_IMS_e**2.)**1.5)  
        #-----------------------------------     
    #---------------------------------------
    
    #---------------------------------------
    #convert to cross section:
    #---------------------------------------
    IMS_GWinsp_123_crosssecAU2[ac,0]    = cs_wfac*IMS_GWinsp_123_nrc[ac,0]
    IMS_GWinsp_123_crosssecAU2[ac,1]    = cs_wfac*IMS_GWinsp_123_nrc[ac,1]
    IMS_GWinsp_123_crosssecAU2[ac,2]    = cs_wfac*IMS_GWinsp_123_nrc[ac,2]
    IMS_GWinsp_123_crosssecAU2[ac,3]    = cs_wfac*IMS_GWinsp_123_nrc[ac,3]
    IMS_GWinsp_123_crosssecAU2[ac,4]    = cs_wfac*IMS_GWinsp_123_nrc[ac,4]
    #---------------------------------------
    
    
#----------------------------------------------------------   
#----------------------------------------------------------










fig, ax1 = plt.subplots(figsize=(6, 5))

#IMS:
fGW_IMS_e_arr   = IMS_GWinsp_123_aeINFO[:,1]
cs_IMS_wfac_arr = IMS_GWinsp_123_aeINFO[:,2]
input_arr       = np.log10(fGW_IMS_e_arr)
warr            = cs_IMS_wfac_arr
ax1.hist(input_arr, bins=250, range = [-8, 0.1], weights = warr, normed=False, linewidth=2.0, histtype='step', linestyle = '-', color = 'black', label='GW inspiral mergers')

#BS:
fGW_BS_e_arr    = IMS_BSinsp_123_aeINFO[:,1]
cs_BS_wfac_arr  = IMS_BSinsp_123_aeINFO[:,2]
input_arr       = np.log10(fGW_BS_e_arr)
warr            = cs_BS_wfac_arr
ax1.hist(input_arr, bins=250, range = [-8, 0.1], weights = warr, normed=False, linewidth=0.5, histtype='step', linestyle = '-', color = 'black', label = 'Post-interaction GW mergers')

#open CR data and overplot:
tf = open('CRdata.txt', "r")
CRdata_x = []
CRdata_y = []
filecsv_read  = csv.reader(tf)
for row in filecsv_read:
    CRdata_x.append(float(row[0]))
    CRdata_y.append(float(row[1]))
tf.close()
CRdata_x    = np.array(CRdata_x)
CRdata_y    = np.array(CRdata_y)
#sort, norm and plot:
pos_sort = np.argsort(CRdata_x)
CRdata_x = CRdata_x[pos_sort]
CRdata_y = CRdata_y[pos_sort]
CRdata_y    = 34.0*CRdata_y/(1.*max(CRdata_y))
ax1.plot(CRdata_x[:], CRdata_y[:], color='black', linestyle = ':', label = 'Rodriguez et al. 2016')

ax1.set_yticklabels([])
ax1.set_ylim([-0.1, 35.0])
ax1.set_xlabel(r'log $e$ [10 Hz]')
ax1.set_ylabel(r'histogram counts')
ax1.set_title(r'BBH Eccentricity Distribution at 10Hz')
plt.legend(loc='upper right', numpoints = 1, fontsize = 10.0, frameon = False)


# These are in unitless percentages of the figure size. (0,0 is bottom left)
left, bottom, width, height = [0.55, 0.35, 0.3, 0.3]
ax2 = fig.add_axes([left, bottom, width, height])
#only BS:
input_arr           = fGW_BS_e_arr
warr                = cs_BS_wfac_arr
ax2.hist(input_arr, bins=500, range = [1e-8, 1.0], weights = warr, cumulative=True, normed=True, linewidth=0.5, linestyle = '-', histtype='step', color = 'black')
#BS+IMS:
input_arr           = np.append(fGW_IMS_e_arr, fGW_BS_e_arr)
warr                = np.append(cs_IMS_wfac_arr, cs_BS_wfac_arr)
ax2.hist(input_arr, bins=500, range = [1e-8, 1.0], weights = warr, cumulative=True, normed=True, linewidth=2.0, linestyle = '-', histtype='step', color = 'black')
#axis settings:
ax2.set_xlabel(r'$e$ [10 Hz]', fontsize=10)
ax2.set_title(r'cumulative histogram', fontsize=10)
ax2.set_xlim([0.00, 1.0])
ax2.set_ylim([0.98, 1.0])
ax2.tick_params(labelsize=10)


plt.savefig('BBH_eccdist_BSGWinsp.eps', bbox_inches='tight')        
plt.show()



exit()














fig, ax1 = plt.subplots(figsize=(5, 4))

#IMS:
fGW_IMS_e_arr       = IMS_GWinsp_123_aeINFO[:,1]
cs_IMS_wfac_arr     = IMS_GWinsp_123_aeINFO[:,2]
#BS:
fGW_BS_e_arr        = IMS_BSinsp_123_aeINFO[:,1]
cs_BS_wfac_arr      = IMS_BSinsp_123_aeINFO[:,2]

input_arr       = fGW_IMS_e_arr
warr            = cs_IMS_wfac_arr
ax1.hist(input_arr, bins=50, range = [0.0, 1.0], weights = warr, normed=False, linewidth=0.5, histtype='step')

input_arr       = fGW_BS_e_arr
warr            = cs_BS_wfac_arr
ax1.hist(input_arr, bins=50, range = [0.0, 1.0], weights = warr, normed=False, linewidth=0.5, histtype='step')

plt.yscale('log')

plt.show()










fig, ax1 = plt.subplots(figsize=(5, 4))


iniSI_IMS_e_arr     = IMS_GWinsp_123_aeINFO[:,0]
fGW_IMS_e_arr       = IMS_GWinsp_123_aeINFO[:,1]
cs_IMS_wfac_arr     = IMS_GWinsp_123_aeINFO[:,2]

iniSI_BS_e_arr      = IMS_BSinsp_123_aeINFO[:,0]
fGW_BS_e_arr        = IMS_BSinsp_123_aeINFO[:,1]
cs_BS_wfac_arr      = IMS_BSinsp_123_aeINFO[:,2]

print type(fGW_IMS_e_arr)

input_arr           = np.append(fGW_IMS_e_arr, fGW_BS_e_arr)
warr                = np.append(cs_IMS_wfac_arr, cs_BS_wfac_arr)
ax1.hist(np.log10(input_arr), bins=50, range = [-8.0, 0.0], weights = warr, cumulative=True, normed=True, linewidth=0.5, histtype='step')




#plt.xscale('log')



plt.show()





fig = plt.figure(figsize=(12, 8))

xa  = MCvalues_SMA_arr[:]*(R_sun_SI/AU_SI)

ycs = IMS_GWinsp_123_crosssecAU2[:,0]
fig.add_subplot(111).plot(xa, ycs, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

ycs = IMS_GWinsp_123_crosssecAU2[:,1]
fig.add_subplot(111).plot(xa, ycs, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

ycs = IMS_GWinsp_123_crosssecAU2[:,2]
fig.add_subplot(111).plot(xa, ycs, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

ycs = IMS_GWinsp_123_crosssecAU2[:,3]
fig.add_subplot(111).plot(xa, ycs, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

ycs = IMS_GWinsp_123_crosssecAU2[:,4]
fig.add_subplot(111).plot(xa, ycs, marker='o', linestyle=':', markersize=5.0, linewidth=1.5, alpha=0.75)

plt.xscale('log')
plt.yscale('log')

plt.show()


#CUT ON:
#GW insp time
#v_kick
#





































