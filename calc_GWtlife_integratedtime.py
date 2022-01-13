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

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)
 
 
#print ((3./2.)*G_new_SI*(20.*M_sun_SI)/(10000.**2))/AU_SI
#exit()
 
#----------------------------------------------------------
#Read in data from files:
#----------------------------------------------------------
#data name:

data_name = 'SIM1_BH202020_NT_WGR_001_1_10kmsec_T2500'#'SIM1_BH101010_NT_WGR_0001_10_10kmsec_T2500'#'TEST2_IMSae'

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

#INPUT: equal mass case:
mSI                 = 20.*M_sun_SI

#calc:
totnr               = nr_SMA*nr_scatterings_per_paramcomb

#read info:
all_iniSMA          = RS_avs_outlist_MC_info_REAL[:,0,:,0].reshape(totnr)
all_bmax            = RS_avs_outlist_MC_info_REAL[:,0,:,2].reshape(totnr)
all_csweight        = all_bmax**2.
all_IDs             = RS_avs_output_Nbody_endstate_INT[:,0,:,0].reshape(totnr)
all_aRsun           = RS_avs_output_Nbody_endstate_REAL[:,0,:,3].reshape(totnr)
all_ecc             = RS_avs_output_Nbody_endstate_REAL[:,0,:,4].reshape(totnr)

#calc inspiral times:
all_aSI                 = all_aRsun*R_sun_SI
all_tbin_yr_circular    = (all_aSI**4.)/(4.*(64./5.)*(G_new_SI**3.)*mSI*mSI*(mSI+mSI)/(c_SI**5.))
all_tbin_yr             = ((768./425.)*all_tbin_yr_circular*((1.-all_ecc**2.)**(7./2.)))/sec_year            

#print info:
print ((MCvalues_SMA_arr[:]*R_sun_SI)**4.)/(4.*(64./5.)*(G_new_SI**3.)*mSI*mSI*(mSI+mSI)/(c_SI**5.))/sec_year

#choose the ini SMA we consider:
useSMA      = MCvalues_SMA_arr[4]
pos_useSMA  = np.where(all_iniSMA[:]  ==  useSMA)[0]
print 'ONLY INCLUDES (AU):', useSMA*R_sun_SI/AU_SI





#derive time merger distributions:
pos_ID3     = np.where(all_IDs[:] == 3)[0]
pos_ID3     = list(set(pos_ID3).intersection(pos_useSMA))
tbin_yr_id3 = all_tbin_yr[pos_ID3]
warr_id3    = all_csweight[pos_ID3]

pos_ID5     = np.where(all_IDs[:] == 5)[0]
pos_ID5     = list(set(pos_ID5).intersection(pos_useSMA))
tbin_yr_id5 = 1e-5 + 0.0*all_tbin_yr[pos_ID5]
warr_id5    = all_csweight[pos_ID5]

nrbins      = 500
logtmin     = 8   #in log years
logtmax     = 10  #in log years
dbin        = (logtmax-logtmin)/(1.*nrbins)
bin_edges   = np.linspace(logtmin, logtmax, nrbins+1)

#allocate arrays:
totint_tbin_yr_id3 = []
totint_warr_id3    = []

totint_tbin_yr_id5 = []
totint_warr_id5    = []

#integrate forward in time:
for bc in range(0,nrbins):  #loop over time bins
    tmin_bc     = 10.**(bin_edges[bc])
    tmax_bc     = 10.**(bin_edges[bc+1])
    tcenter_bc  = 10.**(bin_edges[bc] + 0.5*dbin) 
    
    dt_bc   = tmax_bc-tmin_bc   #delta bin in years
    
    totint_tbin_yr_id3.append(tcenter_bc + tbin_yr_id3)
    totint_warr_id3.append(dt_bc*warr_id3)

    totint_tbin_yr_id5.append((tcenter_bc + tbin_yr_id5))
    totint_warr_id5.append(dt_bc*warr_id5)    










#ratio of cumulative rate:
fig = plt.figure(figsize=(5, 4))

input_arr   = np.log10(totint_tbin_yr_id3)
warr        = totint_warr_id3
hist_id3    = np.histogram(input_arr, bins = bin_edges, weights = warr, normed=False)[0]
cumhist_id3 = np.cumsum(hist_id3)

input_arr   = np.log10(totint_tbin_yr_id5)
warr        = totint_warr_id5
hist_id5    = np.histogram(input_arr, bins = bin_edges, weights = warr, normed=False)[0]
cumhist_id5 = np.cumsum(hist_id5)

hist_ratio      = hist_id5/hist_id3
cumhist_ratio   = cumhist_id5/cumhist_id3

fig.add_subplot(111).plot(bin_edges[0:nrbins]+0.5*dbin, cumhist_ratio[:], label = r'numerical sol. (steady-state)', color = 'red', linewidth=1.5)
#fig.add_subplot(111).plot(bin_edges[0:nrbins]+0.5*dbin, hist_ratio[:], label = r'$\Gamma_{\rm I}(t)/\Gamma_{\rm M}(t)$', color = 'black', linewidth=1.5)

t_yr = (10.**(np.linspace(8.0, 10.0, num=100)))
frac_NI_NM = (9./7.)*(3.3*10.**(-2.))*((20./1.)**(-1./7.))*((1./1.)**(3./7.))*((t_yr/(10.**10))**(-2./7.))
fig.add_subplot(111).plot(np.log10(t_yr), frac_NI_NM, linewidth=1.0, linestyle = '--', color='red', label='analytical sol. (steady-state)')

t_yr = (10.**(np.linspace(0.0, 10.0, num=100)))
frac_NI_NM = (3.3*10.**(-2.))*((20./1.)**(-1./7.))*((1./1.)**(3./7.))*((t_yr/(10.**10))**(-2./7.))
fig.add_subplot(111).plot(np.log10(t_yr), frac_NI_NM, linewidth=1.0, linestyle = '--', color='grey', label='analytical sol. (early-burst)')

#add vertical guideline:
fig.add_subplot(111).plot([8.,8.], [-1e0,1e10], linewidth=1.0, linestyle = ':', color='black')

#axis settings:
fig.add_subplot(111).set_xlim([3,10])
fig.add_subplot(111).set_ylim([10.**(-2.), 10.**(0.)])
plt.yscale('log')
plt.xlabel(r'log $t$ [years]')
plt.ylabel(r'$N_{\rm I}/N_{\rm M}$')
#plt.title(r'Cumulative BBH Mergers')
plt.legend(loc='lower left', numpoints = 1, fontsize = 10.0, frameon = False)

#Save figure:
plt.savefig('ratio_BHmergers.eps', bbox_inches='tight')

plt.show()

exit()










#rate:
fig = plt.figure(figsize=(5, 4))

input_arr   = np.log10(totint_tbin_yr_id3)
warr        = totint_warr_id3
hist        = np.histogram(input_arr, bins = bin_edges, weights = warr, normed=False)[0]
fig.add_subplot(111).plot(bin_edges[0:nrbins]+0.5*dbin, hist[:])

input_arr   = np.log10(totint_tbin_yr_id5)
warr        = totint_warr_id5
hist        = np.histogram(input_arr, bins = bin_edges, weights = warr, normed=False)[0]
fig.add_subplot(111).plot(bin_edges[0:nrbins]+0.5*dbin, hist[:])

plt.yscale('log')

plt.show()




#cumulative rate:
fig = plt.figure(figsize=(5, 4))

input_arr   = np.log10(totint_tbin_yr_id3)
warr        = totint_warr_id3
cumhist     = np.cumsum(np.histogram(input_arr, bins = bin_edges, weights = warr, normed=False)[0])
normfac     = 1.*max(cumhist)
fig.add_subplot(111).plot(bin_edges[0:nrbins]+0.5*dbin, cumhist[:]/normfac, label ='Post-Interaction Mergers', color = 'black', linewidth=1.5)

input_arr   = np.log10(totint_tbin_yr_id5)
warr        = totint_warr_id5
cumhist     = np.cumsum(np.histogram(input_arr, bins = bin_edges, weights = warr, normed=False)[0])
fig.add_subplot(111).plot(bin_edges[0:nrbins]+0.5*dbin, cumhist[:]/normfac, label ='Resonance-Interaction Mergers', color = 'red', linewidth=1.5)

#axis settings:
plt.yscale('log')
plt.title(r'[BH($20 M_{\odot}$), BH($20 M_{\odot}$)] - BH($20 M_{\odot}$)')
plt.xlabel(r'log$_{10}$ time [years]')
plt.ylabel(r'cumulative merger rate')
plt.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)

#Save figure:
plt.savefig('cumu_BHmergers.eps', bbox_inches='tight')

plt.show()

















