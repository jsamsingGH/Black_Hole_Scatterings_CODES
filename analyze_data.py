import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from scipy import integrate as integrate
import matplotlib as mpl
from scipy.optimize import fsolve
from numpy.random import choice


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
#----------------------------------------------------------

#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

#--------------------------------------------------------------
#Sample and Analyze:
#--------------------------------------------------------------
#--------------------------------------
#INPUT:
#--------------------------------------
#input numbers:
obs_time_yrs    = 5.0       #input in YEARS
Ns_tot          = 10000000  #nr of samplings
#input data file:
filename = raw_input('open filename: ')
saveoutput_arr      = np.load(filename + '_saveoutput_arr' + '.npz')['arr_0']
obs_INFO_W_samp_arr = saveoutput_arr
tf = open(filename + '_info_sim_arr.txt', "r")
info_sim_arr        = np.loadtxt(tf, dtype=float)
tf.close()
#define:
[Nsamp, v_kms, vesc_kms, m_Msun, n_nr_pc3, delta_bs, a_HB_SI, a_ej_SI, a_max_SI, a_min_SI] = info_sim_arr[:]
nr_obsbins          = len(obs_INFO_W_samp_arr[:,0])
#--------------------------------------
#--------------------------------------
#incl. obs. time in weight array:
#--------------------------------------
pos_BBHform     = np.where(obs_INFO_W_samp_arr[:,2] == 10)[0]
obs_INFO_W_samp_arr[pos_BBHform,4]  = obs_INFO_W_samp_arr[pos_BBHform,4] + obs_time_yrs
obs_INFO_W_samp_arr[:,4]            = obs_INFO_W_samp_arr[:,4]/sum(obs_INFO_W_samp_arr[:,4])    #normalize 
#--------------------------------------
#SAMPLE:
#--------------------------------------    
pos_obsbins_arr = np.arange(nr_obsbins)
sampweight_arr  = obs_INFO_W_samp_arr[:,4]
pos_weightsampl = np.random.choice(pos_obsbins_arr, Ns_tot, p=sampweight_arr) 
#output arrays:
tinsp_WS    = obs_INFO_W_samp_arr[pos_weightsampl,0]    #in yrs
ecc_WS      = obs_INFO_W_samp_arr[pos_weightsampl,1]
fGW_WS      = obs_INFO_W_samp_arr[pos_weightsampl,3]    #in Hz
mID_WS      = obs_INFO_W_samp_arr[pos_weightsampl,5]    #1 = (2-body merger), 2 = (3-body merger).
#--------------------------------------
#test PLOTs:
#--------------------------------------
plt.figure(figsize=(8,8))

pos_mID1    = np.where(mID_WS == 1)[0]
pos_mID2    = np.where(mID_WS == 2)[0]

plt.subplot(221)
#2-body mergers:
dataplot    = fGW_WS[pos_mID1]
plt.hist(np.log10(dataplot), range=(-3,3), bins=50, histtype='step', stacked=True, fill=False, color = 'black')
#3-body mergers:
dataplot    = fGW_WS[pos_mID2]
plt.hist(np.log10(dataplot), range=(-3,3), bins=50, histtype='step', stacked=True, fill=False, color = 'red')
plt.yscale('log')

plt.subplot(222)
#2-body mergers:
dataplot    = 1.-ecc_WS[pos_mID1]**2
plt.hist(np.log10(dataplot), range=(-8,0), bins=50, histtype='step', stacked=True, fill=False, color = 'black')
#3-body mergers:
dataplot    =  1.-ecc_WS[pos_mID2]**2
plt.hist(np.log10(dataplot), range=(-8,0), bins=50, histtype='step', stacked=True, fill=False, color = 'red')
plt.yscale('log')

plt.show()
#--------------------------------------
exit()
#--------------------------------------------------------------
#--------------------------------------------------------------





















        

