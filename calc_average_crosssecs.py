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

#calc average cs or err of the following datasets:

#calc using new rptid = 1 threshold:

#name_dataset    = ['SIM1_BH102020_NT_WGR_0001_10_1kmsec_T2500', 'SIM2_BH102020_NT_WGR_0001_10_1kmsec_T2500', 'SIM3_BH102020_NT_WGR_0001_10_1kmsec_T2500']
#data_name       = 'average_BH102020_NT_WGR_0001_10_1kmsec_T2500'

#name_dataset    = ['Test_NIC_NS14NS14NS14_2', 'Test_NIC_NS14NS14NS14_3', 'Test_NIC_NS14NS14NS14_4']
#data_name       = 'average_Test_NIC_NS14NS14NS14'

#name_dataset    = ['Test_NIC_NS14NS14NS14_A1', 'Test_NIC_NS14NS14NS14_A2', 'Test_NIC_NS14NS14NS14_A3', 'Test_NIC_NS14NS14NS14_A4']
#data_name       = 'average_Test_NIC_NS14NS14NS14_A1234'

#name_dataset    = ['Test_NIC_NS14NS14NS14_C1', 'Test_NIC_NS14NS14NS14_C2', 'Test_NIC_NS14NS14NS14_C3']
#data_name       = 'average_Test_NIC_NS14NS14NS14_C123'

name_dataset    = ['NIC_WD06NS14NS14_A1', 'NIC_WD06NS14NS14_A2', 'NIC_WD06NS14NS14_A3', 'NIC_WD06NS14NS14_A4', 'NIC_WD06NS14NS14_A5']
data_name       = 'average_NIC_WD06NS14NS14_A12345'


#name_dataset    = ['SIM1_BH102020_NT_WGR_0001_1_1kmsec_T2500', 'SIM2_BH102020_NT_WGR_0001_1_1kmsec_T2500']
#data_name       = 'average_BH102020_NT_WGR_0001_1_1kmsec_T2500'


#name_dataset    = ['SIM2_WD10_NSNS_WT_WGR_0001_01_1kmsec_T1000', 'SIM3_WD10_NSNS_WT_WGR_0001_01_1kmsec_T1000']
#data_name       = 'average_WD10_NSNS_WT_WGR_0001_01_1kmsec_T1000'


#name_dataset    = ['SIM1_WD10_NSNS_NT_WGR_0001_01_1kmsec', 'SIM2_WD10_NSNS_NT_WGR_0001_01_1kmsec']
#data_name       = 'average_WD10_NSNS_NT_WGR_0001_01_1kmsec'



#define:
nr_dataset      = len(name_dataset) 
cs_arr_dataset       = []
cs_err_arr_dataset   = []

#open data:
for nd in range(0,nr_dataset):
    #open data:
    tf = open('cs_data_'        + name_dataset[nd], "r")
    cs_arr_data        = np.loadtxt(tf, dtype=float)
    tf.close()
    tf = open('cs_err_data_'    + name_dataset[nd], "r")
    cs_err_arr_data    = np.loadtxt(tf, dtype=float)
    tf.close()
    #save data:
    cs_arr_dataset.append(cs_arr_data)
    cs_err_arr_dataset.append(cs_err_arr_data)
      
#calc average:
nr_a        = len(cs_arr_dataset[0][0,:])
nr_cs       = 20    #dont change. We just always keep space for ...
save_arr_cs     = np.zeros((nr_cs,nr_a), dtype=np.float64)
save_arr_cs_err = np.zeros((nr_cs,nr_a), dtype=np.float64)
for nd in range(0,nr_dataset):
    save_arr_cs[:,:]        = save_arr_cs[:,:]      + (cs_arr_dataset[nd][:,:])
    save_arr_cs_err[:,:]    = save_arr_cs_err[:,:]  + (cs_err_arr_dataset[nd][:,:])**(2.)
    
save_arr_cs     = (1./(1.*nr_dataset))*save_arr_cs
save_arr_cs_err = (1./(1.*nr_dataset))*np.sqrt(save_arr_cs_err)

#Save average cross section data:
#save cs:
tf = open('cs_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs,   fmt='%5f')
tf.close()
#save cs err:
tf = open('cs_err_data_' + data_name, "w")
np.savetxt(tf, save_arr_cs_err,   fmt='%5f')
tf.close()
