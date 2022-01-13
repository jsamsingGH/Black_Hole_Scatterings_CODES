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
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI 
#----------------------------------------------------------


#----------------------------------------------------------
#SPLIT data:
#----------------------------------------------------------
#data path
path_data           = '/Users/jsamsing/Desktop/TIDES_PROJ/MOCCA_ICs/MC_1/'

#input data name:
name_input_data     = path_data + 'inter-binsin-bhs-13gyr-m100.dat'

#number of lines input data:
f_input     = open(name_input_data, "r")  #OPEN input FILE.
nr_lines    = len(f_input.readlines())

#split file:
name_output_data    = path_data + 'MT1000_inter-binsin-bhs-13gyr-m100.dat'  #'MOCCATEST_inter-binsin-bhs-13gyr-m100.dat'
f_input     = open(name_input_data,  "r")  #OPEN input FILE.
f_output    = open(name_output_data, "w")   
for nl in range(0, 100):
    f_line = f_input.readline()
    f_output.write(f_line)
#close output file:    
f_output.close()

exit()
#----------------------------------------------------------











#----------------------------------------------------------
#TEST open file lines:
#----------------------------------------------------------
name_data   = path_data + 'TESTset1_inter-binsin-bhs-13gyr-m100.dat'
#get total number of lines:
f_data      = open(name_data, "r")
nr_lines    = len(f_data.readlines())
#analyze line per line:
f_data      = open(name_data, "r")
for nl in range(0, nr_lines):
    line        = f_data.readline()
    line_split  = re.split(r'\s{1,}', line)
    #get info from line:
    print float(line_split[110])*float(line_split[110])

exit()
#----------------------------------------------------------




























