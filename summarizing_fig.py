import matplotlib
from matplotlib import pyplot as plt
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.patches import Ellipse, Polygon
import matplotlib as mpl

#------------------------------------------------------------------------------------------
#Units and conversions:
#------------------------------------------------------------------------------------------
#code units: Rsun, Msun, G=1, ...
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
#------------------------------------------------------------------------------------------

#input val for when the cs for inspirals equals cs for colls:
log_a0oR_coll_eq_insp   = 3.5   #read of from fig in paper.
a_to_R                  = 10**(-log_a0oR_coll_eq_insp) 

#M array for plotting:
M_arr_UMsun = 10**(np.arange(-10.0, 10.0, 0.1))

#to calc SMA HB limit we need to input v_infty:
v_HB_arr    = [1., 10., 100.]   #km/sec
nr_v_HB     = len(v_HB_arr)

#------------------------------------------------------------------------------------------
#PLOT AND MAKE FIGURE:
#------------------------------------------------------------------------------------------
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)

fig = plt.figure(figsize=(6, 5))

#PLOT: Hard-Binary (HB) limits:
for i in range(0,nr_v_HB):
    #read v,m:
    v_HB        = v_HB_arr[i]       #km/sec
    M_arr       = M_arr_UMsun[:]    #mass in M_sun
    #calc HB sma a:
    a_HB_URsun_arr  = AU_U*((36.5**2.)*(M_arr/v_HB**2.))  #HB limit SMA in units of Rsun
    print AU_U*((36.5**2.))
    #convert a_HB to RA_HB (using the insp,coll cross over val):
    R_HB_URsun_arr  = a_to_R*a_HB_URsun_arr 
    #plot M,R HB relation:
    fig.add_subplot(111).plot(M_arr,R_HB_URsun_arr, color='grey', linestyle='-', linewidth = 0.5)  
    
    
#PLOT: WD mass, radius relation:
M_WD_UMsun_arr  = np.arange(0.1,1.4,0.001)
R_WD_URsun_arr  = 0.013*((1.43/M_WD_UMsun_arr)**(1./3.))*((1.-M_WD_UMsun_arr/1.43)**(0.447))
#plot M_WD, R_WD relation:
fig.add_subplot(111).plot(M_WD_UMsun_arr, R_WD_URsun_arr, color='black', linewidth=4.0, linestyle = '-', label=r'White Dwarfs')


#PLOT: MS mass, radius relation:
M_MS_UMsun_arr  = np.arange(0.1,100.,0.001)
R_MS_URsun_arr  = M_MS_UMsun_arr**(0.8)
#plot M_MS, R_MS relation:
fig.add_subplot(111).plot(M_MS_UMsun_arr, R_MS_URsun_arr, color='black', linewidth=4.0, linestyle = '--', label=r'Main-Sequence Stars')


##PLOT: PLANET mass, radius relation:
#M_PL_UMsun_arr  = np.arange(1e-7,1e-3,0.00001)
#R_PL_URsun_arr  = 0.6*M_PL_UMsun_arr**(0.25)
##plot M_PL, R_PL relation:
#fig.add_subplot(111).plot(M_PL_UMsun_arr, R_PL_URsun_arr, color='black', linewidth=4.0, linestyle = ':', label=r'Planets')

#plot settings:
fig.add_subplot(111).set_xlim(1e-2,1e2)
fig.add_subplot(111).set_ylim(1e-3,1e1)

#legends and labels:
fig.add_subplot(111).legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = True)

plt.title(r'Tidal Mergers vs Collisions')
plt.xlabel(r'Mass $m/M$')
plt.ylabel(r'Radius $R/R$')

plt.yscale('log')
plt.xscale('log')

plt.savefig('fig_overview_tidseqcoll.pdf', bbox_inches='tight')

plt.show()
#------------------------------------------------------------------------------------------


























