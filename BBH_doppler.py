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
from scipy.integrate import solve_ivp
from scipy.signal import find_peaks
from mpl_toolkits.axes_grid.inset_locator import (inset_axes, InsetPosition,
                                                  mark_inset)



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
tU_tSI     = np.sqrt(R_sun_SI**3./(G_new_SI*M_sun_SI))    
Rsch_1Msun_unitRsun = ((2.*G_new_SI*(1.*M_sun_SI))/(c_SI**2.))/R_sun_SI
#----------------------------------------------------------\





m_1     = (10.)*M_sun_SI    # merging bin   1
m_2     = (10.)*M_sun_SI     # merging bin   2
m_3     = (10.)*M_sun_SI    # single        3
m_tot   = m_1 + m_2 + m_3
mr_12   = (m_1*m_2)/(m_1+m_2)
f_GW    = 10
R_dist  = 0.1*AU_SI#100.*(10000000.)*Rsch_1Msun_unitRsun*R_sun_SI  #1.*AU_SI

const_k = (768./425.)*(5./256.)*(2.**(7./2.))*((6.*np.sqrt(2.))/(85.*np.pi))**(1./2.)
const_1 = const_k/(2.*(np.pi)**(9./2.))*(c_SI**(15./2.)/(G_new_SI**2.))

tm_Torb = const_1*(1./(f_GW**(7./2.)))*(1./(m_1*m_2))*((m_tot/mr_12)**(1./2.))*(1./(R_dist**(3./2.)))


print tm_Torb





















 
 
 
 
 
 
 
 
 
 
 


