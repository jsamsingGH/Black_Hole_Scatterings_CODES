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
from scipy.interpolate import make_interp_spline, BSpline


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



m_BH    = 60.*M_sun_SI      #kg
m_NS    = 30.4*M_sun_SI      #kg

v_disp  = 10.*(10.**3.)     #m/s
n_BHs   = 100000.*(1./(m_parsec**3.)) #nr per m^-3
#f_bs = ...
#n_bs    = 

r_cap   = (((85.*np.pi)/(6.*np.sqrt(2.)))**(2./7.))*G_new_SI*(m_BH**(2./7.))*(m_NS**(2./7.))*((m_BH + m_NS)**(3./7.))/((c_SI**(10./7.))*(v_disp**(4./7.)))
r_bin   = 1.*AU_SI 

s_cap   = np.pi*(r_cap**2.)*(1. + (2.*G_new_SI*(m_BH + m_NS))/(r_cap*(v_disp**2.)))
s_bin   = np.pi*(r_bin**2.)*(1. + (2.*G_new_SI*(m_BH + m_NS))/(r_bin*(v_disp**2.)))

 
t_cap   = 1./(n_BHs*s_cap*v_disp)
t_bin   = 1./(n_BHs*s_bin*v_disp)

#SOMETHING WRONG??? LONG TIMESCALES???



print t_cap/sec_year
print t_bin/sec_year


#print (3./2.)*(G_new_SI*m_BH/(v_disp**2.))/AU_SI
#print (1./6.)*(1./(7./9.) - 1.)*(G_new_SI*m_BH/((5.*v_disp)**2.))/AU_SI
































