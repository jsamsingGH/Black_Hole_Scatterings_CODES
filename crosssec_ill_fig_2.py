from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import quad
from sympy.solvers import solve
from sympy import Symbol
import matplotlib as mpl


def intgrand_insp(ap, beta):
    return (4.*ap/(np.sqrt(3.)*(ap**beta)*((ap-1.0)**(3./2.))))**(1./beta)

def intgrand_coll(ap):
    return (1./ap)


#-------------------------------------------------
#Unit conversions:
#-------------------------------------------------
#code units: Rsun, Msun, G=1, ...
M_sun_SI   = 1.989*(10.**30.)   #kg
R_sun_SI   = 695800000.         #m 
AU_SI      = 149597871000.      #m 
G_new_SI   = 6.67*(10.**(-11.))
AU_U       = AU_SI/R_sun_SI                             #from dist AU to code units (U)
kmsec_U    = 1000./np.sqrt(G_new_SI*M_sun_SI/R_sun_SI)  #from vel km/sec to code units (U)
#-------------------------------------------------

#FOR TIDES WE ASSUME beta and epsilon is the same through out the code! 

#-------------------------------------------------
#Estimate normalizations T_i, T_c:
#-------------------------------------------------
#General model params: (USED IN THE REST OF THE CODE)
beta        = 9.0   #tides
beta_GW     = 7./2. #GWs
#Specific sim params: (USED JUST HERE FOR CALIBRATION)
#We take numbers from a real simulation.
#We now use the WD-NS-NS sim since this is probably the most 'clean'
#sim we can have. From that we can also calibrate both tides and GR.
#We ONLY try to calibrate the ASYMPTOTIC SOLUTION.
#The insp cs should therefore be given a higher value than observed due to coll overlab...
R_A         = 0.00608955302883      #R_sun
a_0         = 0.1*AU_U              #R_sun
m_A         = 1.2                   #M_sun
v_inf       = 10.0                  #km/sec
#T_c: (collisions)
nrcomb_fac  = 2.0       #because we use cs insp [12,13] from sim
cs          = 0.135           #AU^2 coll cross section FOR v_inf
T_c         = (1./nrcomb_fac)*cs*((m_A*R_A/v_inf**2.))**(-1.)
#T_i: (inspirals)
insp_asymp_corr_fac = 1.95   #E_A vary very fast with just small changes in this. Thats ok.
nrcomb_fac  = 2.0       #because we use cs insp [12,13] from sim
cs          = 0.135          #AU^2 insp cross section FOR v_inf
T_i         = insp_asymp_corr_fac*(1./nrcomb_fac)*cs*((m_A*R_A/v_inf**2.)*((a_0/R_A)**(1./beta)))**(-1.)
#GW inspirals:
nrcomb_fac  = 1.0       
cs          = 0.01          #AU^2 coll cross section FOR v_inf
R_s_A       = m_A*(2950.0/R_sun_SI) #Schwarzschild radius R_s of m_A in units of R_sun
T_GWi       = (1./nrcomb_fac)*cs*((m_A*R_s_A/v_inf**2.)*((a_0/R_s_A)**(1./beta_GW)))**(-1.)
print T_i, T_c, T_GWi
#NOTES:
#this way of calibrating assumes that epsilon cannot change in the cs below (because epsilon etc. is included in T)
#in other words, we assume the polytropic index is the same for all the objects we here plot. Thats ok.
#-------------------------------------------------

#-------------------------------------------------
#tests:
#-------------------------------------------------
print 5.0*AU_SI/R_sun_SI
print (1.2*2950.0/R_sun_SI)
print 2.*(1.2*0.006/10**2.)*727.*((5.0*AU_SI/R_sun_SI)/0.006)**(1./9.)
print 1.*(1.2*(1.2*2950.0/R_sun_SI)/10**2.)*2095.*((5.0*AU_SI/R_sun_SI)/(1.2*2950.0/R_sun_SI))**(2./7.)
print 2.*(1.2*0.006/10**2.)*924.
print (2./3.)*(1./G_new_SI)*(10.**6)*(10.**(3.5))*(R_sun_SI/M_sun_SI)
print 0.013*((1.43/1.0)**(1./3.))*((1.-1.0/1.43)**(0.447))
print 1.0*(2950.0/R_sun_SI)
#-------------------------------------------------

#-------------------------------------------------
#make cross sections for plotting:
#-------------------------------------------------
#First find value for epsilon_A (E_A):
ll = 1.0
ul = 2.0
I_insp  = quad(intgrand_insp, ll, ul, args=(beta))[0]
I_coll  = quad(intgrand_coll, ll, ul)[0]
E_A     = ((T_i/T_c)*(I_coll/I_insp))**(beta)

print 'E_A: ', E_A

#a/R range arr:
a0_over_RA_arr  = 10.**(np.arange(0.0,20.0,0.1))  #lin sampled in log

#Collisions:
csp_coll            = T_c + 0.0*a0_over_RA_arr  #last term added just to keep the dimensions correct for plotting

#Tidal inspirals:
#CO-TO: (compact object - tidal object)
csp_insp_COTO   = T_i*(a0_over_RA_arr**(1./beta))   #asymptotic solution
#calc corresponding insp-coll area correction:
E_AB            = E_A
nr_a            = len(a0_over_RA_arr)
f_insp_corr_arr = np.zeros(nr_a, dtype=np.float64)
for i in range(0,nr_a):
    print i, nr_a-1
    #define:
    a0_over_RA = a0_over_RA_arr[i]
    x = Symbol('x')
    #solve for when the e_i,e_c cross each other:
    e_i_x   = (E_AB**(1./beta))*(a0_over_RA**((1./beta)-1.0))*((4.*x/(np.sqrt(3.)*(x**beta)*((x-1.0)**(3./2.))))**(1./beta))
    e_c_x   = (a0_over_RA**(-1.))*(1./x)
    x_sol   = solve(e_i_x - e_c_x, x)[0]
    if (x_sol > ul): x_sol = ul #formalism only valid for ap up to around 2 (ul).
    #fractional change in nr. inspirals:
    ul_x        = float(x_sol)
    insp_A_fac  = (E_AB**(1./beta))*(a0_over_RA**((1./beta)-1.0))
    coll_A_fac  = (a0_over_RA**(-1.))
    A_insp_x    = insp_A_fac*quad(intgrand_insp, ll, ul_x, args=(beta))[0]  #area insp from 1 to ul_x
    A_coll_x    = coll_A_fac*quad(intgrand_coll, ll, ul_x)[0]               #area coll from 1 to ul_x
    A_insp_full = insp_A_fac*quad(intgrand_insp, ll, ul, args=(beta))[0]    #area insp from 1 to ul (2): full range. Used for normalization
    frac_insp_corr  = (A_insp_x-A_coll_x)/A_insp_full    #the fraction the asymptotic solution should be corrected with.
    f_insp_corr_arr[i] = frac_insp_corr    
#define corrected insp cross section:
csp_insp_COTO_CORR  = csp_insp_COTO*f_insp_corr_arr

#GW inspirals:
csp_insp_GW         = T_GWi*(a0_over_RA_arr**(1./beta_GW))

#Save txt:
savearr = np.zeros((4,len(a0_over_RA_arr)), dtype=np.float64)
savearr[0,:] = a0_over_RA_arr
savearr[1,:] = csp_insp_COTO_CORR
savearr[2,:] = csp_insp_GW
savearr[3,:] = csp_coll
tf = open('CS_analytical_test.txt',  "w")
np.savetxt(tf, savearr,   fmt='%5f')
tf.close()
#-------------------------------------------------








#-------------------------------------------------
#Plot:
#-------------------------------------------------
#Setting general font:
font = {'family' : 'serif'}
mpl.rc('font', **font)


fig, ax1 = plt.subplots(figsize=(5,4))
#ax2 = ax1.twinx()

#y limits on ax1,ax2 must match:
ax1_ylim    = [5e0,1e7]
#ax2_ylim    = ax1_ylim/csp_coll[0]
ax1.set_ylim(ax1_ylim)
#ax2.set_ylim(ax2_ylim)

#set x limits:
ax1_xlim    = [0,10]
ax1.set_xlim(ax1_xlim)

#oplot guide lines:
#to calc SMA HB limit we need to input mass, v_infty:
v_HB = 10.0 #km/sec
m_HB = 1.0  #mass
a_HB_limit = AU_U*((36.5**2.)*(m_HB/v_HB**2.))  #HB limit SMA in units of Rsun
#objs for illustration:
R_sun   = 1.0                                                           # m = 1
R_WD    = 0.013*((1.43/m_HB)**(1./3.))*((1.-m_HB/1.43)**(0.447))        # m = 1 MUST MATCH with m_HB above.     
R_S     = (2950.0/R_sun_SI)                                             # m = 1 MUST MATCH with m_HB above: Schwarzschild radius of an object with mass 1Msun    
print R_sun, R_WD, R_S
#plot:
xval    = np.log10(a_HB_limit/R_sun)
ax1.plot([xval,xval], ax1_ylim, color='black', linestyle = ":")
ax1.text(xval+0.2, 60, r'HB-limit (sun)', size = 10,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=-90)

xval    = np.log10(a_HB_limit/R_WD)
ax1.plot([xval,xval], ax1_ylim, color='black', linestyle = ":")
ax1.text(xval+0.2, 60, r'HB-limit (WD)', size = 10,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=-90)

xval    = np.log10(a_HB_limit/R_S)
ax1.plot([xval,xval], ax1_ylim, color='red', linestyle = ":")
ax1.text(xval+0.2, 60, r'HB-limit (CO)', size = 10,
        horizontalalignment='center',
        verticalalignment='center',
        rotation=-90)


#x-axis plot:
px = np.log10(a0_over_RA_arr)

y   = csp_insp_COTO_CORR
py  = (y) 
ax1.plot(px, py,    linewidth=1.5, color='black', linestyle = "-", label=r'Tidal Inspirals (coll corr)')

y   = csp_insp_COTO
py  = (y) 
ax1.plot(px, py,    linewidth=0.5, color='black', linestyle = "-", label=r'Tidal Inspirals (asym sol)')

y   = csp_coll
py  = (y) 
ax1.plot(px, py,    linewidth=2.5, color='black', linestyle = "--", label=r'Collisions')

y   = csp_insp_COTO_CORR+csp_coll
py  = (y) 
ax1.plot(px, py,    linewidth=1.5, color='black', linestyle = "--", label=r'Collisions + Tidal Inspirals')

y   = csp_insp_GW
py  = (y) 
ax1.plot(px, py,    linewidth=2.5, color='red', linestyle = "-.", label=r'GW Inspirals (asym sol)')
#print np.log10(a_0/R_s_A)

#oplot legend:
ax1.legend(loc='upper left', numpoints = 1, fontsize = 10.0, frameon = False)

#Title/labels:
ax1.set_title(r'Inspiral and Collision Cross Sections')
ax1.set_xlabel(r'log  $(a_{0}/\mathcal{R})$')
ax1.set_ylabel(r'$\sigma$ $[$AU$^2]$ $[(m/M_{\odot})(\mathcal{R}/R_{\odot})/(v_{\infty}/$km$\ $s$^{-1})^2]$')
#ax2.set_ylabel(r'$\sigma/\sigma_{coll}$')

#dytext  = 0.07
#ytext   = 0.275
#xtext   = 0.70
#ax1.text(xtext, ytext,r'Tides:',
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax1.transAxes)
#ax1.text(xtext, ytext-(1.-0.3)*dytext,r'$\beta=3$, $n=3$',
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax1.transAxes)
#ax1.text(xtext, ytext-2*dytext,r'GR:',
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax1.transAxes)
#ax1.text(xtext, ytext-(3.-0.3)*dytext,r'$\beta=7/2$',
#     horizontalalignment='left',
#     verticalalignment='center',
#     transform = ax1.transAxes)

#ax1.text(0.02, 0.72, r'HB-limit: $3M_{\odot}$, $10$ $km\ s^{-1}$',
#     horizontalalignment='left',
#     verticalalignment='center',
#     fontsize = 8.0,
#     transform = ax1.transAxes)

#ax1.text(0.025, 0.67, r'$m_{tot}    =3m_{\odot}$',
#     horizontalalignment='left',
#     verticalalignment='center',
#     fontsize = 8,
#     transform = ax1.transAxes)
#ax1.text(0.025, 0.64, r'$v_{\infty} =10$ km/sec',
#     horizontalalignment='left',
#     verticalalignment='center',
#     fontsize = 8,
#     transform = ax1.transAxes)
     
#set log:
ax1.set_yscale('log')
#ax2.set_yscale('log')

#save fig:
plt.savefig('crossec_summaryfig.eps', bbox_inches='tight')
plt.show()
#-------------------------------------------------




exit()





