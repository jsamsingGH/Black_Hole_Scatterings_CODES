#!/usr/bin/python

from matplotlib import pyplot as plt
import numpy as np
from scipy.integrate import odeint
import sys
from itertools import combinations
import itertools
from scipy import linalg
from sympy import mpmath



#------------------------------------------------------------------------------------------
#FUNCTIONS:
#------------------------------------------------------------------------------------------



def func_internal_properties(const_in, q_in, qdot_in):   
    #------------------------
    # Define:
    #------------------------
    #UNPACK const_in:
    mass        = const_in[0]
    radius      = const_in[1]
    gas_n       = const_in[2]
    gas_gamma   = const_in[3]
    Mqp_sph     = const_in[4] 
    #q, qdot:
    q_MAT       = np.matrix(q_in)            
    q_MAT_T     = q_MAT.T 
    q_absdet    = np.abs(linalg.det(q_MAT))
    qdot_MAT    = np.matrix(qdot_in)            
    S_MAT       = q_MAT*q_MAT_T
    S_absdet    = np.abs(linalg.det(S_MAT))
    #------------------------
    # Calc:
    #------------------------
    #kinetic energy:
    E_kinetic   = (1./2.)*Mqp_sph*np.trace(qdot_MAT*qdot_MAT.T)
    #Gravitational self energy:
    OM_sph      = - (3./(5.-gas_n))*((mass**2.)/radius)
    A_MAT       = func_selfgrav_A_matrix(q_MAT)  
    OM_MAT      = (1./2.)*OM_sph*(S_absdet**(-1./2.))*(A_MAT*S_MAT)
    E_selfgrav  = np.trace(OM_MAT)
    #Gas energy
    E_gas       = - (OM_sph/(3.*(gas_gamma-1.)))*(q_absdet**(1-gas_gamma))
    #------------------------
    #Return:
    #------------------------
    return [E_kinetic, E_selfgrav, E_gas]
    #------------------------
    




def func_selfgrav_A_matrix(q_in):   
    #------------------------
    # Define:
    #------------------------
    q_MAT       = np.matrix(q_in)            
    q_MAT_T     = q_MAT.T 
    S_MAT       = q_MAT*q_MAT_T
    S_absdet    = np.abs(linalg.det(S_MAT))
    #------------------------
    # Calc:
    #------------------------
    Eigen_vals_S, Eigen_vecs_S = np.linalg.eigh(S_MAT)
    A_MAT_D     = np.matrix(np.zeros((3,3)))
    R_MAT       = np.matrix(Eigen_vecs_S)
    R_MAT_inv   = R_MAT.I
    s1 = Eigen_vals_S[0]
    s2 = Eigen_vals_S[1]
    s3 = Eigen_vals_S[2]
    for sc in range(0,3):
        si = Eigen_vals_S[sc]
        CSEI_s1_s2_s3_si = mpmath.elliprj(s1, s2, s3, si) #Carlson symmetric elliptic integral (CSEI) of the third kind.
        A_Di = (S_absdet**(1./2.))*((2./3.)*CSEI_s1_s2_s3_si)
        A_MAT_D[sc,sc] = A_Di
    A_MAT = R_MAT*A_MAT_D*R_MAT_inv
    #------------------------
    #Return:
    #------------------------    
    return A_MAT
    #------------------------


#------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------
class C_1Body(object):
#------------------------------------------------------------------------------------------    
    def __init__(self, const, posCM, velCM, q, qdot):
          
        # evolution vector:
        vec_posCM       = np.ravel(posCM)
        vec_velCM       = np.ravel(velCM)
        vec_q           = np.ravel(q)
        vec_qdot        = np.ravel(qdot)
        self.vec_Yevo   = np.concatenate([vec_posCM, vec_velCM, vec_q, vec_qdot])     
    
        #constant vector:
        self.vec_const  = np.ravel(const)
#------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------
class C_NBsystem(list):
#------------------------------------------------------------------------------------------
        
    #------------------------------------------------------------------
    def M_evolve_NBsystem(self,t_start,t_end,nr_t_out):
    #------------------------------------------------------------------    
        #make Yevo IC for Nbody system:
        Yevo_NBsystem_IC        = self.M_make_Yevo_NBsystem
        
        #make vec with constants for Nbody system:
        const_NBsystem          = self.M_Make_const_NBsystem
        
        #make time array for analysis output:
        time_arr_output         = np.linspace(t_start,t_end,nr_t_out)
        
        #Evolve NBsystem:
        func_Y_in_Ydot_out      = self.M_Ydot
        Evolve_NBsystem_result  = odeint(func_Y_in_Ydot_out, Yevo_NBsystem_IC, time_arr_output, args=(const_NBsystem,), full_output=0, atol=1e-12, rtol=1e-12)

        #Return result:
        return Evolve_NBsystem_result
    #------------------------------------------------------------------   
    
    #------------------------------------------------------------------    
    #make evolution vector for 'Nbody system':
    #------------------------------------------------------------------
    @property
    def M_make_Yevo_NBsystem(self):
        Yevo_NBsystem = []
        for body_n in self:
            Yevo_NBsystem.append(body_n.vec_Yevo)
        return np.ravel(np.array(Yevo_NBsystem))
    #------------------------------------------------------------------
    
    #------------------------------------------------------------------
    #make vector with constants for 'Nbody system':
    #------------------------------------------------------------------
    @property
    def M_Make_const_NBsystem(self):
        const_NBsystem = []
        for body_n in self:
            const_NBsystem.append(body_n.vec_const)
        return np.ravel(np.array(const_NBsystem))
    #------------------------------------------------------------------     
    
    #------------------------------------------------------------------    
    #Make Ydot_out from Y_in (Y is the evo vector for the NBsystem)
    #------------------------------------------------------------------
    def M_Ydot(self, Y_in, t, const_NBsystem):
    
        #define:
        #--------------------------------
        nr_body                 = len(self)        
        length_Yevo_per_body    = len(Y_in)/nr_body
        length_const_per_body   = len(const_NBsystem)/nr_body
        
        Ydot_out_body_n_all     = np.zeros((nr_body,length_Yevo_per_body))
        I_MAT                   = np.matrix(np.identity(3)) #3x3 identity matrix.
        #--------------------------------
        
        #UNPACK:
        #--------------------------------
        #unpack Y_in ([vec_posCM(3),vec_velCM(3),vec_q(3x3),vec_qdot(3x3),....] per body):    
        Y_in_body_n_all         = np.reshape(Y_in,(nr_body,length_Yevo_per_body))
        Y_in_body_n_pos         = np.reshape(Y_in_body_n_all[:,0:3],   (nr_body,3))
        Y_in_body_n_vel         = np.reshape(Y_in_body_n_all[:,3:6],   (nr_body,3))
        Y_in_body_n_q           = np.reshape(Y_in_body_n_all[:,6:15],  (nr_body,3,3))
        Y_in_body_n_qdot        = np.reshape(Y_in_body_n_all[:,15:24], (nr_body,3,3))
        #unpack const_NBsystem ([mass,radius,n,gamma,Mqp,...] per body):
        const_body_n_all            = np.reshape(const_NBsystem,(nr_body,length_const_per_body))
        const_body_n_mass           = const_body_n_all[:,0]
        const_body_n_radius         = const_body_n_all[:,1]
        const_body_n_gas_n          = const_body_n_all[:,2]
        const_body_n_gas_gamma      = const_body_n_all[:,3]
        const_body_n_Mqp_sph        = const_body_n_all[:,4]
        const_body_n_evoTides_yesno = const_body_n_all[:,5]
        const_body_n_RigidSph_yesno = const_body_n_all[:,6]
        #-------------------------------- 
        
        #Defin arrays:
        #--------------------------------
        acc_CM_pointmass_NT_ij          = np.zeros((3))
        acc_CM_tidalfield_ij            = np.zeros((3))                    
        body_i_acc_CM_sum_j             = np.zeros((3))
        body_i_acc_tidal_tensor_C_sum_j = np.zeros((3,3))                
        D_MAT_xk                        = np.zeros((3,3,3)) # the d(XX)/dx_k matrix
        #--------------------------------
        

        #-begin loop i (body i):
        #-------------------------------------------------------------        
        for i in range(0,nr_body):
                                               
            #INFO body i:         
            #--------------------------------
            #intrinsic params (constants):
            mass_i              = const_body_n_mass[i]
            radius_i            = const_body_n_radius[i] 
            gas_n_i             = const_body_n_gas_n[i]
            gas_gamma_i         = const_body_n_gas_gamma[i]
            Mqp_sph_i           = const_body_n_Mqp_sph[i]            
            evoTides_yesno_i    = const_body_n_evoTides_yesno[i]    #yes(1)/no(0)
            RigidSph_yesno_i    = const_body_n_RigidSph_yesno[i]    #yes(1)/no(0)
            #time dependent (pos,vel,q,qdot,...):
            pos_i       = Y_in_body_n_pos[i,:]  #pos
            vel_i       = Y_in_body_n_vel[i,:]  #vel
            q_MAT_i     = np.matrix(Y_in_body_n_q[i,:,:])
            qdot_MAT_i  = np.matrix(Y_in_body_n_qdot[i,:,:])            
            q_MAT_i_T   = q_MAT_i.T 
            S_MAT_i     = q_MAT_i*q_MAT_i_T
            #--------------------------------
                        
            #--------------------------------
            # EOM for Center-of-Mass of body i:
            #--------------------------------
            #initialize:
            body_i_acc_CM_sum_j[:]  = 0.
            #-begin loop j (body j):
            #--------------------------------
            for j in range(0,nr_body):         
                if (i != j):                                       
                    
                    #------------------------
                    #INFO body j:         
                    #------------------------
                    #intrinsic params (constants):
                    mass_j              = const_body_n_mass[j]
                    Mqp_sph_j           = const_body_n_Mqp_sph[j]
                    evoTides_yesno_j    = const_body_n_evoTides_yesno[j]    #yes(1)/no(0)
                    RigidSph_yesno_j    = const_body_n_RigidSph_yesno[j]    #yes(1)/no(0)
                    #time dependent (pos,vel,q,qdot,...):
                    pos_j       = Y_in_body_n_pos[j,:]  #pos  
                    vel_j       = Y_in_body_n_vel[j,:]  #vel
                    q_MAT_j     = np.matrix(Y_in_body_n_q[j,:,:])
                    q_MAT_j_T   = q_MAT_j.T 
                    S_MAT_j     = q_MAT_j*q_MAT_j_T
                    #------------------------
                   
                    #------------------------
                    #Define:
                    #------------------------
                    pos_vec_ij  = (pos_i - pos_j) 
                    vel_vec_ij  = (vel_i - vel_j)
                    r_ij        = np.sqrt(pos_vec_ij.dot(pos_vec_ij))
                    v_ij        = np.sqrt(vel_vec_ij.dot(vel_vec_ij))
                    #------------------------
                                   
                    #------------------------
                    # Center-of-mass acc ij:
                    #------------------------                   
                    #------------------
                    #Point Particle contributions:
                    #------------------
                    acc_CM_pointmass_NT_ij[:]  = - (mass_j/(r_ij**2.))*(pos_vec_ij[:]/r_ij)                                    
                    #------------------
                    #Tidal field contribution:
                    #--------------------------------------
                    if (RigidSph_yesno_i == 1 and RigidSph_yesno_j == 1):   #both objects (i,j) are rigid and spherical.
                    #--------------------------------------
                        acc_CM_tidalfield_ij[:] = 0.0
                    #--------------------------------------
                    if (RigidSph_yesno_i == 0 or RigidSph_yesno_j == 0 ):   #at least one of the objects (i,j) are not rigid and spherical.
                    #--------------------------------------    
                        #D_MAT_xk = d(xx)/dx_k:       
                        D_MAT_xk[0,:,:] = [[2*pos_vec_ij[0],pos_vec_ij[1],pos_vec_ij[2]],[pos_vec_ij[1],0,0],[pos_vec_ij[2],0,0]]
                        D_MAT_xk[1,:,:] = [[0,pos_vec_ij[0],0],[pos_vec_ij[0],2*pos_vec_ij[1],pos_vec_ij[2]],[0,pos_vec_ij[2],0]]
                        D_MAT_xk[2,:,:] = [[0,0,pos_vec_ij[0]],[0,0,pos_vec_ij[1]],[pos_vec_ij[0],pos_vec_ij[1],2*pos_vec_ij[2]]]
                        #M_MAT = (m_1M_2S_2 + m_2M_1S_1):
                        M_MAT   = (mass_i*Mqp_sph_j*S_MAT_j + mass_j*Mqp_sph_i*S_MAT_i)                    
                        #XX = (xx^T):
                        X_MAT   = np.matrix(pos_vec_ij)
                        XX_MAT  = X_MAT.T*X_MAT
                        #sum over each pos coordinate (x,y,z):
                        for k in range(0,3):
                            xk          = pos_vec_ij[k]
                            D_MAT       = np.matrix(D_MAT_xk[k,:,:])
                            Cp_MAT_xk   = 3.*(1./(r_ij**5.))*(D_MAT - 5.*xk*(XX_MAT/(r_ij**2.)) + xk*I_MAT)                     #Tidal Deviation Tensor (Cijk)
                            Cp_M_sum_xk = np.trace(Cp_MAT_xk*M_MAT)
                            #final CM-acc from tides:  
                            acc_CM_tidalfield_ij[k] = (1./mass_i)*(1./2.)*Cp_M_sum_xk 
                    #--------------------------------------                       
                    #------------------                               
                    #calc total CM acc:
                    #------------------
                    acc_tot_CM_ij           = acc_CM_pointmass_NT_ij + acc_CM_tidalfield_ij
                    #save and sum terms:
                    #------------------
                    body_i_acc_CM_sum_j[:]  = body_i_acc_CM_sum_j[:] + acc_tot_CM_ij[:]  
                    #------------------------                                           
            
            #-end loop j.
            #--------------------------------
            
            
            #--------------------------------
            # EOM for TIDES (body i):
            #--------------------------------
            if (evoTides_yesno_i == 1):
            #--------------------------------               
                #------------------------
                # External terms
                #------------------------
                #initialize:
                body_i_acc_tidal_tensor_C_sum_j[:,:]    = 0.                
                #-begin loop j (body j):
                for j in range(0,nr_body):         
                    if (i != j):                      
                        #INFO body j:         
                        pos_j       = Y_in_body_n_pos[j,:]
                        mass_j      = const_body_n_mass[j]
                        #Define:
                        pos_vec_ij  = (pos_i - pos_j) 
                        r_ij        = np.sqrt(pos_vec_ij.dot(pos_vec_ij))
                        # For q:
                        X_MAT = np.matrix(pos_vec_ij)
                        Tidal_Tensor_C_body_ij                  = (1./(r_ij**3.))*(3.*(X_MAT.T*X_MAT)/(r_ij**2.) - I_MAT)   #Tidal Tensor (Cij)
                        acc_Tidal_Tensor_C_body_ij              = mass_j*Tidal_Tensor_C_body_ij
                        body_i_acc_tidal_tensor_C_sum_j[:,:]    = body_i_acc_tidal_tensor_C_sum_j[:,:] + acc_Tidal_Tensor_C_body_ij[:,:] 
                #------------------------                
                #------------------------
                # Internal terms
                #------------------------         
                q_MAT_i_inv     = q_MAT_i.I                
                q_MAT_i_T_inv   = q_MAT_i_T.I
                q_i_absdet      = np.abs(linalg.det(q_MAT_i))
                S_i_absdet      = np.abs(linalg.det(S_MAT_i))
                #------------------------
                #calc self gravity tensor:
                #------------------------
                selfgrav_OM_sph = - (3./(5.-gas_n_i))*((mass_i**2.)/radius_i)
                A_MAT           = func_selfgrav_A_matrix(q_MAT_i)  
                selfgrav_OM_MAT = (1./2.)*selfgrav_OM_sph*(S_i_absdet**(-1./2.))*(A_MAT*S_MAT_i)
                #------------------------
                #calc pressure integral:
                #------------------------
                pressure_int_PI_sph = - selfgrav_OM_sph/3. #from virial theorem.
                pressure_int_PI     = pressure_int_PI_sph*(q_i_absdet**(1.-gas_gamma_i)) 
                #------------------------              
                #------------------------
                # FINAL qdotdot: EOM for q_i:
                #------------------------
                qdotdot_MAT_body_i  = (1./Mqp_sph_i)*(pressure_int_PI)*q_MAT_i_T_inv + (1./Mqp_sph_i)*(selfgrav_OM_MAT)*q_MAT_i_T_inv + (body_i_acc_tidal_tensor_C_sum_j)*q_MAT_i
                #------------------------           
            #--------------------------------
        
            
            #--------------------------------
            # if: DO NOT incl TIDES: (body i)
            #--------------------------------
            if (evoTides_yesno_i == 0):
                qdotdot_MAT_body_i  = np.zeros((3,3))
            #--------------------------------
         
                             
            #--------------------------------
            #PACK and SAVE Ydot OUT:
            #--------------------------------
            Ydot_out_body_n_all[i,0:3]    = np.ravel(vel_i)                 #xdot
            Ydot_out_body_n_all[i,3:6]    = np.ravel(body_i_acc_CM_sum_j)   #xdotdot
            Ydot_out_body_n_all[i,6:15]   = np.ravel(qdot_MAT_i)            #qdot
            Ydot_out_body_n_all[i,15:24]  = np.ravel(qdotdot_MAT_body_i)    #qdotdot
            #--------------------------------
                
        #-end loop i.    
        #-------------------------------------------------------------
             
        #-------------------------------------------------------------
        #make Ydot_out vec and RETURN:
        #-------------------------------------------------------------
        Ydot_out = np.ravel(np.array(Ydot_out_body_n_all))
        return Ydot_out
        #-------------------------------------------------------------
        
    #------------------------------------------------------------------          
                 
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------









#------------------------------------------------------------------------------------------




#------------------------------------------------------------------------------------------
#OPEN IC FILE:
#------------------------------------------------------------------------------------------
IC_file = open('input_Nbody_and_Tides.txt', 'r')
#---------------------------
#nr particles:
#---------------------------
fline = IC_file.readline()
nr_body = int(fline[0])
#---------------------------
#Nbody solver params:
#---------------------------
#params list 1:
fline_split = IC_file.readline().split()
Nbody_solver_params_1_INT = np.array([int(i) for i in fline_split[0:10]])
#params list 2:
fline_split = IC_file.readline().split()
Nbody_solver_params_2_REAL = np.array([float(i) for i in fline_split[0:10]])
#---------------------------
#define arrays:
#---------------------------
IC_body_i_const_vec = np.zeros((nr_body,10))
IC_body_i_pos_vec = np.zeros((nr_body,3))
IC_body_i_vel_vec = np.zeros((nr_body,3))
IC_body_i_q_3x3 = np.zeros((nr_body,3,3))
IC_body_i_qdot_3x3 = np.zeros((nr_body,3,3))
#---------------------------
#read info per body_i:
#---------------------------
for nbs in range(0,nr_body):
    #constants:
    fline_split = IC_file.readline().split()
    IC_body_i_const_vec[nbs, :] = np.array([float(i) for i in fline_split[0:10]])
    #pos:
    fline_split = IC_file.readline().split()
    IC_body_i_pos_vec[nbs, :] = np.array([float(i) for i in fline_split[0:3]])
    #vel:
    fline_split = IC_file.readline().split()
    IC_body_i_vel_vec[nbs, :] = np.array([float(i) for i in fline_split[0:3]])
    #q matrix:
    for rc in range(0,3):
        fline_split = IC_file.readline().split()
        IC_body_i_q_3x3[nbs,rc,:] = np.array([float(i) for i in fline_split[0:3]])
    #qdot matrix:
    for rc in range(0,3):
        fline_split = IC_file.readline().split()
        IC_body_i_qdot_3x3[nbs,rc,:] = np.array([float(i) for i in fline_split[0:3]])
#---------------------------
IC_file.close()
#------------------------------------------------------------------------------------------





#------------------------------------------------------------------------------------------
#SIMULATE with input ICs:
#------------------------------------------------------------------------------------------
#---------------------------
#make list of body-objects:
#---------------------------
list_Nbody = []
for nbs in range(0,nr_body):
    #define:
    const_vec   = IC_body_i_const_vec[nbs,:]
    pos_vec     = IC_body_i_pos_vec[nbs,:] 
    vel_vec     = IC_body_i_vel_vec[nbs,:]
    q_3x3       = IC_body_i_q_3x3[nbs,:,:]
    qdot_3x3    = IC_body_i_qdot_3x3[nbs,:,:]
    #create body-obj and save:
    body_i = C_1Body(const_vec, pos_vec, vel_vec, q_3x3, qdot_3x3)    #call class 'C_1Body'
    list_Nbody.append(body_i)
#---------------------------    
#Create Nbody system:
#---------------------------
NBsystem    = C_NBsystem(list_Nbody)
#---------------------------
#Evolve system:
#---------------------------
#odeint (or odepack) does not offer event handling,
#we therefore have to look for events in 'time-steps':


nrtimes_output  = 1000
print 'integration start...'
int_result_NBsystem = NBsystem.M_evolve_NBsystem(0,25,nrtimes_output)
print 'integration stop...'




#------------------------------------------------------------------------------------------









#------------------------------------------------------------------------------------------
#Analyze:
#------------------------------------------------------------------------------------------
len_vec_Yevo                    = len(list_Nbody[0].vec_Yevo)  #same length for all objects.
len_vec_const                   = 10#len(const_1)      #same length for all objects.
int_result_NBsystem_unpacked    = np.reshape(int_result_NBsystem,(nrtimes_output, nr_body, len_vec_Yevo))
NBsystem_body_n_pos             = np.reshape(int_result_NBsystem_unpacked[:,:,0:3],     (nrtimes_output,nr_body,3))
NBsystem_body_n_vel             = np.reshape(int_result_NBsystem_unpacked[:,:,3:6],     (nrtimes_output,nr_body,3))
NBsystem_body_n_q               = np.reshape(int_result_NBsystem_unpacked[:,:,6:15],    (nrtimes_output,nr_body,3,3))
NBsystem_body_n_qdot            = np.reshape(int_result_NBsystem_unpacked[:,:,15:24],   (nrtimes_output,nr_body,3,3))

#allocate space/define:
t_a123_principalaxes_arr     = np.zeros((nrtimes_output,3))
time_arr_for_plotting        = np.linspace(0,100,nrtimes_output)
const_NBsystem_arr           = np.reshape((NBsystem.M_Make_const_NBsystem), (nr_body,len_vec_const))
E_internal                   = np.zeros((nrtimes_output))

color_arr_body_n = [50,250]

for i in range(0,nr_body):
#-------------------------------------------------------------
    
    #POSITION:
    #--------------------------------
    pos_x   = NBsystem_body_n_pos[:,i,0]
    pos_y   = NBsystem_body_n_pos[:,i,1]
    
    #PLOT 1:    
    plt.subplot(2,3,1)
    plt.xlabel('pos x')
    plt.ylabel('pos y')
    plt.title('Orbits')
    for t in range(0,nrtimes_output-1):
        plt.plot(pos_x[t:t+2],pos_y[t:t+2], color=plt.cm.jet(color_arr_body_n[i]))
    #--------------------------------
    
    #TIDES: 
    #--------------------------------  
    const_i             = const_NBsystem_arr[i,:]
    evoTides_yesno_i    = const_i[5]

    if (evoTides_yesno_i == 1):    
        #analyze tidal deformations/internal energy:
        for t in range(0,nrtimes_output):
            q_t     = np.matrix(NBsystem_body_n_q[t,i,:,:])
            qdot_t  = np.matrix(NBsystem_body_n_qdot[t,i,:,:])
            S_t     = q_t*q_t.T        
            Eigen_vals, Eigen_vecs = np.linalg.eigh(S_t)
            a1 = Eigen_vals[0]**(1./2.)
            a2 = Eigen_vals[1]**(1./2.)
            a3 = Eigen_vals[2]**(1./2.)      
            t_a123_principalaxes_arr[t,:] = [a1,a2,a3]
            get_Internal_Info   = func_internal_properties(const_i, q_t, qdot_t) 
            E_internal[t]       = np.sum(get_Internal_Info)   
            
            print q_t
            print qdot_t
            print t  
              
              
              
                
        #PLOT 2:
        plt.subplot(2,3,2)
        plt.xlabel('Time')
        plt.ylabel('ai/R0')
        plt.title('Tidal Deformation (ai)')
        tt = time_arr_for_plotting
        for t in range(0,nrtimes_output-1):
            plt.plot(tt[t:t+2],t_a123_principalaxes_arr[t:t+2,0], color=plt.cm.jet(color_arr_body_n[i]))
            plt.plot(tt[t:t+2],t_a123_principalaxes_arr[t:t+2,1], color=plt.cm.jet(color_arr_body_n[i]))
            plt.plot(tt[t:t+2],t_a123_principalaxes_arr[t:t+2,2], color=plt.cm.jet(color_arr_body_n[i]))
        #plt.axis([0, 100, 0.7, 1.5])
    
        #PLOT 3:
        plt.subplot(2,3,3)
        plt.xlabel('Time')
        plt.ylabel('E_star(t) - E_star(0)')
        E_internal_initial  = E_internal[0]
        delta_E             = E_internal - E_internal_initial
        for t in range(0,nrtimes_output-1):
            plt.plot(tt[t:t+2],delta_E[t:t+2], color=plt.cm.jet(color_arr_body_n[i]))
    #--------------------------------
    
#-------------------------------------------------------------    


#CM Orbital energy and parameters (a,rp):
E_pot_pair_ij_t         = np.zeros((nr_body,nr_body,nrtimes_output))
E_kin_i_t               = np.zeros((nr_body,nrtimes_output))
E_tot_NBsystem_COM_t    = np.zeros((nrtimes_output))

dist_rij_t              = np.zeros((nr_body,nr_body,nrtimes_output))
rpericenter_ij_t        = np.zeros((nr_body,nr_body,nrtimes_output))
asemimajora_ij_t        = np.zeros((nr_body,nr_body,nrtimes_output))

FtidFbin_ij_t           = np.zeros((nr_body,nr_body,nrtimes_output))

for i in range(0,nr_body):
#-------------------------------------------------------------
    const_i = const_NBsystem_arr[i,:]
    pos_i   = np.reshape(NBsystem_body_n_pos[:,i,:], (nrtimes_output,3))
    vel_i   = np.reshape(NBsystem_body_n_vel[:,i,:], (nrtimes_output,3))
    mass_i  = const_i[0]
    Rad_i   = const_i[1]
    
    E_kin_t = (1./2.)*mass_i*(vel_i[:,0]**2 + vel_i[:,1]**2 + vel_i[:,2]**2) 
    E_kin_i_t[i,:]  = E_kin_t  

    for j in range(0,nr_body):
        if (i != j):                                       
            #-------------------------------------------------    
            const_j = const_NBsystem_arr[j,:]
            pos_j   = np.reshape(NBsystem_body_n_pos[:,j,:], (nrtimes_output,3))
            vel_j   = np.reshape(NBsystem_body_n_vel[:,j,:], (nrtimes_output,3))
            mass_j  = const_j[0]
    
            r_vec_ij    = pos_i - pos_j
            r_ij        = np.sqrt(r_vec_ij[:,0]**2 + r_vec_ij[:,1]**2 + r_vec_ij[:,2]**2)
            #energy:
            E_pot_t = mass_i*mass_j/r_ij
            E_pot_pair_ij_t[i,j,:]  = E_pot_t

            Mtot_ij     = mass_i + mass_j
            Mred_ij     = mass_i*mass_j/Mtot_ij
            v_vec_ij    = vel_i - vel_j
            E_kin_ij_t  = (1./2.)*Mred_ij*(v_vec_ij[:,0]**2 + v_vec_ij[:,1]**2 + v_vec_ij[:,2]**2)
            E_tot_ij_t  = E_kin_ij_t - E_pot_t
            a_ij_t      = mass_i*mass_j/(2.*np.abs(E_tot_ij_t))
            
            L_ij_t  = np.zeros((nrtimes_output))
            for t in range(0,nrtimes_output):
                L_vec       = Mred_ij*np.cross(r_vec_ij[t,:],v_vec_ij[t,:])
                L_ij_t[t]   = np.sqrt(L_vec.dot(L_vec))
            e2m1    = 2.*E_tot_ij_t*(L_ij_t**2)/(Mred_ij*((Mred_ij*Mtot_ij)**2))
            ecc     = np.sqrt(1.+e2m1)
            rp_ij_t = a_ij_t*(1.-ecc)
            
            #orbit:
            dist_rij_t[i,j,:]           = r_ij
            asemimajora_ij_t[i,j,:]     = a_ij_t
            rpericenter_ij_t[i,j,:]     = rp_ij_t
            
            #tidal forces (Ftid) relative to binding force (Fbin):
            FtidFbin_ij             = (mass_j/mass_i)*((Rad_i/r_ij)**3.)
            FtidFbin_ij_t[i,j,:]    = FtidFbin_ij 
            #-------------------------------------------------
#-------------------------------------------------------------

#-------------------------------------------------------------
for t in range(0,nrtimes_output):
#-------------------------------------------------------------    
    E_kin_tot_t = np.sum(E_kin_i_t[:,t])
    E_pot_tot_t = np.sum(E_pot_pair_ij_t[:,:,t])/2.
    E_tot_NBsystem_COM_t[t]  = E_kin_tot_t - E_pot_tot_t 
#-------------------------------------------------------------

#PLOT 4:
plt.subplot(2,3,4)
plt.xlabel('Time')
plt.ylabel('CM E_tot NBsystem')
tt = time_arr_for_plotting
for t in range(0,nrtimes_output-1):
    plt.plot(tt[t:t+2],E_tot_NBsystem_COM_t[t:t+2], color=plt.cm.jet(0))

#PLOT 5:
plt.subplot(2,3,5)
plt.xlabel('Time')
plt.ylabel('Distance r (ij)')
plt.axis([0, 100, 1e-5, 20])       
tt = time_arr_for_plotting
for t in range(0,nrtimes_output-1):
    plt.plot(tt[t:t+2],(dist_rij_t[0,1,t:t+2]),      color=plt.cm.jet(0))
    plt.plot(tt[t:t+2],(FtidFbin_ij_t[1,0,t:t+2]),   color=plt.cm.jet(50))
plt.yscale('log')

print FtidFbin_ij_t[1,0,:]

#PLOT 6:
plt.subplot(2,3,6)
plt.xlabel('Time')
plt.ylabel('orbital SMA and rp (ij)')
plt.axis([0, 100, 0, 20])       
tt = time_arr_for_plotting
for t in range(0,nrtimes_output-1):
    plt.plot(tt[t:t+2],asemimajora_ij_t[0,1,t:t+2], color=plt.cm.jet(0))
    plt.plot(tt[t:t+2],rpericenter_ij_t[0,1,t:t+2], color=plt.cm.jet(0))



plt.subplots_adjust(hspace=0.4)
plt.savefig('Encounter.pdf')    
plt.show()
exit()
#------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------


















































