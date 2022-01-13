#!/usr/bin/env python
"""
Following MPI template is from:
https://github.com/jbornschein/mpi4py-examples/blob/master/09-task-pull.py
Demonstrate the task-pull paradigm for high-throughput computing
using mpi4py. Task pull is an efficient way to perform a large number of
independent tasks when there are more tasks than processors, especially
when the run times vary for each task. 
"""

from mpi4py import MPI
import testname
import numpy as np
import sys
from subprocess import call
import subprocess
import time
import itertools
from sympy import *
from sympy.solvers import solve
from sympy import Symbol
from sympy.abc import x



#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def func_make_MC_IC_3body_binsin():
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
    
    #------------------------------------------------------------------------------------------
    #read data and define:
    #------------------------------------------------------------------------------------------
    #data name and folder:    
    data_folder     = '/Users/jsamsing/Desktop/TIDES_PROJ/MOCCA_ICs/MC_1/'
    #data_folder     = '/scratch/network/jsamsing/Stellar_Tides/' 
    data_name       = 'MT1000_case1_'   #CHOOSE ALSO THIS FOR OUTPUT FILE NAME!

    #read data:
    tf = open(data_folder+data_name+'nbody_params_INT.txt', "r")
    nbody_params_INT_all        = np.loadtxt(tf, dtype=int)
    tf.close()
    tf = open(data_folder+data_name+'nbody_params_REAL.txt', "r")
    nbody_params_REAL_all       = np.loadtxt(tf, dtype=float)
    tf.close()
    tf = open(data_folder+data_name+'obj_info.txt', "r")
    obj_info_all                = np.loadtxt(tf, dtype=float)
    tf.close()
    tf = open(data_folder+data_name+'pos_vel.txt', "r")
    pos_vel_all                 = np.loadtxt(tf, dtype=float)
    tf.close()
   
    #define:
    MCinfo_ICNbody_list = []
    nr_tot_scatterings  = len(nbody_params_INT_all)
    #------------------------------------------------------------------------------------------
    
    
    #------------------------------------------------------------------------------------------
    #make MC IC info list:
    #------------------------------------------------------------------------------------------
    #allocate:
    MC_settings_list_INT    = np.zeros(5,dtype='i')
    MC_settings_list_REAL   = np.zeros(5,dtype='d')
    #insert info:
    MC_settings_list_INT[:]     = [0, 0, 0, nr_tot_scatterings, 0]
    MC_settings_list_REAL[:]    = [0.0, 0.0, 0.0, 0.0, 0.0]
    #put into one list:
    MC_settings_list = [MC_settings_list_INT, MC_settings_list_REAL]
    #------------------------------------------------------------------------------------------

    
    #------------------------------------------------------------------------------------------
    #Make MC ICs:
    #------------------------------------------------------------------------------------------
    #initilize IC ID counter:
    icidc = int(0)
    for tc in range(0,nr_tot_scatterings):                                        
            
        #---------------------------------------
        #get IC info from input file:
        #---------------------------------------    
        #Nbody params:
        nbody_params_arr_1_INT      = nbody_params_INT_all[tc,:]
        nbody_params_arr_2_REAL     = nbody_params_REAL_all[tc,:]
        #obj 1:
        b1_const_arr    = obj_info_all[tc, 0:10]
        b1_posxyz_CM    = pos_vel_all[tc,  0:3]
        b1_velxyz_CM    = pos_vel_all[tc,  3:6]
        b1_q            = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
        b1_qdot         = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
        #obj 2:
        b2_const_arr    = obj_info_all[tc, 10:20]
        b2_posxyz_CM    = pos_vel_all[tc,  6:9]
        b2_velxyz_CM    = pos_vel_all[tc,  9:12]
        b2_q            = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
        b2_qdot         = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
        #obj 3:
        b3_const_arr    = obj_info_all[tc, 20:30]
        b3_posxyz_CM    = pos_vel_all[tc,  12:15]
        b3_velxyz_CM    = pos_vel_all[tc,  15:18]
        b3_q            = np.array([[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]], dtype='d')
        b3_qdot         = np.array([[0.,0.,0.],[0.,0.,0.],[0.,0.,0.]], dtype='d')
        #---------------------------------------
        
        #---------------------------------------
        #IC list:
        #---------------------------------------
        IN_IC_code_version                          = 2       #MUST = 2 for parallel version !!!!
        IN_dimlen_IC_nbody                          = 3*9     #(for n_particles=3, = 27)
        #define arrays with a specific format:
        IN_IC_simparams_INT                         = np.zeros(10,dtype='i')
        IN_IC_simparams_REAL                        = np.zeros(10,dtype='d')
        IN_IC_nbody_const_posvel_qqdot_etc_arr      = np.zeros((IN_dimlen_IC_nbody,10),dtype='d')
        #sim parameters:
        IN_IC_simparams_INT[:]                      = nbody_params_arr_1_INT[:]
        IN_IC_simparams_REAL[:]                     = nbody_params_arr_2_REAL[:]
        #obj 1:
        IN_IC_nbody_const_posvel_qqdot_etc_arr[0,:]     = b1_const_arr[:]   #[1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[1,0:3]   = b1_posxyz_CM[:]   #[-0.538861, 93.044669, 64.479373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[2,0:3]   = b1_velxyz_CM[:]   #[0.004368, 0.117192, -0.027313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[3,0:3]   = b1_q[0,:]         #[1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[4,0:3]   = b1_q[1,:]         #[0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[5,0:3]   = b1_q[2,:]         #[0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[6,0:3]   = b1_qdot[0,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[7,0:3]   = b1_qdot[1,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[8,0:3]   = b1_qdot[2,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #obj 2:
        IN_IC_nbody_const_posvel_qqdot_etc_arr[9,:]     = b2_const_arr[:]   #[1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[10,0:3]  = b2_posxyz_CM[:]   #[-22.038986, 93.044669, 64.479373, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[11,0:3]  = b2_velxyz_CM[:]   #[0.004368, -0.187805, -0.027313, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[12,0:3]  = b2_q[0,:]         #[1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[13,0:3]  = b2_q[1,:]         #[0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[14,0:3]  = b2_q[2,:]         #[0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[15,0:3]  = b2_qdot[0,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[16,0:3]  = b2_qdot[1,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[17,0:3]  = b2_qdot[2,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #obj 3:
        IN_IC_nbody_const_posvel_qqdot_etc_arr[18,:]    = b3_const_arr[:]   #[1.0,1.0, 1.500000, 1.666667, 0.102300, 1.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[19,0:3]  = b3_posxyz_CM[:]   #[22.577848, -186.089337, -128.958746, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[20,0:3]  = b3_velxyz_CM[:]   #[-0.008735, 0.070613, 0.054626, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[21,0:3]  = b3_q[0,:]         #[1.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[22,0:3]  = b3_q[1,:]         #[0.0,   1.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[23,0:3]  = b3_q[2,:]         #[0.0,   0.0,   1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[24,0:3]  = b3_qdot[0,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[25,0:3]  = b3_qdot[1,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        IN_IC_nbody_const_posvel_qqdot_etc_arr[26,0:3]  = b3_qdot[2,:]      #[0.0,   0.0,   0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        #Put into one list:
        outlist_IC_Nbody            = [IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL, IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr] 
        #---------------------------------------
        #MC param info list:
        #---------------------------------------
        #allocate:
        MC_info_INT     = np.zeros(10,dtype='i')
        MC_info_REAL    = np.zeros(10,dtype='d')
        #input info:
        MC_info_INT[:]  = [icidc, 0, 0, 0, 0, 0, 0, 0, 0, 0]
        MC_info_REAL[:] = [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
        outlist_MC_info_INT_REAL    = [MC_info_INT,MC_info_REAL]
        #---------------------------------------
        #append MC,IC list to final output list:
        #---------------------------------------
        MCinfo_ICNbody = [outlist_MC_info_INT_REAL, outlist_IC_Nbody]
        MCinfo_ICNbody_list.append(MCinfo_ICNbody)
        #---------------------------------------
        #update ic id counter
        #---------------------------------------
        icidc = icidc + 1
        #---------------------------------------
             
    #------------------------------------------------------------------------------------------
    #------------------------------------------------------------------------------------------
    #Return final info from function:
    #------------------------------------------------------------------------------------------
    final_MC_output_list = [MC_settings_list, MCinfo_ICNbody_list]
    return final_MC_output_list
    #------------------------------------------------------------------------------------------
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def func_run_Nbody(ICnbody):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    #input format: IC_inputNbody = ICnbody = [IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL, IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr]
    #unpack input list:
    IN_IC_code_version                      = ICnbody[0]
    IN_IC_simparams_INT                     = ICnbody[1]
    IN_IC_simparams_REAL                    = ICnbody[2]
    IN_dimlen_IC_nbody                      = ICnbody[3]
    IN_IC_nbody_const_posvel_qqdot_etc_arr  = ICnbody[4]
    #call Nbody code with these input IC:
    OUT_Nbody = testname.interfacesub(IN_IC_code_version, IN_IC_simparams_INT, IN_IC_simparams_REAL, IN_dimlen_IC_nbody, IN_IC_nbody_const_posvel_qqdot_etc_arr)
    return OUT_Nbody    #OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
def enum(*sequential, **named):
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    """
    This function is also from the MPI template example.
    Handy way to fake an enumerated type in Python
    http://stackoverflow.com/questions/36932/how-can-i-represent-an-enum-in-python
    """
    enums = dict(zip(sequential, range(len(sequential))), **named)
    return type('Enum', (), enums)
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------







#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
#Main program: make IC list, Run Nbody using MPI and write output to .txt files
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------


"""
Basics: This full code will be send to and executed by each processor.
At the time the processor receives the code, the processor gets a unique 
'rank' number. The list of these number goes from 0 to N, where N is the total
number of processors. These numbers (or ids) will not change during the computation.
We use this rank number to pass info between the processors.
Normally you set processor 0 (rank '0') to be the 'Master' which sends and receives
info from the 'Worker' processers (rank >0). What part of the code a specific processor
runs is controlled by a simple 'if rank ==...' statement.
"""


#------------------------------------------------------------------------------------------
#Define and initialize:
#------------------------------------------------------------------------------------------
# Define MPI message tags
tags = enum('READY', 'DONE', 'EXIT', 'START', 'WAKEUP')

# Initializations and preliminaries
comm = MPI.COMM_WORLD   # get MPI communicator object
size = comm.size        # total number of processes
rank = comm.rank        # rank of this process
status = MPI.Status()   # get MPI status object
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#WORKER process executes code below
#------------------------------------------------------------------------------------------
if rank != 0:
#------------------------------------------------------------------------------------------       
   
    #---------------------------------------
    #If Worker is active:
    #---------------------------------------
    while True:
        #-------------------------
        #send info: READY to Master (dest=0):
        #-------------------------
        time.sleep(1.0*float(np.random.rand(1)))    #SET TIME BUFFER
        comm.send(None, dest=0, tag=tags.READY)
        #-------------------------
        #Receive info from Master on what to do:
        #-------------------------
        data    = comm.recv(source=0, tag=MPI.ANY_TAG, status=status)
        tag     = status.Get_tag()
        #-------------------------
        #If Master says 'Start':
        #-------------------------
        if tag == tags.START:
            # split data file into info and ICs:
            MC_info_INT_REAL    = data[0] 
            IC_Nbody            = data[1]
            #Run Nbody code with 'IC_Nbody' as input:
            output_Nbody    = func_run_Nbody(IC_Nbody)
            result          = [MC_info_INT_REAL, output_Nbody]
            comm.send(result, dest=0, tag=tags.DONE)
        #-------------------------    
        #If Master says 'Exit' break loop:
        #-------------------------
        if tag == tags.EXIT:
            break
    #---------------------------------------
    
    #---------------------------------------
    #If Worker has been stopped (Exit):
    #---------------------------------------        
    #Send info to Master that Worker has now stopped (Exit): 
    comm.send(None, dest=0, tag=tags.EXIT)  
    #---------------------------------------
    
#------------------------------------------------------------------------------------------
#if rank != 0:    
#------------------------------------------------------------------------------------------


#------------------------------------------------------------------------------------------
#MASTER process executes code below
#------------------------------------------------------------------------------------------
if rank == 0:
#------------------------------------------------------------------------------------------ 
    
    #---------------------------------------
    #wake up all processers: (necessary for e.g. SLURM)
    #---------------------------------------
    for nps in range(1,size):   #size = total number of processes
        comm.send(None, dest=nps, tag=tags.WAKEUP)
        print 'wake up processor:    ', nps
    #---------------------------------------
    
    #---------------------------------------
    #Input to Master process:
    #---------------------------------------
    #-------------------------
    #save data:
    #-------------------------
    data_name   = 'OUT_MT1000_case1_'
    data_folder = '/Users/jsamsing/Desktop/TIDES_PROJ/MC_OUTPUT/'
    #data_folder = '/scratch/network/jsamsing/Stellar_Tides/'
    #-------------------------
    #max sim time in secs:
    #-------------------------
    #when using a scheduler (SLURM): PUT SAME TIME AS GIVEN TO SLURM!!
    sim_time_maxupper   = 1.0*3600.                             #in secs
    sim_time_buffer     = 0.02*sim_time_maxupper + 2.*60        #in secs
    max_sim_secs        = sim_time_maxupper - sim_time_buffer   
    '''
    NOTES on time: We give each task a physical time limit (N*initial orbital times)
    and a wall clock time limit which is the time we have left of the total sim time.
    If a task does not finish before the total sim time is over, the task will
    end with endstate id 12 (wall clock limit). IMPORTANT: Its important that all tasks have been
    distributed and have start running before the sim time is over. The code still works
    if this is not the case, but its really hard to use the dataset. The main idea
    with the wall clock limit is simply to stop the (often small) fraction of interactions
    that just keeps going forever.    
    '''
    #-------------------------
    #---------------------------------------
    
    #---------------------------------------
    #Generate bin-sin ICs: (all settings are set in the function)
    #---------------------------------------
    out = func_make_MC_IC_3body_binsin()
    #split into different lists:
    MC_settings_list_INT    = out[0][0] # MC_settings_list_INT[:]     = [nr_SMA, nr_vinf, nr_scatterings_per_paramcomb, nr_tot_scatterings, 0]
    MC_settings_list_REAL   = out[0][1] # MC_settings_list_REAL[:]    = [0.0, 0.0, 0.0, 0.0, 0.0]
    MC_output_list          = out[1]    # [outlist_MC_info_INT_REAL, outlist_IC_Nbody], where: outlist_MC_info_INT_REAL = [[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0],[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]] 
    #define:
    nr_tot_scatterings      = MC_settings_list_INT[3]
    all_data                = [] 
    #---------------------------------------    

    #---------------------------------------
    #initialize:
    #---------------------------------------
    taskc               = 0         #task counter - counting the number of tasks (or IC scatterings) that have been send to Workers.
    nr_tasks_completed  = 0         #info counter.
    num_workers         = size - 1
    closed_workers      = 0
    time_start          = time.time()
    #---------------------------------------
    
    #---------------------------------------
    #SIMULATE - if not all workers are closed:
    #---------------------------------------
    while closed_workers < num_workers:        
        #-------------------------
        #receive info from Worker:
        #-------------------------
        data    = comm.recv(source=MPI.ANY_SOURCE, tag=MPI.ANY_TAG, status=status)
        source  = status.Get_source()
        tag     = status.Get_tag()
        #-------------------------
        #if Worker is Ready:
        #-------------------------
        if tag == tags.READY:
            #-------------------------
            #if we have more tasks, send Worker a task:
            #-------------------------
            if taskc < nr_tot_scatterings:
                #calc time left of sim:
                time_now            = time.time() - time_start
                time_left           = max_sim_secs - time_now
                #set max wall clock time for Nbody code:
                data_for_worker             = MC_output_list[taskc]
                data_for_worker[1][2][4]    = time_left + (sim_time_buffer/2.)*float(np.random.rand(1))   #1000.*float(np.random.rand(1)) 
                #send data to Worker:
                data_send = data_for_worker
                comm.send(data_send, dest=source, tag=tags.START)
                #update taskc by +1 because we now succesfully send a task:
                taskc += 1
                #TEST: print info:
                #print 'sent:    ', taskc, nr_tot_scatterings
                #print time_now, time_left, max_sim_secs
            #-------------------------
            #if we do not have more tasks, tell Worker to stop:
            #-------------------------
            else:
                comm.send(None, dest=source, tag=tags.EXIT)
            #-------------------------    
        #-------------------------
        #if Worker is Done:
        #-------------------------
        if tag == tags.DONE:
            result = data
            #collect data:
            all_data.append(result)
            #TEST: print info:
            nr_tasks_completed +=1  #nr completed tasks
            print 'comp:    ', nr_tasks_completed, nr_tot_scatterings
            #print result[1][0]
            #print (time.time() - time_start)
        #-------------------------    
        #if Worker is closed down:
        #-------------------------
        if tag == tags.EXIT:
            closed_workers += 1
	    print 'closed_workers:  ', closed_workers 
    #--------------------------------------
    
    #---------------------------------------
    #organize output results:
    #---------------------------------------
    print 'SAVING DATA 1:', (time.time() - time_start)
    MC_output_Nbody     = [None]*nr_tot_scatterings 
    for i in range(0,nr_tot_scatterings):
        data_i              = all_data[i]                       # [outlist_MC_info_INT_REAL, output_Nbody]
        #where: outlist_MC_info_INT_REAL = [[icidc,ac,vc,sc,0],[SMA_bin,vinf_sin,b_max,0.0,0.0]] and output_Nbody: [OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL]
        data_i_ori_id       = data_i[0][0][0]                   # =icidc which is the ori IC ID. This runs from 0->nr_tot_scatterings in the (correct) order we made the ICs in.
        data_i_output_Nbody = data_i[1]                         # output_Nbody for this data_i: [OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL]
        MC_output_Nbody[data_i_ori_id] = data_i_output_Nbody    # Here we use 'icidc' (data_i_ori_id) to re-shuffle our data to the initial order.
    #Now the list: MC_output_Nbody has the ori and same order as the IC was created in.
    #---------------------------------------
    #Unpack Nbody data for saving:
    #---------------------------------------
    print 'SAVING DATA 2:', (time.time() - time_start)    
    output_Nbody_endstate_INT       = [] #1x10 per one MC scattering
    output_Nbody_endstate_REAL      = [] #1x10 per one MC scattering
    output_Nbody_xtra_info_INT      = [] #1x10 per one MC scattering
    output_Nbody_xtra_info_REAL     = [] #1x10 per one MC scattering
    output_Nbody_xtra_2_info_REAL   = [] #1x10 per one MC scattering
    for i in range(0,nr_tot_scatterings):
        #MC_output_Nbody[i][:] = [OUT_out_endstate_INT, OUT_out_endstate_REAL, OUT_out_xtra_info_INT, OUT_out_xtra_info_REAL, OUT_out_xtra_2_info_REAL]
        output_Nbody_endstate_INT.append(MC_output_Nbody[i][0])
        output_Nbody_endstate_REAL.append(MC_output_Nbody[i][1])
        output_Nbody_xtra_info_INT.append(MC_output_Nbody[i][2])
        output_Nbody_xtra_info_REAL.append(MC_output_Nbody[i][3])
        output_Nbody_xtra_2_info_REAL.append(MC_output_Nbody[i][4])
    #Unpack info from 'MC_output_list':
    #MC_output_list          = out[1]    # [outlist_MC_info_INT_REAL, outlist_IC_Nbody], where: outlist_MC_info_INT_REAL = [[icidc, ac, vc, sc, 0, 0, 0, 0, 0, 0],[SMA_bin, vinf_sin, b_max, b_sampsurf, E_tot_binsinsystem, L_tot_binsinsystem, 0.0, 0.0, 0.0, 0.0]] 
    print 'SAVING DATA 3:', (time.time() - time_start)
    outlist_MC_info_INT     = []
    outlist_MC_info_REAL    = []
    for i in range(0,nr_tot_scatterings):
        outlist_MC_info_INT.append(MC_output_list[i][0][0])
        outlist_MC_info_REAL.append(MC_output_list[i][0][1])
    #---------------------------------------          
    #save to .txt files:
    #---------------------------------------
    print 'SAVING DATA 4:', (time.time() - time_start)
    #write data to files in folder:
    tf = open(data_folder+data_name+'MC_settings_list_INT.txt', "w")
    np.savetxt(tf, MC_settings_list_INT[None],   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'MC_settings_list_REAL.txt', "w")
    np.savetxt(tf, MC_settings_list_REAL[None],   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'outlist_MC_info_INT.txt', "w")
    np.savetxt(tf, outlist_MC_info_INT,   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'outlist_MC_info_REAL.txt', "w")
    np.savetxt(tf, outlist_MC_info_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_endstate_INT.txt', "w")
    np.savetxt(tf, output_Nbody_endstate_INT,   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_endstate_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_endstate_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_info_INT.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_info_INT,   fmt='%5i')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_info_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_info_REAL,   fmt='%5f')
    tf.close()
    tf = open(data_folder+data_name+'output_Nbody_xtra_2_info_REAL.txt', "w")
    np.savetxt(tf, output_Nbody_xtra_2_info_REAL,   fmt='%5f')
    tf.close()
    #---------------------------------------
    
#------------------------------------------------------------------------------------------    
#if rank == 0:   
#------------------------------------------------------------------------------------------   


#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------   
#------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
