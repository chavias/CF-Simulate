#!/usr/bin/env python
from __future__ import print_function
from subprocess import call
from subprocess import check_call
import os
import time
import numpy as np

def calc_sim_param(time_mod, n_max, nu_r, seq = 'sequence'):
    tr = 1./nu_r
    dt = time_mod/n_max
    steps_list = range(1000,8001);
    n_sample_list = range(2,201);
    epsilon = [[0 for x in range(0,len(n_sample_list))] \
               for y in range(0,len(steps_list))] ;
    steps_per_cycle = float(tr)/dt;
    value = 1
    for index_steps, steps in enumerate(steps_list):
        for index_n_sample, n_sample in enumerate(n_sample_list):
            n_sample_per_steps = float(steps)/n_sample
            # print('n_sample_per_steps: ', abs(n_sample_per_steps))
            if abs(n_sample_per_steps-steps_per_cycle) < value:
                sol1,sol2 = index_steps, index_n_sample
                value = abs(n_sample_per_steps-steps_per_cycle)
            else: pass
    print('[+] simulation parameter calculated')
    print('[*] epsilon = ' +str(min(min(epsilon))))
    print('[*] value = ', value)
    step = steps_list[sol1]
    n_sampl = n_sample_list[sol2]
    # run simulations in gamma
    arg1 = seq
    arg2 = str(int(step))
    arg3 = str(int(n_max))
    arg4 = str(int(n_sampl))
    arg5 = str(int(nu_r))
    return arg1,arg2,arg3,arg4,arg5

def call_gamma_sim(arg1,arg2,arg3,arg4,arg5,CS):
    check_call(['./gamma_time.csh',arg1,arg2,arg3,arg4,arg5,str(int(CS))])

if __name__=='__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__));
    directory = dir_path.split('/')
    try: os.mkdir('/cluster/scratch/chavezm/' + directory[-1])
    except: pass
    print('---------------------------------------------------')
    print('------------- Numerical Simulation ----------------')
    print('---------------------- of -------------------------')
    print('------------ Increasing repetitions ---------------')
    print('---------------------------------------------------')
    # read parameter from infos for first iteration
    time_mod = 0.004 # sec determind by shapefile
    n_1rep =  20160 # points determinde by shapefile
    nu_r = 10e3
    chemical_shifts = np.linspace(-70e3,70e3,1001) # Hz
    reps = 1 # 1 means pulse scheme is applied once 
    for index, CS in enumerate(chemical_shifts):
        arg1,arg2,arg3,arg4,arg5 = calc_sim_param(\
                reps*time_mod,reps*n_1rep,nu_r,seq='C7_r20')
        call_gamma_sim(arg1,arg2,arg3,arg4,arg5,CS)
        print('---------------------------------------------------')
