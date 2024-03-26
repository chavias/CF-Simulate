#!/usr/bin/env python
from __future__ import print_function
from subprocess import call
from subprocess import check_call
import os
import time
import numpy as np

def calc_sim_param(time_mod, n_max, nu_r, seq = 'sequence'):
    tr = 1./int(nu_r)
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
    seq = seq
    step = str(int(step))
    n_max = str(int(n_max))
    n_sample = str(int(n_sampl))
    nu_r = str(int(nu_r))
    return seq,step,n_max,n_sample,nu_r

def call_gamma_sim(seq,step,n_max,n_sample,nu_r,item):
    check_call(['./gamma_time.csh',seq,step,n_max,n_sample,nu_r,str(int(item))])

if __name__=='__main__':
    dir_path = os.path.dirname(os.path.realpath(__file__));
    directory = dir_path.split('/')
    try: os.mkdir('/cluster/scratch/chavezm/' + directory[-1])
    except: pass
    # read parameter from infos for first iteration
    time_mod = 2e-3 # sec determind by shapefile
    n_1rep =  1000 # points determinde by shapefile
    # nu_r = np.linspace(20e3,150e3,1001) # Hz
    nu_r = 100e3
    sweep = np.linspace(0,300e3,2001)
    reps = [1] # 1 means pulse scheme is applied once
    for rep in reps:
        print(f'============================== rep = {rep} ==============================')
        for item in sweep:
            seq,step,n_max,n_sample,nu_r = calc_sim_param(\
                                                      rep*time_mod,rep*n_1rep,nu_r,seq='CW-1000')
            call_gamma_sim(seq,step,n_max,n_sample,nu_r,item)
            print('---------------------------------------------------')
