#!/usr/bin/env python
"""
Example code that runs the whole workflow. Instead of submitting
jobs to a queue, this just submits them as subprocesses and then waits.

It should be pretty simple just to look at
"""
#############################################################################
# Imports
#############################################################################

import os
import subprocess
import time
#############################################################################
# GLOBALS
#############################################################################

def _print_round_info(round, n_engines):
    print("\n\n")
    print("======================================")
    print("Round %d" % round)
    print("Using %d engines" % n_engines)
    print("======================================")
    print("\n\n")     
    
def run_round(n_engines):

                
    # start your engines!
    pids = set()  # keep track of the pids of all of the engines we start
    for j in range(n_engines):
        # each of these engines will run a single trajectory
        proc = subprocess.Popen(['accelerator', 'OpenMM',
            '--device_index=' + str(j)])
        pids.add(proc.pid)
        time.sleep(1)

    # wait on the engines to finish
    while pids:
        pid, retval = os.wait()
        print('{p} finished'.format(p=pid))
        pids.remove(pid)

    # build MSM
    proc = subprocess.Popen(['accelerator', 'model'])
    proc.wait()      

def exponential_scheduler(out_dir, n_rounds):
    for i in range(n_rounds):
        n_engines = 2 ** i
        _print_round_info(i, n_engines)
        run_round(n_engines)

def constant_scheduler(out_dir, n_rounds, n_engines):
    
    for i in range(n_rounds):
        _print_round_info(i, n_engines)
        run_round(n_engines)
    

    
def run_accel(out_dir, run_type, n_rounds, n_engines):
    # Change to working directory
    os.chdir(os.path.abspath(out_dir))    
    
    # start the server independently
    server = subprocess.Popen(['accelerator', 'serve'])
    
    if run_type == 'constant':
        constant_scheduler(n_rounds, n_engines)
    elif run_type == 'exponential':
        exponential_scheduler(n_rounds)
    else:
        print("Please specify a valid type of scheduler")
    
    
    # we're all done, so lets shut down the server
    server.kill()
