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
import subprocess
import time
#############################################################################
# GLOBALS
#############################################################################

N_ROUNDS = 1
N_ENGINES = 4

#############################################################################
# Script
#############################################################################

subprocess.call("./clean.sh", shell=True)

# start the server independently
server = subprocess.Popen(['accelerator', 'serve'])

# we're going to run N_ROUNDS of adaptive sampling
for i in range(N_ROUNDS):
    # start your engines!
    pids = set()  # keep track of the pids of all of the engines we start
    for j in range(N_ENGINES):
        int_file = "integrator%d.xml" % j
        # each of these engines will run a single trajectory
        proc = subprocess.Popen(['accelerator', 'simulate',
            '--device_index=' + str(j), '--integrator_xml=%s' % int_file])
        pids.add(proc.pid)
        time.sleep(1)

    # wait on the engines to finish
    while pids:
        pid, retval = os.wait()
        print('{p} finished'.format(p=pid))
        pids.remove(pid)

    # build MSM #TODO: implement rmsd in 2d
    # proc = subprocess.Popen(['accelerator', 'model'])
    # proc.wait()

# we're all done, so lets shut down the server
server.kill()
