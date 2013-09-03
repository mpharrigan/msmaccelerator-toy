
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
import numpy as np
#############################################################################
# GLOBALS
#############################################################################

def logistic_beta(i_round, switch_point, tension):
    """A logistic curve for scheduling beta."""
    return 1.0 / (1.0 + np.exp(-(i_round - switch_point) / tension))

def _print_round_info(i_round, n_engines, beta):
    """Print debugging info for this round."""
    print("\n\n")
    print("======================================")
    print("Round %d" % i_round)
    print("Using %d engines" % n_engines)
    print("Beta = %f" % beta)
    print("======================================")
    print("\n\n")

def run_round(n_engines, beta):
    """Run a single round of simulation.

    This function is called by the various schedulers, which can vary the
    parameters of this function over time (i.e. over round number).
    """
    proc = subprocess.Popen(['accelerator', 'interact', '--set_beta=%f' % beta])
    proc.wait()

    # start your engines!
    pids = set()  # keep track of the pids of all of the engines we start
    for j in range(n_engines):
        # each of these engines will run a single trajectory
        proc = subprocess.Popen(['accelerator', 'OpenMM',
            '--device_index=%d' % j ])
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

def exponential_scheduler(n_rounds, beta):
    """Use 2^n engines, where n is the round index."""
    for i in range(n_rounds):
        n_engines = 2 ** i
        _print_round_info(i, n_engines, beta)
        run_round(n_engines, beta)

def constant_scheduler(n_rounds, n_engines, beta):
    """Use a constant number of engines and a constant beta."""
    for i in range(n_rounds):
        _print_round_info(i, n_engines, beta)
        run_round(n_engines, beta)

def vary_beta_scheduler(n_rounds, n_engines, beta_center=None, beta_tension=None):
    """Use a constant number of engines but vary beta.

    Beta is varied according to a logistic curve parameterized by the 'center',
    or change in concavity and the 'tension' which describes how abruptly one
    switches from low to high.
    """
    if beta_center is None:
        beta_center = n_rounds / 2
    if beta_tension is None:
        beta_tension = n_rounds / 15

    print("======================================")
    print("Beta setup.")
    print("Using center = %f" % beta_center)
    print("Using tension = %f" % beta_tension)
    print("======================================")

    for i in range(n_rounds):
        beta = logistic_beta(i, beta_center, beta_tension)
        _print_round_info(i, n_engines, beta)
        run_round(n_engines, beta)



def run_accel(args):
    """Run an accelerated simulation.

    args contains different arguments depending on args.runtype.

    This function sets up the msmaccelerator server and then calls
    the appropriate function to schedule the rounds of sampling.
    """
    run_type = args.runtype

    # Parse arguments
    out_dir = args.out_dir
    n_rounds = args.n_rounds

    # Change to working directory
    os.chdir(os.path.abspath(out_dir))

    # start the server independently
    server = subprocess.Popen(['accelerator', 'serve'])

    if run_type == 'constant':
        n_engines = args.n_engines
        beta = args.beta
        constant_scheduler(n_rounds, n_engines, beta)
    elif run_type == 'exponential':
        beta = args.beta
        exponential_scheduler(n_rounds, beta)
    elif run_type == 'varybeta':
        n_engines = args.n_engines
        vary_beta_scheduler(n_rounds, n_engines)
    else:
        print("Please specify a valid type of scheduler")


    # we're all done, so lets shut down the server
    server.kill()

