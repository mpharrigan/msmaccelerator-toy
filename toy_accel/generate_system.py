#!/usr/bin/env python

from simtk import openmm, unit
from sys import argv
import mdtraj
import mullerforce as mf
import numpy as np
import random
import sys

# Global parameters
nParticles = 100
mass = 12.0 * unit.dalton
temperature = 750 * unit.kelvin
friction = 100 / unit.picosecond
timestep = 10.0 * unit.femtosecond

def generate_openmm(sys_fn, int_fn):

    # Prepare the system
    system = openmm.System()
    mullerforce = mf.MullerForce()
    system.addParticle(mass)
    mullerforce.addParticle(0, [])
    system.addForce(mullerforce)

    print("Writing System file: %s" % sys_fn)
    with open(sys_fn, 'w') as f:
        f.write(openmm.XmlSerializer.serialize(system))

    # Prepare integrator
    print("Writing Integrator file: %s" % int_fn)
    integrator = openmm.LangevinIntegrator(temperature, friction, timestep)
    with open(int_fn, 'w') as f:
        f.write(openmm.XmlSerializer.serialize(integrator))

def generate_starting_structures_fixed(top_fn, x, y, out_fn):
    pdb = mdtraj.load(top_fn)
    xyz = np.array([
            [[x, y, 0.0]]
            ])
    traj = mdtraj.Trajectory(xyz, pdb.topology)
    print("Writing fixed seed structure: %s" % out_fn)
    traj.save(out_fn)

def _randx():
    return random.random() * (1.2 + 1.5) - 1.5
def _randy():
    return random.random() * (2.0 + 0.2) - 0.2

def generate_starting_structures_random(top_fn, num, out_fn):
    pdb = mdtraj.load(top_fn)
    xyz = list()
    for i in xrange(num):
        xyz.append([[_randx(), _randy(), 0.0]])
    traj = mdtraj.Trajectory(np.array(xyz), pdb.topology)
    print("Writing random seed structures: %s" % out_fn)
    traj.save(out_fn)

def generate_config_file(rep_int, n_steps, rel_top_fn, seed_rel_fn, beta, config_fn):
    
    # Prepare config
    configs = list()
    configs.append("c = get_config()")
    configs.append("c.Modeler.use_custom_metric = True")
    configs.append("c.Modeler.custom_metric_path = 'metric.pickl'")
    configs.append("c.Modeler.clusterer = 'hybrid'")
    configs.append("c.OpenMMSimulator.report_interval = %d" % rep_int)
    configs.append("c.OpenMMSimulator.number_of_steps = %d" % n_steps)
    configs.append("c.OpenMMSimulator.minimize = False")
    configs.append("c.AdaptiveServer.topology_pdb = '%s'" % rel_top_fn)
    configs.append("c.BaseSampler.seed_structures = '%s'" % seed_rel_fn)
    configs.append("c.CountsSampler.beta = %d" % beta)
    
    # Write config
    with open(config_fn, 'w') as config_file:
        for line in configs:
            config_file.write(line + '\n')
        


# if __name__ == "__main__":
#     error = "Specify what structures to generate: random, fixed, none"
#     if len(argv) < 2:
#         print(error)
#         sys.exit(1)
#     if argv[1] == 'random':
#         generate_starting_structures_random()
#     elif argv[1] == 'fixed':
#         generate_starting_structures_fixed()
#     elif argv[1] == 'none':
#         pass
#     else:
#         print(error)
#         sys.exit(1)
#     generate_openmm() #TODO: Update these function calls or remove this code


