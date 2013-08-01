#!/usr/bin/env python

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from sys import stdout, argv
import mdtraj
import mullerforce as mf
import numpy as np
import random

# Global parameters
nParticles = 100
mass = 12.0 * dalton
temperature = 750 * kelvin
friction = 100 / picosecond
timestep = 10.0 * femtosecond

def generate_openmm():

    # Prepare the system
    system = System()
    mullerforce = mf.MullerForce()
    system.addParticle(mass)
    mullerforce.addParticle(0, [])
    system.addForce(mullerforce)

    print("Writing system.xml file")
    with open('system.xml', 'w') as f:
        f.write(XmlSerializer.serialize(system))

    # Prepare integrator

    filename = "integrator.xml"
    print("Writing %s file" % filename)
    integrator = LangevinIntegrator(temperature, friction, timestep)
    with open(filename, 'w') as f:
        f.write(XmlSerializer.serialize(integrator))

def generate_starting_structures_fixed():
    pdb = mdtraj.load("single.pdb")
    xyz = np.array([
            [[0.0, 0.0, 0.0]]
            ])
    traj = mdtraj.Trajectory(xyz, pdb.topology)
    print("Writing seed_structures.h5")
    traj.save("seed_structures.h5")

def randx():
    return random.random() * (1.2 + 1.5) - 1.5
def randy():
    return random.random() * (2.0 + 0.2) - 0.2

def generate_starting_structures_random(num=10):
    pdb = mdtraj.load("single.pdb")
    xyz = list()
    for i in xrange(num):
        xyz.append([[randx(), randy(), 0.0]])
    traj = mdtraj.Trajectory(np.array(xyz), pdb.topology)
    print("Writing seed_structures.h5")
    traj.save("seed_structures.h5")


if __name__ == "__main__":
    error = "Specify what structures to generate: random, fixed, none"
    if len(argv) < 2:
        print(error)
        sys.exit(1)
    if argv[1] == 'random':
        generate_starting_structures_random()
    elif argv[1] == 'fixed':
        generate_starting_structures_fixed()
    elif argv[1] == 'none':
        pass
    else:
        print(error)
        sys.exit(1)
    generate_openmm()


