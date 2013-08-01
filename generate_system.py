#!/usr/bin/env python

from simtk.openmm import *
from simtk.openmm.app import *
from simtk.unit import *
from sys import stdout, argv
import mullerforce as mf
import numpy as np
import random
import mdtraj

# Global parameters
nParticles = 100
mass = 1.0 * dalton
temperature = 750 * kelvin
friction = 100 / picosecond
timestep = 10.0 * femtosecond

def generate(num=4):
    pdb = PDBFile('single.pdb')
    forcefield = ForceField('./no-ff.xml')

    system = forcefield.createSystem(pdb.topology,
            nonbondedMethod=CutoffNonPeriodic, removeCMMotion=False)
    # force = CustomExternalForce("0.5*kx*x^2+0.5*ky*y^2; kx=1; ky=10")
    force = CustomExternalForce("kx*(x * (x-2)*(x-4))^2 + ky*(y)^2; kx=0.5; ky=0.1")
    system.addForce(force)
    for i in range(system.getNumParticles()):
        force.addParticle(i, [])

    print("Writing system.xml file")
    with open('system.xml', 'w') as f:
        f.write(XmlSerializer.serialize(system))

    for i in xrange(num):
        filename = "integrator%d.xml" % i
        print("Writing %s file" % filename)
        integrator = LangevinIntegrator(300 * kelvin, 1 / picosecond, 0.002 * picoseconds)
        integrator.setRandomNumberSeed(random.randint(0, 100000))
        with open(filename, 'w') as f:
            f.write(XmlSerializer.serialize(integrator))
            
            
def generate3(num=4):
    
    system = System()
    mullerforce = mf.MullerForce()
    system.addParticle(mass)
    mullerforce.addParticle(0, [])
    system.addForce(mullerforce)

    print("Writing system.xml file")
    with open('system.xml', 'w') as f:
        f.write(XmlSerializer.serialize(system))

    for i in xrange(num):
        filename = "integrator%d.xml" % i
        print("Writing %s file" % filename)
        integrator = LangevinIntegrator(temperature, friction, timestep)
        integrator.setRandomNumberSeed(random.randint(0, 100000))
        with open(filename, 'w') as f:
            f.write(XmlSerializer.serialize(integrator))
            
def generate_starting_structures_traj():
    pdb = mdtraj.load("single.pdb")
    xyz = np.array([
            [[0.0,0.0,0.0]],
            [[-1.0,1.0,0.0]],
            [[0.0,0.5,0.0]]
            ])
    traj = mdtraj.Trajectory(xyz, pdb.topology)
    print("Writing seed_structures.h5")
    traj.save("seed_structures.h5")

# simulation = Simulation(pdb.topology, system, integrator)
# simulation.context.setPositions(pdb.positions)
# simulation.minimizeEnergy()
# simulation.reporters.append(DCDReporter('single-out.dcd', 100))
# simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
# simulation.step(10000)





def generate2():

    # Choose starting conformations uniform on the grid between (-1.5, -0.2) and (1.2, 2)
    startingPositions = (np.random.rand(nParticles, 3) * np.array([2.7, 1.8, 1])) + np.array([-1.5, -0.2, 0])

    system = System()
    mullerforce = mf.MullerForce()
    for i in range(nParticles):
        system.addParticle(mass)
        mullerforce.addParticle(i, [])
    system.addForce(mullerforce)

    print("Writing system.xml file")
    with open('system.xml', 'w') as f:
        f.write(XmlSerializer.serialize(system))

    integrator = LangevinIntegrator(temperature, friction, timestep)

    print("Writing integrator0.xml file")
    with open('integrator0.xml', 'w') as f:
        f.write(XmlSerializer.serialize(integrator))

#     context = Context(system, integrator)
#     context.setPositions(startingPositions)
#     context.setVelocitiesToTemperature(temperature)
# 
#     return context


    # MullerForce.plot(ax=pp.gca())

#     for i in range(1000):
#         x = context.getState(getPositions=True).getPositions(asNumpy=True).value_in_unit(nanometer)
#         pp.scatter(x[:, 0], x[:, 1], edgecolor='none', facecolor='k')
#         integrator.step(100)
#
#     pp.show()

if __name__ == "__main__":
    if len(argv) < 2:
        print("Please specify which method to run. Options: 1 or 2")
        exit(1)
    generate_starting_structures_traj()
    if argv[1] == '1':
        if len(argv) < 3:
            generate()
        else:
            generate(int(argv[2]))
    elif argv[1] == '2':
        generate2()
    elif argv[1]=='3':
        if len(argv) < 3:
            generate3()
        else:
            generate3(int(argv[2]))
    else:
        print("Invalid option.")
