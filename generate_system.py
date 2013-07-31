#!/usr/bin/env python

from simtk.openmm.app import *
from simtk.openmm import *
from simtk.unit import *
from sys import stdout
import random

def generate(num=4):
    pdb = PDBFile('single.pdb')
    forcefield = ForceField('./no-ff.xml')

    system = forcefield.createSystem(pdb.topology,
            nonbondedMethod=CutoffNonPeriodic, removeCMMotion=False)
    #force = CustomExternalForce("0.5*kx*x^2+0.5*ky*y^2; kx=1; ky=10")
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
        integrator.setRandomNumberSeed(random.randint(0,100000))
        with open(filename, 'w') as f:
            f.write(XmlSerializer.serialize(integrator))

if __name__ == "__main__":
    generate()

# simulation = Simulation(pdb.topology, system, integrator)
# simulation.context.setPositions(pdb.positions)
# simulation.minimizeEnergy()
# simulation.reporters.append(DCDReporter('single-out.dcd', 100))
# simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))
# simulation.step(10000)
