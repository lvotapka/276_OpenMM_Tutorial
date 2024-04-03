from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Load initial coordinates of atoms of the Villin protein in a water box.
pdb = PDBFile('input.pdb')

# Load the AMBER14 forcefield for proteins, and the TIP3P forcefield for water.
forcefield = ForceField('amber14-all.xml', 'amber14/tip3pfb.xml')

# Create the system, defining periodic boundary conditions, a nonbonded cutoff, and constrain bonds with hydrogen to remain constant-length.
system = forcefield.createSystem(pdb.topology, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Create a Verlet integrator for Newtonian mechanics.
integrator = VerletIntegrator(0.002*picoseconds)

# Uncomment the following line to use a Langevin integrator for maintaining constant temperature.
#integrator = LangevinMiddleIntegrator(300*kelvin, 1/picosecond, 0.004*picoseconds)

# Uncomment the following lines to create a barostat to maintain constant pressure. NOTE: this barostat should be used with the Langevin integrator, not the Verlet.
#barostat = MonteCarloBarostat(1.0*bar, 300.0*kelvin, 25)
#system.addForce(barostat) 

# Create a simulation object.
simulation = Simulation(pdb.topology, system, integrator)

# Uncomment the following lines to create a simulation that uses one of the GPUs.
# If you do that, you'll want to delete the 'simulation = ...' line above here.
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
#simulation = Simulation(pdb.topology, system, integrator, platform, properties)

# Set the system's atomic positions.
simulation.context.setPositions(pdb.positions)

# Set the system's atomic velocities to a random distribution defined by a Maxwell-Boltzmann distribution.
simulation.context.setVelocitiesToTemperature(300*kelvin)

# Perform gradient descent to find a local energy minimum.
simulation.minimizeEnergy()

# Create a reporter to write the trajectory to a file.
simulation.reporters.append(PDBReporter('output.pdb', 1000))

# Create another reporter to display system info as the simulation runs.
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# Advance time by 10000 timesteps (20 picoseconds for Verlet integrator).
simulation.step(10000)
