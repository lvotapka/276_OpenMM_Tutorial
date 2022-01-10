from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Read the PSF
psf = CharmmPsfFile('charmm-gui-4183445964/step3_pbcsetup.psf')
psf.setBox(64.0*angstroms, 64.0*angstroms, 64.0*angstroms)

# Get the coordinates from the PDB
pdb = PDBFile('charmm-gui-4183445964/step3_pbcsetup.pdb')

# Load the parameter set.
params = CharmmParameterSet(
    'charmm-gui-4183445964/toppar/par_all36m_prot.prm', 
    'charmm-gui-4183445964/toppar/top_all36_prot.rtf', 
    'charmm-gui-4183445964/toppar/toppar_water_ions.str')

# Create the system, defining periodic boundary conditions, a nonbonded cutoff, and constrain bonds with hydrogen to remain constant-length.
system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Create a Langevin integrator for constant temperature.
integrator = LangevinMiddleIntegrator(300.0*kelvin, 1/picosecond, 0.004*picoseconds)

# Define a barostat to maintain constant pressure for equilibration
barostat = MonteCarloBarostat(1.0*bar, 300.0*kelvin, 25)
system.addForce(barostat) 

# Create a simulation object.
simulation = Simulation(psf.topology, system, integrator)

# Uncomment the following line to create a simulation that uses one of the GPUs.
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
#simulation = Simulation(pdb.topology, system, integrator, platform, properties)

# Set the system's atomic positions.
simulation.context.setPositions(pdb.positions)

# Set the system's atomic velocities to a random distribution defined by a Maxwell-Boltzmann distribution.
simulation.context.setVelocitiesToTemperature(300.0*kelvin)

# Perform gradient descent to find a local energy minimum.
simulation.minimizeEnergy()

# Create a reporter to write the trajectory to a file.
simulation.reporters.append(PDBReporter('equilibration.pdb', 1000))

# Create another reporter to display system info as the simulation runs.
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))

# Advance time by 10000 timesteps (40 picoseconds for Langevin integrator).
simulation.step(10000)

# Save a state file to maintain positions and velocities into the production stage.
state = simulation.context.getState(getPositions=True, getVelocities=True, getParameters=False, getIntegratorParameters=False)
xml = XmlSerializer.serialize(state)
with open('equil_done.xml', 'w') as f:
    f.write(xml)

