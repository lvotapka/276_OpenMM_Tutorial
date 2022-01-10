from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

# Read the PSF
psf = CharmmPsfFile('charmm-gui-4183093546/step5_assembly.psf')
psf.setBox(40.68*angstroms, 40.68*angstroms, 92.50*angstroms)

# Get the coordinates from the PDB
pdb = PDBFile('charmm-gui-4183093546/step5_assembly.pdb')

# Load the parameter set.
params = CharmmParameterSet('charmm-gui-4183093546/toppar/par_all36m_prot.prm', 
    'charmm-gui-4183093546/toppar/top_all36_prot.rtf', 
    'charmm-gui-4183093546/toppar/par_all36_lipid.prm', 
    'charmm-gui-4183093546/toppar/top_all36_lipid.rtf', 
    'charmm-gui-4183093546/toppar/toppar_water_ions.str')
    
# Create the system, defining periodic boundary conditions, a nonbonded cutoff, and constrain bonds with hydrogen to remain constant-length.
system = psf.createSystem(params, nonbondedMethod=PME, nonbondedCutoff=1*nanometer, constraints=HBonds)

# Create a Langevin integrator for constant temperature.
integrator = LangevinMiddleIntegrator(300.0*kelvin, 1/picosecond, 0.004*picoseconds)

# This barostat is designed for membrane simulations.
barostat = MonteCarloMembraneBarostat(1.0*bar, 0.0*bar*nanometers, 310.0*kelvin, MonteCarloMembraneBarostat.XYIsotropic, MonteCarloMembraneBarostat.ZFixed)
system.addForce(barostat) 

# Create a simulation object.
simulation = Simulation(psf.topology, system, integrator)

# Uncomment the following line to create a simulation that uses one of the GPUs.
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
#simulation = Simulation(pdb.topology, system, integrator, platform, properties)

# Read the equilibration state file for system positions/velocities.
with open('equil_done.xml', 'r') as f:
    xml = f.read()

simulation.context.setState(XmlSerializer.deserialize(xml))

# Create a reporter to write the trajectory to a file.
simulation.reporters.append(PDBReporter('production.pdb', 1000))

# Create another reporter to display system info as the simulation runs.
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True))

# Advance time by 10000 timesteps (40 picoseconds for Langevin integrator).
simulation.step(10000)

