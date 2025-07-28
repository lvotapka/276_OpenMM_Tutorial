from openmm.app import *
from openmm import *
from openmm.unit import *
from sys import stdout

prmtop = AmberPrmtopFile('alanine_dipeptide.prmtop')
inpcrd = AmberInpcrdFile('alanine_dipeptide.inpcrd')

# Create the system with an implicit solvent and rigid bonds
system = prmtop.createSystem(implicitSolvent=OBC2, implicitSolventKappa=1.0/nanometer, constraints=HBonds)

output_xml_system_filename = "alanine_dipeptide.xml"
print("writing to file:", output_xml_system_filename)
with open(output_xml_system_filename, "w") as f:
    f.write(openmm.XmlSerializer.serialize(system))


exit()
# Create a Langevin integrator for constant temperature.
integrator = LangevinMiddleIntegrator(300.0*kelvin, 1/picosecond, 0.004*picoseconds)

# Create a simulation object.
simulation = Simulation(prmtop.topology, system, integrator)

# Uncomment the following line to create a simulation that uses one of the GPUs.
#platform = Platform.getPlatformByName('CUDA')
#properties = {'CudaDeviceIndex': '0', 'CudaPrecision': 'mixed'}
#simulation = Simulation(prmtop.topology, system, integrator, platform, properties)

# Assign atomic positions
simulation.context.setPositions(inpcrd.positions)

# Set the system's atomic velocities to a random distribution defined by a Maxwell-Boltzmann distribution.
simulation.context.setVelocitiesToTemperature(300.0*kelvin)
    
# Minimize system energy
simulation.minimizeEnergy()

# Create a reporter to write the trajectory to a file.
simulation.reporters.append(PDBReporter('alanine_dipeptide_traj.pdb', 1000))

# Create another reporter to display system info as the simulation runs.
simulation.reporters.append(StateDataReporter(stdout, 1000, step=True, potentialEnergy=True, temperature=True, volume=True))

# Advance time by 10000 timesteps (40 picoseconds for Langevin integrator).
simulation.step(10000)

