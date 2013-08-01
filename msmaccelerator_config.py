# Configuration file for msmaccelerator.

c = get_config()

#------------------------------------------------------------------------------
# App configuration
#------------------------------------------------------------------------------

# This is an application.

# App will inherit config from: Application

# The Logging format template
# c.App.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.App.log_level = 30

#------------------------------------------------------------------------------
# RootApplication configuration
#------------------------------------------------------------------------------

# This is an application.

# RootApplication will inherit config from: App, Application

# The Logging format template
# c.RootApplication.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.RootApplication.log_level = 30

#------------------------------------------------------------------------------
# Device configuration
#------------------------------------------------------------------------------

# This is an application.

# Device will inherit config from: App, Application

# ZeroMQ port to connect to the server on
# c.Device.zmq_port = 12345

# The Logging format template
# c.Device.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.Device.log_level = 30

# URL to connect to server with
# c.Device.zmq_url = '127.0.0.1'

#------------------------------------------------------------------------------
# Modeler configuration
#------------------------------------------------------------------------------

# This is an application.

# Modeler will inherit config from: Device, App, Application

# Set the log level by value or name.
# c.Modeler.log_level = 30

# URL to connect to server with
# c.Modeler.zmq_url = '127.0.0.1'

# Symmetrization method for constructing the reversibile counts matrix.
# c.Modeler.symmetrize = None

# Subsample data by taking only every stride-th point
# c.Modeler.stride = 1

# Lag time for building the model, in units of the stride. Currently, we are not
# doing the step in MSMBuilder that is refered to as "assignment", where you
# assign the remaining data that was not used during clustering to the cluster
# centers that were identified.
# c.Modeler.lag_time = 1

# File containing the indices of atoms to use in the RMSD computation. Using a
# PDB as input, this file can be created with the MSMBuilder script
# CreateAtomIndices.py
# c.Modeler.rmsd_atom_indices = 'AtomIndices.dat'

# Do ergodic trimming when constructing the Markov state model. This is
# generally a good idea for building MSMs in the high-data regime where you wish
# to prevent transitions that appear nonergodic because they've been
# undersampled from influencing your model, but is inappropriate in the sparse-
# data regime when you're using min-counts sampling, because these are
# precisiely the states that you're most interested in.
# c.Modeler.ergodic_trimming = False

# Distance cutoff for clustering, in nanometers. We will continue to create new
# clusters until each data point is within this cutoff from its cluster center.
# c.Modeler.rmsd_distance_cutoff = 0.2

# ZeroMQ port to connect to the server on
# c.Modeler.zmq_port = 12345

# The Logging format template
# c.Modeler.log_format = '[%(name)s] %(message)s'

#------------------------------------------------------------------------------
# Simulator configuration
#------------------------------------------------------------------------------

# This is an application.

# Simulator will inherit config from: Device, App, Application

# Path to the XML file containing the OpenMM Integrator to use
# c.Simulator.integrator_xml = 'integrator.xml'

# Set the log level by value or name.
# c.Simulator.log_level = 30

# Interval at which to save positions to a disk, in units of steps
c.Simulator.report_interval = 1000

# Path to the XML file containing the OpenMM system to propagate
# c.Simulator.system_xml = 'system.xml'

# Do local energy minimization on the configuration that's passed to me, before
# running dynamics
c.Simulator.minimize = False

# OpenMM device index for CUDA or OpenCL platforms. This is used to select which
# GPU will be used on a multi-gpu system. This option is ignored on reference
# platform
# c.Simulator.device_index = 0

# URL to connect to server with
# c.Simulator.zmq_url = '127.0.0.1'

# The OpenMM platform on which to run the simulation
# c.Simulator.platform = 'CUDA'

# Number of steps of dynamics to do
c.Simulator.number_of_steps = 100000

# ZeroMQ port to connect to the server on
# c.Simulator.zmq_port = 12345

# Choose random initial velocities from the Maxwell-Boltzmann distribution
# c.Simulator.random_initial_velocities = True

# The Logging format template
# c.Simulator.log_format = '[%(name)s] %(message)s'

#------------------------------------------------------------------------------
# Interactor configuration
#------------------------------------------------------------------------------

# This is an application.

# Interactor will inherit config from: Device, App, Application

# Go into interactive shell mode
# c.Interactor.shell = False

# Set the log level by value or name.
# c.Interactor.log_level = 30

# URL to connect to server with
# c.Interactor.zmq_url = '127.0.0.1'

# Set the server's beta parameter
# c.Interactor.set_beta = None

# The Logging format template
# c.Interactor.log_format = '[%(name)s] %(message)s'

# ZeroMQ port to connect to the server on
# c.Interactor.zmq_port = 12345

# Path to the database (sqlite3 file)
# c.Interactor.db_path = 'db.sqlite'

#------------------------------------------------------------------------------
# BaseServer configuration
#------------------------------------------------------------------------------

# This is an application.

# BaseServer will inherit config from: App, Application

# The Logging format template
# c.BaseServer.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.BaseServer.log_level = 30

# ZeroMQ port to serve on
# c.BaseServer.zmq_port = 12345

# Path to the database (sqlite3 file)
# c.BaseServer.db_path = 'db.sqlite'

#------------------------------------------------------------------------------
# AdaptiveServer configuration
#------------------------------------------------------------------------------

# This is an application.

# AdaptiveServer will inherit config from: BaseServer, App, Application

# Set the log level by value or name.
# c.AdaptiveServer.log_level = 30

# Path to the XML file containing the OpenMM system to propagate. This is
# required by the server to properly serialize the starting conformations.
# c.AdaptiveServer.system_xml = 'system.xml'

# Directory on the local filesystem where output trajectories will be saved
# c.AdaptiveServer.traj_outdir = 'trajs/'

# Directory on the local filesystem where MSMs will be saved.
# c.AdaptiveServer.models_outdir = 'models/'

# A PDB used to determine the system's topology. This is sent directly to the
# Simulator. Honestly, I'm not sure exactly why we need it. TODO: ask Peter
# about this.
c.AdaptiveServer.topology_pdb = 'single.pdb'

# ZeroMQ port to serve on
# c.AdaptiveServer.zmq_port = 12345

# 'Directory on the local filesystem where starting structures will be stored
# c.AdaptiveServer.starting_states_outdir = 'starting_states'

# The Logging format template
# c.AdaptiveServer.log_format = '[%(name)s] %(message)s'

# Path to the database (sqlite3 file)
# c.AdaptiveServer.db_path = 'db.sqlite'

#------------------------------------------------------------------------------
# MKProfile configuration
#------------------------------------------------------------------------------

# This is an application.

# MKProfile will inherit config from: App, Application

# The Logging format template
# c.MKProfile.log_format = '[%(name)s] %(message)s'

# Set the log level by value or name.
# c.MKProfile.log_level = 30

# Output directory in which to save the file msmbuilder_config.py
# c.MKProfile.output_dir = '.'

#------------------------------------------------------------------------------
# BaseSampler configuration
#------------------------------------------------------------------------------

# Adaptive sampler that simply returns a given initial structure.
# 
# This can be used for the very first round when you have no data, and also as a
# base class for the other samplers

# Trajectory file giving the initial structures that you want to sample from.
# This should be a single PDB or other type of loadable trajectory file. These
# structures will only be used in the beginning, before we have an actual MSM to
# use.
c.BaseSampler.seed_structures = 'seed_structures.h5'

#------------------------------------------------------------------------------
# CentroidSampler configuration
#------------------------------------------------------------------------------

# Adaptive sampler that uses the centroids/"generators"" of the states to choose
# amongst, with a multinomial distibition.
# 
# This class DOES NOT actually contain a method to *set* the weights. That is
# done by subclasses. See, for example, CountsSampler, that sets the weights
# using the counts.

# CentroidSampler will inherit config from: BaseSampler

# Trajectory file giving the initial structures that you want to sample from.
# This should be a single PDB or other type of loadable trajectory file. These
# structures will only be used in the beginning, before we have an actual MSM to
# use.
# c.CentroidSampler.seed_structures = 'ala5.pdb'

#------------------------------------------------------------------------------
# CountsSampler configuration
#------------------------------------------------------------------------------

# CountsSampler will inherit config from: CentroidSampler, BaseSampler

# Temperature factor that controls the level of exploration vs. refinement. When
# beta = 0 (high temp), we do full exploration, putting simulations where few
# counts have been seen. When beta = 1 (room temp), we do uniform sampling from
# the microstate, with no preference based on their counts. When beta > 1 (low
# temp), we do refinement, such that we focus on microstates with a high number
# of counts. At beta = 2, we choose microstates proportional to our estimate of
# their current equilibrium propbability. The explicit formula used is: Prob(
# choose state i ) ~ \sum_j C_{ij} ^{ beta - 1 }
# c.CountsSampler.beta = 1

# Trajectory file giving the initial structures that you want to sample from.
# This should be a single PDB or other type of loadable trajectory file. These
# structures will only be used in the beginning, before we have an actual MSM to
# use.
# c.CountsSampler.seed_structures = 'ala5.pdb'

