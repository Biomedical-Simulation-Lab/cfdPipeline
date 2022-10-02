__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2017-08-10"
__copyright__ = "Copyright (C) 2017 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"


from dolfin import *
import subprocess
from os import getpid, path
from collections import defaultdict
from numpy import array, maximum, zeros
import platform
from datetime import datetime

# UnitSquareMesh(20, 20) # Just due to MPI bug on Scinet
# try:
#from fenicstools import getMemoryUsage
# except:


master = MPI.rank(MPI.comm_world) == 0

## Welcome Message
if master:
    print ("**"*51)
    print ("Aneurysm CFD v.1.2 - by Mehdi Najafi, Copyright (c) 2018-%d - BioMedical Simulation Lab - University of Toronto"%(datetime.today().year))
    print (">>"*51)
    print ('Running on', platform.platform(), 'using',MPI.size(MPI.comm_world),'CPUs.\n           Python:', platform.python_version(), ' and DOLFIN:', cpp.__version__)
    print ("<<"*51)

def getMemoryUsage(rss=True):
    mypid = str(getpid())
    rss = "rss" if rss else "vsz"
    process = subprocess.Popen(['ps', '-o', rss, mypid],
                                stdout=subprocess.PIPE)
    out, _ = process.communicate()
    mymemory = out.split()[1]
    return eval(mymemory) / 1024


parameters["linear_algebra_backend"] = "PETSc"
parameters["form_compiler"]["optimize"] = True
parameters["form_compiler"]["cpp_optimize"] = True
parameters["form_compiler"]["representation"] = "uflacs" #"quadrature"
#parameters["form_compiler"]["quadrature_degree"] = 4
#parameters["form_compiler"]["cache_dir"] = "instant"
parameters["form_compiler"]["cpp_optimize_flags"] = "-O3"
#parameters["mesh_partitioner"] = "ParMETIS"
#parameters["form_compiler"].add("no_ferari", True)
#set_log_active(False)

# Default parameters for all solvers
NS_parameters = dict(
    nu=0.01,             # Kinematic viscosity
    folder='results',    # Relative folder for storing results
    velocity_degree=2,   # default velocity degree
    pressure_degree=1    # default pressure degree
)

# Default parameters NSfracStep solver
NS_parameters.update(
    # Physical constants and solver parameters
    t=0.0,               # Time
    tstep=0,             # Timestep
    T=1.0,               # End time
    dt=0.01,             # Time interval on each timestep

    # Some discretization options
    # Use Adams Bashforth projection as first estimate for pressure on new timestep
    AB_projection_pressure=False,
    solver="IPCS_ABCN",  # "IPCS_ABCN", "IPCS_ABE", "IPCS", "Chorin", "BDFPC", "BDFPC_Fast"

    # Parameters used to tweek solver
    max_iter=1,                 # Number of inner pressure velocity iterations on timestep
    max_error=1e-6,             # Tolerance for inner iterations (pressure velocity iterations)
    iters_on_first_timestep=2,  # Number of iterations on first timestep
    use_krylov_solvers=True,    # Otherwise use LU-solver
    print_intermediate_info=10,
    print_velocity_pressure_convergence=False,

    # Parameters used to tweek output
    plot_interval=10,
    checkpoint=10,       # Overwrite solution in Checkpoint folder each checkpoint
    save_step=10,        # Store solution each save_step
    restart_folder=None, # If restarting solution, set the folder holding the solution to start from here
    output_timeseries_as_vector=True,  # Store velocity as vector in Timeseries

    # Choose LES model and set default parameters
    # NoModel, Smagorinsky, Wale, DynamicLagrangian, ScaleDepDynamicLagrangian
    les_model='NoModel',

    # LES model parameters
    Smagorinsky=dict(Cs=0.1677),              # Standard Cs, same as OpenFOAM
    Wale=dict(Cw=0.325),
    DynamicSmagorinsky=dict(Cs_comp_step=1),  # Time step interval for Cs to be recomputed
    KineticEnergySGS=dict(Ck=0.08, Ce=1.05),

    # Parameter set when enabling test mode
    testing=False,

    # Solver parameters that will be transferred to dolfins parameters['krylov_solver']
    krylov_solvers=dict(
        monitor_convergence=False,
        report=False,
        error_on_nonconvergence=False,
        nonzero_initial_guess=True,
        maximum_iterations=200,
        relative_tolerance=1e-8,
        absolute_tolerance=1e-8),

    # Velocity update
    velocity_update_solver=dict(
        method='default',  # "lumping", "gradient_matrix"
        solver_type='cg',
        preconditioner_type='jacobi',
        low_memory_version=False),

    velocity_krylov_solver=dict(
        solver_type='bicgstab',
        preconditioner_type='jacobi'),

    pressure_krylov_solver=dict(
        solver_type='gmres',
        preconditioner_type='hypre_amg'),

    scalar_krylov_solver=dict(
        solver_type='bicgstab',
        preconditioner_type='jacobi'),

    nut_krylov_solver=dict(
        method='WeightedAverage',  # Or 'default'
        solver_type='cg',
        preconditioner_type='jacobi'),
)


NS_expressions = {}

constrained_domain = None

# To solve for scalars provide a list like ['scalar1', 'scalar2']
scalar_components = []

# With diffusivities given as a Schmidt number defined by:
#   Schmidt = nu / D (= momentum diffusivity / mass diffusivity)
Schmidt = defaultdict(lambda: 1.)
Schmidt_T = defaultdict(lambda: 0.7)  # Turbulent Schmidt number (LES)

Scalar = defaultdict(lambda: dict(Schmidt=1.0,
                                  family="CG",
                                  degree=1))

# The following helper functions are available in dolfin
# They are redefined here for printing only on process 0.
RED = "\033[1;37;31m%s\033[0m"
BLUE = "\033[1;37;34m%s\033[0m"
GREEN = "\033[1;37;32m%s\033[0m"
GRAY = "\033[1;37;30m%s\033[0m"


def info_blue(s, check=True):
    if master and check:
        print(BLUE % s)

def info_green(s, check=True):
    if master and check:
        print(GREEN % s)

def info_red(s, check=True):
    if master and check:
        print(RED % s)

def info_gray(s, check=True):
    if master and check:
        print ("\033[1;37;30m%s\033[0m"%s)

def info_cyan(s, check=True):
    if master and check:
        print ("\033[1;36;30m%s\033[0m"%s)

def info_yellow(s, check=True):
    if master and check:
        print ("\033[1;37;33m%s\033[0m"%s)

def info_purple(s, check=True):
    if master and check:
        print ("\033[1;37;35m%s\033[0m"%s)

def info_bkred(s, check=True):
    if master and check:
        print ("\033[1;37;41m%s\033[0m"%s)

def info_bkblue(s, check=True):
    if master and check:
        print ("\033[1;37;44m%s\033[0m"%s)


class OasisTimer(Timer):
    def __init__(self, task, verbose=False):
        Timer.__init__(self, task)
        info_blue(task, verbose)


class OasisMemoryUsage:
    def __init__(self, s):
        self.memory = 0
        self.memory_vm = 0
        self(s)

    def __call__(self, s, verbose=False):
        self.prev = self.memory
        self.prev_vm = self.memory_vm
        self.memory = MPI.sum(MPI.comm_world, getMemoryUsage())
        self.memory_vm = MPI.sum(MPI.comm_world, getMemoryUsage(False))
        if MPI.rank(MPI.comm_world) == 0 and verbose:
            info_blue('{0:26s}  {1:10d} MB {2:10d} MB {3:10d} MB {4:10d} MB'.format(s,
                        int(self.memory - self.prev), int(self.memory),
                        int(self.memory_vm - self.prev_vm), int(self.memory_vm)))


# Print memory use up til now
initial_memory_use = getMemoryUsage()
oasis_memory = OasisMemoryUsage('Start')


# Convenience functions
def strain(u):
    return 0.5 * (grad(u) + grad(u).T)


def omega(u):
    return 0.5 * (grad(u) - grad(u).T)


def Omega(u):
    return inner(omega(u), omega(u))


def Strain(u):
    return inner(strain(u), strain(u))


def QC(u):
    return Omega(u) - Strain(u)


def recursive_update(dst, src):
    """Update dict dst with items from src deeply ("deep update")."""
    for key, val in src.items():
        if key in dst and isinstance(val, dict) and isinstance(dst[key], dict):
            dst[key] = recursive_update(dst[key], val)
        else:
            dst[key] = val
    return dst


def add_function_to_tstepfiles(function, newfolder, tstepfiles, tstep):
    name = function.name()
    tstepfolder = path.join(newfolder, "Timeseries")
    tstepfiles[name] = XDMFFile(MPI.comm_world,
                                path.join(tstepfolder,
                                          '{}_from_tstep_{}.xdmf'.format(name, tstep)))
    tstepfiles[name].function = function
    tstepfiles[name].parameters["rewrite_function_mesh"] = False


def body_force(mesh, **NS_namespace):
    """Specify body force"""
    return Constant((0,) * mesh.geometry().dim())


def body_force_gravity(mesh, **NS_namespace):
    """Specify body force"""
    # The unit is m/(mili s.s)
    if master:
         print("Adding constant gravity of g=9.81E-3 m/(mili s.s).\n")
    return Constant((0,0,9.81E-3));# * mesh.geometry().dim())


def initialize(**NS_namespace):
    """Initialize solution."""
    pass


def create_bcs(sys_comp, **NS_namespace):
    """Return dictionary of Dirichlet boundary conditions."""
    return dict((ui, []) for ui in sys_comp)


def scalar_hook(**NS_namespace):
    """Called prior to scalar solve."""
    pass


def scalar_source(scalar_components, **NS_namespace):
    """Return a dictionary of scalar sources."""
    return dict((ci, Constant(0)) for ci in scalar_components)


def pre_solve_hook(**NS_namespace):
    """Called just prior to entering time-loop. Must return a dictionary."""
    return {}


def theend_hook(**NS_namespace):
    """Called at the very end."""
    pass


def problem_parameters(**NS_namespace):
    """Updates problem specific parameters, and handles restart"""
    pass


def post_import_problem(NS_parameters, mesh, commandline_kwargs,
                        NS_expressions, **NS_namespace):
    """Called after importing from problem."""

    # Update NS_parameters with all parameters modified through command line
    for key, val in commandline_kwargs.items():
        if isinstance(val, dict):
            NS_parameters[key].update(val)
        else:
            NS_parameters[key] = val

    # If the mesh is a callable function, then create the mesh here.
    if callable(mesh):
        mesh, dS, fd, nors, dim, inout_area = mesh(**NS_parameters)

    assert(isinstance(mesh, Mesh))

    # Returned dictionary to be updated in the NS namespace
    d = dict(mesh=mesh, dS=dS, subdomain_data=fd, normals=nors, dim=dim, inout_area=inout_area)
    d.update(NS_parameters)
    d.update(NS_expressions)
    return d

def velocity_tentative_hook(**NS_namespace):
    """Called just prior to solving for tentative velocity."""
    pass


def pressure_hook(**NS_namespace):
    """Called prior to pressure solve."""
    pass


def start_timestep_hook(**NS_namespace):
    """Called at start of new timestep"""
    pass


def temporal_hook(**NS_namespace):
    """Called at end of a timestep."""
    pass

