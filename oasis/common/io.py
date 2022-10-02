__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-07-26"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__ = "GNU Lesser GPL version 3 or any later version"

from os import makedirs, getcwd, listdir, remove, system, path
import pickle
from dolfin import (MPI, Function, XDMFFile, HDF5File, ##info_red,
    VectorFunctionSpace, FunctionAssigner)

__all__ = ["create_initial_folders", "save_solution", "save_tstep_solution_h5",
           "save_checkpoint_solution_h5", "check_if_kill", "check_if_reset_statistics",
           "init_from_restart"]

master = MPI.rank(MPI.comm_world) == 0


def create_initial_folders(folder, restart_folder, sys_comp, tstep, info_red,
                           scalar_components, output_timeseries_as_vector,
                           **NS_namespace):
    """Create necessary folders."""
    #info_red("Creating initial folders")
    # To avoid writing over old data create a new folder for each run
    if master:
        try:
            makedirs(folder)
        except OSError:
            pass

    MPI.barrier(MPI.comm_world)
    newfolder = path.join(folder, 'data')
    if restart_folder:
        #newfolder = path.join(newfolder, restart_folder.split('/')[-2])
        previous_folder = listdir(newfolder)
        #print ('restart_folder=',restart_folder,'\n newfolder=',newfolder,'\nprevious_folder=',previous_folder)
        previous_folder = max(map(eval, previous_folder)) if previous_folder else 0
        #print ('previous_folder=',previous_folder)
        newfolder = path.join(newfolder, str(previous_folder + 1))
    else:
        if not path.exists(newfolder):
            newfolder = path.join(newfolder, '1')
        else:
            previous_folder = listdir(newfolder)
            previous_folder = max(map(eval, previous_folder)) if previous_folder else 0
            newfolder = path.join(newfolder, str(previous_folder + 1))
    # remove terminate_solver
    if master and path.exists(path.join(newfolder,'terminate_solver')):
        remove(path.join(newfolder,'terminate_solver'))

    MPI.barrier(MPI.comm_world)
    if master:
        if not restart_folder:
            #makedirs(path.join(newfolder, "Timeseries"))
            makedirs(path.join(newfolder, "Checkpoint"))
        else:
            makedirs(path.join(newfolder, "Checkpoint"))

    tstepfolder = path.join(newfolder, "Timeseries")
    tstepfiles = {}
    comps = sys_comp
    if output_timeseries_as_vector:
        comps = ['p', 'u'] + scalar_components

    for ui in comps:
        tstepfiles[ui] = XDMFFile(MPI.comm_world, path.join(
            tstepfolder, ui + '_from_tstep_{}.xdmf'.format(tstep)))
        tstepfiles[ui].parameters["rewrite_function_mesh"] = False
        tstepfiles[ui].parameters["flush_output"] = True

    return newfolder, tstepfiles


def save_solution(tstep, t, q_, q_1, folder, newfolder, save_step, checkpoint,
                  NS_parameters, tstepfiles, u_, u_components, scalar_components,
                  output_timeseries_as_vector, constrained_domain,
                  AssignedVectorFunction, **NS_namespace):
    """Called at end of timestep. Check for kill and save solution if required."""
    NS_parameters.update(t=t, tstep=tstep)
    if tstep % save_step == 0:
        save_tstep_solution_h5(tstep, q_, u_, newfolder, tstepfiles, constrained_domain,
                               output_timeseries_as_vector, u_components, AssignedVectorFunction,
                               scalar_components, NS_parameters)

    terminate_solver = check_if_kill(newfolder)
    if tstep % checkpoint == 0 or terminate_solver:
        save_checkpoint_solution_h5(tstep, q_, q_1, newfolder, u_components,
                                    NS_parameters)

    return terminate_solver


def save_tstep_solution_h5(tstep, q_, u_, newfolder, tstepfiles, constrained_domain,
                           output_timeseries_as_vector, u_components, AssignedVectorFunction,
                           scalar_components, NS_parameters):
    """Store solution on current timestep to XDMF file."""
    timefolder = path.join(newfolder, 'Timeseries')
    if output_timeseries_as_vector:
        # project or store velocity to vector function space
        for comp, tstepfile in tstepfiles.items():
            if comp == "u":
                # Create vector function and assigners
                uv = AssignedVectorFunction(u_)

                # Assign solution to vector
                uv()

                # Store solution vector
                tstepfile.write(uv, float(tstep))

            elif comp in q_:
                tstepfile.write(q_[comp], float(tstep))

            else:
                tstepfile.write(tstepfile.function, float(tstep))

    else:
        for comp, tstepfile in tstepfiles.items():
            tstepfile << (q_[comp], float(tstep))

    if master:
        if not path.exists(path.join(timefolder, "params.dat")):
            f = open(path.join(timefolder, 'params.dat'), 'wb')
            pickle.dump(NS_parameters,  f)
            f.close()


def save_checkpoint_solution_h5(tstep, q_, q_1, newfolder, u_components,
                                NS_parameters):
    """Overwrite solution in Checkpoint folder.

    For safety reasons, in case the solver is interrupted, take backup of
    solution first.

    Must be restarted using the same mesh-partitioning. This will be fixed
    soon. (MM)

    """
    checkpointfolder = path.join(newfolder, "Checkpoint")
    NS_parameters["num_processes"] = MPI.size(MPI.comm_world)
    if master:
        if path.exists(path.join(checkpointfolder, "params.dat")):
            system('cp {0} {1}'.format(path.join(checkpointfolder, "params.dat"),
                                       path.join(checkpointfolder, "params_old.dat")))
        f = open(path.join(checkpointfolder, "params.dat"), 'wb')
        pickle.dump(NS_parameters,  f)
        f.close()

    MPI.barrier(MPI.comm_world)
    for ui in q_:
        h5file = path.join(checkpointfolder, ui + '.h5')
        oldfile = path.join(checkpointfolder, ui + '_old.h5')
        # For safety reasons...
        if path.exists(h5file):
            if master:
                system('cp {0} {1}'.format(h5file, oldfile))
        MPI.barrier(MPI.comm_world)
        ###
        newfile = HDF5File(MPI.comm_world, h5file, 'w')
        newfile.flush()
        newfile.write(q_[ui].vector(), '/current')
        #if ui in u_components:
        newfile.write(q_1[ui].vector(), '/previous')
        if path.exists(oldfile):
            if master:
                system('rm {0}'.format(oldfile))
        MPI.barrier(MPI.comm_world)
    if master and path.exists(path.join(checkpointfolder, "params_old.dat")):
        system('rm {0}'.format(path.join(checkpointfolder, "params_old.dat")))


def check_if_kill(folder):
    """Check if user has put a file named killoasis in folder."""
    found = 0
    if 'terminate_solver' in listdir(folder):
        found = 1
    collective = MPI.sum(MPI.comm_world, found)
    if collective > 0:
        if master:
            remove(path.join(folder, 'terminate_solver'))
            #info_red('terminate_solver Found! Stopping simulations cleanly...')
        return True
    else:
        return False


def check_if_reset_statistics(folder):
    """Check if user has put a file named resetoasis in folder."""
    found = 0
    if 'resetoasis' in listdir(folder):
        found = 1
    collective = MPI.sum(MPI.comm_world, found)
    if collective > 0:
        if master:
            remove(path.join(folder, 'resetoasis'))
            #info_red('resetoasis Found!')
        return True
    else:
        return False


def init_from_restart(restart_folder, sys_comp, uc_comp, u_components,
                      q_, q_1, q_2, **NS_namespace):
    """Initialize solution from checkpoint files """
    if restart_folder:
        for ui in sys_comp:
            filename = path.join(restart_folder, ui + '.h5')
            hdf5_file = HDF5File(MPI.comm_world, filename, "r")
            hdf5_file.read(q_[ui].vector(), "/current", False)
            q_[ui].vector().apply('insert')
            # Check for the solution at a previous timestep as well
            if ui in uc_comp:
                q_1[ui].vector().zero()
                q_1[ui].vector().axpy(1., q_[ui].vector())
                q_1[ui].vector().apply('insert')
                if ui in u_components:
                    hdf5_file.read(q_2[ui].vector(), "/previous", False)
                    q_2[ui].vector().apply('insert')
            else:
                if hdf5_file.has_dataset("/previous"):
                    hdf5_file.read(q_1[ui].vector(), "/previous", False)
                    q_1[ui].vector().apply('insert')


