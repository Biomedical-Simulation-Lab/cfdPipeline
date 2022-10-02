#!/bin/python
# ///////////////////////////////////////////////////////////////
#  HPC Job Command Wrapper specific to Niagara
#  Author: Mehdi Najafi, mnuoft@gmail.com
#  Date: 2018-03-14
#
#  This library is not intended for distributions in any form.
# ///////////////////////////////////////////////////////////////

from __future__ import print_function

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-01-08"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private; to be obtained directly from the author."

import inspect, sys, subprocess, time, os
srcpath = inspect.getsourcefile(lambda:0)
srcpath = srcpath[:srcpath.rfind('/')]
sys.path.append(srcpath+'/../')
import naming

# retrive the case name using the naming convention
def get_case_name(result_folder):
     if result_folder.find('art_') >= 0:
         txt = result_folder[result_folder.rfind("art_")+4:]
         pos = txt.rfind('_I')
         name = ''
         while ( pos > 0 ):
             #print (txt, pos, txt[:pos], txt[pos+2:pos+3])
             #print (result_folder[:result_folder.rfind('art_')] + '../data/' + txt[:pos] + '.h5')
             if txt[pos+2:pos+3].isdigit():
                 #if os.path.exists( result_folder[:result_folder.rfind('art_')] + '../data/' + txt[:pos] + '.h5' ):
                     name = txt[:pos]
             pos = txt.rfind('_I',0,pos)
         if len(name) > 0: return name
         print ("<!> Cannot find a valid case name from its given result folder name:", result_folder)
         return None
     elif result_folder.find('pipe') >= 0:
         return result_folder.split("_ipcs_ab_cn_")[1].split("_constant")[0]
     return result_folder.split('_')[1]

# retrive the case mesh filename using the naming convention and given results folder
def get_case_mesh_filename(result_folder):
     cn = get_case_name(result_folder)
     txt = result_folder[:result_folder.rfind("results/")]
     main_folder = txt if len(txt) > 0 else './'
     return main_folder + '/data/'+cn+'.h5', main_folder, cn

# retrive the period by which the case was simulated in seconds NOT in miliseconds
def get_period(result_folder):
    return naming.get_period_from_info(get_case_name(result_folder)) / 1000.0

# provide a sulrm job script
def get_job_script_header(job_name, req_time=15, nodes=1, run_on_debug=False):
    partition = 'debug' if run_on_debug else 'compute'
    cmdtxt = """#!/bin/bash
#SBATCH --partition=%s
#SBATCH --time=%02d:%02d:00
#SBATCH --mail-type=NONE
#SBATCH --nodes=%d
#SBATCH --ntasks-per-node=80
#SBATCH --job-name %s
#SBATCH --output=hpclog/%s_stdout_%%j.txt
#SBATCH --error=hpclog/%s_stderr_%%j.txt

cd $SLURM_SUBMIT_DIR

module purge
# EXECUTION COMMANDs
\n"""%(partition, int(req_time/60), req_time%60, nodes, job_name, job_name, job_name)

    return cmdtxt

# submit a job and return the job id
def submit_a_job(job_type_name, job_cmd_text, depend):
    temp_file = 'run_me_%s_%s.sh'%(job_type_name, str(time.time()) )
    with open(temp_file, 'w') as job_file:
        job_file.write(job_cmd_text)
    cmd = ['sbatch', temp_file]
    if depend:
        cmd.insert(1, '--dependency=afterok:%s'%depend)
    process = subprocess.Popen(cmd, stdout=subprocess.PIPE)
    out, err = process.communicate()
    print (out)
    os.remove(temp_file)
    if err:
        return -1
    if len(out.split()) == 0: return -1
    return out.split()[-1]

#  provide the requested environment configuration
def get_env_settings(target_name):
    if target_name == 'PARALLEL':
        return "module load gnu-parallel"

    if target_name == 'HDF5':
        # return "module load intel/2018.2 hdf5/1.8.20"
        return "module load intel/2019u4 intelmpi/2019u4 hdf5/1.8.21"+'\n'# python/3.8.5"

    if target_name == 'SOLVER_HOME':
        # return 'SOLVER_HOME=/home/s/steinman/mnajafi/Solver'
        return 'SOLVER_HOME=~/..//mnajafi/.conda/envs/.test'

    if target_name == 'SOLVER_POST_HOME':
        # return 'SOLVER_POST_HOME=/home/s/steinman/mnajafi/Solver/Post'
        return 'SOLVER_POST_HOME=~/../mnajafi/.conda/envs/.test/Post'

    if target_name == 'FEniCS_HOME':
        return 'FEniCS_HOME=/home/s/steinman/mnajafi/.conda/envs/fenics201810'

    if target_name == 'SOLVER':
        return  get_env_settings('FEniCS_HOME') + '\n' + \
                get_env_settings('SOLVER_HOME') + '\n' + \
"""
export DIJITSO_CACHE_DIR=$SCRATCH/.local/dijitso-cache/
export INSTANT_CACHE_DIR=$SCRATCH/.local/instant-cache/
export INSTANT_ERROR_DIR=$SCRATCH/.local/instant-error/

# Turn off implicit threading
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1

module load NiaEnv/2018a
module load anaconda3/5.1.0 intel/2018.2 intelmpi/2018.2 boost/1.66.0 petsc/3.8.4 hdf5-mpi/1.8.20 netcdf-mpi/4.6.1 trilinos/12.12.1 fenics/2018.1.0
source activate $FEniCS_HOME
"""

    if target_name == 'SOLVER_ENV':
        return  get_env_settings('FEniCS_HOME') + '\n' + \
                get_env_settings('SOLVER_HOME') + '\n' + \
                get_env_settings('SOLVER_POST_HOME') + '\n' + \
"""
# Turn off implicit threading
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

module load NiaEnv/2018a
module load anaconda3/5.1.0 intel/2018.2 intelmpi/2018.2 boost/1.66.0 hdf5-mpi/1.8.20
source activate $FEniCS_HOME
"""

    if target_name == 'SOLVER_ENV_MT':
        return  get_env_settings('FEniCS_HOME') + '\n' + \
                get_env_settings('SOLVER_HOME') + '\n' + \
                get_env_settings('SOLVER_POST_HOME') + '\n' + \
"""
# Turn off implicit threading
#export OMP_NUM_THREADS=1
#export OPENBLAS_NUM_THREADS=1
#export MKL_NUM_THREADS=1

module load NiaEnv/2018a
module load anaconda3/5.1.0 intel/2018.2 intelmpi/2018.2 boost/1.66.0 hdf5-mpi/1.8.20
source activate $FEniCS_HOME
"""

    return None


def split_into_groups(a,n):
    k, m = divmod(len(a), n)
    return [a[i * k + min(i, m):(i + 1) * k + min(i + 1, m)] for i in range(n)]

