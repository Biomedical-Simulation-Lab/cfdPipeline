#!/bin/python
# ///////////////////////////////////////////////////////////////
#  Job to Calculate Wall Shear Stress
#  Author: Mehdi Najafi, mnuoft@gmail.com
#  Date: 2017-10-01
#
#  This library is not intended for distributions in any form and
#  distribution of it is not allowed in any form.
# ///////////////////////////////////////////////////////////////

from __future__ import print_function

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2017-10-01"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private; to be obtained directly from the author."

import argparse, multiprocessing
import os, glob
import job_utils

def get_job_env_setting():
    return '# Turn off implicit threading\n export OMP_NUM_THREADS=1 \n module load NiaEnv/2019b' + \
            job_utils.get_env_settings('PARALLEL')+'\n' +\
            job_utils.get_env_settings('HDF5')+'\n'+ \
            job_utils.get_env_settings('SOLVER_POST_HOME') + '\n'

def get_wss_job_cmd(input_folder, nproc, mu, redo=False, force_recalc_matrices=False):

    print ('Viscoity is set to mu = ', mu)

    mesh_filename, case_folder, mesh_name = job_utils.get_case_mesh_filename(input_folder)

    h5_files = glob.glob(os.path.join(input_folder,"*_up.h5"))
    if not h5_files:
        print ("There are no h5 files in the results folders. Exiting...")
        exit(1)

    print ('Found ', len(h5_files),'files to process.')

    if not redo:
       # remove those already computed
       wss_files = glob.glob(os.path.join(input_folder,"wss_files","*_wss.h5"))
       print ('Number of previoulsy computed WSS files:', len(wss_files))
       for wf in wss_files:
          match = wf.replace('wss_files/','').replace('_wss.h5','_up.h5')
          #print ('checking ', match, wf, h5_files[0])
          if match in h5_files:
             # print ('removing ', match) 
             h5_files.remove(match)

    if len(h5_files) == 0:
       print ('OK, No files left to compute WSS.')
       exit(0)

    if not os.path.exists(input_folder+'/wss_files'):
        print ('Creating folder:', input_folder+'/wss_files')
        os.mkdir(input_folder+'/wss_files')
    print ('The results will be stored in:', input_folder+'/wss_files')

    print ('Grouping ', len(h5_files),'files ... ', end='')

    h5_files = sorted(h5_files, key=lambda N: int(N.split("ts=")[1].split("_up.h5")[0]))
    h5_files_group = job_utils.split_into_groups(h5_files, nproc)

    print('into', len(h5_files_group), 'groups of ~%d files.'%len(h5_files_group[0]))

    cmdtxt = '# WSS Calculations\n'

    if not (os.path.exists(mesh_filename[:-3]+'_normals.mat') and os.path.exists(mesh_filename[:-3]+'_gradient.mat') or force_recalc_matrices ):
        cmdtxt += """# Initiate
$SOLVER_POST_HOME/wsscalc %s
""" % (mesh_filename)

    cmdtxt += """
# EXECUTION COMMAND
parallel -j %d <<EOF
""" % (nproc)

    for i, grp in enumerate(h5_files_group):
        cmdtxt += '  $SOLVER_POST_HOME/wsscalc %g '%(mu) + mesh_filename+' '+ ' '.join(grp) + ' > %s/logs/wss_%s_%d.log \n'%(case_folder,mesh_name,i)
    cmdtxt += 'EOF'

    os.remove(mesh_filename[:-3]+'_normals.mat')
    os.remove(mesh_filename[:-3]+'_gradient.mat')

    return cmdtxt

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='WallShearStress Calculation Job - August 2018')
    parser.add_argument('result_folder', 
                        help='the folder where velocity h5 files have been stored.')
    parser.add_argument('-n', '--ncores', dest='ncores', type=int, default=multiprocessing.cpu_count()*2,
                        help='the number of cores to be used.')
    parser.add_argument('-p', '--depend', dest='depend', type=int, default=0,
                        help='the job id to start the job after.')
    parser.add_argument('-t', '--time', dest='required_time', type=int, default=60,
                        help='the required time to finish the job')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False,
                        help='submit the job to a debug node. Note that only ONE debug job is allowed to be running or in queue.')
    parser.add_argument('-m', '--mu', dest='mu', type=float, default=3.70,
                        help='dynamic viscosity by which the velocity gradients are multiplied to compute shear stress. Unit is g/ms.')
    parser.add_argument('-c', '--clean', dest='force_recalc_matrices', action='store_true', default=False,
                        help='recalculate the gradients and normal matrices.')
    parser.add_argument('-redo', '--redo', dest='redo', action='store_true', default=False,
                        help='ignore any previous attempts and recalculate.')
    args = parser.parse_args()

    partition = 'debug' if args.debug else 'compute'
    print ('Calculations will be run on %s partition.'%partition)

    case_name = job_utils.get_case_name(args.result_folder)

    job_cmd_text = job_utils.get_job_script_header('wss_'+case_name, req_time=args.required_time, nodes=1, run_on_debug=args.debug)
    job_cmd_text += get_job_env_setting()
    job_cmd_text += get_wss_job_cmd(args.result_folder, nproc=args.ncores, mu=args.mu, redo=args.redo, force_recalc_matrices=args.force_recalc_matrices)

    job_utils.submit_a_job('wss', job_cmd_text, depend=str(args.depend) if args.depend>0 else '')
