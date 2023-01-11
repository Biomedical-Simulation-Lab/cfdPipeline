#!/bin/python
# ///////////////////////////////////////////////////////////////
#  Job to Calculate Hemodynamics: TAWSS, OSI, SPI, RRT, ...
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
import job_utils

def get_job_env_setting():
    return job_utils.get_env_settings('SOLVER_ENV')

def get_hemo_job_cmd(input_folders, interval, ncores):
    cmdtxt = '# Hemodynamics Calculations\n'
    for input_folder in input_folders:
        mesh_filename, case_folder, case_name = job_utils.get_case_mesh_filename(input_folder)
        period = job_utils.get_period(input_folder)
        cmdtxt += '  python $SOLVER_POST_HOME/hemodynamic_indices.py %s %f %d %d > %s/logs/hemodynamics_%s.log\n'%(input_folder,period,interval,ncores,case_folder,case_name)
    return cmdtxt

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Hemodynamics Calculation Job - August 2018')
    parser.add_argument('result_folder', metavar='folder', nargs='+',
                        help='the folder where wss h5 files have been stored.')
    parser.add_argument('-n', '--ncores', dest='ncores', type=int, default=multiprocessing.cpu_count()*2,
                        help='the number of cores to be used.')
    parser.add_argument('-p', '--depend', dest='depend', type=int, default=0,
                        help='the job id to start the job after.')
    parser.add_argument('-t', '--time', dest='required_time', type=int, default=40,
                        help='the required time to finish the job')
    parser.add_argument('-d', '--debug', dest='debug', action='store_true', default=False,
                        help='submit the job to a debug node. Note that only ONE debug job is allowed to be running or in queue.')
    parser.add_argument('-i', '--interval', dest='interval', type=int, default=1,
                        help='the file increment to use to do modal analysis.')
    parser.add_argument('-m', '--mu', dest='mu', type=float, default=3.70,
                        help='dynamic viscosity by which the velocity gradients are multiplied to compute shear stress. Unit is g/ms.')
    parser.add_argument('-redo', '--redo', dest='redo', action='store_true', default=False,
                        help='ignore any previous attempts and recalculate.')
    args = parser.parse_args()

    partition = 'debug' if args.debug else 'compute'

    case_names = [job_utils.get_case_name(r) for r in args.result_folder]

    print ('Calculations for %d cases will be run on %s partition.'%(len(case_names),partition) )

    job_cmd_text = job_utils.get_job_script_header('hem_'+'+'.join(case_names), req_time=args.required_time, nodes=1, run_on_debug=args.debug)
    job_cmd_text += get_job_env_setting()
    job_cmd_text += get_hemo_job_cmd(args.result_folder, interval=args.interval, ncores=args.ncores)

    job_utils.submit_a_job('hemo', job_cmd_text, depend=str(args.depend) if args.depend>0 else '')
