# ///////////////////////////////////////////////////////////////
#  Courant Number Calculation
#  Author: Mehdi Najafi, mnajafi@sharif.edu
#  Date: 2008-02-11
#
#  This library is not intended for distributions in any form and
#  distribution of it is not allowed in any form.
# ///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnajafi@sharif.edu>"
__date__ = "2008-02-11"
__copyright__ = "Copyright (C) 2008 " + __author__
__license__  = "Private; to be obtained directly from the author."

import sys, os
# os.environ['OPENBLAS_NUM_THREADS'] = '1'
# os.environ['MKL_NUM_THREADS'] = '1'

import h5py, glob
import numpy
import gc
import job_utils
from multiprocessing import sharedctypes
import multiprocessing

################################################################
shared_vars_pool={}
def create_shared_array(names, shape, dtype='d'):
    name_list = names.split(',')
    sz = int(numpy.prod(shape))
    print ('Allocating %dx %5.2f MB = %5.2f MB of memory ... '%(len(name_list),sz*8/(1024*1024),len(name_list)*sz*8/(1024*1024)) , end='', flush=True)
    for name in name_list:
        # shared_vars_pool[name+'_var'] = numpy.ctypeslib.as_ctypes(numpy.zeros(shape))
        # shared_vars_pool[name] = multiprocessing.sharedctypes.RawArray(shared_vars_pool[name+'_var']._type_, shared_vars_pool[name+'_var'])
        shared_vars_pool[name] = multiprocessing.sharedctypes.RawArray(dtype, sz)
        shared_vars_pool[name+'_shape'] = shape
    print('done.', flush=True)
def get_shared_var(name):
    #return numpy.ctypeslib.as_array(shared_vars_pool[name+'_var'])
    return numpy.frombuffer(shared_vars_pool[name]).reshape(shared_vars_pool[name+'_shape'])
################################################################
# get the element dimensions
#  dt*(u/dx+v/dy+w/dz)<=1
#  CFL = dt * (u/min(dx) + v/min(dy) + w/min(dz))
################################################################
# get the element dimensions
def get_mesh_deltas(mesh_h5_filename):
    mesh = h5py.File(mesh_h5_filename, 'r')
    cells = numpy.asarray(mesh['/Mesh/topology'])

    # look for previously generated one and load it
    if 'elem_dx' in mesh['/Mesh'].keys():
        return numpy.asarray(mesh['/Mesh/elem_dx']), cells

    points = numpy.asarray(mesh['/Mesh/coordinates'])
    number_of_cells = cells.shape[0]
    mesh_dx = numpy.zeros((number_of_cells,4), dtype=numpy.float64)
    for k in range(number_of_cells):
        dx = [1E20,1E20,1E20,0]
        point_ids = cells[k]
        for i in range(point_ids.shape[0]):
            for j in range(i+1,point_ids.shape[0]):
                ddx = numpy.fabs(points[point_ids[i]][0] - points[point_ids[j]][0])
                ddy = numpy.fabs(points[point_ids[i]][1] - points[point_ids[j]][1])
                ddz = numpy.fabs(points[point_ids[i]][2] - points[point_ids[j]][2])
                if ddx > 1E-9: dx[0] = min(dx[0], ddx)
                if ddy > 1E-9: dx[1] = min(dx[1], ddy)
                if ddz > 1E-9: dx[2] = min(dx[2], ddz)
        dx[3] = numpy.sqrt(numpy.mean(numpy.square(dx[:3])))
        mesh_dx[k] = numpy.asarray(dx)
    mesh.close()
    h5py.File(mesh_h5_filename,'a').create_dataset("Mesh/elem_dx", dtype=numpy.float64, data=mesh_dx, compression="gzip")
    return mesh_dx, cells

################################################################
# get the element dimensions
def compute_local_u_dx(ids, h5_files, mesh_dx, cells, write_each_timestep=0):
    if len(ids)==1:
        print ('    reading', len(ids), 'file:', ids, h5_files[ids[0]], flush=True)
    else:
        print ('    reading', len(ids), 'files:', ids, h5_files[ids[0]], ' ... ',  h5_files[ids[-1]], flush=True)

    number_of_cells = cells.shape[0]

    # get the shared variable
    u_dx_max = get_shared_var('u_dx_max')

    # loop over files
    for i in ids:
        hw = h5py.File(h5_files[i], 'r')
        lu = numpy.asarray(hw['/Solution/u'])

        u_dx = numpy.zeros((number_of_cells,4), dtype=numpy.float64)
        for j in range(number_of_cells):
            ids = cells[j]
            u_dx[j] = numpy.mean(lu[ids,0]) / mesh_dx[j][0] \
                    + numpy.mean(lu[ids,1]) / mesh_dx[j][1] \
                    + numpy.mean(lu[ids,2]) / mesh_dx[j][2]
        hw.close()

        u_dx_max[i] = numpy.max(u_dx)

        if write_each_timestep:
            # write to h5
            output_filename = h5_files[i].replace('_up.h5', '_u_dx.h5')
            print ('Writing data to %s ...'%(output_filename), end='')
            hf = h5py.File(output_filename, 'w')
            hf.create_dataset("Computed/u_dx", dtype=numpy.float64, data=u_dx, compression="gzip")
            hf.close()
            print (' done.')

################################################################
def compute_cfl(input_folder, interval, nproc, period, detailed):
    pos = -2 if (input_folder[-1] == '/') else -1
    folder_itself = input_folder.split('/')[pos]
    timesteps = int(folder_itself.split("_ts")[-1].split("_cy")[0])

    mesh_h5_filename, case_folder, case_name = job_utils.get_case_mesh_filename(input_folder)

    if not os.path.exists(mesh_h5_filename):
        print ('No mesh file found: %s \nYou may running this script from an incorrect folder.\n'%mesh_filename)
        exit(1)
    print ('Looking for mesh file:', mesh_h5_filename, ' and loading volume mesh (/Mesh/coordiantes, /Mesh/topology).')
    print ('Computing mesh element length scales.')
    mesh_dx, mesh_cells = get_mesh_deltas(mesh_h5_filename)

    print ('Looking inside', input_folder, '...')

    h5_files = glob.glob(input_folder+"/*_up.h5")
    if not h5_files:
        print ('No h5 file found in', input_folder+'. Run the solver and run this script after it finished.\n')
        exit(1)
    # sort the files according to the simulation time
    h5_files = sorted(h5_files, key=lambda N: int(N.split("ts=")[1].split("_up.h5")[0]))
    h5_files = h5_files[0:len(h5_files):interval]
    file_count = len(h5_files)

    print ('   found', file_count, 'up files for', mesh_cells.shape[0], 'elements.')

    # determine the time increment based on the number of samples and the corresponding FFT frequencies
    dt = period/timesteps

    # outfile = open(input_folder+'/%s_mesh_dx.txt'%(folder_itself), 'w')
    # for j in range(mesh_cells.shape[0]):
    #     outfile.write('%f\n'%(mesh_dx[j,]))
    # outfile.close()

     # write the header to output in tecplot format
    output_filename = input_folder+'/'+'%s_max_cfl.tec'%(folder_itself)
    print ('Writing data to %s ...'%(output_filename))
    outfile = open(output_filename, 'w')
    outfile.write('VARIABLES=timestep,t,CFL\nZONE I=%d'%(file_count))
    outfile.close()
    
    # create shared arrays
    create_shared_array('u_dx_max', (file_count))
    # make group and divide the procedure
    step = max(int(file_count / nproc), 1)
    rng = list(range(0,file_count))
    groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

    print ('Reading', len(h5_files), 'files in', len(groups), 'of', step, '.')

    p_list=[]
    for i,g in enumerate(groups):
        pr = multiprocessing.Process(target=compute_local_u_dx, name='Process'+str(i), args=(g,h5_files,mesh_dx,mesh_cells,detailed))
        p_list.append(pr)
    for pr in p_list: pr.start()
    # Wait for all the processes to finish
    for pr in p_list: pr.join()
    gc.collect()

    print (' done.', flush=True)

    u_dx_max = get_shared_var('u_dx_max')

    # Now, all points are done. The output file is to be closed.
    outfile = open(output_filename, 'a+')
    for i in range(file_count):
        fname = h5_files[i][h5_files[i].rfind('/'):]
        ts = int(fname.split("_ts=")[1].split("_")[0])
        t = float(fname.split("_t=")[1].split("_")[0])
        outfile.write('\n%8d %14.5f % 18.12f'%(ts,t,dt*u_dx_max[i]))
    outfile.close()
    print (' done.')

if __name__ == '__main__':
    nargs = len(sys.argv)
    #print('Arguments',nargs, sys.argv)
    # get the number of CPU cores
    if nargs > 1:
        ncore = multiprocessing.cpu_count()
        interval = 1
        period = 0.951
        detailed = 0
        if nargs > 2:
            period = float(sys.argv[2])
        if nargs > 3:
            interval = int(sys.argv[3])
        if nargs > 4:
            ncore = int(sys.argv[4])
        if nargs > 5:
            detailed = int(sys.argv[5])
        print ( 'Performing CFL computation on %d core%s and interval of %d.'%(ncore,'s' if ncore>1 else '',interval) )
        compute_cfl(sys.argv[1], interval, ncore, period, detailed)
