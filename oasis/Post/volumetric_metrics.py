# ///////////////////////////////////////////////////////////////
#  Volumetric Hemodynamic Indices
#  Author: Mehdi Najafi, mnuoft@gmail.com
#  Date: 2018-07-18
#
#  This library is not intended for distributions in any form and
#  distribution of it is not allowed in any form.
# ///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-07-18"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private; to be obtained directly from the author."

import sys, os
# os.environ['OPENBLAS_NUM_THREADS'] = '1'
# os.environ['MKL_NUM_THREADS'] = '1'

import h5py, glob
import numpy
import gc
from numpy.fft import fftfreq, fft, ifft
# from scipy.fftpack import fftfreq, fft, ifft
# from scipy.fftpack import rfftfreq, rfft
import job_utils
from multiprocessing import sharedctypes
import multiprocessing

# use_half_fft_signal = True  # old codes mistake
use_half_fft_signal = False

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
# read h5 files in chunks
def read_h5_files(ids, rng, h5_files):
    if len(ids)==1:
        print ('    reading', len(ids), 'file:', ids, h5_files[ids[0]], flush=True)
    else:
        print ('    reading', len(ids), 'files:', ids, h5_files[ids[0]], ' ... ',  h5_files[ids[-1]], flush=True)
    # get the variables
    u_1 = get_shared_var('u_1')
    u_2 = get_shared_var('u_2')
    u_3 = get_shared_var('u_3')
    u_abs = get_shared_var('u_abs')
    #p = get_shared_var('p')
    # loop over files
    for i in ids:
        hw = h5py.File(h5_files[i], 'r')
        hwu = numpy.asarray(hw['/Solution/u'])
        #hwp = numpy.asarray(hw['/Solution/p'])
        lu = hwu[rng,:]
        u_1[:,i] = lu[:,0]
        u_2[:,i] = lu[:,1]
        u_3[:,i] = lu[:,2]
        u_abs[:,i] = numpy.sqrt(u_1[:,i]*u_1[:,i]+u_2[:,i]*u_2[:,i]+u_3[:,i]*u_3[:,i])
        #p[:,i] = hwp[rng,0]
        hw.close()
################################################################
################################################################
def filter_SPI(U, W_cut):    
    # "withmean":
    # U_fft = fft(U)
    # "withoutmean":
    # U_fft = fft(U-numpy.mean(U, axis=0), axis=1)
    U_fft = fft(U, axis=1)

    # cut off the upper frequency to 25Hz
    U_fft_25Hz = U_fft.copy()
    U_fft[ :, W_cut[0] ] = 0
    U_fft_25Hz[ :, W_cut[1] ] = 0

    # compute the absolute value
    m0 = numpy.absolute(U_fft)
    m1 = numpy.absolute(U_fft_25Hz)

    power_25Hz = numpy.sum ( m1 * m1 , axis=1)
    power_0Hz  = numpy.sum ( m0 * m0 , axis=1)
    spi = power_25Hz

    for i in range(power_0Hz.shape[0]):
        if power_0Hz[i] < 1e-5:
            spi[i] = 0
        else:
            spi[i] /= power_0Hz[i]
    return spi
################################################################
def filter_TKE(U, W_cut):
    # swith to frequency domain
    U_fft = fft(U, axis=1)
    # print ('U_fft.shape', U_fft.shape)
    # Cut off the upper frequency to low_cut (default: 25Hz)
    U_fft[ :, W_cut ] = 0
    # return back to time domain
    return numpy.absolute( ifft(U_fft, axis=1) )
################################################################
def compute_hemo(ids, W_low_cut):
    print ('    working on', len(ids), 'points:', ids[0], ' to ', ids[-1], flush=True)
    # get the variables
    u_1 = get_shared_var('u_1')
    u_2 = get_shared_var('u_2')
    u_3 = get_shared_var('u_3')
    u_abs = get_shared_var('u_abs')
    TAVEL = get_shared_var('TAVEL')
    TKE = get_shared_var('TKE')
    SPI = get_shared_var('SPI')
    SPI_a = get_shared_var('SPI_a')
    SPI_p = get_shared_var('SPI_p')

    # Time Averaged Velocity
    TAVEL[ids] = numpy.mean(u_abs[ids], axis=1)
    
    # TKE
    u1_ms = filter_TKE(u_1[ids], W_low_cut[1])
    u2_ms = filter_TKE(u_2[ids], W_low_cut[1])
    u3_ms = filter_TKE(u_3[ids], W_low_cut[1])
    # print ('u1_ms.shape', u1_ms.shape)
    TKE[ids] = 0.5 * numpy.mean(u1_ms*u1_ms + u2_ms*u2_ms + u3_ms*u3_ms, axis=1)

    # SPI
    u1_spi = filter_SPI(u_1[ids], W_low_cut)
    u2_spi = filter_SPI(u_2[ids], W_low_cut)
    u3_spi = filter_SPI(u_3[ids], W_low_cut)
    SPI[ids] = numpy.sqrt(u1_spi*u1_spi + u2_spi*u2_spi + u3_spi*u3_spi) / numpy.sqrt(3)

    # SPI based on magnitude
    SPI_a[ids] = filter_SPI(u_abs[ids], W_low_cut)

    # SPI based on projection of on time averaged vector
    u1_avg = numpy.mean(u_1[ids], axis=1).reshape( u_1[ids].shape[0], 1)
    u2_avg = numpy.mean(u_2[ids], axis=1).reshape( u_2[ids].shape[0], 1)
    u3_avg = numpy.mean(u_3[ids], axis=1).reshape( u_3[ids].shape[0], 1)
    
    u_projected_on_uav = (u1_avg*u_1[ids] + u2_avg*u_2[ids] + u3_avg*u_3[ids]) #/ TAVEL[ids].reshape(TAVEL[ids].shape[0],1)
    SPI_p[ids] = filter_SPI(u_projected_on_uav, W_low_cut)

################################################################
def volumetric_indices(input_folder, interval, nproc, period, low_cut=25):
    pos = -2 if (input_folder[-1] == '/') else -1
    folder_itself = input_folder.split('/')[pos]
    timesteps = int(folder_itself.split("_ts")[-1].split("_cy")[0])

    mesh_filename, case_folder, case_name = job_utils.get_case_mesh_filename(input_folder)
    if not os.path.exists(mesh_filename):
        print ('No mesh file found: %s \nYou may running this script from an incorrect folder.\n'%mesh_filename)
        exit(1)
    print ('Looking for mesh file:', mesh_filename, ' and loading volume mesh (/Mesh/coordiantes, /Mesh/topology).')
    mesh = h5py.File(mesh_filename, 'r')
    number_of_points = mesh['/Mesh/coordinates'].shape[0]

    print ('Looking inside', input_folder, '...')

    h5_files = glob.glob(input_folder+"/*_up.h5")
    if not h5_files:
        print ('No h5 file found in', input_folder+'. Run the solver and run this script after it finished.\n')
        exit(1)
    # sort the files according to the simulation time
    h5_files = sorted(h5_files, key=lambda N: int(N.split("ts=")[1].split("_up.h5")[0]))
    h5_files = h5_files[0:len(h5_files):interval]
    file_count = len(h5_files)

    print ('   found', file_count, 'up files for', number_of_points, 'points.')

    # determine the time increment based on the number of samples and the corresponding FFT frequencies
    dt = period/(file_count-1)
    W = fftfreq(file_count, d=dt)
    print ('period=',period, '  samples or file count=',file_count, '  dt=', dt, '  frequencies=',W)
    print ('max frequency=',numpy.max(W))
    print ('min frequency=',numpy.min(W))

    # W_low_cut = numpy.where( W < low_cut ) # wrong assumption in old codes
    W_low_cut = numpy.where( numpy.abs(W) == 0 ) + numpy.where( numpy.abs(W) < low_cut )
    print ( 'w_low_cut:', W_low_cut )

    # ! NOTE !
    # Due to memory limitations it may not be feasible to load the whole dataset. Instead it load chunks.
    chunk_size = int( number_of_points / max( int(number_of_points*file_count / 1.05E9), 1 ) )
    point_rng = list(range(0,number_of_points))
    point_groups = [point_rng[i:i+chunk_size] for i  in range(point_rng[0], point_rng[-1]+1, chunk_size)]
    
    print ('Number of point groups:', len(point_groups), flush=True)

     # write the header to output in tecplot format
    output_filename = input_folder+'/'+'%s_volumetric_metrics%s.tec'%(folder_itself, '' if use_half_fft_signal else '_full')
    print ('Writing data to %s ...'%(output_filename))
    coord = numpy.asarray(mesh['Mesh/coordinates'])
    elems = numpy.asarray(mesh['Mesh/topology'])
    outfile = open(output_filename, 'w')
    vars = 'VARIABLES = X,Y,Z,TAVEL,TKE,SPI,SPI_mag,SPI_prj'
    outfile.write(vars + '\nZONE N=%d,E=%d,F=FEPOINT,ET=TETRAHEDRON\n'%(coord.shape[0], elems.shape[0]))
    outfile.close()
    
    # create shared arrays
    num_pt = len(point_groups[0])
    create_shared_array('u_1,u_2,u_3,u_abs,p', (num_pt,file_count))
    create_shared_array('TAVEL,TKE,SPI,SPI_a,SPI_p', (num_pt))
    for k, point_list in enumerate(point_groups):
        num_pt = len(point_list)
        
        print('Working on %d/%d - %d points out of %d points from %d to %d:'%(k+1,len(point_groups),num_pt,number_of_points, point_list[0], point_list[-1]))

        print ('Reading', file_count, 'files into 4 arrays of shape', [num_pt, file_count],' ...')#, end='')

        # make group and divide the procedure
        step = max(int(file_count / nproc), 1)
        rng = list(range(0,file_count))
        groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

        p_list=[]
        for i,g in enumerate(groups):
            pr = multiprocessing.Process(target=read_h5_files, name='Process'+str(i), args=(g,point_list,h5_files))
            p_list.append(pr)
        for pr in p_list: pr.start()
        # Wait for all the processes to finish
        for pr in p_list: pr.join()
        gc.collect()

        print (' done.', flush=True)

        # now compute desired variables
        print ('Now, computing TAVEL, TKE, SPI and so on for %d/%d - %d points out of %d points from %d to %d:'%(k+1,len(point_groups),num_pt,number_of_points, point_list[0], point_list[-1]))


        # make group and divide the procedure
        ##step = max(int(len(point_list) / nproc), 1)
        ##rng = point_list
        ##groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

        # p_list=[]
        # for i,g in enumerate(groups):
        #     pr = multiprocessing.Process(target=compute_hemo, name='Process'+str(i), args=(g,u_1,u_2,u_3,u_abs,TransWSS,TAWSS,OSI,SPI,W,file_count))
        #     p_list.append(pr)
        # for pr in p_list: pr.start()
        # # Wait for all the processes to finish
        # for pr in p_list: pr.join()
        ##for i,g in enumerate(point_list):
        ##    compute_hemo(g,u_1,u_2,u_3,u_abs,TAVEL,TKE,SPI,SPI_a,SPI_p, W_low_cut,file_count)

        compute_hemo(list(range(0,num_pt)), W_low_cut)

        # gc.collect()
        print ('done.', flush=True)

        TAVEL = get_shared_var('TAVEL')
        TKE = get_shared_var('TKE')
        SPI = get_shared_var('SPI')
        SPI_a = get_shared_var('SPI_a')
        SPI_p = get_shared_var('SPI_p')

        # write the results to output in tecplot format
        print ('Writing data to %s ...'%(output_filename), end='')
        outfile = open(output_filename, 'a+')
        for i in range(num_pt):
            xyz = coord[i+point_list[0]]
            outfile.write('% 16.12f % 16.12f % 16.12f  % 16.12f % 16.12f % 16.12f % 16.12f % 16.12f\n'%
                            (xyz[0],xyz[1],xyz[2],TAVEL[i],TKE[i],SPI[i],SPI_a[i],SPI_p[i]))
        outfile.close()
        print (' done.')
        gc.collect()

    # Now, all points are done. The output file is to be closed.
    outfile = open(output_filename, 'a+')
    for i in range(elems.shape[0]):
        c = elems[i]
        outfile.write('\n%d %d %d %d'%(c[0]+1,c[1]+1,c[2]+1,c[3]+1))
    outfile.close()
    mesh.close()
    print (' done.')

if __name__ == '__main__':
    nargs = len(sys.argv)
    #print('Arguments',nargs, sys.argv)
    # get the number of CPU cores
    if nargs > 1:
        ncore = multiprocessing.cpu_count()
        interval = 1
        period = 0.951
        if nargs > 2:
            period = float(sys.argv[2])
        if nargs > 3:
            interval = int(sys.argv[3])
        if nargs > 4:
            ncore = int(sys.argv[4])
        print ( 'Performing volumetric hemodynamics computation on %d core%s and interval of %d.'%(ncore,'s' if ncore>1 else '',interval) )
        volumetric_indices(sys.argv[1], interval, ncore, period)
