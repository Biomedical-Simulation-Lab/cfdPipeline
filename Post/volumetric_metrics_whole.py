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
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

import h5py, glob
import numpy
import gc
from scipy.fftpack import fftfreq, fft, ifft
from scipy.fftpack import rfftfreq, rfft
#from numpy.fft import fftfreq, fft, ifft

import multiprocessing
from multiprocessing import sharedctypes

# use_half_fft_signal = True  # old codes mistake
use_half_fft_signal = False

################################################################
def create_shared_var(size, dtype=numpy.float64):
    a = numpy.ctypeslib.as_ctypes( numpy.zeros(size, dtype=dtype) )
    return sharedctypes.Array(a._type_, a,  lock=False)

def get_shared_var(a):
    return sharedctypes.Array(a._type_, a,  lock=False)

def get_array(shared_var):
    return numpy.ctypeslib.as_array(shared_var)

################################################################
# read all h5 files
def read_h5_files(ids, h5_files, _u_1, _u_2, _u_3, _u_abs, _p):
    u_1 = get_array(_u_1)
    u_2 = get_array(_u_2)
    u_3 = get_array(_u_3)
    u_abs = get_array(_u_abs)
    p = get_array(_p)
    if len(ids)==1:
        print ('    reading', len(ids), 'file:', ids, h5_files[ids[0]], flush=True)
    else:
        print ('    reading', len(ids), 'files:', ids, h5_files[ids[0]], ' ... ',  h5_files[ids[-1]], flush=True)
    for i in ids:
        hw = h5py.File(h5_files[i], 'r')
        lu = hw['/Solution/u']
        u_1[:,i] = lu[:,0]
        u_2[:,i] = lu[:,1]
        u_3[:,i] = lu[:,2]
        u_abs[:,i] = numpy.sqrt(u_1[:,i]*u_1[:,i]+u_2[:,i]*u_2[:,i]+u_3[:,i]*u_3[:,i])
        p[:,i] = hw['/Solution/p'][:,0]
        hw.close()
################################################################
################################################################
def filter_SPI(U, W_cut):    
    # "withmean":
    # U_fft = fft(U)
    # "withoutmean":
    U_fft = fft(U-numpy.mean(U, axis=0), axis=1)
    
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
def compute_hemo(ids,_u_1,_u_2,_u_3,_u_abs,_TAVEL,_TKE,_SPI, W_low_cut, number_of_steps):
    u_1 = get_array(_u_1)
    u_2 = get_array(_u_2)
    u_3 = get_array(_u_3)
    u_abs = get_array(_u_abs)

    TAVEL = get_array(_TAVEL)
    TKE = get_array(_TKE)
    SPI = get_array(_SPI)

    print ('    working on', len(ids), 'points:', ids, flush=True)

    # print (u_abs.shape)
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
    SPI[ids] = numpy.sqrt(u1_spi*u1_spi + u2_spi*u2_spi + u3_spi*u3_spi)

    # SPI[ids] = filter_SPI(u_abs[ids], W_low_cut)

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

    h5_files = glob.glob(input_folder+"*_up.h5")
    if not h5_files:
        print ('No h5 file found in', input_folder+'. Run the solver and run this script after it finished.\n')
        exit(1)
    # sort the files according to the simulation time
    h5_files = sorted(h5_files, key=lambda N: int(N.split("ts=")[1].split("_up.h5")[0]))

    h5_files = h5_files[0:len(h5_files):interval]

    file_count = len(h5_files)

    print ('   found', file_count, 'up files for', number_of_points, 'points.')


    print ('Allocating %5.2f'%(number_of_points*file_count*10/(1024*1024))+'MB of memory ...', end='', flush=True)

    u_1 = create_shared_var([number_of_points, file_count])
    u_2 = create_shared_var([number_of_points, file_count])
    u_3 = create_shared_var([number_of_points, file_count])
    u_abs = create_shared_var([number_of_points, file_count])
    p = create_shared_var([number_of_points, file_count])

    print ('done.', flush=True)

    print ('Reading', len(h5_files), 'files into 4 arrays of shape', [number_of_points, file_count],' ...')#, end='')

    # make group and divide the procedure
    step = max(int(file_count / nproc), 1)
    rng = list(range(0,file_count))
    groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

    p_list=[]
    for i,g in enumerate(groups):
        pr = multiprocessing.Process(target=read_h5_files, name='Process'+str(i), args=(g,h5_files,u_1,u_2,u_3,u_abs,p))
        p_list.append(pr)
    for pr in p_list: pr.start()
    # Wait for all the processes to finish
    for pr in p_list: pr.join()
    gc.collect()

    # u_1 = get_array(u_1)
    # u_2 = get_array(u_2)
    # u_3 = get_array(u_3)
    # u_abs = get_array(u_abs)

    print (' done.', flush=True)

    # now compute desired variables
    print ('Now, computing  TAVEL, TKE, SPI and ...')

    TAVEL = create_shared_var([number_of_points])
    TKE = create_shared_var([number_of_points])
    SPI = create_shared_var([number_of_points])

    #time_line = numpy.linspace(0,0.951,file_count)
    #dt = time_line[1] - time_line[0]
    dt = period/(file_count-1)
    W = fftfreq(file_count, d=dt)

    # W_low_cut = numpy.where( W < low_cut ) # wrong assumption in old codes
    W_low_cut = numpy.where( numpy.abs(W) == 0 ) + numpy.where( numpy.abs(W) < low_cut )

    # make group and divide the procedure
    step = max(int(number_of_points / nproc), 1)
    rng = list(range(0,number_of_points))
    groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

    # p_list=[]
    # for i,g in enumerate(groups):
    #     pr = multiprocessing.Process(target=compute_hemo, name='Process'+str(i), args=(g,u_1,u_2,u_3,u_abs,TransWSS,TAWSS,OSI,SPI,W,file_count))
    #     p_list.append(pr)
    # for pr in p_list: pr.start()
    # # Wait for all the processes to finish
    # for pr in p_list: pr.join()
    for i,g in enumerate(groups):
        compute_hemo(g,u_1,u_2,u_3,u_abs,TAVEL,TKE,SPI, W_low_cut,file_count)
    gc.collect()

    TAVEL = get_array(TAVEL)
    TKE = get_array(TKE)
    SPI = get_array(SPI)

    print ('done.', flush=True)

     # write the results to output in tecplot format
    output_filename = input_folder+'/'+'%s_volumetric_metrics%s.tec'%(folder_itself, '' if use_half_fft_signal else '_full')
    print ('Writing data to %s ...'%(output_filename), end='')
    coord = mesh['Mesh/coordinates']
    elems = mesh['Mesh/topology']
    outfile = open(output_filename, 'w')
    vars = 'VARIABLES = X,Y,Z,TAVEL,TKE,SPI'
    outfile.write(vars + '\nZONE N=%d,E=%d,F=FEPOINT,ET=TETRAHEDRON\n'%(coord.shape[0], elems.shape[0]))
    # number_of_points == coord.shape[0]
    for i in range(number_of_points):
        x = coord[i][0]; y = coord[i][1]; z = coord[i][2];
        outfile.write('% 16.12f % 16.12f % 16.12f  % 16.12f % 16.12f % 16.12f\n'%
                        (x,y,z,TAVEL[i],TKE[i],SPI[i]))
    for c in elems:
        outfile.write('%d %d %d %d\n'%(c[0]+1,c[1]+1,c[2]+1,c[3]+1))
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
        print ( 'Performing volumetric indices computation on %d core%s and interval of %d.'%(ncore,'s' if ncore>1 else '',interval) )
        volumetric_indices(sys.argv[1], interval, ncore, period)

