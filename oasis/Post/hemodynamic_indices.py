# ///////////////////////////////////////////////////////////////
#  Hemodynamic Indices Script
#  Author: Mehdi Najafi, mnuoft@gmail.com
#  Date: 2018-07-16, updated: 2020-03-23
#
#  This library is not intended for distributions in any form and
#  distribution of it is not allowed in any form.
# ///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-07-16"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private; to be obtained directly from the author."

import sys, os
os.environ['OPENBLAS_NUM_THREADS'] = '1'
os.environ['MKL_NUM_THREADS'] = '1'

import h5py, glob
import numpy
import gc
# from scipy.fftpack import fftfreq, fft, ifft
# from scipy.fftpack import rfftfreq, rfft
from numpy.fft import fftfreq, fft, ifft
import job_utils
import multiprocessing
from multiprocessing import sharedctypes
import xdmf_tools

# to render older spi calculations, change to True, otherwise use False
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
def read_h5_files(ids, h5_files, _wss_1, _wss_2, _wss_3, _wss_abs):
    wss_1 = get_array(_wss_1)
    wss_2 = get_array(_wss_2)
    wss_3 = get_array(_wss_3)
    wss_abs = get_array(_wss_abs)
    if len(ids)==1:
        print ('    reading', len(ids), 'file:', ids, h5_files[ids[0]])
    else:
        print ('    reading', len(ids), 'files:', ids, h5_files[ids[0]], ' ... ',  h5_files[ids[-1]])
    for i in ids:
        hw = h5py.File(h5_files[i], 'r')
        ls = numpy.asarray(hw['Computed/wss'])
        la = numpy.asarray(hw['Computed/wss_abs'])
        for j in range(ls.shape[0]):
            wss_1[j][i] = ls[j][0]
            wss_2[j][i] = ls[j][1]
            wss_3[j][i] = ls[j][2]
            wss_abs[j][i] = la[j]
        hw.close()
################################################################
# read all ASCII files
def read_ascii_files(ids, ascii_files, _wss_1, _wss_2, _wss_3, _wss_abs):
    wss_1 = get_array(_wss_1)
    wss_2 = get_array(_wss_2)
    wss_3 = get_array(_wss_3)
    wss_abs = get_array(_wss_abs)
    if len(ids)==1:
        print ('    reading', len(ids), 'file:', ids, ascii_files[ids[0]])
    else:
        print ('    reading', len(ids), 'files:', ids, ascii_files[ids[0]], ' ... ',  ascii_files[ids[-1]])

    for i in ids:
        fw = open(ascii_files[i], 'r')
        j = 0
        for line in fw.readlines():
            v = line.split()
            wss_1[j][i] = float(v[0])
            wss_2[j][i] = float(v[1])
            wss_3[j][i] = float(v[2])
            wss_abs[j][i] = numpy.sqrt( wss_1[j][i]*wss_1[j][i] + wss_2[j][i]*wss_2[j][i] + wss_3[j][i]*wss_3[j][i] )
            j = j+1
        fw.close()
################################################################
def filter_SPI_bands(U,W_low_cut,tag, bands):
    if tag=="withmean":
        U_fft = fft(U)
    else:
        U_fft = fft(U-numpy.mean(U))

    Power_xHz = []
    for i in range(len(bands)): Power_xHz.append(0.0)

    # filter any amplitude corresponding frequency equal to 0Hz
    U_fft_0Hz = U_fft.copy()
    U_fft_0Hz[W_low_cut[0]] = 0

    # Compute the absolute value
    Power_0Hz  = numpy.sum ( numpy.power(numpy.absolute(U_fft_0Hz),2) )
    if Power_0Hz<1e-5:
        return Power_xHz

    for i, band in enumerate(bands):
        U_fft_xHz = U_fft.copy()
        # filter any amplitude corresponding to those frequencies
        U_fft_xHz [W_low_cut[2+i]] = 0
        pw = numpy.sum ( numpy.power(numpy.absolute(U_fft_xHz),2) )
        Power_xHz[i] = pw / Power_0Hz

    return Power_xHz
################################################################
def filter_SPI(U, W_low_cut, tag):
    #for HI
    if tag=="withmean":
        U_fft = fft(U)
    else:
        U_fft = fft(U-numpy.mean(U))
    # filter any amplitude corresponding frequency equal to 0Hz
    U_fft[W_low_cut[0]] = 0
    # filter any amplitude corresponding frequency lower to 25Hz
    U_fft_25Hz = U_fft.copy()
    U_fft_25Hz[W_low_cut[1]] = 0
    #Compute the absolute value
    Power_25Hz = numpy.sum ( numpy.power( numpy.absolute(U_fft_25Hz),2))
    Power_0Hz  = numpy.sum ( numpy.power( numpy.absolute(U_fft     ),2))
    if Power_0Hz < 1e-5:
        return 0
    return Power_25Hz/Power_0Hz
################################################################
def compute_hemo(ids,_wss_1,_wss_2,_wss_3,_wss_abs,_TransWSS,_TAWSS,_OSI,_SPI,_SPI_a,_SPI_p, _SPIw, bands, W_low_cut, normals,number_of_steps):
    wss_1 = get_array(_wss_1)
    wss_2 = get_array(_wss_2)
    wss_3 = get_array(_wss_3)
    wss_abs = get_array(_wss_abs)

    TransWSS = get_array(_TransWSS)
    TAWSS = get_array(_TAWSS)
    OSI = get_array(_OSI)
    SPI = get_array(_SPI)
    SPIw = get_array(_SPIw)
    SPI_a = get_array(_SPI_a)
    SPI_p = get_array(_SPI_p)

    print ('    working on', len(ids), 'points:', ids)

    for j in ids:
        TAWSS[j] = numpy.mean(wss_abs[j])
        TAWSS0 = numpy.mean(wss_1[j])
        TAWSS1 = numpy.mean(wss_2[j])
        TAWSS2 = numpy.mean(wss_3[j])
        TAWSSVM = numpy.sqrt(TAWSS0*TAWSS0+TAWSS1*TAWSS1+TAWSS2*TAWSS2)
        #print ('mean values:',TAWSS0,TAWSS1,TAWSS2,TAWSSVM, normals[j])
        #Normalize the WSS vectors
        TAWSSV = numpy.array([TAWSS0,TAWSS1,TAWSS2])/TAWSSVM
        TransWSS_term1 = numpy.cross(normals[j],TAWSSV)
        TransWSS_temp = 0
         #Loop over all steps
        for k in range(number_of_steps): 
            TransWSS_temp += numpy.absolute(numpy.dot([wss_1[j][k],wss_2[j][k],wss_3[j][k]],TransWSS_term1))
        #Compute TransWSS
        TransWSS[j] = TransWSS_temp/number_of_steps
        OSI[j] = 0.5*(1-TAWSSVM/TAWSS[j])
        a1 = filter_SPI(wss_1[j],W_low_cut,"withoutmean")
        a2 = filter_SPI(wss_2[j],W_low_cut,"withoutmean")
        a3 = filter_SPI(wss_3[j],W_low_cut,"withoutmean")
        SPI[j] = numpy.sqrt(a1*a1+a2*a2+a3*a3)/numpy.sqrt(3.0)
        SPI_a[j] = filter_SPI(wss_abs[j],W_low_cut,"withoutmean")

        wss_projected_on_tawss = ( TAWSS0*wss_1[j] + TAWSS1*wss_2[j] + TAWSS2*wss_3[j] ) / TAWSSVM

        SPI_p[j] = filter_SPI(wss_projected_on_tawss, W_low_cut,"withoutmean") 

        spiw = filter_SPI_bands(wss_abs[j],W_low_cut,"withoutmean", bands)
        for k,bnd in enumerate(bands):
            SPIw[j][k] = spiw[k]
    #print (SPI)
    #print (numpy.max(SPI))

################################################################
def hemodynamics(input_folder, interval, nproc, period, ascii=False):
    pos = -2 if (input_folder[-1] == '/') else -1
    folder_itself = input_folder.split('/')[pos]
    timesteps = int(folder_itself.split("_ts")[-1].split("_cy")[0])

    mesh_filename, case_folder, case_name = job_utils.get_case_mesh_filename(input_folder)
    print ('Looking for mesh file:', mesh_filename, ' and loading wall normal vectors (/Mesh/Wall/normal).')
    mesh = h5py.File(mesh_filename, 'r')
    normals = numpy.asarray(mesh['/Mesh/Wall/normal'])
    number_of_points = normals.shape[0]

    print ('Looking at', input_folder+'/wss_files', '...')

    if ascii:
        h5_files = glob.glob(input_folder+"/wss_files/*.dat")
    else:
        h5_files = glob.glob(input_folder+"/wss_files/*_wss.h5")

    if not h5_files:
        print ('No h5 file found in', input_folder+'/wss_files. Re-compute wss and run this script again.\n')
        exit(1)
    
    if ascii:
        h5_files = sorted(h5_files,key=lambda TS: int(TS.split("_000000_wss.")[0].split("=")[1]))
        print (h5_files)
    else:
        h5_files = sorted(h5_files, key=lambda N: int(N.split("ts=")[1].split("_wss.h5")[0]))

    h5_files = h5_files[0:len(h5_files):interval]

    file_count = len(h5_files)

    print ('   found', file_count, 'wss files for', number_of_points, 'points.')


    print ('Allocating %5.2f'%(number_of_points*file_count*8/(1024*1024))+'MB of memory ...', end='')

    wss_1 = create_shared_var([number_of_points, file_count])
    wss_2 = create_shared_var([number_of_points, file_count])
    wss_3 = create_shared_var([number_of_points, file_count])
    wss_abs = create_shared_var([number_of_points, file_count])

    print ('done.', flush=True)

    print ('Reading', len(h5_files), 'files into 4 arrays of shape', [number_of_points, file_count],' ...')#, end='')

    # make group and divide the procedure
    step = max(int(file_count / nproc), 1)
    rng = list(range(0,file_count))
    groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

    p_list=[]
    if ascii:
        for i,g in enumerate(groups):
            p = multiprocessing.Process(target=read_ascii_files, name='Process'+str(i), args=(g,h5_files,wss_1,wss_2,wss_3,wss_abs))
            p_list.append(p)
    else:
        for i,g in enumerate(groups):
            p = multiprocessing.Process(target=read_h5_files, name='Process'+str(i), args=(g,h5_files,wss_1,wss_2,wss_3,wss_abs))
            p_list.append(p)
    for p in p_list: p.start()
    # Wait for all the processes to finish
    for p in p_list: p.join()
    gc.collect()

    # wss_1 = get_array(wss_1)
    # wss_2 = get_array(wss_2)
    # wss_3 = get_array(wss_3)
    # wss_abs = get_array(wss_abs)

    print (' done.', flush=True)
 

    # time_line = numpy.linspace(0,0.951,file_count)
    # dt = time_line[1] - time_line[0]
    dt = period/float(file_count-1)
    W = fftfreq(file_count, d=dt)
    
    if use_half_fft_signal:
        # wrong assumption in old codes
        W_low_cut = numpy.where( W < 0 ) + numpy.where( W < 25.0 )
    else:
        W_low_cut = numpy.where( numpy.abs(W) == 0 ) + numpy.where( numpy.abs(W) < 25.0 )

    print ('period=',period, '  file_count=',file_count, '  dt=', dt, '  frequencies=',W)
    print ('max frequency=',numpy.max(W))
    print ('min frequency=',numpy.min(W))


    # bands
    width = 25
    high = 500
    bands = [[i,i+width] for i  in range(0, high, width)]
    bands.append([high, 1E10])
    bands.append([50,100])
    bands.append([100,200])
    bands.append([200,400])
    bands.append([400,1E10])
    nbands = len(bands)

    print ('Looking at', nbands,'frequency bands:', bands)

    if use_half_fft_signal:
        for bnd in bands:
            # wrong assumption in old codes
            W_low_cut += numpy.where( W < bnd[0] and W > bnd[1] )
    else:
        for bnd in bands:
            W_low_cut += numpy.where( numpy.abs(W) < bnd[0] ) and numpy.where( numpy.abs(W) > bnd[1] )

    numpy.set_printoptions(threshold=sys.maxsize)
    print ( 'w_low_cut:', W_low_cut )


    # now compute desired variables
    print ('Now, computing  TAWSS, OSI, SPI and ...')

    TransWSS = create_shared_var([number_of_points])
    TAWSS = create_shared_var([number_of_points])
    OSI = create_shared_var([number_of_points])
    SPI = create_shared_var([number_of_points])
    SPI_a = create_shared_var([number_of_points])
    SPI_p = create_shared_var([number_of_points])

    SPIw = create_shared_var([number_of_points, nbands])


    # make group and divide the procedure
    step = max(int(number_of_points / nproc), 1)
    rng = list(range(0,number_of_points))
    groups = [rng[i:i+step] for i  in range(rng[0], rng[-1]+1, step)]

    # p_list=[]
    # for i,g in enumerate(groups):
    #     p = multiprocessing.Process(target=compute_hemo, name='Process'+str(i), args=(g,wss_1,wss_2,wss_3,wss_abs,TransWSS,TAWSS,OSI,SPI,SPIw,bands,W_low_cut,normals,file_count))
    #     p_list.append(p)
    # for p in p_list: p.start()
    # # Wait for all the processes to finish
    # for p in p_list: p.join()
    for i,g in enumerate(groups):
        compute_hemo(g,wss_1,wss_2,wss_3,wss_abs,TransWSS,TAWSS,OSI,SPI, SPI_a,SPI_p, SPIw, bands, W_low_cut, normals,file_count)
    gc.collect()

    TransWSS = get_array(TransWSS)
    TAWSS = get_array(TAWSS)
    OSI = get_array(OSI)
    SPI = get_array(SPI)
    SPIw = get_array(SPIw)
    SPI_a = get_array(SPI_a)
    SPI_p = get_array(SPI_p)

    nTransWSS = TransWSS/TAWSS
    #print ('TAWSS:\n',TAWSS, '\n\nOSI:\n',OSI)
    RRT = 1.0 / ( (1.0 - 2.0 * OSI) * TAWSS )
    print ('done.', flush=True)

    # write the results to output in tecplot format
    output_filename = input_folder+'/'+'%s_hemodynamics_w%s.tec'%(folder_itself, '' if use_half_fft_signal else '_full')
    print ('Writing data to %s ...'%(output_filename), end='')
    coord = numpy.asarray(mesh['Mesh/Wall/coordinates'])
    elems = numpy.asarray(mesh['Mesh/Wall/topology'])
    outfile = open(output_filename, 'w')
    vars = 'VARIABLES = X,Y,Z,TAWSS,OSI,SPI,SPI_prj,SPI_mehdi,TransWSS,nTransWSS,RRT'
    for bnd in bands: vars += ',SPI_%d_%d'%(bnd[0],bnd[1])
    outfile.write(vars + '\nZONE N=%d,E=%d,F=FEPOINT,ET=TRIANGLE\n'%(coord.shape[0], elems.shape[0]))
    # number_of_points == coord.shape[0]
    for i in range(number_of_points):
        outfile.write('% 16.12f % 16.12f % 16.12f  % 16.12f % 16.12f % 16.12f % 16.12f % 16.12f % 16.12f % 16.12f % 16.12f'%
                        (coord[i,0],coord[i,1],coord[i,2],TAWSS[i],OSI[i],SPI_a[i],SPI_p[i],SPI[i],TransWSS[i],nTransWSS[i],RRT[i]))
        for j in range(nbands):
            outfile.write(' % 16.12f'%SPIw[i][j])
        outfile.write('\n')
    for i in range(elems.shape[0]):
        c = elems[i]
        outfile.write('\n%d %d %d'%(c[0]+1,c[1]+1,c[2]+1))
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
        print ( 'Performing hemodynamics computation on %d core%s and interval of %d.'%(ncore,'s' if ncore>1 else '',interval) )
        hemodynamics(sys.argv[1], interval, ncore, period)
        xdmf_tools.make_xdmf_files(sys.argv[1])

