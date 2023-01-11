#!/usr/bin/python

################################################################################
##   Project:   Data Format Conversion Tools
##   Date:      2010/10/02 10:15:34
##   Version:   1.0
##   Author:    Mehdi Najafi, (mnajafi@sharif.edu).
##              All rights reserved.
################################################################################

from __future__ import print_function

import numpy, h5py, vtk
import glob, sys, os, gzip
import xml.etree.ElementTree as ET

def read_vtu_dataset(filename, array_name=None):
    reader = vtk.vtkXMLUnstructuredGridReader()
    reader.SetFileName(filename)
    reader.Update()
    polydata = reader.GetOutput()
    return vtk.util.numpy_support.vtk_to_numpy( polydata.GetPointData().GetArray(0) )

def read_xml(filename, cells, num_points = None):
    # Use iterparse() to avoid loading the entire file via parse(). iterparse()
    # allows to discard elements (via clear()) after they have been processed.
    # See <https://stackoverflow.com/a/326541/353337>.
    for event, elem in ET.iterparse(filename, events=("start", "end")):
        if event == "end":
            continue

        if elem.tag == "dolfin":
            # Don't be too strict with the assertion. Some mesh files don't have the
            # proper tags.
            # assert elem.attrib['nsmap'] \
            #     == '{\'dolfin\': \'https://fenicsproject.org/\'}'
            pass
        elif elem.tag == "function_data":
            data_size = int(elem.attrib["size"])
            if num_points is None: num_points = data_size
            data = numpy.zeros((num_points,), dtype=numpy.float64)
        elif elem.tag == "dof":
            lid = int(elem.attrib['cell_dof_index'])
            if lid < 4:
                data[cells[int(elem.attrib['cell_index']),lid]] = eval(elem.attrib["value"])
        else:
            print("Unknown entry %s. Ignoring.", elem.tag)

        elem.clear()

    return data

def read_xml_gz(filename, cells, num_points=None, tmp_folder='./'):
    print ('- loading [xml] ', filename)
    extracted_file = False
    if filename[-3:]=='.gz':
        ofilename = os.path.join( tmp_folder , filename[:-3].replace('/','__') )
        with gzip.open(filename, 'rb') as f:
            open(ofilename,'wb').write(f.read())
        filename = ofilename
        extracted_file = True
    data = read_xml(filename, cells, num_points)
    if extracted_file: os.remove(filename)
    return data

def read_npz(filename):
    print ('- loading [npz] ', filename)
    return numpy.load(filename)

def convert_xml_npz(u0_xml_file,u1_xml_file,u2_xml_file,u_npz_file, cells, num_points, redo=False):
    # make the h5 name
    if u_npz_file:
        h5_filename = u_npz_file.replace('_u_','_')
        h5_filename = h5_filename.replace('.npz','_up.h5')
    else:
        h5_filename = u0_xml_file.replace('_u0_','_')
        h5_filename = h5_filename.replace('.xml.gz','_up.h5')
    h5_filename = h5_filename.replace(input_folder_name, output_folder_name)
    
    # check if the h5 file exists
    if not redo and os.path.exists(h5_filename):
        return 1

    # read velocity
    if u0_xml_file:
        u0 = read_xml_gz(u0_xml_file, cells,num_points, output_folder_name)
        u1 = read_xml_gz(u1_xml_file, cells,num_points, output_folder_name)
        u2 = read_xml_gz(u2_xml_file, cells,num_points, output_folder_name)
        u = numpy.column_stack((u0,u1,u2))
    else:
        filename = u_npz_file.replace('.npz','_000000.vtu')
        pos = filename.rfind('/')
        filename = filename[:pos] + '/vtu_files'+filename[pos:]
        u = read_vtu_dataset(filename)

    # find and read pressure
    if u0_xml_file:
        p_xml_filename = u0_xml_file.replace('_u0_','_p_')
    else:
        p_xml_filename = u_npz_file.replace('_u_','_p_')
        p_xml_filename = p_xml_filename.replace('.npz','.xml.gz')
    p_npz_filename = p_xml_filename.replace('.xml.gz','.npz')
    if os.path.exists(p_xml_filename):
        p = read_xml_gz(p_xml_filename, cells,num_points, output_folder_name)
    elif os.path.exists(p_npz_filename):
        p = read_npz(p_npz_filename)['p']
    else:
        print ('No pressure data found:',p_xml_filename, ' or ',p_npz_filename )
        return 0
    
    # write u and p to h5
    print ('- writing [h5] ', h5_filename)
    # print(u.shape, p.shape)
    hf = h5py.File(h5_filename, 'w')
    hf.create_dataset('/Solution/u', dtype=numpy.float64, data=u, compression="gzip")
    hf.create_dataset('/Solution/p', dtype=numpy.float64, data=p, compression="gzip")
    hf.close()


if __name__ == '__main__':
    print ('~'*79)
    print (' XML/NPZ/H5 Data Format Conversion Tool, by Mehdi Najafi, Copyright 2010-2020. ')
    print ('~'*79)
    print (' Usage:\n  python dolfin_data_to_h5.py input_folder output_folder mesh_file_name')
    print ('~'*79)

    redo = False
    if len(sys.argv) > 1:
        input_folder_name = sys.argv[1]
        output_folder_name = input_folder_name
        if len(sys.argv) > 2:
            output_folder_name = sys.argv[2]
            if len(sys.argv) > 3:
                mesh_h5_filename = sys.argv[3]
                if len(sys.argv) > 4:
                    redo = True
    else:
        print ('<!> No input folder name given.')
        exit(1)

    # retrive the mesh cells and point ids
    print('- reading mesh topology:',mesh_h5_filename)
    hf = h5py.File(mesh_h5_filename, 'r')
    cells = numpy.asarray( hf['Mesh/topology'] )
    num_points = numpy.asarray( hf['/Mesh/coordinates'] ).shape[0]
    print ('Number of Mesh Points:', num_points, '\nNumber of Mesh Cells:', cells.shape[0])
    hf.close()

    # retrive the list of files in the given folder
    u0_xml_files = glob.glob(input_folder_name + "/*_u0_*.xml.gz")
    u1_xml_files = glob.glob(input_folder_name + "/*_u1_*.xml.gz")
    u2_xml_files = glob.glob(input_folder_name + "/*_u2_*.xml.gz")
    p_xml_files  = glob.glob(input_folder_name + "/*_p_*.xml.gz" )
    u_npz_files = glob.glob(input_folder_name + "/*_u_*.npz")
    p_npz_files = glob.glob(input_folder_name + "/*_p_*.npz")

    # sort the file names 
    u0_xml_files = sorted(u0_xml_files,key=lambda N:int(N.split(".xml.gz")[0].split("=")[1]))
    u1_xml_files = sorted(u1_xml_files,key=lambda N:int(N.split(".xml.gz")[0].split("=")[1]))
    u2_xml_files = sorted(u2_xml_files,key=lambda N:int(N.split(".xml.gz")[0].split("=")[1]))
    p_xml_files  = sorted(p_xml_files, key=lambda N:int(N.split(".xml.gz")[0].split("=")[1]))
    u_npz_files = sorted(u_npz_files,key=lambda N:int(N.split(".npz")[0].split("=")[1]))
    p_npz_files = sorted(p_npz_files,key=lambda N:int(N.split(".npz")[0].split("=")[1]))

    print ('Number of [xml] files:  ','u0:',len(u0_xml_files), 'u1:', len(u1_xml_files), 'u2:', len(u2_xml_files), 'p:', len(p_xml_files))
    print ('Number of [npz] files:  ','u:',len(u_npz_files),'p:',len(p_npz_files))

    # check number of files
    if len(u0_xml_files) != len(u1_xml_files) or len(u1_xml_files) != len(u2_xml_files) or len(u0_xml_files) != len(u2_xml_files):
        print ('<!> Incomplete data files for velocity.')
        exit(1)

    if len(u_npz_files) + len(u0_xml_files) != len(p_npz_files) + len(p_xml_files):
        print ('<!> Incomplete data files for velocity and pressure.')
        exit(1)

    file_count_xml = len(u0_xml_files)
    file_count_npz = len(u_npz_files)
    
    for i in range(file_count_xml):
        convert_xml_npz(u0_xml_files[i],u1_xml_files[i],u2_xml_files[i],None,cells,num_points,redo)

    for i in range(file_count_npz):
        convert_xml_npz(None,None,None,u_npz_files[i],cells,num_points,redo)


# Bullshit stuff to write down
# sampleData/p302/results/art_p302_edited_cl_pr_I2_FC_MCA_10_Q285_Per951_Newt370_ts9600_cy2_uO1/p302_edited_cl_pr.h5