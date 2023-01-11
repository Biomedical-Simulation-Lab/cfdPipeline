#!/bin/python
# ///////////////////////////////////////////////////////////////
#  XDMF Tools
#  Author: Mehdi Najafi, mnajafi@sharif.edu
#  Date: 2013-10-21
#
#  This library is not intended for distributions in any form and
#  distribution of it is not allowed in any form.
# ///////////////////////////////////////////////////////////////

from __future__ import print_function
import os, glob, sys, h5py
import job_utils

# ---------------------------------------------------------------------
step_name_leading_zeros = 0
time_name_leading_zeros = 0

def get_timestep_time_string(file_name):
    step_stamp = file_name.split('_ts=')[1].split('_')[0][step_name_leading_zeros:]
    time_stamp = file_name.split('_t=')[1].split('_')[0][time_name_leading_zeros:]
    return step_stamp, time_stamp
# ---------------------------------------------------------------------
def get_mesh_filename(input_folder, move=True):
    mesh_filename, case_folder, case_name = job_utils.get_case_mesh_filename(input_folder)
    if move:
      os.system('mv %s/%s.h5 %s/%s_simple.h5; cp %s %s/'%(input_folder,case_name,input_folder,case_name,mesh_filename,input_folder))
    return case_name + '.h5', mesh_filename

# ---------------------------------------------------------------------

header = """
<?xml version="1.0"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>
    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
"""

ending = """
    </Grid>
</Domain>"""

# ---------------------------------------------------------------------
def make_xdmf_files(input_folder, interval=1):

    if input_folder[-1] == '/': input_folder = input_folder[:-1]
    
    mesh_filename, mesh_pathfilename = get_mesh_filename(input_folder)

    # check the mesh file and get the number of points and elements
    hw = h5py.File(mesh_pathfilename, 'r')
    num_points = hw['/Mesh/coordinates'].shape[0]
    num_elems = hw['/Mesh/topology'].shape[0]
    num_wall_points = hw['/Mesh/Wall/coordinates'].shape[0]
    num_wall_elems = hw['/Mesh/Wall/topology'].shape[0]
    hw.close()

    # prepare the headers
    header_surface_mesh = """
        <Grid Name="mesh" GridType="Uniform">
            <Topology NumberOfElements="%d" TopologyType="Triangle" NodesPerElement="3">
                <DataItem Dimensions="%d 3" NumberType="UInt" Format="HDF">%s:/Mesh/Wall/topology</DataItem>
            </Topology>
            <Geometry GeometryType="XYZ">
                <DataItem Dimensions="%d 3" Format="HDF">%s:/Mesh/Wall/coordinates</DataItem>
            </Geometry>
        </Grid>
    """%(num_wall_elems, num_wall_elems, mesh_filename, num_wall_points, mesh_filename)

    header_volume_mesh = """
        <Grid Name="mesh" GridType="Uniform">
            <Topology NumberOfElements="%d" TopologyType="Tetrahedron" NodesPerElement="4">
                <DataItem Dimensions="%d 4" NumberType="UInt" Format="HDF">%s:/Mesh/topology</DataItem>
            </Topology>
            <Geometry GeometryType="XYZ">
                <DataItem Dimensions="%d 3" Format="HDF">%s:/Mesh/coordinates</DataItem>
            </Geometry>
        </Grid>
    """%(num_elems, num_elems, mesh_filename, num_points, mesh_filename)

    item_up = """
          <Grid Name="%%s">
            <xi:include xpointer="xpointer(//Grid[@Name=&quot;TimeSeries&quot;]/Grid[1]/*[self::Topology or self::Geometry])" />
            <Time Value="%%s" />
            <Attribute Name="u" AttributeType="Vector" Center="Node">
              <DataItem Dimensions="%d 3" Format="HDF">%%s:/Solution/u</DataItem>
            </Attribute>
            <Attribute Name="p" AttributeType="Scalar" Center="Node">
              <DataItem Dimensions="%d 1" Format="HDF">%%s:/Solution/p</DataItem>
            </Attribute>
          </Grid>
    """%(num_points, num_points)

    item_q = """
          <Grid Name="%%s">
            <xi:include xpointer="xpointer(//Grid[@Name=&quot;TimeSeries&quot;]/Grid[1]/*[self::Topology or self::Geometry])" />
            <Time Value="%%s" />
            <Attribute Name="Q" AttributeType="Scalar" Center="Node">
              <DataItem Dimensions="%d 1" Format="HDF">%%s:/Computed/qcriterion</DataItem>
            </Attribute>
          </Grid>
    """%(num_points)

    item_wss = """
          <Grid Name="%%s">
            <xi:include xpointer="xpointer(//Grid[@Name=&quot;TimeSeries&quot;]/Grid[1]/*[self::Topology or self::Geometry])" />
            <Time Value="%%s" />
            <Attribute Name="wss" AttributeType="Vector" Center="Node">
              <DataItem Dimensions="%d 3" Format="HDF">%%s:/Computed/wss</DataItem>
            </Attribute>
            <Attribute Name="wss_mag" AttributeType="Scalar" Center="Node">
              <DataItem Dimensions="%d 1" Format="HDF">%%s:/Computed/wss_abs</DataItem>
            </Attribute>
          </Grid>
    """%(num_wall_points, num_wall_points)

    up_h5_files = glob.glob(input_folder+"/*_up.h5")
    wss_h5_files = glob.glob(input_folder+"/wss_files/*_wss.h5")
    if not up_h5_files or not wss_h5_files:
        print ("There are no h5 files in the results folders. Exiting...")
        exit(1)

    up_h5_files = ['/'.join(fname.split('/')[2:]) for fname in up_h5_files]
    wss_h5_files = ['/'.join(fname.split('/')[2:]) for fname in wss_h5_files]

    up_h5_files = sorted(up_h5_files, key=lambda N: int(N.split("ts=")[1].split("_up.h5")[0]))
    up_h5_files = up_h5_files[0:len(up_h5_files):interval]
    up_file_count = len(up_h5_files)

    wss_h5_files = sorted(wss_h5_files, key=lambda N: int(N.split("ts=")[1].split("_wss.h5")[0]))
    wss_h5_files = wss_h5_files[0:len(wss_h5_files):interval]
    wss_file_count = len(wss_h5_files)

    global step_name_leading_zeros, time_name_leading_zeros
    step_stamp,time_stamp = get_timestep_time_string(up_h5_files[-1])
    while step_stamp[step_name_leading_zeros]=='0': step_name_leading_zeros += 1
    while time_stamp[time_name_leading_zeros]=='0': time_name_leading_zeros += 1

    case_full_name = os.path.basename(input_folder)

    print("- Making ", case_full_name+"_up.xdmf")
    out = open(os.path.join(input_folder, case_full_name+'_up.xdmf'), 'w')
    out.write(header+header_volume_mesh)
    for name in up_h5_files: out.write( item_up%(*get_timestep_time_string(name), name, name) )
    out.write(ending)
    out.close()

    print("- Making ", case_full_name+"_wss.xdmf")
    out = open(os.path.join(input_folder, case_full_name+'_wss.xdmf'), 'w')
    out.write(header+header_surface_mesh)
    for name in wss_h5_files: out.write( item_wss%(*get_timestep_time_string(name), name, name) )
    out.write(ending)
    out.close()

    print("- Making ", case_full_name+"_q.xdmf")
    out = open(os.path.join(input_folder, case_full_name+'_q.xdmf'), 'w')
    out.write(header+header_volume_mesh)
    for name in wss_h5_files: out.write( item_q%(*get_timestep_time_string(name), name) )
    out.write(ending)
    out.close()

if __name__ == '__main__':
    input_result_folder = sys.argv[1]
    interval = 1
    if len(sys.argv)==3: interval = int(sys.argv[2])

    make_xdmf_files (input_result_folder, interval)
  
