#///////////////////////////////////////////////////////////////
#// h5io.py
#// Copyright (C) 2018 Mehdi Najafi (mnuoft@gmail.com)
#//
#// Distribution of this library is not allowed in any form.
#///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-07-30"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private; to be obtained directly from the author."


import os, sys
import numpy as np
from dolfin import compile_cpp_code, Function, FunctionSpace, Mesh, HDF5File, XDMFFile, MPI
from ufl import *

#///////////////////////////////////////////////////////////////
# MPI node identification and size
# mpi_size = MPI.COMM_WORLD.size
# mpi_rank = MPI.COMM_WORLD.rank
# mpi_comm_world = MPI.COMM_WORLD

# master = mpi_rank == 0

#///////////////////////////////////////////////////////////////
# Standard IO class
class HDF5StdIO:
    def init(self, output_folder, filename_time_template):
        self.timeseries_count = 1
        self.t = 0.0
        self.timestep = 0
        self.output_folder = output_folder
        self.filename_time_template = filename_time_template
        self.module = self.init_code()

    def init_code(self):
        return compile_cpp_code('''
            #include <array>
            #include <map>
            #include <numeric>
            #include <set>
            #include <utility>
            #include <unordered_map>
            #include <typeinfo>
            #include <string>
            #include <utility>
            #include <vector>
            #include <iostream>
            #include <pybind11/pybind11.h>
            #include <pybind11/embed.h>
            #include <pybind11/stl.h>
            #include <pybind11/numpy.h>
            #include <pybind11/cast.h>
            #include <pybind11/pytypes.h>
            #include <dolfin.h>
            #include <dolfin/common/MPI.h>
            #include <dolfin/common/Variable.h>
            #include <dolfin/io/HDF5Interface.h>
            #include <dolfin/function/Function.h>
            #include <dolfin/function/FunctionSpace.h>
            #include <dolfin/mesh/Mesh.h>
            #include <dolfin/mesh/DistributedMeshTools.h>

            using namespace dolfin;

            template<typename T> std::string to_string(const T& t)
            { std::ostringstream sstr; sstr << t; return sstr.str(); }

            std::int64_t get_padded_width(const dolfin::Function& u)
            {
                std::int64_t width = u.value_size();
                std::int64_t rank = u.value_rank();
                if (rank == 1 and width == 2)
                    return 3;
                else if (rank == 2 and width == 4)
                    return 9;
                return width;
            }

            std::vector<double> get_point_data_values(const Function& u)
            {
                const auto mesh = u.function_space()->mesh();

                std::vector<double> data_values;
                //dolfin_assert(mesh->geometry().degree() == 1);

                u.compute_vertex_values(data_values, *mesh);

                std::int64_t width = get_padded_width(u);

                if (u.value_rank() > 0)
                {
                    // Transpose vector/tensor data arrays
                    const std::size_t num_local_vertices = mesh->num_entities(0);
                    const std::size_t value_size = u.value_size();
                    std::vector<double> _data_values(width*num_local_vertices, 0.0);
                    for (std::size_t i = 0; i < num_local_vertices; i++)
                    {
                    for (std::size_t j = 0; j < value_size; j++)
                    {
                        std::size_t tensor_2d_offset = (j > 1 && value_size == 4) ? 1 : 0;
                        _data_values[i*width + j + tensor_2d_offset]
                        = data_values[i + j*num_local_vertices];
                    }
                    }
                    data_values = _data_values;
                }

                // Remove duplicates for vertex-based data in parallel
                if (MPI::size(mesh->mpi_comm()) > 1)
                {
                    DistributedMeshTools::reorder_values_by_global_indices(*mesh,
                                                                        data_values, width);
                }

                return data_values;
            }

            void write_vars(const std::string filename, const std::string file_mode,
                        const std::string attributes, pybind11::dict explicit_attributes, 
                        const std::string dataset_name,
                        dolfin::Function& uf0, dolfin::Function& uf1, dolfin::Function& uf2,
                        dolfin::Function& pf, const long int mpi_comm)
            {
                std::size_t i,j;
                std::vector<double> uc[3];
                uc[0] = get_point_data_values(uf0);
                uc[1] = get_point_data_values(uf1);
                uc[2] = get_point_data_values(uf2);
                std::vector<double> p = get_point_data_values(pf);
                std::vector<double> u(3*uc[0].size());
                for (j = i = 0; i < uc[0].size(); i++)
                {
                    u[j++] = uc[0][i];
                    u[j++] = uc[1][i];
                    u[j++] = uc[2][i];
                }

                std::int64_t num_values = pf.function_space()->mesh()->num_entities_global(0);

                std::vector<std::int64_t> shape_u = {num_values, 3};
                std::vector<std::int64_t> shape_p = {num_values, 1};

                // std::int64_t num_items_total_u = 1, num_items_total_p = 1;
                // for (auto n : shape_u) num_items_total_u *= n;
                // for (auto n : shape_p) num_items_total_p *= n;
                // dolfin_assert(num_items_total_u == (std::int64_t) MPI::sum(mpi_comm, u.size()));
                // dolfin_assert(num_items_total_p == (std::int64_t) MPI::sum(mpi_comm, p.size()));

                // Compute data offset and range of values
                std::int64_t local_shape_u = u.size();
                for (i = 1; i < shape_u.size(); ++i) local_shape_u /= shape_u[i];
                std::int64_t local_shape_p = p.size();
                for (i = 1; i < shape_p.size(); ++i) local_shape_p /= shape_p[i];

                const std::int64_t offset_u = MPI::global_offset(mpi_comm, local_shape_u, true);
                const std::int64_t offset_p = MPI::global_offset(mpi_comm, local_shape_p, true);
                const std::pair<std::int64_t, std::int64_t> local_range_u = {offset_u, offset_u + local_shape_u};
                const std::pair<std::int64_t, std::int64_t> local_range_p = {offset_p, offset_p + local_shape_p};

                // Write data to HDF5 file
                const bool use_mpi_io = true; // (MPI::size(comm) > 1);
                const bool chunking = false;

                hid_t hdf5_file_id = HDF5Interface::open_file(mpi_comm, filename, file_mode, use_mpi_io);
                HDF5Interface::write_dataset(hdf5_file_id, dataset_name+"/u", u,
                                            local_range_u, shape_u, use_mpi_io, chunking);
                HDF5Interface::write_dataset(hdf5_file_id, dataset_name+"/p", p,
                                            local_range_p, shape_p, use_mpi_io, chunking);
                for (auto item: explicit_attributes)
                    HDF5Interface::add_attribute(hdf5_file_id, dataset_name, to_string(item.first), to_string(item.second));
                HDF5Interface::add_attribute(hdf5_file_id, dataset_name, "Parameters", attributes);
                HDF5Interface::add_attribute(hdf5_file_id, dataset_name, "About", std::string("Aneurysm CFD v.1.2 - Copyright (c) Mehdi Najafi - BioMedical Simulation Lab - University of Toronto.") );
                HDF5Interface::add_attribute(hdf5_file_id, dataset_name, "Notes", std::string("The results you are looking at is out of a code that I have been developing for some years. It is expected to invite me for collaborations in your study and send me your feedbacks."));
                HDF5Interface::close_file(hdf5_file_id);
            }

            PYBIND11_MODULE(SIGNATURE, m)
            {
                m.def("write_vars", &write_vars);
            }
        ''')

    def SetTime(self, t, timestep):
        self.t = t
        self.timestep = timestep

    def Save(self, cycle, t, timestep, Q_ins, Q_outs, parameters, name, q, mpi_comm_world): #=mpi_comm_world):
        self.SetTime(t, timestep)
        filename = self.filename_time_template%(cycle, self.t, self.timestep)
        txt = self.xml_node_grid_vector_scalar_tmp%(name, self.t, filename, filename)
        self.XDMFText = self.XDMFText[:self.insert_pos] + txt + self.XDMFText[self.insert_pos:]
        self.insert_pos += len(txt)
        iflux = str(Q_ins)
        oflux = str(Q_outs)
        self.module.write_vars( os.path.join(self.output_folder, filename), 'w', str(parameters), 
                {'t':t, 'timestep':timestep, 'timesteps':parameters['time_steps'], 'cycle':cycle, 'cycles':parameters['no_of_cycles'], 'Qin(mL/s)':iflux, 'Qout(mL/s)':oflux},
                '/Solution', q['u0'].cpp_object(), q['u1'].cpp_object(), q['u2'].cpp_object(), q['p'].cpp_object(), mpi_comm_world)
        self.timeseries_count += 1

    def SaveXDMF(self, filename):
        f = open(filename, 'w')
        f.write(self.XDMFText)
        f.close()

    def SetMeshInfo(self, mesh_h5_filepathname, mesh_h5_filename, num_cells, num_points, u_degree=1, p_degree=1):
        self.num_cells   = num_cells
        self.num_points  = num_points
        self.num_u_points  = num_points * u_degree
        self.num_p_points  = num_points * p_degree
        self.XDMFText = '''<?xml version="1.0"?>
<!DOCTYPE Xdmf SYSTEM "Xdmf.dtd" []>
<Xdmf Version="3.0" xmlns:xi="http://www.w3.org/2001/XInclude">
<Domain>
    <Grid Name="TimeSeries" GridType="Collection" CollectionType="Temporal">
        <Grid Name="mesh" GridType="Uniform">
        <Topology NumberOfElements="%d" TopologyType="Tetrahedron" NodesPerElement="4">
            <DataItem Dimensions="%d 4" NumberType="UInt" Format="HDF">%s:/Mesh/topology</DataItem>
        </Topology>
        <Geometry GeometryType="XYZ">
            <DataItem Dimensions="%d 3" Format="HDF">%s:/Mesh/coordinates</DataItem>
        </Geometry>
        </Grid>
    </Grid>
</Domain>
</Xdmf>
'''%(self.num_cells, self.num_cells, mesh_h5_filename, self.num_points, mesh_h5_filename)
        self.insert_pos = self.XDMFText.rfind('</Grid>')
        self.insert_pos = self.XDMFText.rfind('\n', 0, self.insert_pos)+1
        self.xml_node_grid_vector_scalar_tmp = '''      <Grid Name="%%s">
	<xi:include xpointer="xpointer(//Grid[@Name=&quot;TimeSeries&quot;]/Grid[1]/*[self::Topology or self::Geometry])" />
        <Time Value="%%g" />
        <Attribute Name="u" AttributeType="Vector" Center="Node">
          <DataItem Dimensions="%d 3" Format="HDF">%%s:/Solution/u</DataItem>
        </Attribute>
        <Attribute Name="p" AttributeType="Scalar" Center="Node">
          <DataItem Dimensions="%d 1" Format="HDF">%%s:/Solution/p</DataItem>
        </Attribute>
      </Grid>
'''%(self.num_u_points, self.num_p_points)
