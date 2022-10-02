#///////////////////////////////////////////////////////////////
#// Artery.py
#// Copyright (C) 2018 Mehdi Najafi (mnuoft@gmail.com)
#//
#// Distribution of this library is not allowed in any form.
#///////////////////////////////////////////////////////////////

__author__ = "Mehdi Najafi <mnuoft@gmail.com>"
__date__ = "2018-07-12"
__copyright__ = "Copyright (C) 2018 " + __author__
__license__  = "Private"

#///////////////////////////////////////////////////////////////

import numpy as np
from probe import *
import sys, time
from os import path, makedirs, getcwd
import pickle
from Womersley import *
import CustomFunction
from problem_init import *
import naming
import os, h5io

#///////////////////////////////////////////////////////////////
# MPI node identification and size
mpi_size = MPI.size(MPI.comm_world)
mpi_rank = MPI.rank(MPI.comm_world)
master = mpi_rank == 0
if master:
    from dolfin import Timer
    g_t = Timer()
    initial_wall_time = g_t.elapsed()[0]

h5stdio = h5io.HDF5StdIO()

max_wall_time_before_being_killed = (23.7*60*60)
#///////////////////////////////////////////////////////////////
def print_section_header(title, l=100):
    if master:
        print ("-"*l)
        print(title)
        sys.stdout.flush()
def print_section_footer(l=100):
    if master:
        print ("-"*l)
        sys.stdout.flush()
def tuple2str(t, fmt='%12.10f'):
    return ','.join([fmt]*len(t))%tuple(t)
#///////////////////////////////////////////////////////////////
# Get the zero leading string for given time step No.
def step_str(i, l=10):
    a = str(i)
    la = l - len(a)
    return '0'*la + a

#///////////////////////////////////////////////////////////////
# Type check parser
def get_cmdarg(cmdline, key, default_value = None):
    if key in cmdline.keys():
        value = cmdline[key]
        if default_value:
            if type(default_value) is int:
                return int(value)
            if type(default_value) is float:
                return float(value)
            if type(default_value) is bool:
                return bool(eval(value))
        return value
    return default_value

#///////////////////////////////////////////////////////////////
# Retrive the boundary information from the info file
def read_mesh_info(mesh_info_path, key):
    # Extract inflow rate and outflow split ratios
    # Sample
    # p162
    # <INLETS>
    # 3 ICA_V27:FC_MCA_10 (3.4278309480,13.32740523503,-28.2355071983) (-0.2155297545,0.12005896011,-0.9690886291) 1.7649351905 9.7860492613   A*0.27
    #
    # <OUTLETS>
    # 1  None  (-16.8228963362,-1.42906111694,17.7845745237)  (-0.9847740285,-0.16502447805,-0.0546537689)   0.8568220486   2.3063814691   0.3318956234191912
    # 2  None    (7.9704219814,-9.22385217254,15.0763376349)   (0.7251682706,-0.49373752256,-0.4799523290)   1.3398999124   5.6402011155   0.6681043765808088

    #info = open(path.splitext(mesh_path)[0]+'.info', 'r').read()
    info = open(mesh_info_path, 'r').read()

    # Looking for the given key
    p1 = info.find(key)
    if p1<0:
        return [], [], []
    p1 += len(key)
    p2 = info.find('<', p1)
    if p2<0:
        buf = info[p1:]
    else:
        buf = info[p1:p2-1]
    lines = buf.split('\n')

    ids = []
    idfr = []
    ida = []
    idr = []
    fcs = []
    # Reaing at the key values
    for line in lines:
        ls = line.split()
        if len(ls) > 1:
            # id
            ids.append(int(ls[0]))# eval(ls[ 0]))
            # radius
            idr.append(eval(ls[-3]))
            # area
            ida.append(eval(ls[-2]))
            # flowrate or arearatio
            s = ls[-1].replace('A[','a[').replace('R[','r[')
            s = s.replace('A',ls[-2]).replace('R',ls[-3])
            idfr.append(s)
            # wave form
            fcs.append(ls[1])
    # evaluate all the expressions in the flowrates and outflow ratios
    for i,expr in enumerate(idfr):
        for j,k in enumerate(ids):
             expr = expr.replace( 'r[%d]'%k, str(idr[j])).replace( 'a[%d]'%k, str(ida[j]))
        idfr[i] = eval(expr)

    # sum of area ratio correction:
    if key == '<OUTLETS>':
        idfr[-1] = 1.0 - sum(idfr[:-1])

    # print the summary
    for i,s in enumerate(idfr):
        if master and key == '<INLETS>': print ('Inlet  id:', ids[i], ' flowrate:', s, ' mL/s')
        if master and key == '<OUTLETS>': print ('Outlet id:', ids[i], ' flowrate ratio:', s)


    return ids, idfr, ida, fcs

# check if the period is mentioned in the fc waveform file
def _not_used_get_period_from_fcs(fcs):
    periods = [951.0 for f in fcs]
    for i,f in enumerate(fcs):
        fcs_i_filename = f.split(':')[-1]
        if path.exists( path.join('./data', fcs_i_filename) ):
            fcs_ifname = path.join('./data', fcs_i_filename)
        else:
            fcs_ifname = path.join(path.dirname(path.abspath(__file__)), 'data', fcs_i_filename)
            if not path.exists( fcs_ifname ):
                print ('<!> Cannot find the waveform:', fcs_i_filename)
        for line in open(fcs_ifname,'r').readlines():
            if line.strip()[0] in ['#','!','/']:
                p = line.find('period_ms')
                if p > 0:
                    periods[i] = float(''.join((ch if ch in '0123456789.-e' else ' ') for ch in line[p+9:]).strip().split(' ')[0])
    return periods

#/////////////////////////////////////////////////////////////////////////
# Override some problem specific parameters: (Parameters are in mm and ms)
def problem_parameters(commandline_kwargs, NS_parameters, **NS_namespace):
    restart_folder = get_cmdarg(commandline_kwargs, 'restart_folder')
    if restart_folder:
        f = open(path.join(restart_folder, 'params.dat'), 'rb')
        NS_parameters.update(pickle.load(f))
        f.close()
        NS_parameters['restart_folder'] = restart_folder
        #Assign case name: pipe_ipcs_ab_cn_VWI_03k_LICA_constant_ts10000_cycles3_uOrder1
        case_name = NS_parameters['case_name']
        case_fullname = NS_parameters['case_fullname']


        mesh_path = path.join('./data', get_cmdarg(commandline_kwargs, 'mesh_name')+'.xml.gz')
        mesh_info_path = path.join('./data', get_cmdarg(commandline_kwargs, 'mesh_name')+'.info')
        if mesh_path == None:
            print('<!> Unable to run without a mesh file.')

        # Extract inflow rate and outflow split ratios
        if master: print('Reading mesh information:', mesh_info_path)
        id_in, Q_means, inlet_area, fcs = read_mesh_info(mesh_info_path, '<INLETS>')
        id_out, area_ratio, outlet_area, _ = read_mesh_info(mesh_info_path, '<OUTLETS>')

        NS_parameters.update( max_allowed_wall_time = get_cmdarg(commandline_kwargs, 'maxwtime', max_wall_time_before_being_killed) )

    else:
        mesh_path = path.join('./data', get_cmdarg(commandline_kwargs, 'mesh_name')+'.xml.gz')
        mesh_info_path = path.join('./data', get_cmdarg(commandline_kwargs, 'mesh_name')+'.info')
        if mesh_path == None:
            print('<!> Unable to run without a mesh file.')

        # Extract inflow rate and outflow split ratios
        if master: print('Reading mesh information:', mesh_info_path)
        id_in, Q_means, inlet_area, fcs = read_mesh_info(mesh_info_path, '<INLETS>')
        id_out, area_ratio, outlet_area, _ = read_mesh_info(mesh_info_path, '<OUTLETS>')

        case_name = get_cmdarg(commandline_kwargs, 'mesh_name') #NS_parameters["mesh_path"].split(".")[0]
        #Assign case full name: pipe_ipcs_ab_cn_VWI_03k_LICA_constant_ts10000_cycles3_uOrder1_curcyc_2_t_1902.1_u_ts=20000_000000
        # case_fullname = "pipe_ipcs_ab_cn_" + get_cmdarg(commandline_kwargs, 'mesh_name') + \
        #             "_constant" + \
        #             "_ts" + str(get_cmdarg(commandline_kwargs, 'timesteps')) + \
        #             "_cycles" + str(get_cmdarg(commandline_kwargs, 'cycles')) + \
        #             "_uOrder" + str(get_cmdarg(commandline_kwargs, 'uOrder'))

        # check if the period is mentioned in the fc waveform file
        if 'period' in commandline_kwargs.keys():
            period = get_cmdarg(commandline_kwargs, 'period', 951.0)
        else:
            periods = naming.get_period_from_fcs(fcs)
            if len(periods) > 1:
                if master:
                    print('Found out multiple periods: %s (ms)'%str(periods))
                    print('<!> The current implementation cannot handle multiple periods. The minimum is considered for all inlets.')
                period = np.min(periods)
            else:
                period = periods[0]
        if master:
            print('Found out period=%s (ms)'%str(period))

        timesteps = get_cmdarg(commandline_kwargs, 'timesteps', 2000)

        txt = ''
        for i, id in enumerate(id_in):
            txt += '_I%d_%s_Q%d'%(id,fcs[i].replace(':','_'),int(Q_means[i]*100))
        txt += '_Per%d'%int(period)

        case_fullname = "art_" + get_cmdarg(commandline_kwargs, 'mesh_name') + \
                    txt + \
                    "_Newt370" + \
                    "_ts" + str(get_cmdarg(commandline_kwargs, 'timesteps')) + \
                    "_cy" + str(get_cmdarg(commandline_kwargs, 'cycles')) + \
                    "_uO" + str(get_cmdarg(commandline_kwargs, 'uOrder'))

        # Parameters are in mm and ms
        NS_parameters.update(
            max_allowed_wall_time = get_cmdarg(commandline_kwargs, 'maxwtime', max_wall_time_before_being_killed),
            case_name = case_name,
            case_fullname = case_fullname,
            nu = get_cmdarg(commandline_kwargs, 'viscosity', 0.0035),
            period = period,
            T  = period * get_cmdarg(commandline_kwargs, 'cycles', 1),
            dt = period / timesteps,
            time_steps = timesteps,
            velocity_degree = get_cmdarg(commandline_kwargs, 'uOrder', 1),
            folder = "./results/" + case_fullname,
            save_step = 100000,
            checkpoint = 500,
            no_of_cycles = get_cmdarg(commandline_kwargs, 'cycles', 2),
            mesh_path = mesh_path, # commandline_kwargs["mesh_path"],
            id_in = id_in,
            id_out = id_out,
            area_ratio = area_ratio,
            not_zero_pressure_outlets = not get_cmdarg(commandline_kwargs, 'zero_pressure_outlets', False),
            include_gravitational_effects = get_cmdarg(commandline_kwargs, 'include_gravitational_effects', False),
            flat_profile_at_intlet_bc = get_cmdarg(commandline_kwargs, 'flat_profile_at_intlet_bc', False),
            dump_stats = 500,
            store_data = get_cmdarg(commandline_kwargs, 'save_frequency', 5),
            compute_flux = 5,
            #plot_interval = 10e10,
            print_intermediate_info = 100,
            use_krylov_solvers = True,
            krylov_solvers = dict(monitor_convergence=False)
        )

    if master: info_gray(str(NS_parameters))

    f = int(np.log10(NS_parameters['T']))+2
    g = int(np.log10(NS_parameters['time_steps'])*NS_parameters['no_of_cycles'])+1
    h5stdio.init(NS_parameters['folder'], case_fullname+'_curcyc_%%d_t=%%0%d.4f_ts=%%0%dd_up.h5'%(f+4,g))

    #caseName = NS_parameters["mesh_path"].split(".")[0]
    #NS_parameters["folder"] = path.join(NS_parameters["folder"], caseName)

#///////////////////////////////////////////////////////////////
# Create a mesh
def mpi_comm():
    return MPI.comm_world
def mesh(mesh_path, **NS_namespace):
    print_section_header('Loading mesh file: ' + mesh_path)

    mesh_folder = mesh_path #path.join(path.dirname(path.abspath(__file__)), mesh_path)

    m =  Mesh(mesh_folder)
    m.mpi_comm = mpi_comm

    # m = h5stdio.init_mesh(mesh_folder, int(MPI.comm_world.underlying_comm()))
    # print (m, type(m))
    # m.mpi_comm = mpi_comm

    num_points = Function(FunctionSpace(m, "CG", 1)).vector().size()
    vol = MPI.sum(MPI.comm_world, assemble(Constant(1)*dx(m)))
    cell_dia= [Cell(m,i).circumradius() for i in range (m.num_cells())]
    avg_cell_dia = sum(cell_dia) / len(cell_dia)
    num_cells   = int( MPI.sum(MPI.comm_world, m.num_cells()) )
    # num_points   = int( MPI.sum(MPI.comm_world, m.num_vertices()) ) // shared points?
    hmin        = MPI.min(MPI.comm_world, m.hmin())
    hmax        = MPI.max(MPI.comm_world, m.hmax())
    num_facets = int( MPI.sum(MPI.comm_world, m.num_facets()) )
    pss = mesh_folder.rfind('/')
    pss = 0 if pss < 0 else pss+1
    pos = mesh_folder.rfind('.xml.gz')
    if pos < 0: pos = mesh_folder.rfind('.')
    mesh_h5_filename = mesh_folder[pss:pos]+'.h5'
    mesh_h5_filepathname = os.path.join( NS_namespace['folder'], mesh_h5_filename)

    if master:
        #print ("-"*100)
        print ("Mesh Name:              ", mesh_folder)
        print ("Number of cells:        ", num_cells)
        print ("Number of points:       ", num_points)
        print ("Number of facets:       ", num_facets)
        print ("Mesh Volume:            ", vol)
        print ("Min cell diameter:      ", hmin)
        print ("Max cell diameter:      ", hmax)
        print ("Average cell diameter:  ", avg_cell_dia)
        sys.stdout.flush()
        #info(m, False)

    fd = MeshFunction("size_t", m, m.geometry().dim() - 1, m.domains())

    inout_area = {}
    dS = {}
    for id in NS_namespace['id_in']:
        dS[id] = ds(id, domain=m, subdomain_data=fd)
        inout_area[id] = abs( assemble(1.0*dS[id]) )
    for id in NS_namespace['id_out']:
        dS[id] = ds(id, domain=m, subdomain_data=fd)
        inout_area[id] = abs( assemble(1.0*dS[id]) )
    NS_namespace['inout_area'] = inout_area

    normals = FacetNormal(m)

    if master:
        print('writing ', mesh_h5_filepathname)
        sys.stdout.flush()

    ## Output the Mesh file into HDF5 format
    Hdf = HDF5File(m.mpi_comm(), mesh_h5_filepathname, "w")
    Hdf.write(m, '/Mesh')
    Hdf.close()

    h5stdio.SetMeshInfo(mesh_h5_filepathname, mesh_h5_filename, num_cells, num_points)

    print_section_footer()

    return m, dS, fd, normals, m.geometry().dim(), inout_area

#///////////////////////////////////////////////////////////////
# Read Inflow wave form and return the flow rate at all times
def flow_waveform(Qmean, cycles, period, time_steps, FC):
    omega = (2.0 * np.pi / period) #* cycles
    an = []
    bn = []

    #Load the Fourier Coefficients
    infile_FC = open( path.join(path.dirname(path.abspath(__file__)), 'data', FC), 'r').readlines()
    for line in infile_FC:
        abn = line.split()
        an.append(float(abn[0]))
        bn.append(float(abn[1]))

    t_values = np.linspace(0, period*cycles, num=time_steps)
    Q_values = []
    for t in t_values:
        Qn = 0 + 0j
        t1 = t / cycles
        for i in range (len(an)):
            Qn = Qn + (an[i]-bn[i]*1j)*np.exp(1j*i*omega*t1)
        Qn = abs(Qn)
        Q_values.append( Qmean * Qn )
        #print (t, Qn)
    return t_values, Q_values


#///////////////////////////////////////////////////////////////
# Boundary conditions
def create_bcs(u_, p_, p_1, t, NS_expressions, V, Q, area_ratio, mesh, subdomain_data, 
               dS, normals, folder, mesh_path, nu,
               id_in, id_out, velocity_degree, pressure_degree, no_of_cycles,
               T, not_zero_pressure_outlets, flat_profile_at_intlet_bc, **NS_namespace):

    print_section_header('Inspecting boundaries and making boundary conditions:')

    # Mesh function / boundaries
    fd = subdomain_data
    #fd = MeshFunction("size_t", mesh, mesh.geometry().dim() - 1, mesh.domains())

    #print('domains:')
    #f = dolfin.cpp.MeshEntity(mesh, 0, 0)
    #for e in cpp.entities(f, 1):
    #    print (e.index())

    # Extract inflow rate and outflow split ratios
    mesh_info_path = path.join('./data', NS_namespace["mesh_name"]+'.info')
    id_in, Q_means, inlet_area, fcs = read_mesh_info(mesh_info_path, '<INLETS>')
    id_out, area_ratio, outlet_area, _ = read_mesh_info(mesh_info_path, '<OUTLETS>')

    # Noslip condition at wall
    # Create Boundary conditions for the velocity
    wall = Constant(0.0)
    bc_wall = DirichletBC(V, wall, fd, 0) # wall is always with the id zero
    bc_wall_len = len(bc_wall.get_boundary_values())
    if master:
        print( 'Wall BC on ' + str(bc_wall_len) , 'cells')

    # Womersley boundary condition at inlet
    id_in_count = len(id_in)
    if master:
        print('Inlet', 'BCs' if id_in_count > 1 else 'BC', 'on boundaries:' if id_in_count > 1 else 'on boundary', id_in)
        firststr = '    %8s    %-12s    %10s    %15s    %6s'%('inlet_id','wave_form','period(ms)','flowrate(mL/s)','cells')
        secondstr = 'Inlets & Outlets Information\n'+'  id   %-45s  %-45s   %-12s   %-12s'%('center','normal','radius','area')

    # Create inlet boundary conditions
    inlets = []
    inout_area = {}
    bc_inlet_u = [[],[],[]]
    for i in range(id_in_count):
        fcs_i_filename = fcs[i].split(':')[-1]
        if path.exists( path.join('./data', fcs_i_filename) ):
            fcs_ifname = path.join('./data', fcs_i_filename)
        else:
            fcs_ifname = path.join(path.dirname(path.abspath(__file__)), 'data', fcs_i_filename)
            if not path.exists( fcs_ifname ):
                print ('<!> Cannot find the waveform:', fcs_i_filename)

        tmp_a, tmp_c, tmp_r, tmp_n = compute_boundary_geometry_acrn(mesh, dS[id_in[i]], normals)
        # print ('- Wave Form:', fcs_ifname)
        fcs_ifname_base = fcs_ifname[fcs_ifname.rfind('/')+1:]
        if fcs_ifname_base[0:3] == 'FC_':
            if master:
                print ('- loading inflow wave form:', fcs_ifname)
            inlet_i = make_womersley_bcs_2(NS_namespace["period"], Q_means[i], fcs_ifname, mesh, nu, tmp_a, tmp_c, tmp_r, tmp_n, velocity_degree, flat_profile_at_intlet_bc)
        else:
            if master:
                print ('- loading custom inflowrate function:', fcs_ifname)
            inlet_i = CustomFunction.make_custom_function_bcs(NS_namespace["period"], Q_means[i], fcs_ifname, mesh, nu, tmp_a, tmp_c, tmp_r, tmp_n, velocity_degree, flat_profile_at_intlet_bc)
        inlets.append(inlet_i)
        bci = [DirichletBC(V, ilt, fd, id_in[i]) for ilt in inlet_i]
        for j in range(3): bc_inlet_u[j].append(bci[j])

        count = len( bci[0].get_boundary_values() )
        inout_area[id_in[i]] = tmp_a
        
        if master:
            # print (dir(dS[id_in[i]]))
            firststr += "\n    %8d    %-12s    %10g    %15.8g    %6d"%(id_in[i], fcs_i_filename, NS_namespace["period"], Q_means[i], count)
            secondstr += "\nI %2d   %-45s  %-45s   %-12.10f   %-12.10f"%( id_in[i], tuple2str(tmp_c), tuple2str(tmp_n), tmp_r, tmp_a)

    NS_expressions["inlet"] = inlets

    # reset the time in boundary condition expressions
    for inlet in NS_expressions["inlet"]:
        for uc in inlet:
            uc.set_t(t)

    if master: print(firststr)

    #bc_inlet_u = [[DirichletBC(V, inlets[j][i], str(id_in[j])) for j in range(id_in_count)] for i in range(3)]
    #bc_inlet_u = [[DirichletBC(V, inlets[j][i], boundaries, id_in[j]) for j in range(id_in_count)] for i in range(3)]

    for bc_u in bc_inlet_u:
        bc_u.append(bc_wall)

    # Create outlet boundary conditions
    id_out_count = len(id_out)
    bc_p = []
    if master:
        print('Outlet', 'BCs' if id_out_count > 1 else 'BC', 'on boundaries:' if id_out_count > 1 else 'on boundary', id_out)
        print("    outlet_id    mass_flow_ratio      cells")
    for i, ind in enumerate(id_out):
        tmp_a, tmp_c, tmp_r, tmp_n = compute_boundary_geometry_acrn(mesh, dS[id_out[i]], normals)
        inout_area[ind] = tmp_a
        if master:
            secondstr += "\nO %2d   %-45s  %-45s   %-12.10f   %-12.10f"%( ind, tuple2str(tmp_c), tuple2str(tmp_n), tmp_r, tmp_a)
        if not_zero_pressure_outlets:
            if NS_parameters['restart_folder']:
                p_initial = assemble(p_*dS[ind]) / inout_area[ind]
            else:
                p_initial = area_ratio[i]
        else:
            p_initial = 0
        outflow = Expression("p", p=p_initial, degree=pressure_degree)
        bc = DirichletBC(Q, outflow, fd, ind)
        bc_p.append(bc)
        NS_expressions[ind] = outflow
        count  = len(bc.get_boundary_values())
        if master: print(' '*8, '%4d    %14.12f    %8d'%(ind, p_initial, count))

    if master: print(secondstr)

    #vars().update( dict(inout_area=inout_area) )


    # pressure_out = {}
    # for id in id_out:
    #     #inout_area[id] = abs( assemble(1.0*ds(id, domain=mesh, subdomain_data=fd)) )
    #     pressure_out[id] = assemble(p_*dS[id]) / inout_area[id]
    #     print ('%d  %1.15f'%(id, pressure_out[id]))


    print_section_footer(132)

    # Return boundary conditions in dictionary
    return dict(u0=bc_inlet_u[0], u1=bc_inlet_u[1], u2=bc_inlet_u[2], p=bc_p)

#///////////////////////////////////////////////////////////////
def get_file_paths(folder):
    if master:
        counter = 1
        to_check = path.join(folder, "data", "%s")
        while path.isdir(to_check % str(counter)):
            counter += 1

        if counter > 1:
            counter -= 1
        # if not path.exists(path.join(to_check % str(counter), "VTK")):
        #     makedirs(path.join(to_check % str(counter), "VTK"))
    else:
        counter = 0

    counter = MPI.max(MPI.comm_world, counter)

    common_path = path.join(folder, "data", str(counter), "VTK")
    file_u = [path.join(common_path, "u%d.h5" % i) for i in range(3)]
    file_p = path.join(common_path, "p.h5")
    file_nu = path.join(common_path, "nut.h5")
    file_u_mean = [path.join(common_path, "u%d_mean.h5" % i) for i in range(3)]
    files = {"u": file_u, "p": file_p, "u_mean": file_u_mean, "nut": file_nu}

    return files


#///////////////////////////////////////////function////////////////////
def pre_solve_hook(mesh, V, Q, newfolder, folder, u_, mesh_path,
                   restart_folder, velocity_degree, **NS_namespace):

    # Create point for evaluation
    # fd = MeshFunction("size_t", mesh, mesh.geometry().dim()-1, mesh.domains())

    # dS = {}
    # inout_area = {}
    # for id in NS_namespace['id_in']:
    #     dS[id] = ds(id, domain=mesh, subdomain_data=fd)
    #     inout_area[id] = abs( assemble(dS[id]) )
    # for id in NS_namespace['id_out']:
    #     dS[id] = ds(id, domain=mesh, subdomain_data=fd)
    #     inout_area[id] = abs( assemble(dS[id]) )

    # n = FacetNormal(mesh)

    eval_dict = {}
    # rel_path = path.join(path.dirname(path.abspath(__file__)), mesh_path.split(".")[0],
    #                     "{}_probe_point".format(mesh_path.split(".")[0]))
    # probe_points = np.load(rel_path)

    # # Store points file in checkpoint
    # if master:
    #     probe_points.dump(path.join(newfolder, "Checkpoint", "points"))

    # eval_dict["centerline_u_x_probes"] = Probes(probe_points.flatten(), V)
    # eval_dict["centerline_u_y_probes"] = Probes(probe_points.flatten(), V)
    # eval_dict["centerline_u_z_probes"] = Probes(probe_points.flatten(), V)
    # eval_dict["centerline_p_probes"] = Probes(probe_points.flatten(), Q)

    # Link for io
    #hdf5_link = HDF5Link().link

    if restart_folder is None:
        # Get files to store results
        files = get_file_paths(folder)
        NS_parameters.update(dict(files=files))
    else:
        files = NS_namespace["files"]

    return dict(eval_dict=eval_dict, hdf5_link=h5stdio,
                files=files, #inout_area=NS_parameters['inout_area'],
                final_time=NS_namespace['T'], current_cycle=0, 
                timesteps=NS_namespace['time_steps'], total_cycles=NS_namespace['no_of_cycles'],
                timestep_cpu_time=0, current_time=time.time(), cpu_time=0)

#///////////////////////////////////////////////////////////////
def beta(err, p):
    if p < 0:
        if err >= 0.1:
            return 0.5
        else:
            return 1.0 - 5*err**2
    else:
        if err >= 0.1:
            return 1.5
        else:
            return 1.0  + 5*err**2

#///////////////////////////////////////////////////////////////
def w(P):
    return 1.0 / ( 1.0 + 20.0*abs(P))
#///////////////////////////////////////////////////////////////
def temporal_hook(u_, p_, p, q_, mesh, tstep, compute_flux,
                  dump_stats, eval_dict, newfolder, id_in, files, id_out, inout_area,
                  normals, store_data, hdf5_link, NS_expressions, current_cycle,
                  total_cycles, area_ratio, t, dS, timestep_cpu_time, current_time, 
                  cpu_time, final_time, timesteps, hh, not_zero_pressure_outlets, **NS_namespace):

    # update the current cycles
    current_cycle = int(tstep / timesteps)
    if master:
        print ('cycle:', current_cycle, 'tstep', tstep , ' timesteps', timesteps)

    # Update boundary condition
    for inlet in NS_expressions["inlet"]:
        for uc in inlet:
            uc.set_t(t)

    timestep_cpu_time = time.time() - current_time
    current_time = time.time()
    cpu_time += timestep_cpu_time

    # Do not proceed if the time step is less than 3
    if tstep < 3: return

    Q_ideals = {}
    # In-Going Flux & pressure
    flux_in = {}
    Q_ins = {}
    pressure_in = {}
    for id in id_in:
        #inout_area[id] = abs( assemble(1.0*ds(id, domain=mesh, subdomain_data=fd)) )
        pressure_in[id] = -assemble(p_*dS[id]) / inout_area[id]
        flux_in[id] = assemble(dot(u_, normals)*dS[id])
        Q_ins[id] = abs(flux_in[id])
    Q_ins_sum = sum(Q_ins.values())
    if master: print ('Q_ins:', Q_ins_sum, Q_ins)

    # Out-Going Flux
    flux_out = {}
    Q_outs =  {}
    pressure_out = {}
    for id in id_out:
        #inout_area[id] = abs( assemble(1.0*ds(id, domain=mesh, subdomain_data=fd)) )
        pressure_out[id] = assemble(p_*dS[id]) / inout_area[id]
        flux_out[id] = assemble(dot(u_, normals)*dS[id])
        Q_outs[id] = abs(flux_out[id])
    Q_outs_sum = sum(Q_outs.values())

    # Compute flux and update pressure condition
    if not_zero_pressure_outlets:
      for i, out_id in enumerate(id_out):
        Q_ideals[i] = area_ratio[i]*Q_ins_sum
        p_old = NS_expressions[out_id].p

        # Gin and Steinman et al., A Dual-Pressure Boundary Condition
        # for use in Simulations of Bifurcating Conduits
        R_optimal = area_ratio[i]
        R_actual = Q_outs[out_id] / Q_ins_sum

        M_err = abs(R_optimal / R_actual)
        R_err = abs(R_optimal - R_actual)

        if p_old < 0:
            E = 1 + R_err / R_optimal
        else:
            E = -1 * ( 1 + R_err / R_optimal )

        # 1) Linear update to converge first 100 tsteps of the first cycle
        delta = (R_optimal - R_actual) / R_optimal
        if tstep < 100:
            h = 0.1
            if p_old > 1 and delta < 0:
                NS_expressions[out_id].p  = p_old
            else:
                NS_expressions[out_id].p  = p_old * ( 1 - delta*h)

        # 2) Dual pressure BC
        else:
            if p_old > 2 and delta < 0:
                NS_expressions[out_id].p  = p_old
            else:
                NS_expressions[out_id].p  = p_old * beta(R_err,p_old) * M_err ** E
        

    #Print the flow rates, fluxes, pressure
    if master:
        # #Flux In/Out Flux, Velocity, Pressure
        print ("~"*88)
        print ("Flow Rate or Flux Error is: ", 100.*(abs(sum(flux_in.values()))-abs(sum(flux_out.values())))/abs(sum(flux_in.values())),"%")
        print ("~"*88)
        print ("%3s  %2s          %-15s  %-18s  %-18s  %-18s  %-18s"%('I/O','id','Flux', 'Ideal_Flux', 'Velocity', 'Pressure', 'New Pressure'))
        for id in id_in:
            print ("%-3s  %2d  % 16.15f  % 16.15f  % 16.15f  % 16.15f"%('In', id, flux_in[id], flux_in[id], flux_in[id]/inout_area[id], pressure_in[id]))
        for i, id in enumerate(id_out):
            print ("%-3s  %2d  % 16.15f  % 16.15f  % 16.15f  % 16.15f  % 16.15f"%('Out', id, flux_out[id], area_ratio[i]*Q_ins_sum, flux_out[id]/inout_area[id], pressure_out[id], NS_expressions[id].p))
        print ("~"*88)

        sys.stdout.flush()

        elapsed_wall_time = g_t.elapsed()[0] - initial_wall_time

        if elapsed_wall_time > NS_parameters["max_allowed_wall_time"]:
            open(path.join(newfolder,'terminate_solver'),'w').close()


        # # Print progress
        # ss = cpu_time * (final_time/t-1)
        # hh, ss = divmod(ss, 60*60)
        # mm, ss = divmod(ss, 60)
        # print ()
        # s = "Time step %d finished in %.2f seconds, %.1f%% done (t=%.3g, T=%g; %02d:%02d:%02d remaining)." \
        #     % (tstep, timestep_cpu_time, 100.0*t/final_time, t, final_time, hh, mm, ss)
        # print ("-"*len(s))
        # print (s)
        # print ("-"*len(s))


    # # Sample velocity in points
    # eval_dict["centerline_u_x_probes"](u_[0])
    # eval_dict["centerline_u_y_probes"](u_[1])
    # eval_dict["centerline_u_z_probes"](u_[2])
    # eval_dict["centerline_p_probes"](p_)

    # # Store sampled velocity
    # if tstep % dump_stats == 0:
    #     filepath = path.join(newfolder, "Stats")
    #     if master:
    #         if not path.exists(filepath):
    #             makedirs(filepath)

    #     arr_u_x = eval_dict["centerline_u_x_probes"].array()
    #     arr_u_y = eval_dict["centerline_u_y_probes"].array()
    #     arr_u_z = eval_dict["centerline_u_z_probes"].array()
    #     arr_p = eval_dict["centerline_p_probes"].array()

    #     # Dump stats
    #     if master:
    #         num = eval_dict["centerline_u_x_probes"].number_of_evaluations()
    #         pp = (path.join(filepath, "u_x_%s.probes" % str(tstep)))
    #         arr_u_x.dump(path.join(filepath, "u_x_%s.probes" % str(tstep)))
    #         arr_u_y.dump(path.join(filepath, "u_y_%s.probes" % str(tstep)))
    #         arr_u_z.dump(path.join(filepath, "u_z_%s.probes" % str(tstep)))
    #         arr_p.dump(path.join(filepath, "p_%s.probes" % str(tstep)))

    #     # Clear stats
    #     MPI.barrier(mpi_comm_world
    #     eval_dict["centerline_u_x_probes"].clear()
    #     eval_dict["centerline_u_y_probes"].clear()
    #     eval_dict["centerline_u_z_probes"].clear()
    #     eval_dict["centerline_p_probes"].clear()

    # Save velocity and pressure
    if current_cycle == total_cycles-1:
        if tstep % store_data == 0:
            h5stdio.Save( current_cycle, t, tstep, Q_ins, Q_outs, NS_parameters, 'Step-%06d'%tstep, q_, int(MPI.comm_world.underlying_comm()) )
            if master:
                h5stdio.SaveXDMF( os.path.join(NS_parameters['folder'], NS_parameters['case_fullname']+'.xdmf') )
            # # Evaluate points
            # u_[0].rename("u0", "velocity-x")
            # u_[1].rename("u1", "velocity-y")
            # u_[2].rename("u2", "velocity-z")
            # p_.rename("p", "pressure")

            # # Store files 
            # components = {"u0": u_[0], "u1": u_[1], "u2": u_[2], "p": p_}

            # for key in components.keys():
            #     field_name = "velocity" if "u" in key else "pressure"
            #     if "u" in key and key != "nut":
            #         f = files["u"][int(key[-1])]
            #     else:
            #         f = files[key]
            #     save_hdf5(f, field_name, components[key], tstep, hdf5_link)

#/////////////////////////////////////////////////////////////////////////////////
def theend_hook(stop, newfolder, hh):
    if master:
        if stop:
            if path.exists(path.join(newfolder,'complete')):
                os.remove(path.join(newfolder,'complete'))
            last_lines = open(path.join(newfolder,'incomplete'),'w')
        else:
            if path.exists(path.join(newfolder,'incomplete')):
                os.remove(path.join(newfolder,'incomplete'))
            last_lines = open(path.join(newfolder,'complete'),'w')
        last_lines.write('\nTry: ' + newfolder.split('/')[-1])
        last_lines.write('\nCheckpoint: ' + newfolder)
        last_lines.write('\nRemaining Time(hr): %d'%hh)
        last_lines.close()
    print ('Process %d/%d Terminated.'%(mpi_rank,mpi_size))

#/////////////////////////////////////////////////////////////////////////////////
