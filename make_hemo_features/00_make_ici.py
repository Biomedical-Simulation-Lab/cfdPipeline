""" Calculate inflow concentration index (ICI)
ICI = (Q_in / Q_out) / (A_in / A_out)

Load sim data
Load extracted neck
"""

from pathlib import Path 
import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = '/scratch/s/steinman/macdo708/.config/matplotlib'

import matplotlib.pyplot as plt 

import pyvista as pv
from bsl.dataset import Dataset
import numpy as np 
from geometry_tools import common as cc 

import pandas as pd 

# neck_surf_dir = Path('/scratch/s/steinman/macdo708/surge_cfd_analysis/data/neck_surfs_new')

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    data_out_folder = Path(sys.argv[2])
    neck_surf_dir = Path(sys.argv[3])
    meshing_dir = Path(sys.argv[4])
    parent_slice_dir = Path(sys.argv[5])

    data_out_file = data_out_folder / (folder.stem + '.csv')
    df = pd.DataFrame()

    if not data_out_file.exists():
        dd = Dataset(folder)
        case_name_short = dd.case_name.split('_cl')[0]

        print('*'*100)
        print(case_name_short)

        neck_surf_file = neck_surf_dir / (dd.case_name + '.vtp')
        neck_surf = pv.read(neck_surf_file)
        neck_surf = neck_surf.decimate(0.75)

        # meshing_dir = Path('/scratch/s/steinman/macdo708/surge_cfd/mesh/')
        vtu_file = meshing_dir / case_name_short / '03_mesh' / (dd.case_name + '.vtu')
        
        dd = dd.assemble_mesh().assemble_surface(mesh_file=vtu_file)
        
        # Surface meshing file
        surf_dir = meshing_dir / case_name_short / '02_processed'
        
        surf = pv.read(surf_dir / (dd.case_name + '.vtp'))
        sac_mask = surf.point_arrays['Mask'] == 1
        sac = surf.extract_points(sac_mask, include_cells=True, adjacent_cells=True)
        sac = pv.PolyData(sac.points, faces=sac.cells)

        dd.create_sac_mask(sac=sac)

        # Dev ICI
        indices = np.arange(0,len(dd),8)  

        # FIRST, need Q for parent.
        # There is def a nicer way to do this (identify this region during meshing)
        # but using the arrays that are already generated.
        origins = cc.get_parent_slice_location_from_sac_zones(surf)

        # Returns slices of the parent where we'll calc flow rate.
        # Note, there's two slices, calc flow rate for both,
        # should be same.
        slices = [cc.get_nearest_slice(dd.mesh, o, n) for o, n in zip(origins.points, origins.point_arrays['Normals'])]
        slices = [pv.PolyData(sl.points, sl.cells) for sl in slices]
        slices = [sl.triangulate() for sl in slices]
        slices = [sl.decimate(0.75) for sl in slices]
        parent_slice_file = parent_slice_dir / (dd.case_name + '.vtp')
        slices[0].save(parent_slice_file)
        
        ICI, ICI_timeseries, other_data = dd.compute_ICI(neck_surf, slices[0], sac, indices=indices)
        ICI2, ICI2_timeseries, _ = dd.compute_ICI(neck_surf, slices[1], sac, indices=indices)

        print(dd.case_name_short, ICI, ICI2)
    
        df.at[case_name_short, 'ICI'] = np.mean([ICI, ICI2])
        
        df.at[case_name_short, 'u_in_mean'] = other_data['u_in_mean']
        df.at[case_name_short, 'u_out_mean'] = other_data['u_out_mean']

        df.at[case_name_short, 'u_n_in_mean'] = other_data['u_n_in_mean']
        df.at[case_name_short, 'u_n_out_mean'] = other_data['u_n_out_mean']

        df.at[case_name_short, 'Q_i_mean'] = other_data['Q_i_mean']
        df.at[case_name_short, 'Q_v_mean'] = other_data['Q_v_mean']

        df.at[case_name_short, 'A_i_mean'] = other_data['A_i_mean']
        df.at[case_name_short, 'A_o'] = other_data['A_o']

        df.to_csv(data_out_file)