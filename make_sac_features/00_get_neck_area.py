""" Get the neck plane and calculate area.
"""

""" Generate spectrograms using BSL tools.

Note, surface meshing file is uploaded, so extracting sac from that.
"""

import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = '/scratch/s/steinman/macdo708/.config/matplotlib'

import matplotlib.pyplot as plt 

import pyvista as pv
from bsl.dataset import Dataset
import numpy as np 
import pandas as pd 
from scipy.spatial import cKDTree as KDTree
from scipy.spatial import distance

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    neck_dir = Path(sys.argv[2])
    csv_dir = Path(sys.argv[3])
    meshing_dir = Path(sys.argv[4])

    dd = Dataset(folder)
    # case_name_short = dd.case_name.split('_mesh')[0]
    case_name_short = dd.case_name.split('_cl')[0]
    neck_surf_out = neck_dir / (dd.case_name + '.vtp')
    csv_out = csv_dir / (dd.case_name + '.csv')
    sac_out = neck_surf_out.parents[0] / 'sac_vols' / (dd.case_name + '.vtu')

    if not sac_out.parents[0].exists():
        sac_out.parents[0].mkdir(parents=True, exist_ok=True)

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

    neck = dd.extract_neck()
    neck = neck.compute_cell_sizes(area=True)
    neck.save(neck_surf_out)

    # Create SDF for sacs
    mask = dd.mesh.point_arrays['SacMask'] == 1
    sac_vol = dd.mesh.extract_points(mask)
    sac_vol = sac_vol.compute_cell_sizes()

    sac_surf = sac_vol.extract_surface()
    sac_surf = sac_surf.extract_points(sac_surf['SurfaceSacMask'] == 1, adjacent_cells=False)
    sac_surf = sac_surf.compute_cell_sizes()
    sac_surf_dist = distance.cdist(sac_surf.points, sac_surf.points)

    # meshsurf = dd.mesh.extract_surface()
    # tree = KDTree(meshsurf.points)
    # dist, idx = tree.query(sac_vol.points)

    # sac_vol.point_arrays['distance'] = dist
    sac_vol.save(sac_out)

    df = pd.DataFrame()
    df.at[dd.case_name, 'neck_area'] = neck.area #cell_arrays['Area'].sum()
    df.at[dd.case_name, 'surf_area'] = sac_surf.area #['Area'].sum() 
    df.at[dd.case_name, 'sac_vol'] = sac_vol.volume #['Volume'].sum()
    df.at[dd.case_name, 'sac_max_dist'] = np.percentile(sac_surf_dist, 95)

    # df.at[dd.case_name, 'sdf_max'] = dist.max()
    # df.at[dd.case_name, 'sdf_99per'] = np.percentile(dist, 99)
    # df.at[dd.case_name, 'sdf_neck_max'] = dist[sac_vol['SacMask'] == 0].max()
    # df.at[dd.case_name, 'sdf_neck_99per'] = np.percentile(dist[sac_vol['SacMask'] == 0], 99)
    
    df.to_csv(csv_out)