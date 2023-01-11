"""
Compute TAWSS, SCI, LSA.
"""

import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = '/scratch/s/steinman/macdo708/.config/matplotlib'

import matplotlib.pyplot as plt 

import pyvista as pv
from bsl.dataset import Dataset
import numpy as np 
from scipy.spatial import cKDTree as KDTree 
import pandas as pd 

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    out_dir = Path(sys.argv[2])

    case_wss_file = out_dir / (folder.stem + '.csv')

    dd = Dataset(folder)
    case_name_short = dd.case_name.split('_mesh')[0]

    meshing_dir = Path('/scratch/s/steinman/macdo708/surge_cfd/mesh/')
    vtu_file = meshing_dir / case_name_short / '03_mesh' / (dd.case_name + '.vtu')
    
    dd = dd.assemble_mesh().assemble_surface(mesh_file=vtu_file)
    
    # Surface meshing file
    surf_dir = meshing_dir / case_name_short / '02b_processed'
    
    surf = pv.read(surf_dir / (dd.case_name + '.vtp'))
    sac_mask = surf.point_arrays['Mask'] == 1
    sac = surf.extract_points(sac_mask, include_cells=True, adjacent_cells=True)
    sac = pv.PolyData(sac.points, faces=sac.cells)

    dd.create_sac_mask(sac=sac)

    sac_mask = dd.mesh.point_arrays['SacMask'] == 1

    # Want to get SCI and sac indices
    tree = KDTree(surf.points)
    _, ii = tree.query(dd.surf.points)
    dd.surf.point_arrays['sac_zones'] = surf.point_arrays['sac_zones'][ii]
    dd.surf.point_arrays['SacMask'] = surf.point_arrays['Mask'][ii]

    dd.surf.point_arrays['sac_zones'][dd.surf.point_arrays['SacMask'] == 1] = 0

    dd.compute_wss_parameters(mask_key='sac_zones', sac_mask_key='SacMask', indices=np.arange(960)[::int(32)]) # was 32

    df = pd.DataFrame()
    df.at[dd.case_name, 'SCI'] = dd.SCI['SacMask']
    df.at[dd.case_name, 'LSA'] = dd.LSA['SacMask']
    # df.at[dd.case_name, 'TAWSS_sac'] = dd.TAWSS_sac['SacMask']
    # df.at[dd.case_name, 'TAWSS_parent'] = dd.TAWSS_parent['SacMask']

    df.to_csv(case_wss_file)