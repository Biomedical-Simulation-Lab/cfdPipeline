""" Extract sac-average OSI, SPI, TAWSS from hemodynamics.tec files.
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

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    df_out_dir = Path(sys.argv[2])
    meshing_dir = Path(sys.argv[3])

    df_out = df_out_dir / (folder.stem + '.csv')

    df = pd.DataFrame()
    case_name = folder.stem.split('_')[1]
    
    dd = Dataset(folder)
    case_name_short = dd.case_name.split('_cl_pr')[0]

    vtu_file = meshing_dir / case_name_short / '03_mesh' / (dd.case_name + '.vtu')
    
    dd = dd.assemble_mesh().assemble_surface(mesh_file=vtu_file)
    
    # Surface meshing file
    surf_dir = meshing_dir / case_name_short / '02_processed'
    
    surf = pv.read(surf_dir / (dd.case_name + '.vtp'))
    sac_mask = surf.point_arrays['Mask'] == 1
    sac = surf.extract_points(sac_mask, include_cells=True, adjacent_cells=True)
    sac = pv.PolyData(sac.points, faces=sac.cells)

    dd.create_sac_mask(sac=sac)

    sac_mask = dd.surf.point_arrays['SurfaceSacMask'] == 1
    sac_vol_mask = dd.mesh.point_arrays['SacMask'] == 1

    tec_file = sorted(folder.glob('*hemodynamics_w_full_all.tec'))[0]
    tec = pv.read(tec_file)

    fields = ['OSI', 'SPI', 'TAWSS', 'RRT']

    for f in fields:
        df.at[case_name, f] = tec.point_arrays[f][sac_mask].mean()

    tec_file_vol = sorted(folder.glob('*volumetric_metrics_full.tec'))[0]
    tec_vol = pv.read(tec_file_vol)

    fields = ['SPI', 'TAVEL']

    for f in fields:
        df.at[case_name, f] = tec_vol.point_arrays[f][sac_vol_mask].mean()

    # Lastly, get SCI and LSA
    tree = KDTree(surf.points)
    _, ii = tree.query(dd.surf.points)
    dd.surf.point_arrays['sac_zones'] = surf.point_arrays['sac_zones'][ii]
    dd.surf.point_arrays['SacMask'] = surf.point_arrays['Mask'][ii]

    dd.surf.point_arrays['sac_zones'][dd.surf.point_arrays['SacMask'] == 1] = 0

    dd.compute_wss_parameters(mask_key='sac_zones', sac_mask_key='SacMask', indices=np.arange(960)[::int(8)]) # was 32

    df.at[case_name, 'SCI'] = dd.SCI['SacMask']
    df.at[case_name, 'LSA'] = dd.LSA['SacMask']

    df.to_csv(df_out)