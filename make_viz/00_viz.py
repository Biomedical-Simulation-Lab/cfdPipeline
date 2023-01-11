""" Generate spectrograms using BSL tools.

Note, surface meshing file is uploaded, so extracting sac from that.
"""

import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = '/scratch/s/steinman/macdo708/.config/matplotlib'

import pyvista as pv
from bsl.dataset import Dataset
import numpy as np 

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    viz_dir = Path(sys.argv[2])
    viz_type = str(sys.argv[3])
    step = int(sys.argv[4])
    meshing_dir = Path(sys.argv[5])
    cpos_file = Path(sys.argv[6])

    if not viz_dir.exists():
        viz_dir.mkdir(exist_ok=True, parents=True)

    dd = Dataset(folder)
    case_name_short = dd.case_name.split('_cl')[0]

    vtu_file = meshing_dir / case_name_short / '03_mesh' / (dd.case_name + '.vtu')
    # cpos_file = cpos_dir / (case_name_short + '.npz')
    cpos_data = np.load(cpos_file)

    dd = dd.assemble_mesh().assemble_surface(mesh_file=vtu_file)

    if viz_type == 'u_mag':
        dd.viz_contour(
            array_name='u_mag', 
            contour_val=0.7, 
            output_folder=viz_dir, 
            indices=np.arange(960)[::step], # was 8
            cpos=cpos_data['cpos'],
            parallel_scale=cpos_data['parallel_scale'],
            # clip_surf=clip_box,
            )

    if viz_type == 'wss_mag':
        dd.viz_surface(
            array_name='wss_mag',
            output_folder=viz_dir,
            indices=np.arange(960)[::step], # was 8,
            cpos=cpos_data['cpos'],
            parallel_scale=cpos_data['parallel_scale'],
            clim=[1,20],
            log_scale=True,
            cmap='rainbow',
        )

    if viz_type == 'qcriterion_nd':
        dd.viz_contour(
            array_name='qcriterion_nd', 
            contour_val=1.5, 
            output_folder=viz_dir, 
            indices=np.arange(960)[::step],
            cpos=cpos_data['cpos'],
            parallel_scale=cpos_data['parallel_scale'],
            color=None,
            scalars='u',
            # clip_surf=clip_box,
            )