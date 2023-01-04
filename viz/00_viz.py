""" Generate contour viz.
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

    cpos_d = Path('/scratch/s/steinman/macdo708/pereira_rupture_study_analysis/assets/cpos')

    if not viz_dir.exists():
        viz_dir.mkdir(exist_ok=True)

    dd = Dataset(folder)
    cpos_f = cpos_d / (dd.case_name + '.npz')
    cpos_data = np.load(cpos_f)
    mesh_file = folder.parents[1] / 'data' / (dd.case_name + '.h5')

    dd = dd.assemble_mesh().assemble_surface(mesh_file=mesh_file)

    if viz_type == 'u_mag':
        dd.viz_contour(
            array_name='u_mag', 
            contour_val=0.7, 
            output_folder=viz_dir, 
            indices=np.arange(len(dd))[::step], # was 8
            cpos=cpos_data['cpos'],
            parallel_scale=cpos_data['parallel_scale'],
            window_size=[768,768],
            )

    if viz_type == 'wss_mag':
        dd.viz_surface(
            array_name='wss_mag',
            output_folder=viz_dir,
            indices=np.arange(len(dd))[::step], # was 8,
            cpos=cpos_data['cpos'],
            parallel_scale=cpos_data['parallel_scale'],
            clim=[0.5,10],
            window_size=[768,768],
        )

    if viz_type == 'qcriterion_nd':
        dd.viz_contour(
            array_name='qcriterion_nd', 
            contour_val=1.5, 
            output_folder=viz_dir, 
            indices=np.arange(len(dd))[::step],
            cpos=cpos_data['cpos'],
            parallel_scale=cpos_data['parallel_scale'],
            color=None,
            scalars='u',
            window_size=[768,768],
            )