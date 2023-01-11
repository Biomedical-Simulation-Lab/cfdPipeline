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

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    spec_data_out_folder = Path(sys.argv[2])
    spec_img_out_folder = Path(sys.argv[3])
    meshing_dir = Path(sys.argv[4])

    spec_out_file = spec_data_out_folder / (folder.stem + '.npz')

    for ff in [spec_img_out_folder, spec_data_out_folder]:
        if not ff.exists():
            ff.mkdir(parents=True, exist_ok=True)

    if not spec_out_file.exists():
        dd = Dataset(folder)
        case_name_short = dd.case_name.split('_cl')[0]

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

        sac_mask = dd.mesh.point_arrays['SacMask'] == 1

        dd.assemble_matrix(array_key='u_mag', quantity='u_mag', mask='SacMask')

        # Spectrograms
        dd.spectrogram(array_key='u_mag')
        spec_data = dd.spectrogram_data['u_mag']
        np.savez(spec_out_file, **dd.spectrogram_data['u_mag'])

    else:
        spec_data = np.load(spec_out_file)

    # bins = spec_data['bins']
    # freqs = spec_data['freqs']
    # S = spec_data['S']

    # S[S < -20] = -20 

    # # Plotting
    # plt.pcolormesh(bins, freqs, S, shading='gouraud')
    # plt.xlabel('Time (s)')
    # plt.ylabel('Freq (Hz)')

    # # if 'standard' in folder.stem:
    # #     neck_type = 'std'
    # # elif 'surge' in folder.stem:
    # #     neck_type = 'srg'
        
    # title = [folder.stem.split('_')[1], '{:.2f}, {:.2f}'.format(S.min(), S.max())]
    # title = ', '.join(title)

    # plt.title(title)
    # plt.savefig(spec_img_out_folder / (folder.stem + '.png'))
