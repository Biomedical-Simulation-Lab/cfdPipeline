""" Generate spectrograms using BSL tools.
"""

import sys
import os
from pathlib import Path
os.environ['MPLCONFIGDIR'] = '/scratch/s/steinman/macdo708/.config/matplotlib'

import matplotlib.pyplot as plt 

import pyvista as pv
from bsl.dataset import Dataset
import numpy as np 

size = 12
plt.rc('font', size=size) #controls default text size
plt.rc('axes', titlesize=size) #fontsize of the title
plt.rc('axes', labelsize=size) #fontsize of the x and y labels
plt.rc('xtick', labelsize=size) #fontsize of the x tick labels
plt.rc('ytick', labelsize=size) #fontsize of the y tick labels
plt.rc('legend', fontsize=size) #fontsize of the legend

if __name__ == "__main__":
    folder = Path(sys.argv[1])
    spec_data_out_folder = Path(sys.argv[2])
    spec_img_out_folder = Path(sys.argv[3])

    spec_out_file = spec_data_out_folder / (folder.stem + '.npz')

    if not spec_out_file.exists():
        dd = Dataset(folder)
        case_name_short = dd.case_name.split('_mesh')[0]

        mesh_file = folder.parents[1] / 'data' / (dd.case_name + '.h5')
        surf_file = folder.parents[1] / 'data' / (dd.case_name + '.vtp')

        dd = dd.assemble_mesh().assemble_surface(mesh_file=mesh_file)
        
        surf = pv.read(surf_file)
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

    bins = spec_data['bins']
    freqs = spec_data['freqs']
    S = spec_data['S']

    S[S < -20] = -20 

    # Plotting
    fig, ax = plt.subplots(1,1, figsize=(4,4))
    ax.pcolormesh(bins, freqs, S, shading='gouraud')
    ax.set_xlabel('Time (s)', labelpad=-5)
    ax.set_ylabel('Freq (Hz)', labelpad=-10)
    ax.set_xticks([0, 0.9])
    ax.set_xticklabels(['0.0', '0.9'])
    ax.set_yticks([0, 500])
    ax.set_yticklabels(['0', '500'])


    if spec_data_out_folder.parent.stem == 'standard':
        flow_type = 'std'
    elif spec_data_out_folder.parent.stem == 'mean_flow':
        flow_type = 'mean'
        
    # title = [folder.stem.split('_cl')[0], flow_type, '{:.2f}, {:.2f}'.format(S.min(), S.max())]
    # title = ', '.join(title)
    title = folder.stem.split('_cl')[0].split('art_')[1]

    ax.set_title(title)
    plt.tight_layout
    plt.savefig(spec_img_out_folder / (folder.stem + '.png'))#, transparent=True)
