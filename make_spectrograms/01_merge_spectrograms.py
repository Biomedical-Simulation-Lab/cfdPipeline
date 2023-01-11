""" Merge each pair into one img.
"""

import numpy as np
from PIL import Image
from pathlib import Path 

if __name__ == "__main__":
    spec_folder = Path("/scratch/s/steinman/macdo708/surge_cfd_analysis/spectrogram_imgs")
    out_folder = spec_folder / 'merged'

    if not out_folder.exists():
        out_folder.mkdir(parents=True)

    case_folders = sorted(spec_folder.glob('art*'))

    assert len(case_folders) % 2 == 0

    # Split list into pairs
    chunk_size = 2
    paired_list = [case_folders[i:i+chunk_size] for i in range(0, len(case_folders), chunk_size)]

    # Double check all pairs match
    for pair in paired_list:
        assert pair[0].stem.split('_')[1] == pair[1].stem.split('_')[1]

    case_names = sorted(list(set([x.stem.split('_')[1] for x in case_folders])))

    # Now grab frames, merge horizontally, then make video
    for case, pair in zip(case_names, paired_list):
        img_file_out = out_folder / (case + '.png')

        img_std = np.array(Image.open(pair[0]))
        img_srg = np.array(Image.open(pair[1]))

        img_concat = np.concatenate([img_std, img_srg], axis=1)
        Image.fromarray(img_concat).save(img_file_out)
