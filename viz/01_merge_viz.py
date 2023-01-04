""" Merge each pair into one video.
"""

import numpy as np
from PIL import Image
from pathlib import Path 
from bsl import viz 
import sys 

if __name__ == "__main__":
    # viz_type = str(sys.argv[1])

    # viz_folder = Path("/scratch/s/steinman/macdo708/surge_cfd_analysis/viz")
    # out_folder = Path("/scratch/s/steinman/macdo708/surge_cfd_analysis/viz_merge")
    # vid_folder = Path("/scratch/s/steinman/macdo708/surge_cfd_analysis/vids")

    read_folder = Path(sys.argv[1])
    write_folder = Path(sys.argv[2])
    vid_folder = Path(sys.argv[3])
    # viz_type = Path(sys.argv[4])

    # viz_folder = viz_folder / viz_type
    # out_folder = out_folder / viz_type
    # vid_folder = vid_folder / viz_type

    for ff in [write_folder, vid_folder]:
        if not ff.exists():
            ff.mkdir(parents=True)

    case_folders = sorted(read_folder.glob('art*'))

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
        assert case in pair[0].stem
        
        case_out = write_folder / case 

        if not case_out.exists():
            case_out.mkdir()

        img_files_std = sorted(pair[0].glob('*.png'))
        img_files_srg = sorted(pair[1].glob('*.png'))

        imgs_out = []
        for im_st, im_sr in zip(img_files_std, img_files_srg):
            img_file_out = case_out / (im_st.stem.split('cy2_')[1] + '.png')
            
            img_std = np.array(Image.open(im_st))
            img_srg = np.array(Image.open(im_sr))

            img_concat = np.concatenate([img_std, img_srg], axis=1)
            Image.fromarray(img_concat).save(img_file_out)

            imgs_out.append(img_file_out)
            
        viz.images_to_movie(imgs_out, vid_folder / (case + '.mp4'), fps=30)