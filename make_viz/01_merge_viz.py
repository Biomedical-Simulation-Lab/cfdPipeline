""" Merge each pair into one video.

NOTE: Not priority for thesis.
"""

import numpy as np
from PIL import Image
from pathlib import Path 
from bsl import viz 
import sys 

def zoom_at(img, x, y, zoom):
    w, h = img.size
    zoom2 = zoom * 2
    img = img.crop((x - w / zoom2, y - h / zoom2, 
                    x + w / zoom2, y + h / zoom2))
    return img.resize((w, h), Image.LANCZOS)

def get_zoom_locations(seq, length):
    coords = np.array([(450,450), (335,390), (350,300)])
    if seq == 'seq_001':
        coords = coords[[0,1]]
    elif seq == 'seq_005':
        coords = coords[[1,2]]
    elif seq == 'seq_006':
        coords = coords[[2,0]]

    span = np.linspace(0,1,length,endpoint=True)
    out = [coords[0] + frac * (coords[1] - coords[0]) for frac in span]

    return out

if __name__ == "__main__":
    read_folder = Path(sys.argv[1])
    write_folder = Path(sys.argv[2])
    vid_folder = Path(sys.argv[3])

    seq = read_folder.parents[1].stem

    for ff in [write_folder, vid_folder]:
        if not ff.exists():
            ff.mkdir(parents=True)

    case_folders = sorted(read_folder.glob('art*'))
    case_folders = case_folders[::2]
    num_imgs = [len(sorted(ff.glob('*.png'))) for ff in case_folders]
    assert len(set(num_imgs)) == 1
    
    img_files = [sorted(ff.glob('*.png')) for ff in case_folders]
    img_files = dict(zip(case_folders, img_files))

    imgs_out = []

    zoom_coords = get_zoom_locations(seq, len(case_folders))

    # u_mag_imgs = [zoom_at(im, z[0], z[1], 1.4) for im, z in zip(u_mag_imgs, zoom_coords)]

    for ldx in range(num_imgs[0]):
        imgs_N = [img_files[ff][ldx] for ff in case_folders]
        imgs = [Image.open(im_f) for im_f in imgs_N]
        imgs = [zoom_at(im, zoom_coords[imdx][0], zoom_coords[imdx][1], 1.4) for imdx, im in enumerate(imgs)]
        imgs = [np.array(im) for im in imgs]


        img_concat = np.concatenate(imgs, axis=1)
        img_file_out = write_folder / ('img_{:04d}.png'.format(ldx))

        Image.fromarray(img_concat).save(img_file_out)

        imgs_out.append(img_file_out)
        
    viz.images_to_movie(imgs_out, vid_folder / (seq + '.mp4'), fps=30)

